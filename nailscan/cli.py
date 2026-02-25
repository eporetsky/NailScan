"""Click CLI: single and batch modes; mirrors nailscan.single.sh / nailscan.batch.sh."""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import click


def _default_hmm() -> Path:
    return Path.cwd() / "data" / "pfam" / "pfam_a.hmm"


def _default_results() -> Path:
    return Path.cwd() / "results"


def _check_nail() -> None:
    if shutil.which("nail") is None:
        click.echo("Error: nail not in PATH. Activate env: conda activate nailscan", err=True)
        sys.exit(1)


def _build_hmm_map(hmm_path: Path) -> dict[str, tuple[str, str]]:
    """Parse HMM file for NAME -> (ACC, DESC). Same logic as make_pfam_map.sh."""
    name, acc, desc = "", "", ""
    out: dict[str, tuple[str, str]] = {}
    with open(hmm_path) as f:
        for line in f:
            if line.startswith("NAME"):
                if name:
                    out[name] = (acc, desc)
                parts = line.split(None, 1)
                name = parts[1].strip() if len(parts) > 1 else ""
                acc, desc = "", ""
            elif line.startswith("ACC"):
                parts = line.split(None, 1)
                acc = parts[1].strip() if len(parts) > 1 else ""
            elif line.startswith("DESC"):
                desc = line[3:].strip()
        if name:
            out[name] = (acc, desc)
    return out


def _run_nail(hmm: Path, fasta: Path, work_dir: Path, threads: int) -> Path:
    tbl = work_dir / "results.tbl"
    cmd = [
        "nail", "search",
        "-t", str(threads),
        "--tbl-out", str(tbl),
        str(hmm),
        str(fasta),
    ]
    subprocess.run(cmd, cwd=work_dir, check=True)
    if not tbl.exists():
        click.echo("Error: nail did not produce results.tbl", err=True)
        sys.exit(1)
    return tbl


HEADER_NO_MAP = "target\tNAME\ttarget_start\ttarget_end\tquery_start\tquery_end\tscore\tbias\tevalue\tcell_frac"
HEADER_MAP = "target\tNAME\tACC\tDESC\ttarget_start\ttarget_end\tquery_start\tquery_end\tscore\tbias\tevalue\tcell_frac"


@click.group()
def cli() -> None:
    """NailScan: fast HMM-profile search (nail) for InterPro-style annotation."""
    pass


@cli.command()
@click.option("-f", "--fasta", "fasta_path", required=True, type=click.Path(path_type=Path, exists=True),
              help="Input FASTA file.")
@click.option("-h", "hmm_path", "--hmm", default=None, type=click.Path(path_type=Path, exists=True),
              help="HMM file path (default: data/pfam/pfam_a.hmm from cwd).")
@click.option("-o", "--output", "output_dir", default=None, type=click.Path(path_type=Path),
              help="Output directory (default: results/).")
@click.option("-t", "--threads", default=None, type=int,
              help="Threads (default: nproc or 16).")
def single(
    fasta_path: Path,
    hmm_path: Path | None,
    output_dir: Path | None,
    threads: int | None,
) -> None:
    """Run nail on one FASTA; write one TSV with optional ACC/DESC from HMM map."""
    _check_nail()
    hmm_path = hmm_path or _default_hmm()
    output_dir = output_dir or _default_results()
    if not hmm_path.exists():
        click.echo(f"Error: HMM not found: {hmm_path}", err=True)
        sys.exit(1)
    try:
        threads = threads or int(subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip() or "16")
    except Exception:
        threads = 16

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = fasta_path.resolve()
    hmm_path = hmm_path.resolve()

    stem = fasta_path.stem
    out_tsv = output_dir / f"{stem}.pfam.tsv"

    with tempfile.TemporaryDirectory() as work:
        work_path = Path(work)
        tbl = _run_nail(hmm_path, fasta_path, work_path, threads)
        map_file = hmm_path.with_suffix(hmm_path.suffix + ".map")
        if not map_file.exists():
            hmm_map = _build_hmm_map(hmm_path)
            with open(map_file, "w") as mf:
                for name, (acc, desc) in hmm_map.items():
                    mf.write(f"{name}\t{acc}\t{desc}\n")
        else:
            hmm_map = {}
            with open(map_file) as mf:
                for line in mf:
                    parts = line.strip().split("\t", 2)
                    if len(parts) >= 3:
                        hmm_map[parts[0]] = (parts[1], parts[2])
                    elif len(parts) == 2:
                        hmm_map[parts[0]] = (parts[1], "")

        with open(out_tsv, "w") as out:
            out.write(HEADER_MAP if hmm_map else HEADER_NO_MAP)
            out.write("\n")
            with open(tbl) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 10:
                        continue
                    target, query = fields[0], fields[1]
                    if hmm_map:
                        acc, desc = hmm_map.get(query, ("", ""))
                        out.write(f"{target}\t{query}\t{acc}\t{desc}\t" + "\t".join(fields[2:10]) + "\n")
                    else:
                        out.write("\t".join(fields[:10]) + "\n")

    click.echo(f"Wrote: {out_tsv}")


@cli.command()
@click.option("-f", "--fasta-dir", "fasta_dir", default=None, type=click.Path(path_type=Path, exists=True, file_okay=False),
              help="Directory of .fa/.faa/.fasta files (default: fasta/).")
@click.option("-h", "hmm_path", "--hmm", default=None, type=click.Path(path_type=Path, exists=True),
              help="HMM file path (default: data/pfam/pfam_a.hmm from cwd).")
@click.option("-o", "--output", "output_dir", default=None, type=click.Path(path_type=Path),
              help="Output directory for per-genome .tsv.gz (default: results/pfam/).")
@click.option("-t", "--threads", default=None, type=int,
              help="Threads (default: nproc or 16).")
def batch(
    fasta_dir: Path | None,
    hmm_path: Path | None,
    output_dir: Path | None,
    threads: int | None,
) -> None:
    """Combine FASTAs in a directory, run one nail search, split to per-genome .tsv.gz."""
    _check_nail()
    fasta_dir = fasta_dir or Path.cwd() / "fasta"
    hmm_path = hmm_path or _default_hmm()
    output_dir = output_dir or Path.cwd() / "results" / "pfam"
    if not fasta_dir.exists():
        click.echo(f"Error: Fasta directory not found: {fasta_dir}", err=True)
        sys.exit(1)
    if not hmm_path.exists():
        click.echo(f"Error: HMM not found: {hmm_path}", err=True)
        sys.exit(1)
    try:
        threads = threads or int(subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip() or "16")
    except Exception:
        threads = 16

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    hmm_path = hmm_path.resolve()
    fasta_dir = fasta_dir.resolve()

    suffixes = (".fa", ".faa", ".fasta")
    fasta_files = []
    for p in sorted(fasta_dir.iterdir()):
        if p.is_file() and (p.suffix == ".fa" or p.suffix == ".faa" or p.name.endswith(".fasta")):
            fasta_files.append(p)

    if not fasta_files:
        click.echo(f"No .fa/.faa/.fasta files in {fasta_dir}", err=True)
        sys.exit(1)

    with tempfile.TemporaryDirectory() as work:
        work_path = Path(work)
        batch_fa = work_path / "batch.fa"
        with open(batch_fa, "w") as bf:
            for fa in fasta_files:
                genome_id = fa.stem
                if fa.name.endswith(".fasta"):
                    genome_id = fa.name[:-7]
                with open(fa) as f:
                    for line in f:
                        if line.startswith(">"):
                            bf.write(f">{genome_id}|{line[1:]}")
                        else:
                            bf.write(line)

        n_seqs = sum(1 for _ in open(batch_fa) if _.startswith(">"))
        click.echo(f"Combined {n_seqs} sequences from {len(fasta_files)} file(s).")
        tbl = _run_nail(hmm_path, batch_fa, work_path, threads)

        map_file = hmm_path.with_suffix(hmm_path.suffix + ".map")
        if not map_file.exists():
            hmm_map = _build_hmm_map(hmm_path)
            with open(map_file, "w") as mf:
                for name, (acc, desc) in hmm_map.items():
                    mf.write(f"{name}\t{acc}\t{desc}\n")
        else:
            hmm_map = {}
            with open(map_file) as mf:
                for line in mf:
                    parts = line.strip().split("\t", 2)
                    if len(parts) >= 3:
                        hmm_map[parts[0]] = (parts[1], parts[2])
                    elif len(parts) == 2:
                        hmm_map[parts[0]] = (parts[1], "")

        per_genome: dict[str, list[str]] = {}
        with open(tbl) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue
                target, query = fields[0], fields[1]
                idx = target.find("|")
                if idx < 0:
                    continue
                genome_id = target[:idx]
                gene_id = target[idx + 1:]
                if hmm_map:
                    acc, desc = hmm_map.get(query, ("", ""))
                    row = f"{gene_id}\t{query}\t{acc}\t{desc}\t" + "\t".join(fields[2:10])
                else:
                    row = f"{gene_id}\t" + "\t".join(fields[1:10])
                if genome_id not in per_genome:
                    per_genome[genome_id] = []
                per_genome[genome_id].append(row)

        for genome_id, rows in per_genome.items():
            tsv_path = output_dir / f"{genome_id}.tsv.gz"
            with gzip.open(tsv_path, "wt") as zf:
                zf.write(HEADER_MAP if hmm_map else HEADER_NO_MAP)
                zf.write("\n")
                zf.write("\n".join(rows))
                if rows:
                    zf.write("\n")

    click.echo(f"Results in {output_dir}/")