#!/usr/bin/env python3
"""
Filter HMMER3 HMM file: remove models with model_length == 1 (LENG  1) and log them.
Nail and some other tools cannot handle single-node profile HMMs; InterProScan/HMMER3 can.
Reads from stdin or first arg, writes to second arg (or overwrites first arg with a temp file).
Logs removed model names to install_log_path (default: data/install.log next to repo root).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Match LENG followed by optional whitespace and digits (capture the number)
LENG_RE = re.compile(r"LENG\s+(\d+)", re.IGNORECASE)


def find_repo_root() -> Path:
    """Script lives in scripts/; repo root is parent of scripts."""
    return Path(__file__).resolve().parent.parent


def parse_record_header(record_lines: list[str]) -> tuple[str | None, int | None]:
    """Extract NAME and LENG from the first part of a HMMER3 record. Returns (name, length)."""
    name = None
    leng = None
    for line in record_lines:
        if line.startswith("NAME"):
            parts = line.split(None, 1)
            name = parts[1].strip().split()[0] if len(parts) > 1 else None
        m = LENG_RE.search(line)
        if m:
            try:
                leng = int(m.group(1))
            except ValueError:
                pass
            break  # LENG comes after NAME
    return name, leng


def filter_hmm_file(
    input_path: Path,
    output_path: Path,
    install_log_path: Path,
    min_length: int = 2,
) -> int:
    """
    Stream through input_path; write records with LENG >= min_length to output_path.
    Append removed model names (LENG < min_length) to install_log_path.
    Returns number of models removed.
    """
    removed: list[str] = []
    buffer: list[str] = []
    seen_end = False

    with open(input_path, "r") as f:
        for line in f:
            if line.rstrip() == "//":
                if buffer:
                    name, leng = parse_record_header(buffer)
                    if leng is not None and leng < min_length:
                        removed.append(name if name else "(no name)")
                    else:
                        # Keep this record: write buffer + "//"
                        pass  # we write below when we have a kept buffer
                    buffer.append(line)
                    # Flush kept records to output (we'll open output after first pass to know removed count, or we can stream)
                    buffer = []
                seen_end = True
                continue
            buffer.append(line)

        if buffer and not (buffer[-1].rstrip() == "//"):
            # Incomplete last record; treat as keep
            pass

    # Second pass: write kept records to output
    buffer = []
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.rstrip() == "//":
                if buffer:
                    buffer.append(line)
                    name, leng = parse_record_header(buffer)
                    if leng is not None and leng < min_length and name:
                        pass  # skip writing
                    else:
                        fout.writelines(buffer)
                    buffer = []
                continue
            buffer.append(line)
        if buffer:
            name, leng = parse_record_header(buffer)
            if not (leng is not None and leng < min_length):
                fout.writelines(buffer)

    if removed:
        install_log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(install_log_path, "a") as log:
            log.write(
                f"[filter_length1_hmms] Removed {len(removed)} HMM model(s) with model_length=1 "
                f"(incompatible with nail) from {input_path.name}: {', '.join(removed)}\n"
            )

    return len(removed)


def main() -> None:
    repo = find_repo_root()
    if len(sys.argv) < 2:
        print("Usage: filter_length1_hmms.py <hmm_file> [output_file] [install.log path]", file=sys.stderr)
        sys.exit(1)
    input_path = Path(sys.argv[1]).resolve()
    out_arg = sys.argv[2].strip() if len(sys.argv) > 2 else ""
    out_arg = Path(out_arg).resolve() if out_arg else None
    install_log = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else repo / "data" / "install.log"

    if not input_path.is_file():
        print(f"Error: not a file: {input_path}", file=sys.stderr)
        sys.exit(1)

    # If writing in-place, use a temp file then replace (avoid read/write same file).
    if not out_arg or out_arg == input_path:
        output_path = input_path.parent / (input_path.name + ".filtered.tmp")
        replace_in_place = True
    else:
        output_path = out_arg
        replace_in_place = False

    n_removed = filter_hmm_file(input_path, output_path, install_log)
    if replace_in_place:
        output_path.replace(input_path)
    if n_removed > 0:
        print(f"Removed {n_removed} HMM model(s) with model_length=1; log appended to {install_log}")
    if not replace_in_place:
        print(f"Wrote filtered HMM to {output_path}")


if __name__ == "__main__":
    main()
