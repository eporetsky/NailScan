# InterPro member databases and NailScan

Nail accepts **HMMER3 profile HMM** (ASCII) only. Below is the status of each database present in InterProScan's `data/` (from your InterProScan 5.77 install and the [EBI iprscan FTP](https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/)).

## Compatible with Nail (single or concatenated HMMER3 HMM)

| Database   | Main HMM file (after unpack) | In config.json | Notes |
|-----------|-------------------------------|----------------|--------|
| **Pfam**  | pfam_a.hmm                    | yes            | Already supported. |
| **NCBIfam** | ncbifam.hmm                 | yes            | Same format as Pfam. |
| **Superfamily** | hmmlib_1.75              | yes            | Single concatenated HMMER3 file. |
| **AntiFam** | AntiFam.hmm                 | yes            | Spurious ORF / contamination models. |
| **HAMAP** | hamap.hmm.lib                | yes            | Single concatenated HMMER3. |
| **SFLD**  | sfld.hmm                     | yes            | Structure-Function Linkage. |
| **PIRSF** | sf_hmm_all                   | yes            | Superfamily HMMs. |
| **PIRSR** | sr_hmm_all                   | yes            | PIRSR HMMs (2025_05 to match InterProScan 5.77). |
| **PANTHER** | binHmm (in famhmm/)       | yes            | Large (~3.3 GB); HMMER3 ASCII. |
| **CATH/Gene3D** | gene3d_main.hmm          | yes            | In FTP as cath; Gene3D models. |

All of the above are in config.json with FTP URLs and hmm_name. Run ./install.sh (or NAILSCAN_DBS="antifam" ./install.sh etc.) to download and unpack into data/<dbname>/.

## Not compatible with Nail (as-is)

| Database   | Reason |
|-----------|--------|
| **Prosite** | Profiles are .prf (Prosite format). InterProScan uses PfsearchV3. No standard prf to HMM conversion; not supported. |
| **SMART**  | Distribution is HMMER2 (smart.HMMs). Nail expects HMMER3. Could be converted with hmmconvert (HMMER2 to 3) and then used; not in config by default. |
| **CDD**    | Uses RPS-BLAST, not profile HMM. |
| **PRINTS** | Fingerprint-based. |
| **Phobius / TMHMM** | Signal peptide / transmembrane predictors; different pipeline. |

## InterPro and GO lookup (--iprlookup, --goterms)

To add InterPro and/or Gene Ontology columns to the TSV output, you need mapping files in `data/`:

1. **interpro2go** — InterPro ID to GO terms (from InterPro current release).
2. **signature2interpro.tsv** — Member DB signature (e.g. Pfam ACC) to InterPro ID, built from **interpro.xml.gz**.

**One-time setup:** run from the repo root:

```bash
./scripts/download_interpro_mappings.sh
```

This downloads `interpro2go` and `interpro.xml.gz` into `data/`, then builds `data/signature2interpro.tsv` from the XML. Use `--force` to re-download and rebuild.

**Sources:**
- interpro2go: https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go
- interpro.xml.gz: https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz

Then run single or batch with `--iprlookup` and/or `--goterms` to add the extra columns.

## FTP base and adding a new DB

- Base: https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/
- Each DB has a directory; inside: <dbname>-<version>.tar.gz and .md5.
- Add an entry to config.json with url, md5_url, and hmm_name (basename of the main HMM file after unpack). If the tarball layout differs, add a case in place_db_files() in install.sh.

## Reference

- InterProScan properties (HMM paths): interproscan-5.77-108.0-64-bit/interproscan-5.77-108.0/interproscan.properties
- EBI iprscan 6.0: https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/6.0/
