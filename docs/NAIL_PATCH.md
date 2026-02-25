# Superfamily HMM format fix (COMPO lines)

## Problem

The Superfamily database (`hmmlib_1.75`) is distributed in **HMMER2 binary format**
(`HMMER3/b`), which `nail` cannot read.  `hmmconvert` converts it to HMMER3/f,
but the resulting file lacks `COMPO` lines — these are optional per the HMMER3
spec, but `nail`'s internal `from_p7hmm` parser requires them.

### Why COMPO matters for nail

`from_p7hmm` uses a small state machine to parse each model's header section:

```
HMM   A C D E F G H I K L M N P Q R S T V W Y        ← column header
      m->m m->i m->d i->m i->i d->m d->d              ← transition header
  COMPO   <20 values>      ← composition (optional per spec, required by nail)
          <20 values>      ← position-0 insert emissions
          <7 values>       ← position-0 transition probabilities
  1   <20 values>          ← model position 1  ...
```

The parser only advances from the **MatchEmissions** state into **InsertEmissions**
when it sees the `COMPO` token.  Without it, the parser stays stuck on every
subsequent line, silently treating them all as no-ops.  The model body (match
scores, insert emissions, transition probabilities) is never read — every model
parses as a zero-emission profile — so mmseqs finds no seeds and nail returns
0 hits.

## Fix

`install.sh` adds a `COMPO` line to every model that lacks one, immediately
after the transition header (`m->m m->i m->d ...`).  The COMPO values are taken
from the first numeric data line that is already present (the position-0 insert
emissions, which represent the standard HMMER3 null-model amino-acid background
frequencies — semantically identical to what COMPO should contain).

This is a **file-level fix only** — no changes to the `nail` source code.

## How it is applied

The relevant section of `install.sh`:

1. Backs up the original HMMER2 file as `hmmlib_1.75.raw`.
2. Converts HMMER2 → HMMER3 using `hmmconvert`.
3. Inserts `COMPO` lines via an embedded Python script.
4. Writes the result back to `hmmlib_1.75` (the path used by `config.json`).

On subsequent `install.sh` runs the script checks for `hmmlib_1.75.raw` **and**
at least one `  COMPO` line in `hmmlib_1.75`; if both are present the step is
skipped.  If the backup exists but COMPO lines are missing (upgrade from an older
install), the conversion + COMPO-insertion is re-run automatically.
