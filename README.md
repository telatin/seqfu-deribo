# deribo 

A fast k-mer/minimizer based rRNA read filter, written in Nim. This is the
**first viable product** milestone toward the full design in

## What it does

`deribo` classifies reads as "matches a reference database" or not, by
counting canonical (strand-invariant) k-mers/minimizers shared between a read
and an indexed set of reference FASTA sequences -- no Smith-Waterman/Myers
alignment yet (see "Scope" below). This is deliberately closer to a
deacon-style fast minimizer classifier than to SortMeRNA's seed+align
pipeline; see the brainstorming discussion that led here.

## Build

```sh
nimble build          # -> bin/deribo
nimble test           # unit tests, including a real-fixture check against
                       # sortmerna/test_data (run from the repo root)
```

Requires Nim >= 2.2, a C compiler, and `gzip` on PATH (used for gzip
input/output instead of a C zlib binding, for v0.1 simplicity). If `pigz` is
on PATH it is used instead, which is much faster for `--gzip-out`.

## Usage

The reference **database argument accepts either a single FASTA file or a
directory of FASTA files** -- every `*.fasta`/`*.fa`/`*.fna`/`*.ffn`/`*.frn`/
`*.fsa` (optionally `.gz`) file in the directory is indexed, in sorted order
for determinism. This is required to work with an arbitrary local checkout
such as `sortmerna-4.3.4/` (SILVA/Rfam default databases) or
`sortmerna/test_data/refs/`.

```sh
# Build an index from a directory of reference FASTA files
bin/deribo index --db INPUTDIR/ --output deribo.idx
# (repeat --db to combine multiple files/directories in one index)

# Inspect what got built
bin/deribo inspect --index deribo.idx

# Single-end filtering
bin/deribo filter --index deribo.idx --reads1 reads.fastq.gz \
  --output out/sample --mappings

# Two-file paired-end filtering (either/both/independent pair policy)
bin/deribo filter --index deribo.idx \
  --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \
  --pair-policy either --output out/sample --threads 8
```

`--threads` (default 4) sets the classifying worker threads. Reads are
processed in batches and written back in input order, so every output file
is byte-identical at any thread count. With `--gzip-out`, each output file
is compressed by its own `pigz`/`gzip` process while filtering runs.

Outputs: `<prefix>_accepted[_1/_2].fast(a|q)`, `<prefix>_rejected[_1/_2].fast(a|q)`,
`<prefix>_summary.json` (machine-readable), and optionally
`<prefix>_mappings.tsv` (read name, hits, seed count, covered bases, longest
covered run, best-matching reference).

Tunables: `-k` (seed length, 1..32, default 21), `-w` (minimizer window;
1 = every k-mer / most sensitive; >1 thins to one minimizer per window,
recommended for large multi-file directories like `sortmerna-4.3.4/`).
A read is accepted only if all of these hold:
`--min-hits` (seeds found in the index, default 2), `--min-coverage`
(fraction of read bases covered by hit seeds, default 0.25) and `--min-run`
(longest contiguous covered stretch in bp, default 0 = off).

The defaults were calibrated on `testset/`. The original rule (k=17, hits >= 2)
flagged 39-56% of pairs as rRNA against SortMeRNA's 25-42%, because one 18 bp
exact match, or two chance 17-mer hits against a ~130 Mbp database, was
enough. With k=21 and coverage >= 0.25, totals are within 0.7% of SortMeRNA
and the Jaccard overlap of rRNA calls is 0.940 / 0.972 (was 0.641 / 0.752).
Passing `--min-coverage 0` with a k=17 index reproduces the old behaviour.

