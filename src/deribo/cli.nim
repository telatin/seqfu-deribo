## Command-line surface, built with cligen. Three subcommands for v0.1:
## `index`, `filter`, `inspect` (map/validate follow later milestones per
## DERIBO_PLAN.md / NIM_SORTMERNA_PLAN.md §3).

import cligen
import ./pipeline
import ./classify

proc cmdIndex(db: seq[string] = @[], k = 21, w = 1, output = "deribo.idx") =
  ## Build a deribo seed index from one or more reference FASTA files or
  ## directories of FASTA files (e.g. `--db sortmerna-4.3.4/` indexes every
  ## *.fasta file found there). k = seed length (1..32; default 21 -- at
  ## k=17 chance matches against a ~100 Mbp rRNA database are common enough
  ## to call ~10% of unrelated read pairs). w = minimizer
  ## window; w<=1 keeps every canonical k-mer (most sensitive, most memory);
  ## w>1 thins to one minimizer per w-window (deacon-style), recommended for
  ## large multi-gigabase reference directories.
  runIndex(db, k, w, output)

proc cmdFilter(index: string, reads1: string, reads2 = "", output = "deribo_out",
                minHits = DefaultMinHits, minCoverage = DefaultMinCoverage,
                minRun = DefaultMinRun,
                pairPolicy = "either", gzipOut = false,
                mappings = false, summaryJson = "", threads = 4) =
  ## Classify reads against a deribo index and split them into accepted
  ## (database-like) and rejected FASTA/FASTQ output. Pass --reads2 for
  ## two-file paired-end input. --pair-policy is one of: either, both,
  ## independent. A read is accepted only if all three hold: at least
  ## --min-hits seeds found, at least --min-coverage of its bases covered by
  ## hit seeds, and a contiguous covered stretch of at least --min-run bp.
  ## --threads sets the classifying worker threads (output is identical at
  ## any count); with --gzip-out each output file is also compressed by its
  ## own pigz (if installed) using up to that many threads.
  runFilter(index, reads1, reads2, output,
            thresholds(minHits, minCoverage, minRun), pairPolicy, gzipOut,
            mappings, summaryJson, threads)

proc cmdInspect(index: string) =
  ## Print metadata (k, w, reference count, seed count, size) about a
  ## deribo index file.
  runInspect(index)

proc runCli*() =
  dispatchMulti(["multi", cmdName = "deribo"],
                 [cmdIndex, cmdName = "index"],
                 [cmdFilter, cmdName = "filter"],
                 [cmdInspect, cmdName = "inspect"])
