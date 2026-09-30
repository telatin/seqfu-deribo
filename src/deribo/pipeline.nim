## Orchestrates the `index`, `filter`, and `inspect` commands.
##
## `filter` reads input in fixed-size batches. Each batch is classified by
## `--threads` worker threads (plain std threads, one contiguous slice each)
## while the main thread writes the previous batch and reads the next one,
## so output order -- and therefore every output file -- is identical at any
## thread count.

import std/[os, json, times, strutils, strformat, typedthreads]
import ./fastx
import ./index
import ./classify

type
  PipelineError* = object of CatchableError

  PairPolicy* = enum
    ppEither, ppBoth, ppIndependent

proc parsePairPolicy*(s: string): PairPolicy =
  case s.toLowerAscii
  of "either": ppEither
  of "both": ppBoth
  of "independent": ppIndependent
  else: raise newException(PipelineError,
    "unknown --pair-policy '" & s & "' (expected either, both, or independent)")

# --- index command -----------------------------------------------------------

proc runIndex*(db: seq[string], k, w: int, output: string) =
  if db.len == 0:
    raise newException(PipelineError, "no reference given: pass one or more --db FILE_OR_DIRECTORY")
  if k < 1 or k > 32:
    raise newException(PipelineError, "-k must be between 1 and 32 (got " & $k & ")")
  let t0 = epochTime()
  let files = collectFastaFiles(db)
  stderr.writeLine &"[deribo index] {files.len} FASTA file(s) resolved from {db.len} input path(s)"
  let idx = buildIndex(files, k, w)
  let buildElapsed = epochTime() - t0
  saveIndex(idx, output)
  let totalElapsed = epochTime() - t0
  let fileSize = getFileSize(output)
  stderr.writeLine &"[deribo index] references: {idx.refs.len}  total bases: {idx.totalBases}"
  stderr.writeLine &"[deribo index] distinct seeds (k={idx.k}, w={idx.w}): {idx.table.len}"
  stderr.writeLine &"[deribo index] build time: {buildElapsed:.2f}s  total time: {totalElapsed:.2f}s"
  stderr.writeLine &"[deribo index] wrote {output} ({fileSize} bytes)"

# --- inspect command ----------------------------------------------------------

proc runInspect*(indexPath: string) =
  let idx = loadIndex(indexPath)
  echo &"index:          {indexPath}"
  echo &"k:              {idx.k}"
  echo &"w (minimizer):  {idx.w}"
  echo &"references:     {idx.refs.len}"
  echo &"total bases:    {idx.totalBases}"
  echo &"distinct seeds: {idx.table.len}"
  echo &"file size:      {getFileSize(indexPath)} bytes"
  let n = min(5, idx.refs.len)
  if n > 0:
    echo "first reference(s):"
    for i in 0 ..< n:
      echo &"  [{i}] {idx.refs[i].name}  ({idx.refs[i].length} bp)"

# --- filter command ------------------------------------------------------------

type
  RunStats = object
    inputUnits: int      ## reads (single-end) or pairs (paired-end)
    acceptedUnits: int
    rejectedUnits: int
    mate1Accepted, mate1Rejected: int
    mate2Accepted, mate2Rejected: int

proc mappingLine(name: string, res: ClassifyResult, idx: DeriboIndex, mate: string): string =
  let refName = if res.bestRef >= 0: idx.refs[res.bestRef].name else: "*"
  let status = if res.accepted: "accepted" else: "rejected"
  &"{name}\t{mate}\t{status}\t{res.hits}\t{res.seedCount}\t{res.coveredBases}\t{res.longestRun}\t{refName}"

const BatchUnits = 16384
  ## Reads (single-end) or pairs (paired-end) per batch.

type
  Batch = object
    ## `recs` holds one record per unit, or mate 1 / mate 2 interleaved
    ## (2*i, 2*i + 1) for pairs; `res` is parallel to it. Both are reused
    ## across batches so record strings keep their capacity.
    recs: seq[FastxRecord]
    res: seq[ClassifyResult]
    n: int                 ## records in use

  ClassifyJob = object
    idx: ptr DeriboIndex
    t: Thresholds
    batch: ptr Batch
    first, last: int

proc classifyRange(job: ClassifyJob) {.thread.} =
  var scratch: ClassifyScratch
  for i in job.first ..< job.last:
    job.batch.res[i] = classifySequence(job.idx[], job.batch.recs[i].sequence,
                                        job.t, scratch)

proc startClassify(workers: var seq[Thread[ClassifyJob]], idx: ptr DeriboIndex,
                   t: Thresholds, batch: ptr Batch) =
  ## Splits `batch` across `workers` (already sized) and starts them.
  batch.res.setLen(batch.n)
  let per = (batch.n + workers.len - 1) div workers.len
  for w in 0 ..< workers.len:
    let first = min(batch.n, w * per)
    let last = min(batch.n, first + per)
    createThread(workers[w], classifyRange,
                 ClassifyJob(idx: idx, t: t, batch: batch, first: first, last: last))

proc readBatch(r1: var FastxReader, r2: var FastxReader, paired: bool,
               reads1, reads2: string, batch: var Batch) =
  let perUnit = if paired: 2 else: 1
  if batch.recs.len < BatchUnits * perUnit:
    batch.recs.setLen(BatchUnits * perUnit)
  batch.n = 0
  while batch.n < BatchUnits * perUnit:
    let got1 = r1.next(batch.recs[batch.n])
    if paired:
      let got2 = r2.next(batch.recs[batch.n + 1])
      if got1 != got2:
        raise newException(PipelineError,
          "paired input files have different numbers of records: " & reads1 & " vs " & reads2)
    if not got1: break
    batch.n += perUnit

proc runFilterBatches(idx: DeriboIndex, reads1, reads2, outPrefix: string, t: Thresholds,
                      policy: PairPolicy, gzipOut: bool, mappingsOut: File,
                      threads: int): RunStats =
  let paired = reads2.len > 0
  var r1 = openFastx(reads1)
  var r2: FastxReader
  if paired: r2 = openFastx(reads2)
  defer:
    close(r1)
    if paired: close(r2)
  let fmt = r1.format
  let ext = defaultExt(fmt)
  # Writers: [accepted_1, rejected_1, accepted_2, rejected_2] (paired) or
  # [accepted, rejected] (single-end).
  var writers: seq[FastxWriter]
  if paired:
    for name in ["_accepted_1", "_rejected_1", "_accepted_2", "_rejected_2"]:
      writers.add(openFastxWriter(outPrefix & name & ext, fmt, gzipOut, threads))
  else:
    for name in ["_accepted", "_rejected"]:
      writers.add(openFastxWriter(outPrefix & name & ext, fmt, gzipOut, threads))
  defer:
    for w in writers.mitems: close(w)

  # A template, not a closure: a closure would copy the captured `idx`.
  template emit(batch: Batch, stats: var RunStats) =
    if paired:
      for i in countup(0, batch.n - 1, 2):
        let (rec1, rec2) = (batch.recs[i], batch.recs[i + 1])
        let (res1, res2) = (batch.res[i], batch.res[i + 1])
        if not mappingsOut.isNil:
          mappingsOut.writeLine(mappingLine(rec1.name, res1, idx, "1"))
          mappingsOut.writeLine(mappingLine(rec2.name, res2, idx, "2"))
        inc stats.inputUnits
        case policy
        of ppEither, ppBoth:
          let pairAccepted = if policy == ppEither: (res1.accepted or res2.accepted)
                              else: (res1.accepted and res2.accepted)
          if pairAccepted:
            inc stats.acceptedUnits
            writers[0].write(rec1); writers[2].write(rec2)
          else:
            inc stats.rejectedUnits
            writers[1].write(rec1); writers[3].write(rec2)
        of ppIndependent:
          if res1.accepted: (inc stats.mate1Accepted; writers[0].write(rec1))
          else: (inc stats.mate1Rejected; writers[1].write(rec1))
          if res2.accepted: (inc stats.mate2Accepted; writers[2].write(rec2))
          else: (inc stats.mate2Rejected; writers[3].write(rec2))
    else:
      for i in 0 ..< batch.n:
        let res = batch.res[i]
        if not mappingsOut.isNil:
          mappingsOut.writeLine(mappingLine(batch.recs[i].name, res, idx, "single"))
        inc stats.inputUnits
        if res.accepted:
          inc stats.acceptedUnits
          writers[0].write(batch.recs[i])
        else:
          inc stats.rejectedUnits
          writers[1].write(batch.recs[i])

  var bufs: array[2, Batch]
  readBatch(r1, r2, paired, reads1, reads2, bufs[0])
  if threads <= 1:
    var scratch: ClassifyScratch
    while bufs[0].n > 0:
      bufs[0].res.setLen(bufs[0].n)
      for i in 0 ..< bufs[0].n:
        bufs[0].res[i] = classifySequence(idx, bufs[0].recs[i].sequence, t, scratch)
      emit(bufs[0], result)
      readBatch(r1, r2, paired, reads1, reads2, bufs[0])
    return

  # Double buffering: workers classify bufs[cur] while this thread writes the
  # previous batch (bufs[1 - cur]) and then refills that buffer.
  var workers = newSeq[Thread[ClassifyJob]](threads)
  var cur = 0
  var havePrev = false
  while bufs[cur].n > 0:
    startClassify(workers, unsafeAddr idx, t, addr bufs[cur])
    if havePrev: emit(bufs[1 - cur], result)
    readBatch(r1, r2, paired, reads1, reads2, bufs[1 - cur])
    joinThreads(workers)
    havePrev = true
    cur = 1 - cur
  if havePrev: emit(bufs[1 - cur], result)

proc runFilter*(indexPath: string, reads1: string, reads2: string, outPrefix: string,
                 t: Thresholds, pairPolicyStr: string, gzipOut: bool, mappings: bool,
                 summaryJson: string, threads = 1) =
  if reads1.len == 0:
    raise newException(PipelineError, "no input reads given: pass --reads FILE")
  if threads < 1:
    raise newException(PipelineError, "--threads must be >= 1 (got " & $threads & ")")
  if t.minCoverage < 0.0 or t.minCoverage > 1.0:
    raise newException(PipelineError, "--min-coverage must be between 0 and 1 (got " & $t.minCoverage & ")")
  let policy = parsePairPolicy(pairPolicyStr)
  let paired = reads2.len > 0
  let t0 = epochTime()
  let idx = loadIndex(indexPath)
  let loadElapsed = epochTime() - t0

  var mappingsFile: File
  var hasMappings = false
  if mappings:
    if not open(mappingsFile, outPrefix & "_mappings.tsv", fmWrite):
      raise newException(PipelineError, "could not open mappings TSV for writing")
    hasMappings = true
    mappingsFile.writeLine("read_name\tmate\tstatus\thits\tseed_count\tcovered_bases\tlongest_run\tbest_reference")

  let t1 = epochTime()
  let stats = runFilterBatches(idx, reads1, reads2, outPrefix, t, policy, gzipOut,
                               (if hasMappings: mappingsFile else: nil),
                               max(1, threads))
  let classifyElapsed = epochTime() - t1

  if hasMappings:
    mappingsFile.close()

  let unit = if paired: "pairs" else: "reads"
  stderr.writeLine &"[deribo filter] index load time: {loadElapsed:.2f}s"
  stderr.writeLine &"[deribo filter] classified {stats.inputUnits} {unit} in {classifyElapsed:.2f}s"
  if paired and policy == ppIndependent:
    stderr.writeLine &"[deribo filter] mate1: {stats.mate1Accepted} accepted / {stats.mate1Rejected} rejected"
    stderr.writeLine &"[deribo filter] mate2: {stats.mate2Accepted} accepted / {stats.mate2Rejected} rejected"
  else:
    stderr.writeLine &"[deribo filter] accepted: {stats.acceptedUnits}  rejected: {stats.rejectedUnits}"

  var summary = newJObject()
  summary["deribo_version"] = %"0.1.0"
  summary["index"] = %indexPath
  summary["index_k"] = %idx.k
  summary["index_w"] = %idx.w
  summary["reads1"] = %reads1
  summary["reads2"] = (if paired: %reads2 else: newJNull())
  summary["paired"] = %paired
  summary["pair_policy"] = (if paired: %pairPolicyStr else: newJNull())
  summary["min_hits"] = %t.minHits
  summary["min_coverage"] = %t.minCoverage
  summary["min_run"] = %t.minRun
  summary["input_units"] = %stats.inputUnits
  summary["unit"] = %unit
  summary["accepted_units"] = %stats.acceptedUnits
  summary["rejected_units"] = %stats.rejectedUnits
  summary["mate1_accepted"] = %stats.mate1Accepted
  summary["mate1_rejected"] = %stats.mate1Rejected
  summary["mate2_accepted"] = %stats.mate2Accepted
  summary["mate2_rejected"] = %stats.mate2Rejected
  summary["index_load_seconds"] = %loadElapsed
  summary["classify_seconds"] = %classifyElapsed
  summary["threads"] = %max(1, threads)
  let summaryPath = if summaryJson.len > 0: summaryJson else: outPrefix & "_summary.json"
  writeFile(summaryPath, pretty(summary))
  stderr.writeLine &"[deribo filter] wrote summary: {summaryPath}"
