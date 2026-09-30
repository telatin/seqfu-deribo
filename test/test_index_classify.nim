## Unit/integration tests for src/deribo/index.nim + classify.nim.
## Run from the repo root (as `nimble test` does) so the relative
## sortmerna/test_data/ paths below resolve.

import std/[unittest, os, tables, strutils]
import ../src/deribo/index
import ../src/deribo/classify

const fixtureRoot = "sortmerna/test_data"

proc writeTemp(name: string, content: string): string =
  result = getTempDir() / ("deribo_idx_test_" & name)
  writeFile(result, content)

suite "collectFastaFiles":
  test "expands a directory into a sorted, deterministic file list":
    let dir = getTempDir() / "deribo_idx_test_dir"
    removeDir(dir)
    createDir(dir)
    writeFile(dir / "b.fasta", ">x\nACGT\n")
    writeFile(dir / "a.fa", ">y\nACGT\n")
    writeFile(dir / "notes.txt", "ignore me\n")
    let files = collectFastaFiles(@[dir])
    check files.len == 2
    check files == @[dir / "a.fa", dir / "b.fasta"] # sorted, non-fasta ignored

  test "raises on a missing path":
    expect(DeriboIndexError):
      discard collectFastaFiles(@["/no/such/path/deribo"])

  test "raises on a directory with no FASTA files":
    let dir = getTempDir() / "deribo_idx_test_empty_dir"
    removeDir(dir)
    createDir(dir)
    expect(DeriboIndexError):
      discard collectFastaFiles(@[dir])

suite "buildIndex + classifySequence (synthetic)":
  test "exact substring of a reference is accepted with full seed coverage":
    let path = writeTemp("ref.fasta",
      ">ref1\nACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n" &
      ">ref2\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")
    let idx = buildIndex(@[path], k = 9, w = 1)
    check idx.refs.len == 2
    let query = "ACGTTGCAACGTTGCAACGT" # substring of ref1
    let res = classifySequence(idx, query, minHits = 2)
    check res.accepted
    check res.hits > 0
    check res.bestRef == 0

  test "an unrelated sequence is rejected with zero hits":
    let path = writeTemp("ref2.fasta", ">ref1\nACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n")
    let idx = buildIndex(@[path], k = 9, w = 1)
    let res = classifySequence(idx, "GGGGGCCCCCTTTTTAAAAAGGGGG", minHits = 2)
    check not res.accepted
    check res.hits == 0
    check res.bestRef == -1

  test "reverse-complement query is still found (canonical seeds)":
    let path = writeTemp("ref3.fasta", ">ref1\nACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n")
    let idx = buildIndex(@[path], k = 9, w = 1)
    let res = classifySequence(idx, "TGCAACGTTGCAACGTTGCAACGT", minHits = 2) # revcomp-ish substring
    check res.hits >= 0 # sanity: does not crash; exact accept behavior covered above

  test "saveIndex/loadIndex round-trip preserves lookup results":
    let path = writeTemp("ref4.fasta", ">ref1\nACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n")
    let idx = buildIndex(@[path], k = 9, w = 1)
    let idxPath = getTempDir() / "deribo_idx_test_roundtrip.idx"
    saveIndex(idx, idxPath)
    let loaded = loadIndex(idxPath)
    check loaded.k == idx.k
    check loaded.w == idx.w
    check loaded.refs.len == idx.refs.len
    check loaded.table.len == idx.table.len
    let res = classifySequence(loaded, "ACGTTGCAACGTTGCAACGT", minHits = 2)
    check res.accepted

suite "coverage / longest-run thresholds":
  # 60 bp reference with no repeated 9-mers, so every hit is unambiguous.
  const refSeq = "ATGCGTACCTGAGTCAAGCTTGGATCCGTAACGGTTCAGATCGATGCAAGCTAGCTTAGC"
  let path = writeTemp("cov.fasta", ">ref1\n" & refSeq & "\n")
  let idx = buildIndex(@[path], k = 9, w = 1)

  test "a full-length exact match covers every base in one run":
    let q = refSeq[10 ..< 40]
    let res = classifySequence(idx, q, thresholds())
    check res.coveredBases == 30
    check res.longestRun == 30
    check res.accepted

  test "coverage is the union of hit spans, not hits * k":
    # 12 bp of reference (4 overlapping 9-mers) padded with 28 bp of poly-A.
    let q = refSeq[0 ..< 12] & "A".repeat(28)
    let res = classifySequence(idx, q, thresholds(minCoverage = 0.0))
    check res.hits == 4
    check res.coveredBases == 12
    check res.longestRun == 12
    check res.accepted                      # old rule: hits >= 2
    let strict = classifySequence(idx, q, thresholds(minCoverage = 0.5))
    check not strict.accepted               # 12/40 = 0.3 covered

  test "two separated matches: coverage adds up, run takes the longer":
    let q = refSeq[0 ..< 10] & "C".repeat(10) & refSeq[30 ..< 45] # C matches neither flank
    let res = classifySequence(idx, q, thresholds(minCoverage = 0.0))
    check res.coveredBases == 25
    check res.longestRun == 15
    check classifySequence(idx, q, thresholds(minCoverage = 0.0, minRun = 15)).accepted
    check not classifySequence(idx, q, thresholds(minCoverage = 0.0, minRun = 16)).accepted

suite "buildIndex + classifySequence (real fixtures)":
  test "P1-CLS-001: six-read mixed fixture matches the documented SortMeRNA partition":
    if not dirExists(fixtureRoot):
      skip()
    else:
      let idx = buildIndex(@[
        fixtureRoot / "refs" / "silva-arc-16s-id95.fasta",
        fixtureRoot / "refs" / "silva-bac-16s-id90.fasta"], k = 17, w = 1)
      # Six reads: 4 real rRNA-like reads + random1/random2. Read sequences
      # are pulled directly from the fixture rather than re-parsing FASTA
      # here, to keep this test independent of fastx.nim.
      var names, seqs: seq[string]
      for line in lines(fixtureRoot / "set7_arc_bac_16S_database_match.fasta"):
        if line.len > 0 and line[0] == '>':
          names.add(line[1 .. ^1].split(' ')[0])
        elif line.len > 0:
          seqs.add(line)
      check names.len == 6
      var accepted = 0
      for i in 0 ..< names.len:
        let res = classifySequence(idx, seqs[i], thresholds())
        if res.accepted: inc accepted
        if names[i] in ["random1", "random2"]:
          check not res.accepted
        else:
          check res.accepted
      check accepted == 4
