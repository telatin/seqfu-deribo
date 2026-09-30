## Per-read classification against a DeriboIndex.
##
## This plays the role of REIMPLEMENTATION_GUIDE.md's seed-matching + LIS
## + Smith-Waterman pipeline, but reduced to its v0.1 MVP form: find the
## canonical seeds (see kmer.nim) the read shares with the index, and accept
## if they clear three thresholds. There is deliberately no alignment step
## yet (no Myers, no Parasail) -- see DERIBO_PLAN.md's milestone roadmap for
## where that lands. Governing constraint from NIM_SORTMERNA_PLAN.md: this
## is a coarser classifier than the planned one, so its thresholds should be
## treated as tunables, not calibrated defaults.
##
## Why not seed count alone: overlapping seeds from one short exact match
## inflate the count (an 18 bp match gives two 17-mers), and at small k a
## few scattered chance hits per read are expected against a large database.
## Covered bases (union of hit k-mer spans) and the longest contiguous
## covered run measure the actual amount of shared sequence instead.

import std/algorithm
import ./index
import ./kmer

const
  # Calibrated on testset/ (SRR30521510, SRR30521528) against SortMeRNA and
  # Kraken 2 calls with a k=21 index: coverage 0.25 matches SortMeRNA's total
  # rRNA calls to within ~0.6% while cutting extra calls ~16-20x vs the old
  # k=17/hits>=2 rule. minRun is off by default (any hit already implies a
  # run of k bp); it is kept as a tunable for short-k indexes.
  DefaultMinHits* = 2
  DefaultMinCoverage* = 0.25
  DefaultMinRun* = 0

type
  Thresholds* = object
    minHits*: int          ## minimum number of read seeds found in the index
    minCoverage*: float    ## minimum fraction of read bases covered by hit seeds
    minRun*: int           ## minimum longest contiguous covered span, in bp

  ClassifyResult* = object
    hits*: int         ## number of read seeds found in the index
    seedCount*: int    ## total number of seeds extracted from the read
    coveredBases*: int ## read bases covered by at least one hit seed
    longestRun*: int   ## longest contiguous stretch of covered bases
    bestRef*: int      ## index into idx.refs of the most-hit reference, or -1
    accepted*: bool

  ClassifyScratch* = object
    ## Per-thread buffers reused across reads, so classifying allocates
    ## nothing once they have grown to the longest read.
    seeds: seq[Seed]
    slots: seq[int]
    refIds: seq[uint32]

proc prefetch(p: pointer) {.importc: "__builtin_prefetch", nodecl.}

func thresholds*(minHits = DefaultMinHits, minCoverage = DefaultMinCoverage,
                 minRun = DefaultMinRun): Thresholds =
  Thresholds(minHits: minHits, minCoverage: minCoverage, minRun: minRun)

proc bestRefOf(refIds: var seq[uint32]): int =
  ## Most frequent reference among the hits; ties go to the lowest refId.
  if refIds.len == 0: return -1
  refIds.sort()
  var bestCount = 0
  var i = 0
  while i < refIds.len:
    var j = i + 1
    while j < refIds.len and refIds[j] == refIds[i]: inc j
    if j - i > bestCount:
      bestCount = j - i
      result = refIds[i].int
    i = j

proc classifySequence*(idx: DeriboIndex, sequence: string, t: Thresholds,
                       scratch: var ClassifyScratch): ClassifyResult =
  extractSeedsPos(sequence, idx.k, idx.w, scratch.seeds)
  let n = scratch.seeds.len
  # Hash every seed and prefetch its slot before probing any of them: the
  # table is far larger than cache, so the lookups are DRAM-latency bound and
  # issuing them together lets the misses overlap.
  scratch.slots.setLen(n)
  for i in 0 ..< n:
    let slot = idx.table.slotOf(scratch.seeds[i].val)
    scratch.slots[i] = slot
    prefetch(idx.table.keyAddr(slot))
  scratch.refIds.setLen(0)
  var hits = 0
  # Seeds arrive in increasing start order, so the union of hit spans can be
  # swept in one pass: [runStart, runEnd) is the current merged interval.
  var covered = 0
  var longest = 0
  var runStart = -1
  var runEnd = -1
  for i in 0 ..< n:
    let seed = scratch.seeds[i]
    let slot = idx.table.findFrom(seed.val, scratch.slots[i])
    if slot >= 0:
      inc hits
      scratch.refIds.add(idx.table.entryAt(slot).refId)
      let e = seed.pos + idx.k
      if seed.pos > runEnd:
        runStart = seed.pos
        covered += idx.k
      else:
        covered += max(0, e - runEnd)
      runEnd = max(runEnd, e)
      longest = max(longest, runEnd - runStart)
  result.hits = hits
  result.seedCount = n
  result.coveredBases = covered
  result.longestRun = longest
  let coverage = if sequence.len > 0: covered / sequence.len else: 0.0
  result.accepted = hits >= t.minHits and coverage >= t.minCoverage and
                    longest >= t.minRun
  result.bestRef = bestRefOf(scratch.refIds)

proc classifySequence*(idx: DeriboIndex, sequence: string, t: Thresholds): ClassifyResult =
  var scratch: ClassifyScratch
  classifySequence(idx, sequence, t, scratch)

proc classifySequence*(idx: DeriboIndex, sequence: string, minHits: int): ClassifyResult =
  ## Seed-count-only rule (the original v0.1 behaviour).
  classifySequence(idx, sequence, thresholds(minHits = minHits))
