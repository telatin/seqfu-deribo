## Canonical k-mer and windowed-minimizer extraction.
##
## This is deribo's v0.1 "seed" layer (see DERIBO_PLAN.md's seeds/ module
## group): a strand-invariant canonical k-mer stream, optionally thinned to
## one minimizer per window `w` (deacon-style) to keep large reference
## directories tractable. There is no burst trie, MPHF, or Levenshtein
## automaton here -- see REIMPLEMENTATION_GUIDE.md Component 3/4 for that
## design if/when deribo grows an exact SortMeRNA-compatible index.

import std/deques

const InvalidBase*: uint8 = 0xFF

func baseCode*(c: char): uint8 {.inline.} =
  ## 2-bit code for A/C/G/T (U treated as T). Anything else is invalid
  ## (ambiguity codes, gaps, whitespace) and returns InvalidBase.
  case c
  of 'A', 'a': 0'u8
  of 'C', 'c': 1'u8
  of 'G', 'g': 2'u8
  of 'T', 't', 'U', 'u': 3'u8
  else: InvalidBase

type
  Seed* = tuple[pos: int, val: uint64]
    ## A seed value plus the 0-based start of its k-mer in the sequence.

iterator canonicalKmersPos*(sequence: string, k: int): Seed =
  ## Yields the canonical (min of forward / reverse-complement) 2-bit-packed
  ## k-mer for every valid window of `sequence`, with its start offset, in
  ## O(len(sequence)).
  ## Windows containing an ambiguous/invalid base are skipped entirely; the
  ## rolling state resets so no k-mer spans across the gap.
  assert k >= 1 and k <= 32, "k must be in 1..32 to fit a uint64"
  let mask = if k == 32: uint64.high else: (1'u64 shl (2 * k)) - 1
  var fwd: uint64 = 0
  var rc: uint64 = 0
  var validRun = 0
  let shiftTop = 2 * (k - 1)
  for i in 0 ..< sequence.len:
    let code = baseCode(sequence[i])
    if code == InvalidBase:
      validRun = 0
      fwd = 0
      rc = 0
      continue
    fwd = ((fwd shl 2) or code.uint64) and mask
    let comp = 3'u8 - code
    rc = (rc shr 2) or (comp.uint64 shl shiftTop)
    inc validRun
    if validRun >= k:
      yield (pos: i - k + 1, val: (if fwd < rc: fwd else: rc))

iterator canonicalKmers*(sequence: string, k: int): uint64 =
  ## Value-only view of `canonicalKmersPos`.
  for s in canonicalKmersPos(sequence, k):
    yield s.val

proc extractSeedsPos*(sequence: string, k: int, w: int, result: var seq[Seed]) =
  ## Fills `result` (cleared first, capacity reused) with the seed stream
  ## used for both indexing and querying, in increasing start order:
  ## - w <= 1: every canonical k-mer (maximum sensitivity).
  ## - w >  1: one canonical k-mer per window of `w` consecutive k-mer
  ##   positions (the minimum value in the window), deduplicated when the
  ##   window minimum repeats -- a standard minimizer sketch, computed with
  ##   a monotonic deque in O(n).
  result.setLen(0)
  if w <= 1:
    for s in canonicalKmersPos(sequence, k):
      result.add(s)
    return
  var dq: Deque[tuple[ord: int, seed: Seed]] = initDeque[tuple[ord: int, seed: Seed]]()
  var ord = 0
  var hasLastEmitted = false
  var lastEmittedVal = 0'u64
  for s in canonicalKmersPos(sequence, k):
    while dq.len > 0 and dq.peekLast().seed.val >= s.val:
      discard dq.popLast()
    dq.addLast((ord: ord, seed: s))
    while dq.peekFirst().ord <= ord - w:
      discard dq.popFirst()
    if ord >= w - 1:
      let m = dq.peekFirst().seed
      # Dedup by *value*, not position: ties are resolved in favor of the
      # rightmost window, so a run of identical minimizer values (common in
      # homopolymers and tandem repeats) would otherwise re-emit on every
      # single-base shift even though nothing new was found.
      if not hasLastEmitted or m.val != lastEmittedVal:
        result.add(m)
        lastEmittedVal = m.val
        hasLastEmitted = true
    inc ord

proc extractSeedsPos*(sequence: string, k: int, w: int): seq[Seed] =
  extractSeedsPos(sequence, k, w, result)

proc extractSeeds*(sequence: string, k: int, w: int): seq[uint64] =
  ## Value-only view of `extractSeedsPos` (what the index stores).
  for s in extractSeedsPos(sequence, k, w):
    result.add(s.val)
