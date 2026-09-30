## Unit tests for src/deribo/kmer.nim: canonical k-mer extraction and the
## windowed minimizer thinning used by `deribo index -w`.

import std/[unittest, sets, sequtils, tables, strutils]
import ../src/deribo/kmer

suite "baseCode":
  test "encodes ACGT and U, case-insensitively":
    check baseCode('A') == 0'u8
    check baseCode('a') == 0'u8
    check baseCode('C') == 1'u8
    check baseCode('G') == 2'u8
    check baseCode('T') == 3'u8
    check baseCode('U') == 3'u8
    check baseCode('u') == 3'u8

  test "flags ambiguity codes and gaps as invalid":
    check baseCode('N') == InvalidBase
    check baseCode('-') == InvalidBase
    check baseCode('R') == InvalidBase

suite "canonicalKmers":
  test "yields len-k+1 kmers for a clean sequence":
    let s = "ACGTACGTACGT" # 12 bases
    let ks = toSeq(canonicalKmers(s, 4))
    check ks.len == 9

  test "is strand-invariant: seq and its reverse complement share the same canonical multiset":
    proc revcomp(s: string): string =
      const comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.toTable
      result = newString(s.len)
      for i, c in s:
        result[s.len - 1 - i] = comp[c]
    let s = "GATTACAGATTACAGGGCATTAG"
    let fwd = toSeq(canonicalKmers(s, 9)).toHashSet
    let rev = toSeq(canonicalKmers(revcomp(s), 9)).toHashSet
    check fwd == rev

  test "resets across ambiguous bases (no kmer spans an N)":
    let s = "ACGTACGTAC" & "N" & "ACGTACGTAC" # two 10-base runs, k=10
    let ks = toSeq(canonicalKmers(s, 10))
    check ks.len == 2 # exactly one full-length window on each side of the N

  test "handles k=32 (uint64 mask edge case) without crashing":
    let s = "ACGT".repeat(10) # 40 bases
    let ks = toSeq(canonicalKmers(s, 32))
    check ks.len == 9

suite "extractSeeds":
  test "w<=1 returns exactly the canonical kmer stream":
    let s = "ACGTACGTACGTACGTACGT"
    check extractSeeds(s, 6, 1) == toSeq(canonicalKmers(s, 6))
    check extractSeeds(s, 6, 0) == toSeq(canonicalKmers(s, 6))

  test "w>1 thins the stream and every emitted seed is a real canonical kmer":
    let s = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"
    let allValues = toSeq(canonicalKmers(s, 9))
    let full = allValues.toHashSet
    let thinned = extractSeeds(s, 9, 5)
    check thinned.len <= allValues.len # thinning never produces more than one per window
    check thinned.len > 0
    for v in thinned:
      check v in full

  test "w>1 never emits the same minimizer for two consecutive windows":
    let s = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" # homopolymer -> constant kmer
    let thinned = extractSeeds(s, 5, 4)
    # a homopolymer has exactly one distinct canonical kmer value repeated;
    # dedup-on-repeat means it must collapse to a single emitted seed.
    check thinned.len == 1

suite "extractSeedsPos":
  test "positions are k-mer start offsets and skip across ambiguous bases":
    let s = "ACGTAC" & "N" & "GGTACC" # two 6-base runs, k=5
    let seeds = extractSeedsPos(s, 5, 1)
    check seeds.len == 4
    check seeds.mapIt(it.pos) == @[0, 1, 7, 8]
    check seeds.mapIt(it.val) == toSeq(canonicalKmers(s, 5))

  test "minimizer positions are non-decreasing (classifier relies on this)":
    let s = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"
    let seeds = extractSeedsPos(s, 9, 5)
    for i in 1 ..< seeds.len:
      check seeds[i].pos >= seeds[i - 1].pos
    check seeds.mapIt(it.val) == extractSeeds(s, 9, 5)
