## Reference index: build, save, load, and directory-as-database resolution.
##
## v0.1 MVP index (see DERIBO_PLAN.md §"Dependency and tooling decisions" --
## this deliberately skips the mmap'd section-checksummed format and the
## burst trie / MPHF from REIMPLEMENTATION_GUIDE.md Component 3; it is a
## single in-memory open-addressing hash table whose slot arrays are written
## to disk verbatim, so loading is two bulk reads rather than a rehash). The one hard requirement carried over from that plan: the
## reference *database* argument may be a directory of FASTA files, exactly
## like sortmerna-4.3.4/ or sortmerna/test_data/refs/.

import std/[os, streams, strutils, algorithm]
import ./fastx
import ./kmer

const
  IndexMagic = "DRB1"
  IndexVersion = 2'u32
    ## v1: (key, refId, count) records, rebuilt into a table on load.
    ## v2: raw SeedTable slot arrays. v1 files are still readable.
  EmptyKey = uint64.high
    ## Never a canonical k-mer: for k = 32 the all-T forward k-mer has the
    ## all-A (0) reverse complement, so min(fwd, rc) can't be all ones.
  MaxLoad = 0.7
  fastaExts = ["fasta", "fa", "fna", "ffn", "frn", "fsa"]

type
  RefInfo* = object
    name*: string
    length*: int

  SeedEntry* = object
    refId*: uint32
    count*: uint32

  SeedTable* = object
    ## Flat linear-probing hash table, uint64 seed -> SeedEntry. Keys and
    ## values live in separate arrays so a probe sequence (the common miss
    ## case when filtering) usually stays within one cache line of keys.
    keys: seq[uint64]
    vals: seq[SeedEntry]
    mask: int
    count: int

  DeriboIndex* = object
    k*: int
    w*: int
    refs*: seq[RefInfo]
    table*: SeedTable

  DeriboIndexError* = object of CatchableError

# --- SeedTable ------------------------------------------------------------------

func slotOf*(t: SeedTable, key: uint64): int {.inline.} =
  ## Home slot of `key` (Fibonacci hashing: the multiply mixes all key bits
  ## into the high bits, which are the ones kept).
  int((key * 0x9E3779B97F4A7C15'u64) shr 32) and t.mask

proc initSeedTable*(expected = 0): SeedTable =
  var cap = 1024
  while expected.float > cap.float * MaxLoad: cap = cap * 2
  result.keys = newSeq[uint64](cap)
  for k in result.keys.mitems: k = EmptyKey
  result.vals = newSeq[SeedEntry](cap)
  result.mask = cap - 1

func len*(t: SeedTable): int {.inline.} = t.count

func capacity*(t: SeedTable): int {.inline.} = t.keys.len

func findFrom*(t: SeedTable, key: uint64, slot: int): int {.inline.} =
  ## Slot holding `key`, probing from `slot` (normally `slotOf(key)`), or -1.
  var i = slot
  while true:
    let k = t.keys[i]
    if k == key: return i
    if k == EmptyKey: return -1
    i = (i + 1) and t.mask

func find*(t: SeedTable, key: uint64): int {.inline.} =
  t.findFrom(key, t.slotOf(key))

func contains*(t: SeedTable, key: uint64): bool {.inline.} = t.find(key) >= 0

func entryAt*(t: SeedTable, slot: int): SeedEntry {.inline.} = t.vals[slot]

func keyAddr*(t: SeedTable, slot: int): pointer {.inline.} =
  unsafeAddr t.keys[slot]

proc `[]`*(t: SeedTable, key: uint64): SeedEntry =
  let i = t.find(key)
  if i < 0: raise newException(KeyError, "seed not in table: " & $key)
  t.vals[i]

proc grow(t: var SeedTable)

proc slotFor(t: var SeedTable, key: uint64, fresh: var bool): int =
  ## Slot for `key`, inserting it (zeroed value) if absent.
  if (t.count + 1).float > t.keys.len.float * MaxLoad: t.grow()
  var i = t.slotOf(key)
  while true:
    let k = t.keys[i]
    if k == key:
      fresh = false
      return i
    if k == EmptyKey:
      t.keys[i] = key
      inc t.count
      fresh = true
      return i
    i = (i + 1) and t.mask

proc grow(t: var SeedTable) =
  var bigger = initSeedTable(t.keys.len * 2)
  for i in 0 ..< t.keys.len:
    if t.keys[i] != EmptyKey:
      var fresh: bool
      bigger.vals[bigger.slotFor(t.keys[i], fresh)] = t.vals[i]
  t = move bigger

proc `[]=`*(t: var SeedTable, key: uint64, val: SeedEntry) =
  var fresh: bool
  t.vals[t.slotFor(key, fresh)] = val

iterator pairs*(t: SeedTable): (uint64, SeedEntry) =
  for i in 0 ..< t.keys.len:
    if t.keys[i] != EmptyKey: yield (t.keys[i], t.vals[i])

# --- Reference files --------------------------------------------------------------

proc isFastaPath(path: string): bool =
  var name = path.extractFilename.toLowerAscii
  if name.endsWith(".gz"):
    name = name[0 ..< ^3]
  let dot = name.rfind('.')
  if dot < 0: return false
  let ext = name[dot + 1 .. ^1]
  ext in fastaExts

proc collectFastaFiles*(inputs: seq[string]): seq[string] =
  ## Expands directories into their sorted (deterministic) FASTA file
  ## listing; leaves explicit file paths untouched. This is what lets
  ## `deribo index --db sortmerna-4.3.4/` and
  ## `deribo index --db path/to/one.fasta` both work.
  result = @[]
  for input in inputs:
    if dirExists(input):
      var found: seq[string] = @[]
      for kind, path in walkDir(input):
        if kind == pcFile and isFastaPath(path):
          found.add(path)
      if found.len == 0:
        raise newException(DeriboIndexError,
          "no FASTA files (." & fastaExts.join(", .") & ", optionally .gz) found in directory: " & input)
      found.sort()
      result.add(found)
    elif fileExists(input):
      result.add(input)
    else:
      raise newException(DeriboIndexError, "reference path not found: " & input)

proc buildIndex*(paths: seq[string], k: int, w: int): DeriboIndex =
  result.k = k
  result.w = w
  result.refs = @[]
  result.table = initSeedTable()
  for path in paths:
    var reader = openFastx(path)
    if reader.format != sfFasta:
      close(reader)
      raise newException(DeriboIndexError, "reference file is not FASTA: " & path)
    var rec: FastxRecord
    while reader.next(rec):
      let refId = result.refs.len.uint32
      result.refs.add(RefInfo(name: rec.name, length: rec.sequence.len))
      for seed in extractSeeds(rec.sequence, k, w):
        var fresh: bool
        let slot = result.table.slotFor(seed, fresh)
        if fresh: result.table.vals[slot].refId = refId
        inc result.table.vals[slot].count
    close(reader)
  if result.refs.len == 0:
    raise newException(DeriboIndexError, "no reference sequences found in: " & paths.join(", "))

proc saveIndex*(idx: DeriboIndex, path: string) =
  var s = newFileStream(path, fmWrite)
  if s.isNil:
    raise newException(DeriboIndexError, "could not open index file for writing: " & path)
  defer: s.close()
  s.write(IndexMagic)
  s.write(IndexVersion)
  s.write(idx.k.uint32)
  s.write(idx.w.uint32)
  s.write(idx.refs.len.uint32)
  s.write(idx.table.len.uint64)
  for r in idx.refs:
    s.write(r.name.len.uint32)
    if r.name.len > 0:
      s.writeData(unsafeAddr r.name[0], r.name.len)
    s.write(r.length.uint32)
  s.write(idx.table.keys.len.uint64)
  s.writeData(unsafeAddr idx.table.keys[0], idx.table.keys.len * sizeof(uint64))
  s.writeData(unsafeAddr idx.table.vals[0], idx.table.vals.len * sizeof(SeedEntry))

proc readU32(s: Stream): uint32 =
  s.read(result)

proc readU64(s: Stream): uint64 =
  s.read(result)

proc loadIndex*(path: string): DeriboIndex =
  if not fileExists(path):
    raise newException(DeriboIndexError, "index file not found: " & path)
  var s = newFileStream(path, fmRead)
  if s.isNil:
    raise newException(DeriboIndexError, "could not open index file: " & path)
  defer: s.close()
  let magic = s.readStr(4)
  if magic != IndexMagic:
    raise newException(DeriboIndexError, "not a deribo index file (bad magic): " & path)
  let version = readU32(s)
  if version notin [1'u32, IndexVersion]:
    raise newException(DeriboIndexError, "unsupported deribo index version " & $version & " in " & path)
  result.k = readU32(s).int
  result.w = readU32(s).int
  let refCount = readU32(s).int
  let kmerCount = readU64(s).int
  result.refs = newSeq[RefInfo](refCount)
  for i in 0 ..< refCount:
    let nameLen = readU32(s).int
    var name = newString(nameLen)
    if nameLen > 0:
      discard s.readData(addr name[0], nameLen)
    let length = readU32(s).int
    result.refs[i] = RefInfo(name: name, length: length)
  if version == 1:
    result.table = initSeedTable(kmerCount)
    for i in 0 ..< kmerCount:
      let key = readU64(s)
      let refId = readU32(s)
      let count = readU32(s)
      result.table[key] = SeedEntry(refId: refId, count: count)
  else:
    let cap = readU64(s).int
    if cap < 1 or (cap and (cap - 1)) != 0 or kmerCount > cap:
      raise newException(DeriboIndexError, "corrupt deribo index (bad table size): " & path)
    result.table.keys = newSeqUninit[uint64](cap)
    result.table.vals = newSeqUninit[SeedEntry](cap)
    let keyBytes = cap * sizeof(uint64)
    let valBytes = cap * sizeof(SeedEntry)
    if s.readData(addr result.table.keys[0], keyBytes) != keyBytes or
       s.readData(addr result.table.vals[0], valBytes) != valBytes:
      raise newException(DeriboIndexError, "truncated deribo index: " & path)
    result.table.mask = cap - 1
    result.table.count = kmerCount

proc totalBases*(idx: DeriboIndex): int =
  for r in idx.refs: result += r.length
