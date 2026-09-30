## Minimal FASTA/FASTQ reader and writer.
##
## Scope (v0.1 MVP -- see ../../DERIBO_PLAN.md for the target architecture):
## - Auto-detects FASTA vs FASTQ by peeking the first record marker ('>' or '@').
## - Auto-detects gzip input by file extension and decompresses via `pigz -dc`
##   (or `zcat` if pigz is not on PATH); gzip output is streamed through
##   `pigz`/`gzip` as it is written. This avoids a hard C/zlib binding
##   dependency and keeps (de)compression off the classifying threads.
## - Handles multi-line FASTA sequences and strict 4-line FASTQ records.
## - Preserves the original sequence/quality text verbatim; only used for
##   classification internally, never rewritten in place.

import std/[streams, strutils, osproc, os]

type
  SeqFormat* = enum
    sfFasta, sfFastq

  FastxRecord* = object
    name*: string
    comment*: string
    sequence*: string
    quality*: string ## empty for FASTA

  FastxReader* = object
    strm: Stream
    childProc: Process        ## nil unless the input is gzip-compressed
    format*: SeqFormat
    path*: string
    pendingHeader: string  ## FASTA only: header line already consumed
    hasPending: bool
    atEof: bool
    plusLine: string       ## FASTQ scratch for the '+' line

  FastxError* = object of CatchableError

proc isGz(path: string): bool =
  path.toLowerAscii.endsWith(".gz")

proc gzipTool(): string =
  ## "pigz" when available (multi-threaded, faster to decompress too), else "gzip".
  if findExe("pigz").len > 0: "pigz" else: "gzip"

proc openFastx*(path: string): FastxReader =
  ## Opens `path` for reading, auto-detecting format and gzip compression.
  if not fileExists(path):
    raise newException(FastxError, "input file not found: " & path)
  result.path = path
  if isGz(path):
    result.childProc = startProcess(gzipTool(), args = ["-dc", path],
                                 options = {poUsePath, poStdErrToStdOut})
    result.strm = result.childProc.outputStream
  else:
    let fs = newFileStream(path, fmRead)
    if fs.isNil:
      raise newException(FastxError, "could not open input file: " & path)
    result.strm = fs

  # Peek the first non-blank line to decide FASTA vs FASTQ.
  var line: string
  while true:
    if not result.strm.readLine(line):
      raise newException(FastxError, "input file is empty or unreadable: " & path)
    if line.len > 0: break
  if line[0] == '>':
    result.format = sfFasta
    result.pendingHeader = line
    result.hasPending = true
  elif line[0] == '@':
    result.format = sfFastq
    result.pendingHeader = line ## reused as "first header already read"
    result.hasPending = true
  else:
    raise newException(FastxError,
      "unrecognized sequence format in " & path & " (expected '>' or '@', got '" &
      (if line.len > 0: $line[0] else: "<empty>") & "')")

proc close*(r: var FastxReader) =
  if not r.childProc.isNil:
    close(r.childProc)
  elif not r.strm.isNil:
    r.strm.close()

proc parseHeader(line: string, name, comment: var string) =
  ## Splits ">name comment" / "@name comment" at the first space or tab.
  var cut = 1
  while cut < line.len and line[cut] notin {' ', '\t'}: inc cut
  name.setLen(cut - 1)
  if cut > 1: copyMem(addr name[0], unsafeAddr line[1], cut - 1)
  if cut >= line.len:
    comment.setLen(0)
  else:
    comment = line[cut + 1 .. ^1].strip()

proc stripCr(s: string): string =
  if s.len > 0 and s[^1] == '\r': s[0 ..< ^1] else: s

proc stripCrInPlace(s: var string) {.inline.} =
  if s.len > 0 and s[^1] == '\r': s.setLen(s.len - 1)

proc nextFasta(r: var FastxReader, rec: var FastxRecord): bool =
  if not r.hasPending:
    return false
  parseHeader(r.pendingHeader, rec.name, rec.comment)
  rec.sequence.setLen(0)
  rec.quality.setLen(0)
  r.hasPending = false
  var line: string
  while r.strm.readLine(line):
    if line.len > 0 and line[0] == '>':
      r.pendingHeader = line
      r.hasPending = true
      return true
    rec.sequence.add(stripCr(line))
  return true

proc nextFastq(r: var FastxReader, rec: var FastxRecord): bool =
  if r.hasPending:
    r.hasPending = false
  else:
    if r.atEof: return false
    if not r.strm.readLine(r.pendingHeader):
      r.atEof = true
      return false
    if r.pendingHeader.len == 0:
      r.atEof = true
      return false
  parseHeader(r.pendingHeader, rec.name, rec.comment)
  if not r.strm.readLine(rec.sequence):
    raise newException(FastxError, "truncated FASTQ record (missing sequence line) for read '" & rec.name & "' in " & r.path)
  if not r.strm.readLine(r.plusLine):
    raise newException(FastxError, "truncated FASTQ record (missing '+' line) for read '" & rec.name & "' in " & r.path)
  if r.plusLine.len == 0 or r.plusLine[0] != '+':
    raise newException(FastxError, "malformed FASTQ record (expected '+' separator) for read '" & rec.name & "' in " & r.path)
  if not r.strm.readLine(rec.quality):
    raise newException(FastxError, "truncated FASTQ record (missing quality line) for read '" & rec.name & "' in " & r.path)
  stripCrInPlace(rec.sequence)
  stripCrInPlace(rec.quality)
  if rec.quality.len != rec.sequence.len:
    raise newException(FastxError, "quality length (" & $rec.quality.len &
      ") does not match sequence length (" & $rec.sequence.len & ") for read '" &
      rec.name & "' in " & r.path)
  return true

proc next*(r: var FastxReader, rec: var FastxRecord): bool =
  ## Reads the next record into `rec`. Returns false at end of input.
  case r.format
  of sfFasta: nextFasta(r, rec)
  of sfFastq: nextFastq(r, rec)

# --- Writer -----------------------------------------------------------------

type
  FastxWriter* = object
    f: File
    format*: SeqFormat
    path*: string
    gzip: bool

proc popen(cmd, mode: cstring): File {.importc, header: "<stdio.h>".}
proc pclose(f: File): cint {.importc, header: "<stdio.h>".}

proc openFastxWriter*(path: string, format: SeqFormat, gzip: bool,
                      compressThreads = 1): FastxWriter =
  ## With `gzip`, output goes to `path & ".gz"` through a `pigz`/`gzip`
  ## pipe, compressed concurrently with classification.
  result.path = path
  result.format = format
  result.gzip = gzip
  if gzip:
    let tool = gzipTool()
    let threadArg = if tool == "pigz": " -p " & $max(1, compressThreads) else: ""
    let cmd = tool & threadArg & " -c > " & quoteShell(path & ".gz")
    result.f = popen(cmd.cstring, "w")
    if result.f.isNil:
      raise newException(FastxError, "could not start compressor for output file: " & path & ".gz")
  elif not open(result.f, path, fmWrite):
    raise newException(FastxError, "could not open output file for writing: " & path)

proc writeHeader(f: File, marker: char, rec: FastxRecord) {.inline.} =
  f.write(marker)
  f.write(rec.name)
  if rec.comment.len > 0:
    f.write(' ')
    f.write(rec.comment)
  f.write('\n')

proc write*(w: var FastxWriter, rec: FastxRecord) =
  case w.format
  of sfFasta:
    w.f.writeHeader('>', rec)
    w.f.write(rec.sequence)
    w.f.write('\n')
  of sfFastq:
    w.f.writeHeader('@', rec)
    w.f.write(rec.sequence)
    w.f.write("\n+\n")
    w.f.write(rec.quality)
    w.f.write('\n')

proc close*(w: var FastxWriter) =
  if w.gzip:
    let rc = pclose(w.f)
    if rc != 0:
      raise newException(FastxError, "gzip compression failed for output file: " & w.path & ".gz")
  else:
    w.f.close()

proc outputPath*(w: FastxWriter): string =
  if w.gzip: w.path & ".gz" else: w.path

proc defaultExt*(format: SeqFormat): string =
  case format
  of sfFasta: ".fasta"
  of sfFastq: ".fastq"
