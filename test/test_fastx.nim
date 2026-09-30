## Unit tests for src/deribo/fastx.nim.
## Uses small temp files rather than the sortmerna/ fixtures so these run
## standalone; test_index_classify.nim exercises the real fixtures.

import std/[unittest, os]
import ../src/deribo/fastx

proc writeTemp(name: string, content: string): string =
  result = getTempDir() / ("deribo_test_" & name)
  writeFile(result, content)

suite "FASTA parsing":
  test "multi-line sequence with no terminal newline":
    let path = writeTemp("multiline.fasta", ">r1 some comment\nACGT\nACGT\nAC")
    var r = openFastx(path)
    defer: close(r)
    check r.format == sfFasta
    var rec: FastxRecord
    check r.next(rec)
    check rec.name == "r1"
    check rec.comment == "some comment"
    check rec.sequence == "ACGTACGTAC"
    check not r.next(rec)

  test "multiple records":
    let path = writeTemp("multi.fasta", ">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n")
    var r = openFastx(path)
    defer: close(r)
    var names: seq[string]
    var rec: FastxRecord
    while r.next(rec): names.add(rec.name)
    check names == @["a", "b", "c"]

suite "FASTQ parsing":
  test "strict 4-line record":
    let path = writeTemp("ok.fastq", "@read1 desc\nACGTN\n+\nIIIII\n")
    var r = openFastx(path)
    defer: close(r)
    check r.format == sfFastq
    var rec: FastxRecord
    check r.next(rec)
    check rec.name == "read1"
    check rec.comment == "desc"
    check rec.sequence == "ACGTN"
    check rec.quality == "IIIII"

  test "quality/sequence length mismatch raises FastxError":
    let path = writeTemp("bad.fastq", "@read1\nACGTACGT\n+\nIII\n")
    var r = openFastx(path)
    defer: close(r)
    var rec: FastxRecord
    expect(FastxError):
      discard r.next(rec)

  test "truncated record raises FastxError":
    let path = writeTemp("trunc.fastq", "@read1\nACGT\n")
    var r = openFastx(path)
    defer: close(r)
    var rec: FastxRecord
    expect(FastxError):
      discard r.next(rec)

suite "format/compression detection":
  test "unrecognized first character raises FastxError":
    let path = writeTemp("bogus.txt", "not a fastx file\n")
    expect(FastxError):
      discard openFastx(path)

  test "gzip input round-trips through zcat":
    let plain = writeTemp("gz_src.fasta", ">g1\nACGTACGT\n")
    let gz = plain & ".gz"
    doAssert execShellCmd("gzip -kf " & plain) == 0
    var r = openFastx(gz)
    defer: close(r)
    var rec: FastxRecord
    check r.next(rec)
    check rec.name == "g1"
    check rec.sequence == "ACGTACGT"

suite "writer":
  test "writes FASTA and gzip-compresses on close":
    let outPath = getTempDir() / "deribo_test_writer_out.fasta"
    var w = openFastxWriter(outPath, sfFasta, gzip = true)
    w.write(FastxRecord(name: "x", comment: "", sequence: "ACGT"))
    close(w)
    check fileExists(outPath & ".gz")
    check not fileExists(outPath)
