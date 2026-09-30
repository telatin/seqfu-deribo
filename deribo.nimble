# Package

version       = "0.1.0"
author        = "QIB"
description   = "Fast k-mer/minimizer based rRNA read filter (v0.1 MVP toward the deribo design in DERIBO_PLAN.md)"
license       = "MIT"
srcDir        = "src"
bin           = @["deribo"]
binDir        = "bin"

# Dependencies

requires "nim >= 2.2.0"
requires "cligen >= 1.7.0"

task test, "Run the deribo test suite":
  exec "nim c -r --hints:off --mm:orc test/test_kmer.nim"
  exec "nim c -r --hints:off --mm:orc test/test_index_classify.nim"
  exec "nim c -r --hints:off --mm:orc test/test_fastx.nim"
