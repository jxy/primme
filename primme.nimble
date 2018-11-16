version     = "2.1.0"
author      = "Xiao-Yong Jin"
description = "Nim interface for PRIMME: PReconditioned Iterative MultiMethod Eigensolver"
license     = "MIT"
srcDir      = "src"

requires "nim >= 0.19.9"

task test, "Runs the test suite":
  withDir "test":
    exec "nim c -r -d:release t"

before test:
  echo "#"
  echo "#     Modify `test/config.nims' for compiler and link flags for PRIMME"
  echo "#     Test for singular values may fail because of the runtime dynamic method"
  echo "#"

after install:
  echo "#"
  echo "#     To compile and link against PRIMME, define `primmeDir' and `lapackLib' properly."
  echo "#     For example:"
  echo "#        nim c -d:primmeDir='/path/to/primme' \\"
  echo "#              -d:lapackLib='-L/path/to/lapack/lib -llapack -lblas' \\"
  echo "#              -d:release \\"
  echo "#              YourApp"
  echo "#     Check the beginning of the file `primme.nim' for details."
  echo "#"
