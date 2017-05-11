switch("path", "../src/")

const
  primmeDir {.strdefine.} = "$HOME/pkgs/src/primme"
  primmeLib {.strdefine.} = primmeDir&"/lib/libprimme.a"
  lapackLib {.strdefine.} = "$HOME/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran"
switch("define", "primmeDir="&primmeDir)
switch("define", "primmeLib="&primmeLib)
switch("define", "lapackLib="&lapackLib)
