switch("path", "../src/")

const
  primmeDir {.strdefine.} = "$HOME/src/primme"
  primmeLib {.strdefine.} = primmeDir&"/lib/libprimme.a"
  lapackLib {.strdefine.} = "/usr/local/lib/libopenblas.a /usr/local/lib/libopenlibm.a -lgfortran -L/usr/local/lib/gcc8"
switch("define", "primmeDir="&primmeDir)
switch("define", "primmeLib="&primmeLib)
switch("define", "lapackLib="&lapackLib)
