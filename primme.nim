##  Nim wrapper for PRIMME.
##
##  The C library can be found at https://github.com/primme/primme

import os
const
  homeDir = getHomeDir()
  primmeDir = homeDir&"/pkgs/src/primme"
  lapackLib = homeDir&"/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran"
{.passC: "-I"&primmeDir&"/include".}
{.passL: primmeDir&"/lib/libprimme.a "&lapackLib.}

type complex[T] = object
  re,im: T
type PRIMME_INT* = int64
type PRIMME_COMPLEX_FLOAT*{.importc: "float complex", header: "complex.h".} = complex[cfloat]
type PRIMME_COMPLEX_DOUBLE*{.importc: "double complex", header: "complex.h".} = complex[cdouble]

import primme_eigs, primme_svds

proc primme*(evals: var openarray[cfloat]; evecs: var openarray[cfloat]; resNorms: var openarray[cfloat];
             param: var primme_params): int =
  sprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, param.addr)
proc primme*(evals: var openarray[cfloat]; evecs: var openarray[complex[cfloat]]; resNorms: var openarray[cfloat];
             param: var primme_params): int =
  cprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, param.addr)
proc primme*(evals: var openarray[cdouble]; evecs: var openarray[cdouble]; resNorms: var openarray[cdouble];
             param: var primme_params): int =
  dprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, param.addr)
proc primme*(evals: var openarray[cdouble]; evecs: var openarray[complex[cdouble]]; resNorms: var openarray[cdouble];
             param: var primme_params): int =
  zprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, param.addr)
proc primme_initialize*:primme_params =
  result.addr.primme_initialize
  result
proc set_method*(params: var primme_params; preset_method: primme_preset_method): int =
  primme_set_method(preset_method, params.addr).int
template display_params*(params: primme_params) = params.primme_display_params
proc free*(param: var primme_params) = param.addr.primme_free

template asarray*[T](p:pointer):auto =
  type A{.unchecked.} = array[0..0,T]
  cast[ptr A](p)

export primme_eigs, primme_svds

when isMainModule:
  import strutils
  proc laplacianMatVec[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                          params:ptr primme_params, err:ptr cint) {.noconv.} =
    let
      x = asarray[T] x
      dx = ldx[]
      y = asarray[T] y
      dy = ldy[]
    for i in 0..<blocksize[]:
      let
        idx = i*dx
        idy = i*dy
      for row in 0..<params.n:
        let
          nx = row + idx
          ny = row + idy
        y[ny] = 2*x[nx]
        if row>0: y[ny] += -x[nx-1]
        if row+1<params.n: y[ny] += -x[nx+1]
    err[] = 0
  proc laplacianPrecond[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                           params:ptr primme_params, err:ptr cint) {.noconv.} =
    let
      x = asarray[T] x
      dx = ldx[]
      y = asarray[T] y
      dy = ldy[]
    for i in 0..<blocksize[]:
      let
        idx = i*dx
        idy = i*dy
      for row in 0..<params.n:
        y[row + idx] = x[row + idy]/2.0
    err[] = 0
  var param = primme_initialize() # Set default values in primme
  param.matrixMatvec = laplacianMatVec[float] # Function for matrix-vector product A*x for solving A*x=l*x
  param.applyPreconditioner = laplacianPrecond[float] # Optional preconditioner
  param.n = 100                  # Problem dimension
  param.numEvals = 10            # Number of wanted eigenpairs
  param.eps = 1e-9               # ||r|| <= eps * ||matrix||
  param.target = primme_smallest # Want the smallest eigenvalues
  param.correctionParams.precondition = 1 # Use the preconditioner
  # param.maxBasisSize = 14                 # Optional
  # param.minRestartSize = 4                # Optional
  # param.maxBlockSize = 1                  # Optional
  # param.maxMatvecs = 1000                 # Optional
  # var ret = param.primme_set_method PRIMME_DEFAULT_MIN_TIME # Method to solve the problem
  var ret = param.set_method PRIMME_DYNAMIC
  if 0 != ret: # Method to solve the problem
    echo "Error: set_method returned with nonzero exit status: ", ret
    quit QuitFailure
  param.display_params              # Optional display PRIMME configuration struct
  var
    evals = newseq[float](param.numEvals)
    evecs = newseq[float](param.n * param.numEvals)
    rnorms = newseq[float](param.numEvals)
  ret = primme(evals, evecs, rnorms, param)
  if ret != 0:
    echo "Error: primme returned with nonzero exit status: ", ret
    quit QuitFailure
  # Reporting (optional)
  template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
  for i in 0..<param.initSize:
    echo "Eval[",i,"]: ",evals[i].ff," rnorm: ",rnorms[i].ff
  echo " ",param.initSize," eigenpairs converged"
  echo "Tolerance  : ",ff param.aNorm*param.eps
  echo "Iterations : ",param.stats.numOuterIterations
  echo "Restarts   : ",param.stats.numRestarts
  echo "Matvecs    : ",param.stats.numMatvecs
  echo "Preconds   : ",param.stats.numPreconds
  if param.locking != 0 and param.intWork != nil and param.intWork[] == 1:
    echo "\nA locking problem has occurred."
    echo "Some eigenpairs do not have a residual norm less than the tolerance."
    echo "However, the subspace of evecs is accurate to the required tolerance."
  case param.dynamicMethodSwitch:
  of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
  of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
  of -3: echo "Recommended method for next run: DYNAMIC (close call)"
  else: discard
  param.free
