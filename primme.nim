##  Nim wrapper for PRIMME.
##
##  The C library can be found at https://github.com/primme/primme

import os
const
  homeDir = getHomeDir()
  primmeDir {.strdefine.} = homeDir&"/pkgs/src/primme"
  lapackLib {.strdefine.} = homeDir&"/pkg/lib/libopenblas.a -fopenmp -lm -lgfortran"
{.passC: "-I"&primmeDir&"/include".}
{.passL: primmeDir&"/lib/libprimme.a "&lapackLib.}

type complex[T] = tuple[re,im:T]
converter toComplex[T](x:T):complex[T] =
  result.re = x
  result.im = 0
proc `-`*[T](x:complex[T]):complex[T] =
  result.re = -x.re
  result.im = -x.im
proc `+`*[T](x:T, y:complex[T]):complex[T] =
  result.re = x+y.re
  result.im = y.im
proc `+`*[T](x:complex[T], y:T):complex[T] =
  result.re = x.re+y
  result.im = x.im
proc `+`*[T](x:complex[T], y:complex[T]):complex[T] =
  result.re = x.re+y.re
  result.im = x.im+y.im
proc `-`*[T](x:T, y:complex[T]):complex[T] =
  result.re = x-y.re
  result.im = -y.im
proc `-`*[T](x:complex[T], y:T):complex[T] =
  result.re = x.re-y
  result.im = x.im
proc `-`*[T](x:complex[T], y:complex[T]):complex[T] =
  result.re = x.re-y.re
  result.im = x.im-y.im
proc `*`*[T](x:T, y:complex[T]):complex[T] =
  result.re = x*y.re
  result.im = x*y.im
proc `*`*[T](x:complex[T], y:T):complex[T] =
  result.re = x.re*y
  result.im = x.im*y
proc `*`*[T](x:complex[T], y:complex[T]):complex[T] =
  result.re = x.re*y.re - x.im*y.im
  result.im = x.im*y.re + x.re*y.im
proc `/`*[T](x:T, y:complex[T]):complex[T] =
  if y.re.abs < y.im.abs:
    let
      r = y.re / y.im
      den = y.im + r*y.re
    result.re = x*r/den
    result.im = -x/den
  else:
    let
      r = y.im / y.re
      den = y.re + r*y.im
    result.re = x/den
    result.im = -x*r/den
proc `/`*[T](x:complex[T], y:T):complex[T] =
  result.re = x.re/y
  result.im = x.im/y
proc `/`*[T](x:complex[T], y:complex[T]):complex[T] =
  if y.re.abs < y.im.abs:
    let
      r = y.re / y.im
      den = y.im + r*y.re
    result.re = (x.re*r + x.im)/den
    result.im = (x.im*r - x.re)/den
  else:
    let
      r = y.im / y.re
      den = y.re + r*y.im
    result.re = (x.re + r*x.im)/den
    result.im = (x.im - r*x.re)/den
proc `+=`*[T](x:var complex[T], y:T) =
  x.re += y
proc `+=`*[T](x:var complex[T], y:complex[T]) =
  x.re += y.re
  x.im += y.im
proc `-=`*[T](x:var complex[T], y:T) =
  x.re -= y
proc `-=`*[T](x:var complex[T], y:complex[T]) =
  x.re -= y.re
  x.im -= y.im
proc `*=`*[T](x:var complex[T], y:T) =
  x.re *= y
  x.im *= y
proc `*=`*[T](x:var complex[T], y:complex[T]) =
  let im = x.re*y.im + x.im*y.re
  x.re = x.re*y.re - x.im*y.im
  x.im = im
proc `/=`*[T](x:var complex[T], y:T) =
  x.re /= y
  x.im /= y
proc `/=`*[T](x:var complex[T], y:complex[T]) =
  if y.re.abs < y.im.abs:
    let
      r = y.re / y.im
      den = y.im + r*y.re
      im = (x.im*r - x.re)/den
    x.re = (x.re*r + x.im)/den
    x.im = im
  else:
    let
      r = y.im / y.re
      den = y.re + r*y.im
      im = (x.im - r*x.re)/den
    x.re = (x.re + r*x.im)/den
    x.im = im
type PRIMME_INT* = int64
type PRIMME_COMPLEX_FLOAT*{.importc: "float complex", header: "complex.h".} = complex[cfloat]
type PRIMME_COMPLEX_DOUBLE*{.importc: "double complex", header: "complex.h".} = complex[cdouble]

import primme_eigs, primme_svds

proc run*(primme: var primme_params;
          evals: var openarray[cfloat]; evecs: var openarray[cfloat];
          resNorms: var openarray[cfloat]): int =
  sprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cfloat]; evecs: var openarray[complex[cfloat]];
          resNorms: var openarray[cfloat]): int =
  cprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cdouble]; evecs: var openarray[cdouble];
          resNorms: var openarray[cdouble]): int =
  dprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cdouble]; evecs: var openarray[complex[cdouble]];
          resNorms: var openarray[cdouble]): int =
  zprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cfloat]; svecs: var openarray[cfloat];
          resNorms: var openarray[cfloat]): int =
  sprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cfloat]; svecs: var openarray[complex[cfloat]];
          resNorms: var openarray[cfloat]): int =
  cprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cdouble]; svecs: var openarray[cdouble];
          resNorms: var openarray[cdouble]): int =
  dprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cdouble]; svecs: var openarray[complex[cdouble]];
          resNorms: var openarray[cdouble]): int =
  zprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc primme_initialize*:primme_params =
  result.addr.primme_initialize
proc primme_svds_initialize*:primme_svds_params =
  result.addr.primme_svds_initialize
proc set_method*(primme: var primme_params; preset_method: primme_preset_method): int =
  primme_set_method(preset_method, primme.addr).int
proc set_method*(primme: var primme_svds_params; preset_method: primme_svds_preset_method;
                 methodStage1: primme_preset_method; methodStage2: primme_preset_method): int =
  primme_svds_set_method(preset_method, methodStage1, methodStage2, primme.addr).int
template display_params*(primme: primme_params) = primme.primme_display_params
template display_params*(primme: primme_svds_params) = primme.primme_svds_display_params
proc free*(primme: var primme_params) = primme.addr.primme_free
proc free*(primme: var primme_svds_params) = primme.addr.primme_svds_free

template asarray*[T](p:pointer):auto =
  type A{.unchecked.} = array[0..0,T]
  cast[ptr A](p)

export primme_eigs, primme_svds

when isMainModule:
  import strutils, random
  template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
  #----------------------------------------------------------------------
  # Examples of eigenvalue problem.
  #----------------------------------------------------------------------
  proc laplacianMatVec[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                          primme:ptr primme_params, err:ptr cint) {.noconv.} =
    var
      x = asarray[T] x
      dx = ldx[]
      y = asarray[T] y
      dy = ldy[]
    for i in 0..<blocksize[]:
      let
        idx = i*dx
        idy = i*dy
      for row in 0..<primme.n:
        let
          nx = row + idx
          ny = row + idy
        y[ny] = 2.0*x[nx]
        if row>0: y[ny] += -x[nx-1]
        if row+1<primme.n: y[ny] += -x[nx+1]
    err[] = 0
  proc laplacianPrecond[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                           primme:ptr primme_params, err:ptr cint) {.noconv.} =
    var
      x = asarray[T] x
      dx = ldx[]
      y = asarray[T] y
      dy = ldy[]
    for i in 0..<blocksize[]:
      let
        idx = i*dx
        idy = i*dy
      for row in 0..<primme.n:
        y[row + idx] = x[row + idy]/2.0
    err[] = 0
  proc ini:primme_params =
    var primme = primme_initialize() # Set default values in primme_params object
    primme.n = 100                   # Problem dimension
    primme.numEvals = 10            # Number of wanted eigenpairs
    primme.eps = 1e-9               # ||r|| <= eps * ||matrix||
    primme.target = primme_smallest # Want the smallest eigenvalues
    primme.correctionParams.precondition = 1 # Use the preconditioner
    # primme.maxBasisSize = 14                 # Optional
    # primme.minRestartSize = 4                # Optional
    # primme.maxBlockSize = 1                  # Optional
    # primme.maxMatvecs = 1000                 # Optional
    # let ret = primme.primme_set_method PRIMME_DEFAULT_MIN_TIME # Method to solve the problem
    let ret = primme.set_method PRIMME_DYNAMIC
    if 0 != ret: # Method to solve the problem
      echo "Error: set_method returned with nonzero exit status: ", ret
      quit QuitFailure
    primme.display_params              # Optional display PRIMME configuration struct
    return primme
  proc runreport(primme:var auto, evals:var auto, evecs:var auto, rnorms:var auto) =
    let ret = primme.run(evals, evecs, rnorms)
    if ret != 0:
      echo "Error: primme returned with nonzero exit status: ", ret
      quit QuitFailure
    for i in 0..<primme.initSize:
      echo "Eval[",i,"]: ",evals[i].ff," rnorm: ",rnorms[i].ff
    echo " ",primme.initSize," eigenpairs converged"
    echo "Tolerance  : ",ff primme.aNorm*primme.eps
    echo "Iterations : ",primme.stats.numOuterIterations
    echo "Restarts   : ",primme.stats.numRestarts
    echo "Matvecs    : ",primme.stats.numMatvecs
    echo "Preconds   : ",primme.stats.numPreconds
    if primme.locking != 0 and primme.intWork != nil and primme.intWork[] == 1:
      echo "\nA locking problem has occurred."
      echo "Some eigenpairs do not have a residual norm less than the tolerance."
      echo "However, the subspace of evecs is accurate to the required tolerance."
    case primme.dynamicMethodSwitch:
    of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
    of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
    of -3: echo "Recommended method for next run: DYNAMIC (close call)"
    else: discard
  block doublePrecision:
    var primme = ini()
    primme.matrixMatvec = laplacianMatVec[float] # Function for matrix-vector product A*x for solving A*x=l*x
    primme.applyPreconditioner = laplacianPrecond[float] # Optional preconditioner
    var
      evals = newseq[float](primme.numEvals)
      evecs = newseq[float](primme.n * primme.numEvals)
      rnorms = newseq[float](primme.numEvals)
    primme.runreport evals, evecs, rnorms
    primme.free                    # Free required before changing problem (matrix, type, domain)
  block doublePrecisionComplex:
    var primme = ini()
    primme.matrixMatvec = laplacianMatVec[complex[float]] # Function for matrix-vector product A*x for solving A*x=l*x
    primme.applyPreconditioner = laplacianPrecond[complex[float]] # Optional preconditioner
    var
      evals = newseq[float](primme.numEvals)
      evecs = newseq[complex[float]](primme.n * primme.numEvals)
      rnorms = newseq[float](primme.numEvals)
    primme.runreport evals, evecs, rnorms
    # Call primme again with different parameters
    # Find the 5 eigenpairs closest to 0.5
    primme.numTargetShifts = 1
    var targetShifts = 0.5
    primme.targetShifts = targetShifts.addr
    primme.target = primme_closest_abs
    primme.numEvals = 5
    primme.initSize = 0 # initSize may be nonzero after primme; set it to zero to avoid reusing previous eigenvectors
    primme.runreport evals, evecs, rnorms
    # Perturb the approximate eigenvectors in evecs and use them as initial solution.
    for i in 0..<int(primme.n*primme.numEvals): evecs[i] += random(1.0)*1e-4
    primme.initSize = primme.numEvals
    primme.runreport evals, evecs, rnorms
    # Find the next 5 eigenpairs cloest to 0.5
    primme.initSize = 0
    primme.numEvals = 5
    primme.numOrthoConst = 5 # Solver will find solutions orthogonal to the eigenvectors in evecs
    primme.runreport evals, evecs, rnorms
    primme.free
  #----------------------------------------------------------------------
  # Examples of singular value problem.
  #----------------------------------------------------------------------
  #[ lauchli block matrix-vector product, y = a * x (or y = a^t * x), where
   - x, input dense matrix of size primme_svds.n (or primme_svds.m) x blocksize;
   - y, output dense matrix of size primme_svds.m (or primme_svds.n) x blocksize;
   - a, lauchli matrix of dimensions primme_svds.m x (primme_svds.m+1) with this form:
        [ 1  1  1  1  1 ...   1 ],  ei = 1 - (1 - mu)*i/(min(m,n) - 1)
        [e0  0  0  0  0 ...   0 ]
        [ 0 e1  0  0  0 ...   0 ]
         ...
        [ 0  0  0  0  0 ... en-1]
  ]#
  proc lauchliMatvec[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                        transpose:ptr cint, primme:ptr primme_svds_params, err:ptr cint) {.noconv.} =
    let
      x = asarray[T] x
      dx = ldx[]
      y = asarray[T] y
      dy = ldy[]
      minMN = min(primme.m, primme.n)
      mu = cast[ptr cdouble](primme.matrix)[]
    if 0 == transpose[]:        # Do y <- A * x
      for i in 0..<blocksize[]:
        let
          x = asarray[T] x[i*dx].addr
          y = asarray[T] y[i*dy].addr
        y[0] = 0.0
        for j in 0..<primme.n: y[0] += x[j]
        for j in 1..<primme.m: y[j] = (if j-1<primme.n: x[j-1]*(1 - (1-mu)*(j-1).float/(minMN-1).float) else: 0.0)
    else:                       # Do y <- A^t * x
      for i in 0..<blocksize[]:
        let
          x = asarray[T] x[i*dx].addr
          y = asarray[T] y[i*dy].addr
        for j in 0..<primme.n:
          y[j] = x[0]
          if j+1<primme.m: y[j] += x[j+1]*(1 - (1-mu)*j.float/(minMN-1).float)
    err[] = 0
  #[ This performs Y = M^{-1} * X, where
   - X, input dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - Y, output dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - M, preconditioner for A^t*A (or A*A^t or [0 A^t; A 0]), where A is the Lauchli matrix.
  ]#
  proc lauchliPrecond[T](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT, blocksize:ptr cint,
                         mode:ptr cint, primme:ptr primme_svds_params, err:ptr cint) {.noconv.} =
    let
      xa = asarray[T] x
      dx = ldx[]
      ya = asarray[T] y
      dy = ldy[]
      mode = mode[].primme_svds_operator
      mu = cast[ptr cdouble](primme.matrix)[]
      minMN = min(primme.m, primme.n)
    var
      modeAtA = primme_svds_op_AtA.cint
      modeAAt = primme_svds_op_AAt.cint
      notrans:cint = 0
      trans:cint = 1
    if primme_svds_op_AtA == mode:
      # Preconditioner for A^t*A, diag(A^t*A)^{-1}
      for i in 0..<blocksize[]:
        let
          x = asarray[T] xa[i*dx].addr
          y = asarray[T] ya[i*dy].addr
        for j in 0..<primme.n:
          let ei = (if j<primme.m: 1 - (1-mu)*j.float/(minMN-1).float else: 0)
          y[j] = x[j]/(1 + ei*ei)
    elif primme_svds_op_AAt == mode:
      # Preconditioner for A*A^t, diag(A*A^t)^{-1}
      for i in 0..<blocksize[]:
        let
          x = asarray[T] xa[i*dx].addr
          y = asarray[T] ya[i*dy].addr
        y[0] = x[0]/primme.m.float
        for j in 1..<primme.m:
          let ei = (if j<primme.n: 1 - (1-mu)*j.float/(minMN-1).float else: 1)
          y[j] = x[j]/ei/ei
    elif primme_svds_op_augmented == mode:
      # Preconditioner for [0 A^t; A 0], [diag(A^t*A) 0; 0 diag(A*A^t)]^{-1}*[0 A^t; A 0] */
      # [y0; y1] <- [0 A^t; A 0] * [x0; x1]
      var
        daux:PRIMME_INT = primme.n + primme.m
        aux = newseq[T](blocksize[]*daux)
      primme.matrixMatvec(x, ldx, aux[primme.n.int].addr, daux.addr, blocksize, notrans.addr, primme, err)
      primme.matrixMatvec(xa[primme.n].addr, ldx, aux[0].addr, daux.addr, blocksize, trans.addr, primme, err)
      # y0 <- preconditioner for A^t*A * y0
      lauchliPrecond[float](aux[0].addr, daux.addr, y, ldy, blocksize, modeAtA.addr, primme, err)
      # y1 <- preconditioner for A*A^t * y1
      lauchliPrecond[float](aux[primme.n.int].addr, daux.addr, ya[primme.n].addr, ldy, blocksize, modeAAt.addr, primme, err)
    err[] = 0
  block doublePrecisionSVDS:
    var
      mu = 1e-5
      primme = primme_svds_initialize() # Set default values in primme_svds_params object
    primme.matrixMatvec = lauchliMatvec[float] # Function for matrix-vector products A*x and A^t*x
    primme.matrix = mu.addr
    primme.m = 500
    primme.n = 100                # Problem dimension
    primme.numSvals = 4           # Number of wanted singular values
    primme.eps = 1e-12            # ||r|| <= eps * ||matrix||
    primme.target = primme_svds_smallest # Seeking for the smallest singular values
    primme.applyPreconditioner = lauchliPrecond[float] # Preconditioner (optional)
    block:
      # Method to solve the singular value problem and the underneath eigenvalue problem (optional)
      let ret = primme.set_method(primme_svds_default, PRIMME_DEFAULT_METHOD, PRIMME_DEFAULT_METHOD)
      # primme_svds_default: devs choice, now being hybrid, which first solve the normal equation
      # and then the augmented problem.
      if 0 != ret:
        echo "Error: set_method returned with nonzero exit status: ", ret
        quit QuitFailure
    primme.printLevel = 3
    # Configuration for 1st stage */
    #[
    primme.primme.maxBasisSize = 14;
    primme.primme.minRestartSize = 6;
    primme.primme.maxBlockSize = 2;
    ]#
    # Configuration for 2nd stage */
    #[
    primme.primmeStage2.maxBasisSize = 30;
    primme.primmeStage2.minRestartSize = 15;
    primme.primmeStage2.maxBlockSize = 1;
    ]#
    primme.display_params
    var
      svals = newseq[float](primme.numSvals)
      svecs = newseq[float]((primme.n+primme.m)*primme.numSvals)
      rnorms = newseq[float](primme.numSvals)
    block:
      let ret = primme.run(svals, svecs, rnorms)
      if 0 != ret:
        echo "Error: primme returned with nonzero exit status: ", ret
        quit QuitFailure
    for i in 0..<primme.initSize:
      echo "Sval[",i,"]: ",svals[i].ff," rnorm: ",rnorms[i].ff
    echo " ",primme.initSize," singular triplets converged"
    echo "Tolerance  : ",ff primme.aNorm*primme.eps
    echo "Iterations : ",primme.stats.numOuterIterations
    echo "Restarts   : ",primme.stats.numRestarts
    echo "Matvecs    : ",primme.stats.numMatvecs
    echo "Preconds   : ",primme.stats.numPreconds
    if primme.primme.locking != 0 and primme.primme.intWork != nil and primme.primme.intWork[] == 1:
      echo "\nA locking problem has occurred."
      echo "Some triplets do not have a residual norm less than the tolerance."
      echo "However, the subspace of evecs is accurate to the required tolerance."
    primme.free
  block doublePrecisionComplexSVDS:
    var
      mu = 1e-5
      primme = primme_svds_initialize() # Set default values in primme_svds_params object
    primme.matrixMatvec = lauchliMatvec[complex[float]] # Function for matrix-vector products A*x and A^t*x
    primme.matrix = mu.addr
    primme.m = 500
    primme.n = 100              # Problem dimension
    primme.numSvals = 4         # Number of wanted singular values
    primme.eps = 1e-12          # ||r|| <= eps * ||matrix||
    primme.target = primme_svds_smallest # Seeking for the smallest singular values
    primme.applyPreconditioner = lauchliPrecond[complex[float]] # Set preconditioner (optional)
    block:
      # Method to solve the singular value problem and the underneath eigenvalue problem (optional)
      let ret = primme.set_method(primme_svds_hybrid, PRIMME_DYNAMIC, PRIMME_DEFAULT_MIN_TIME)
      if 0 != ret:
        echo "Error: set_method returned with nonzero exit status: ", ret
        quit QuitFailure
    primme.printLevel = 3
    # Set parameters for the underneath eigensolver if you know what you are doing (opitonal)
    primme.primme.locking = 1
    primme.primme.restartingParams.maxPrevRetain = 3
    # Display PRIMME SVDS configuration struct (optional)
    primme.display_params
    var
      svals = newseq[float](primme.numSvals)
      svecs = newseq[complex[float]]((primme.n+primme.m)*primme.numSvals)
      rnorms = newseq[float](primme.numSvals)
    block:
      let ret = primme.run(svals, svecs, rnorms)
      if 0 != ret:
        echo "Error: primme returned with nonzero exit status: ", ret
        quit QuitFailure
    for i in 0..<primme.initSize:
      echo "Sval[",i,"]: ",svals[i].ff," rnorm: ",rnorms[i].ff
    echo " ",primme.initSize," singular triplets converged"
    echo "Tolerance  : ",ff primme.aNorm*primme.eps
    echo "Iterations : ",primme.stats.numOuterIterations
    echo "Restarts   : ",primme.stats.numRestarts
    echo "Matvecs    : ",primme.stats.numMatvecs
    echo "Preconds   : ",primme.stats.numPreconds
    if primme.primme.locking != 0 and primme.primme.intWork != nil and primme.primme.intWork[] == 1:
      echo "\nA locking problem has occurred."
      echo "Some triplets do not have a residual norm less than the tolerance."
      echo "However, the subspace of evecs is accurate to the required tolerance."
    primme.free
