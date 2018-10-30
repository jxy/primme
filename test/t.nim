import primme, primme/ccomplex
import strutils, random, unittest
const REPORT {.intdefine.} = 0
type myComplex = ccomplex[float]
template ff(x:untyped):auto = formatFloat(x,ffScientific,17)
var CT = 1e-10                  # comparison tolerance
proc `=~`(x,y:float):bool = abs(x-y)/max(abs(x),abs(y)) < CT
#----------------------------------------------------------------------
# Examples of eigenvalue problem.
#----------------------------------------------------------------------
#[ 1-D Laplacian block matrix-vector product, Y = A * X, where
   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - A, tridiagonal square matrix of dimension primme.n with this form:
        [ 2 -1  0  0  0 ... ]
        [-1  2 -1  0  0 ... ]
        [ 0 -1  2 -1  0 ... ]
         ...
]#
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
#[ This performs Y = M^{-1} * X, where
   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
]#
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
proc runreport(primme:var primme_params, evals:var auto, evecs:var auto, rnorms:var auto) =
  if 0 != REPORT:
    primme.display_params # Optional display PRIMME configuration struct
  let ret = primme.run(evals, evecs, rnorms)
  if ret != 0:
    echo "Error: primme returned with nonzero exit status: ", ret
    quit QuitFailure
  if 0 == REPORT: return
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
const
  ev = [
    9.674354160238160E-04,
    3.868805732811409E-03,
    8.701304061963075E-03,
    1.546025527344750E-02,
    2.413912051848673E-02,
    3.472950355547286E-02,
    4.722115887278573E-02,
    6.160200160066791E-02,
    7.785811920255102E-02,
    9.597378493454023E-02]
  ev5 = [
    4.903541216934860E-01,
    5.318829424810771E-01,
    4.502857857942200E-01,
    5.748320717049855E-01,
    4.117166983104929E-01]
  ev5n = [
    6.191599588565045E-01,
    3.746841723435005E-01,
    3.392240344704057E-01,
    6.648237195676904E-01,
    3.053705900844449E-01]
suite "Eigenvalues":
  setup:
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
  teardown:
    primme.free
  test "double precision matrix":
    primme.matrixMatvec = laplacianMatVec[float] # Function for matrix-vector product A*x for solving A*x=l*x
    primme.applyPreconditioner = laplacianPrecond[float] # Optional preconditioner
    var
      evals = newseq[float](primme.numEvals)
      evecs = newseq[float](primme.n * primme.numEvals)
      rnorms = newseq[float](primme.numEvals)
    primme.runreport evals, evecs, rnorms
    check:
      ev[0] =~ evals[0]
      ev[1] =~ evals[1]
      ev[2] =~ evals[2]
      ev[3] =~ evals[3]
      ev[4] =~ evals[4]
      ev[5] =~ evals[5]
      ev[6] =~ evals[6]
      ev[7] =~ evals[7]
      ev[8] =~ evals[8]
      ev[9] =~ evals[9]
  test "double precision complex matrix":
    primme.matrixMatvec = laplacianMatVec[myComplex] # Function for matrix-vector product A*x for solving A*x=l*x
    primme.applyPreconditioner = laplacianPrecond[myComplex] # Optional preconditioner
    var
      evals = newseq[float](primme.numEvals)
      evecs = newseq[myComplex](primme.n * primme.numEvals)
      rnorms = newseq[float](primme.numEvals)
    primme.runreport evals, evecs, rnorms
    check:
      ev[0] =~ evals[0]
      ev[1] =~ evals[1]
      ev[2] =~ evals[2]
      ev[3] =~ evals[3]
      ev[4] =~ evals[4]
      ev[5] =~ evals[5]
      ev[6] =~ evals[6]
      ev[7] =~ evals[7]
      ev[8] =~ evals[8]
      ev[9] =~ evals[9]
    # Call primme again with different parameters
    # Find the 5 eigenpairs closest to 0.5
    primme.numTargetShifts = 1
    var targetShifts = 0.5
    primme.targetShifts = targetShifts.addr
    primme.target = primme_closest_abs
    primme.numEvals = 5
    primme.initSize = 0 # initSize may be nonzero after primme; set it to zero to avoid reusing previous eigenvectors
    primme.runreport evals, evecs, rnorms
    check:
      ev5[0] =~ evals[0]
      ev5[1] =~ evals[1]
      ev5[2] =~ evals[2]
      ev5[3] =~ evals[3]
      ev5[4] =~ evals[4]
    # Perturb the approximate eigenvectors in evecs and use them as initial solution.
    for i in 0..<int(primme.n*primme.numEvals): evecs[i] += rand(1.0)*1e-4
    primme.initSize = primme.numEvals
    primme.runreport evals, evecs, rnorms
    check:
      ev5[0] =~ evals[0]
      ev5[1] =~ evals[1]
      ev5[2] =~ evals[2]
      ev5[3] =~ evals[3]
      ev5[4] =~ evals[4]
    # Find the next 5 eigenpairs cloest to 0.5
    primme.initSize = 0
    primme.numEvals = 5
    primme.numOrthoConst = 5 # Solver will find solutions orthogonal to the eigenvectors in evecs
    primme.runreport evals, evecs, rnorms
    check:
      ev5n[0] =~ evals[0]
      ev5n[1] =~ evals[1]
      ev5n[2] =~ evals[2]
      ev5n[3] =~ evals[3]
      ev5n[4] =~ evals[4]
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
proc runreport(primme:var primme_svds_params, svals:var auto, svecs:var auto, rnorms:var auto) =
  if 0 != REPORT:
    primme.display_params # Display PRIMME SVDS configuration struct (optional)
  block:
    let ret = primme.run(svals, svecs, rnorms)
    if 0 != ret:
      echo "Error: primme returned with nonzero exit status: ", ret
      quit QuitFailure
  if 0 == REPORT: return
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
const sv = [
  6.537385514282974E-03,
  1.583473739145508E-02,
  2.571238674965088E-02,
  3.572252806311886E-02]
suite "Singular values":
  setup:
    var
      mu = 1e-5
      primme = primme_svds_initialize() # Set default values in primme_svds_params object
    primme.matrix = mu.addr
    primme.m = 500
    primme.n = 100                # Problem dimension
    primme.numSvals = 4           # Number of wanted singular values
    primme.eps = 1e-12            # ||r|| <= eps * ||matrix||
    primme.target = primme_svds_smallest # Seeking for the smallest singular values
    if 0 != REPORT: primme.printLevel = 3
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
  teardown:
    primme.free
  test "double precision matrix":
    primme.matrixMatvec = lauchliMatvec[float] # Function for matrix-vector products A*x and A^t*x
    primme.applyPreconditioner = lauchliPrecond[float] # Preconditioner (optional)
    block:
      # Method to solve the singular value problem and the underneath eigenvalue problem (optional)
      let ret = primme.set_method(primme_svds_default, PRIMME_DEFAULT_METHOD, PRIMME_DEFAULT_METHOD)
      # primme_svds_default: devs choice, now being hybrid, which first solve the normal equation
      # and then the augmented problem.
      if 0 != ret:
        echo "Error: set_method returned with nonzero exit status: ", ret
        quit QuitFailure
    var
      svals = newseq[float](primme.numSvals)
      svecs = newseq[float]((primme.n+primme.m)*primme.numSvals)
      rnorms = newseq[float](primme.numSvals)
    primme.runreport svals, svecs, rnorms
    check:
      sv[0] =~ svals[0]
      sv[1] =~ svals[1]
      sv[2] =~ svals[2]
      sv[3] =~ svals[3]
  test "double precision complex matrix":
    primme.matrixMatvec = lauchliMatvec[myComplex] # Function for matrix-vector products A*x and A^t*x
    primme.applyPreconditioner = lauchliPrecond[myComplex] # Set preconditioner (optional)
    block:
      # Method to solve the singular value problem and the underneath eigenvalue problem (optional)
      let ret = primme.set_method(primme_svds_hybrid, PRIMME_DYNAMIC, PRIMME_DEFAULT_MIN_TIME)
      if 0 != ret:
        echo "Error: set_method returned with nonzero exit status: ", ret
        quit QuitFailure
    # Set parameters for the underneath eigensolver if you know what you are doing (opitonal)
    primme.primme.locking = 1
    primme.primme.restartingParams.maxPrevRetain = 3
    var
      svals = newseq[float](primme.numSvals)
      svecs = newseq[myComplex]((primme.n+primme.m)*primme.numSvals)
      rnorms = newseq[float](primme.numSvals)
    primme.runreport svals, svecs, rnorms
    check:
      sv[0] =~ svals[0]
      sv[1] =~ svals[1]
      sv[2] =~ svals[2]
      sv[3] =~ svals[3]
