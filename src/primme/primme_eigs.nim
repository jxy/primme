## ******************************************************************************
##  Copyright (c) 2017, Xiao-Yong Jin
##  Copyright (c) 2017, College of William & Mary
##  All rights reserved.
##
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met:
##      * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
##      * Redistributions in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the distribution.
##      * Neither the name of the College of William & Mary nor the
##        names of its contributors may be used to endorse or promote products
##        derived from this software without specific prior written permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
##  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
##  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
##  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
##  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
##  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
##  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
##  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##  PRIMME: https://github.com/primme/primme
##  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
## ******************************************************************************
##  File: primme_eigs.h
##
##  Purpose - Main header with the PRIMME EIGS C interface functions.
##
## ****************************************************************************

import primme

type
  primme_target* {.size: sizeof(cint).} = enum
    primme_smallest,          ##  leftmost eigenvalues
    primme_largest,           ##  rightmost eigenvalues
    primme_closest_geq,       ##  leftmost but greater than the target shift
    primme_closest_leq,       ##  rightmost but less than the target shift
    primme_closest_abs,       ##  the closest to the target shift
    primme_largest_abs        ##  the farthest to the target shift


##  projection methods for extraction

type                          ##  Initially fill up the search subspace with:
  primme_projection* {.size: sizeof(cint).} = enum
    primme_proj_default, primme_proj_RR, ##  Rayleigh-Ritz
    primme_proj_harmonic,     ##  Harmonic Rayleigh-Ritz
    primme_proj_refined       ##  refined with fixed target
  primme_init* {.size: sizeof(cint).} = enum
    primme_init_default, primme_init_krylov, ##  a) Krylov with the last vector provided by the user or random
    primme_init_random,       ##  b) just random vectors
    primme_init_user          ##  c) provided vectors or a single random vector
  primme_restartscheme* {.size: sizeof(cint).} = enum
    primme_thick, primme_dtr
  primme_convergencetest* {.size: sizeof(cint).} = enum
    primme_full_LTolerance, primme_decreasing_LTolerance,
    primme_adaptive_ETolerance, primme_adaptive





##  Identifies the type of event for which monitor is being called

type
  primme_event* {.size: sizeof(cint).} = enum
    primme_event_outer_iteration, ##  report at every outer iteration
    primme_event_inner_iteration, ##  report at every QMR iteration
    primme_event_restart,     ##  report at every basis restart
    primme_event_reset,       ##  event launch if basis reset
    primme_event_converged,   ##  report new pair marked as converged
    primme_event_locked       ##  report new pair marked as locked
  primme_stats* {.importc: "primme_stats", header: "primme.h".} = object
    numOuterIterations* {.importc: "numOuterIterations".}: PRIMME_INT
    numRestarts* {.importc: "numRestarts".}: PRIMME_INT
    numMatvecs* {.importc: "numMatvecs".}: PRIMME_INT
    numPreconds* {.importc: "numPreconds".}: PRIMME_INT
    numGlobalSum* {.importc: "numGlobalSum".}: PRIMME_INT ##  times called globalSumReal
    volumeGlobalSum* {.importc: "volumeGlobalSum".}: PRIMME_INT ##  number of SCALARs reduced by globalSumReal
    numOrthoInnerProds* {.importc: "numOrthoInnerProds".}: cdouble ##  number of inner prods done by Ortho
    elapsedTime* {.importc: "elapsedTime".}: cdouble
    timeMatvec* {.importc: "timeMatvec".}: cdouble ##  time expend by matrixMatvec
    timePrecond* {.importc: "timePrecond".}: cdouble ##  time expend by applyPreconditioner
    timeOrtho* {.importc: "timeOrtho".}: cdouble ##  time expend by ortho
    timeGlobalSum* {.importc: "timeGlobalSum".}: cdouble ##  time expend by globalSumReal
    estimateMinEVal* {.importc: "estimateMinEVal".}: cdouble ##  the leftmost Ritz value seen
    estimateMaxEVal* {.importc: "estimateMaxEVal".}: cdouble ##  the rightmost Ritz value seen
    estimateLargestSVal* {.importc: "estimateLargestSVal".}: cdouble ##  absolute value of the farthest to zero Ritz value seen
    maxConvTol* {.importc: "maxConvTol".}: cdouble ##  largest norm residual of a locked eigenpair
    estimateResidualError* {.importc: "estimateResidualError".}: cdouble ##  accumulated error in V and W

  JD_projectors* {.importc: "JD_projectors", header: "primme.h".} = object
    LeftQ* {.importc: "LeftQ".}: cint
    LeftX* {.importc: "LeftX".}: cint
    RightQ* {.importc: "RightQ".}: cint
    RightX* {.importc: "RightX".}: cint
    SkewQ* {.importc: "SkewQ".}: cint
    SkewX* {.importc: "SkewX".}: cint

  projection_params* {.importc: "projection_params", header: "primme.h".} = object
    projection* {.importc: "projection".}: primme_projection

  correction_params* {.importc: "correction_params", header: "primme.h".} = object
    precondition* {.importc: "precondition".}: cint
    robustShifts* {.importc: "robustShifts".}: cint
    maxInnerIterations* {.importc: "maxInnerIterations".}: cint
    projectors* {.importc: "projectors".}: JD_projectors
    convTest* {.importc: "convTest".}: primme_convergencetest
    relTolBase* {.importc: "relTolBase".}: cdouble

  restarting_params* {.importc: "restarting_params", header: "primme.h".} = object
    scheme* {.importc: "scheme".}: primme_restartscheme
    maxPrevRetain* {.importc: "maxPrevRetain".}: cint



## --------------------------------------------------------------------------

type
  primme_params* {.bycopy, importc: "primme_params", header: "primme.h".} = object
    n* {.importc: "n".}: PRIMME_INT ##  The user must input at least the following two arguments
    matrixMatvec* {.importc: "matrixMatvec".}: proc (x: pointer; ldx: ptr PRIMME_INT;
        y: pointer; ldy: ptr PRIMME_INT; blockSize: ptr cint; primme: ptr primme_params;
        ierr: ptr cint) {.noconv.}        ##  Preconditioner applied on block of vectors (if available)
    applyPreconditioner* {.importc: "applyPreconditioner".}: proc (x: pointer;
        ldx: ptr PRIMME_INT; y: pointer; ldy: ptr PRIMME_INT; blockSize: ptr cint;
        primme: ptr primme_params; ierr: ptr cint) {.noconv.} ##  Matrix times a multivector for mass matrix B for generalized Ax = xBl
    massMatrixMatvec* {.importc: "massMatrixMatvec".}: proc (x: pointer;
        ldx: ptr PRIMME_INT; y: pointer; ldy: ptr PRIMME_INT; blockSize: ptr cint;
        primme: ptr primme_params; ierr: ptr cint) {.noconv.} ##  input for the following is only required for parallel programs
    numProcs* {.importc: "numProcs".}: cint
    procID* {.importc: "procID".}: cint
    nLocal* {.importc: "nLocal".}: PRIMME_INT
    commInfo* {.importc: "commInfo".}: pointer
    globalSumReal* {.importc: "globalSumReal".}: proc (sendBuf: pointer;
        recvBuf: pointer; count: ptr cint; primme: ptr primme_params; ierr: ptr cint) {.noconv.} ## Though primme_initialize will assign defaults, most users will set these
    numEvals* {.importc: "numEvals".}: cint
    target* {.importc: "target".}: primme_target
    numTargetShifts* {.importc: "numTargetShifts".}: cint ##  For targeting interior epairs,
    targetShifts* {.importc: "targetShifts".}: ptr cdouble ##  at least one shift must also be set
                                                      ##  the following will be given default values depending on the method
    dynamicMethodSwitch* {.importc: "dynamicMethodSwitch".}: cint
    locking* {.importc: "locking".}: cint
    initSize* {.importc: "initSize".}: cint
    numOrthoConst* {.importc: "numOrthoConst".}: cint
    maxBasisSize* {.importc: "maxBasisSize".}: cint
    minRestartSize* {.importc: "minRestartSize".}: cint
    maxBlockSize* {.importc: "maxBlockSize".}: cint
    maxMatvecs* {.importc: "maxMatvecs".}: PRIMME_INT
    maxOuterIterations* {.importc: "maxOuterIterations".}: PRIMME_INT
    intWorkSize* {.importc: "intWorkSize".}: cint
    realWorkSize* {.importc: "realWorkSize".}: csize
    iseed* {.importc: "iseed".}: array[4, PRIMME_INT]
    intWork* {.importc: "intWork".}: ptr cint
    realWork* {.importc: "realWork".}: pointer
    aNorm* {.importc: "aNorm".}: cdouble
    eps* {.importc: "eps".}: cdouble
    printLevel* {.importc: "printLevel".}: cint
    outputFile* {.importc: "outputFile".}: ptr FILE
    matrix* {.importc: "matrix".}: pointer
    preconditioner* {.importc: "preconditioner".}: pointer
    ShiftsForPreconditioner* {.importc: "ShiftsForPreconditioner".}: ptr cdouble
    initBasisMode* {.importc: "initBasisMode".}: primme_init
    ldevecs* {.importc: "ldevecs".}: PRIMME_INT
    ldOPs* {.importc: "ldOPs".}: PRIMME_INT
    projectionParams* {.importc: "projectionParams".}: projection_params
    restartingParams* {.importc: "restartingParams".}: restarting_params
    correctionParams* {.importc: "correctionParams".}: correction_params
    stats* {.importc: "stats".}: primme_stats
    convTestFun* {.importc: "convTestFun".}: proc (eval: ptr cdouble; evec: pointer;
        rNorm: ptr cdouble; isconv: ptr cint; primme: ptr primme_params; ierr: ptr cint) {.noconv.}
    convtest* {.importc: "convtest".}: pointer
    monitorFun* {.importc: "monitorFun".}: proc (basisEvals: pointer;
        basisSize: ptr cint; basisFlags: ptr cint; iblock: ptr cint; blockSize: ptr cint;
        basisNorms: pointer; numConverged: ptr cint; lockedEvals: pointer;
        numLocked: ptr cint; lockedFlags: ptr cint; lockedNorms: pointer;
        inner_its: ptr cint; LSRes: pointer; event: ptr primme_event;
        primme: ptr primme_params; err: ptr cint) {.noconv.}
    monitor* {.importc: "monitor".}: pointer


## ---------------------------------------------------------------------------

type
  primme_preset_method* {.size: sizeof(cint).} = enum
    PRIMME_DEFAULT_METHOD, PRIMME_DYNAMIC, PRIMME_DEFAULT_MIN_TIME,
    PRIMME_DEFAULT_MIN_MATVECS, PRIMME_Arnoldi, PRIMME_GD, PRIMME_GD_plusK,
    PRIMME_GD_Olsen_plusK, PRIMME_JD_Olsen_plusK, PRIMME_RQI, PRIMME_JDQR,
    PRIMME_JDQMR, PRIMME_JDQMR_ETol, PRIMME_STEEPEST_DESCENT,
    PRIMME_LOBPCG_OrthoBasis, PRIMME_LOBPCG_OrthoBasis_Window
  primme_type* {.size: sizeof(cint).} = enum
    primme_int, primme_double, primme_pointer
  primme_params_label* {.size: sizeof(cint).} = enum
    PRIMME_n = 0, PRIMME_matrixMatvec = 1, PRIMME_applyPreconditioner = 2,
    PRIMME_numProcs = 3, PRIMME_procID = 4, PRIMME_commInfo = 5, PRIMME_nLocal = 6,
    PRIMME_globalSumReal = 7, PRIMME_numEvals = 8, PRIMME_target = 9,
    PRIMME_numTargetShifts = 10, PRIMME_targetShifts = 11, PRIMME_locking = 12,
    PRIMME_initSize = 13, PRIMME_numOrthoConst = 14, PRIMME_maxBasisSize = 15,
    PRIMME_minRestartSize = 16, PRIMME_maxBlockSize = 17, PRIMME_maxMatvecs = 18,
    PRIMME_maxOuterIterations = 19, PRIMME_intWorkSize = 20, PRIMME_realWorkSize = 21,
    PRIMME_iseed = 22, PRIMME_intWork = 23, PRIMME_realWork = 24, PRIMME_aNorm = 25,
    PRIMME_eps = 26, PRIMME_printLevel = 27, PRIMME_outputFile = 28, PRIMME_matrix = 29,
    PRIMME_preconditioner = 30, PRIMME_restartingParams_scheme = 31,
    PRIMME_restartingParams_maxPrevRetain = 32,
    PRIMME_correctionParams_precondition = 33,
    PRIMME_correctionParams_robustShifts = 34,
    PRIMME_correctionParams_maxInnerIterations = 35,
    PRIMME_correctionParams_projectors_LeftQ = 36,
    PRIMME_correctionParams_projectors_LeftX = 37,
    PRIMME_correctionParams_projectors_RightQ = 38,
    PRIMME_correctionParams_projectors_RightX = 39,
    PRIMME_correctionParams_projectors_SkewQ = 40,
    PRIMME_correctionParams_projectors_SkewX = 41,
    PRIMME_correctionParams_convTest = 42, PRIMME_correctionParams_relTolBase = 43,
    PRIMME_stats_numOuterIterations = 44, PRIMME_stats_numRestarts = 45,
    PRIMME_stats_numMatvecs = 46, PRIMME_stats_numPreconds = 47,
    PRIMME_stats_elapsedTime = 48, PRIMME_dynamicMethodSwitch = 49,
    PRIMME_massMatrixMatvec = 50, PRIMME_convTestFun = 51, PRIMME_ldevecs = 52,
    PRIMME_ldOPs = 53, PRIMME_monitorFun = 54, PRIMME_monitor = 55,
    PRIMME_initBasisMode = 301, PRIMME_projectionParams_projection = 302,
    PRIMME_stats_numGlobalSum = 471, PRIMME_stats_volumeGlobalSum = 472,
    PRIMME_stats_numOrthoInnerProds = 473, PRIMME_stats_estimateMinEVal = 481,
    PRIMME_stats_estimateMaxEVal = 482, PRIMME_stats_estimateLargestSVal = 483,
    PRIMME_stats_maxConvTol = 484, PRIMME_convtest = 510,
    PRIMME_stats_timeMatvec = 4801, PRIMME_stats_timePrecond = 4802,
    PRIMME_stats_timeOrtho = 4803, PRIMME_stats_timeGlobalSum = 4804




proc sprimme*(evals: ptr cfloat; evecs: ptr cfloat; resNorms: ptr cfloat;
             primme: ptr primme_params): cint {.importc: "sprimme",
    header: "primme.h".}
proc cprimme*(evals: ptr cfloat; evecs: ptr PRIMME_COMPLEX_FLOAT; resNorms: ptr cfloat;
             primme: ptr primme_params): cint {.importc: "cprimme",
    header: "primme.h".}
proc dprimme*(evals: ptr cdouble; evecs: ptr cdouble; resNorms: ptr cdouble;
             primme: ptr primme_params): cint {.importc: "dprimme",
    header: "primme.h".}
proc zprimme*(evals: ptr cdouble; evecs: ptr PRIMME_COMPLEX_DOUBLE;
             resNorms: ptr cdouble; primme: ptr primme_params): cint {.
    importc: "zprimme", header: "primme.h".}
proc primme_initialize*(primme: ptr primme_params) {.importc: "primme_initialize",
    header: "primme.h".}
proc primme_set_method*(`method`: primme_preset_method; params: ptr primme_params): cint {.
    importc: "primme_set_method", header: "primme.h".}
proc primme_display_params*(primme: primme_params) {.
    importc: "primme_display_params", header: "primme.h".}
proc primme_free*(primme: ptr primme_params) {.importc: "primme_free",
    header: "primme.h".}
proc primme_get_member*(primme: ptr primme_params; label: primme_params_label;
                       value: pointer): cint {.importc: "primme_get_member",
    header: "primme.h".}
proc primme_set_member*(primme: ptr primme_params; label: primme_params_label;
                       value: pointer): cint {.importc: "primme_set_member",
    header: "primme.h".}
proc primme_member_info*(label: ptr primme_params_label; label_name: cstringArray;
                        `type`: ptr primme_type; arity: ptr cint): cint {.
    importc: "primme_member_info", header: "primme.h".}
proc primme_constant_info*(label_name: cstring; value: ptr cint): cint {.
    importc: "primme_constant_info", header: "primme.h".}
