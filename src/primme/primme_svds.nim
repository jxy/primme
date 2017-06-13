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
## *********************************************************************
##  File: primme_svds.h
##
##  Purpose - Main header with the PRIMME SVDS C interface functions.
##
## ****************************************************************************

import primme, primme_eigs

type
  primme_svds_target* {.size: sizeof(cint).} = enum
    primme_svds_largest, primme_svds_smallest, primme_svds_closest_abs
  primme_svds_preset_method* {.size: sizeof(cint).} = enum
    primme_svds_default, primme_svds_hybrid, primme_svds_normalequations, ##  At*A or A*At
    primme_svds_augmented
  primme_svds_operator* {.size: sizeof(cint).} = enum
    primme_svds_op_none, primme_svds_op_AtA, primme_svds_op_AAt,
    primme_svds_op_augmented
  primme_svds_stats* {.importc: "primme_svds_stats", header: "primme.h".} = object
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

  primme_svds_params* {.bycopy, importc: "primme_svds_params", header: "primme.h".} = object
    primme* {.importc: "primme".}: primme_params ## *** Low interface: configuration for the eigensolver
    ##  Keep it as first field to access primme_svds_params from
    ##                              primme_params
    primmeStage2* {.importc: "primmeStage2".}: primme_params ##  other primme_params, used by hybrid
                                                         ##  Specify the size of the rectangular matrix A
    m* {.importc: "m".}: PRIMME_INT ##  number of rows
    n* {.importc: "n".}: PRIMME_INT ##  number of columns
                                ## **** High interface: these values are transferred to primme and primmeStage2 properly
    matrixMatvec* {.importc: "matrixMatvec".}: proc (x: pointer; ldx: ptr PRIMME_INT;
        y: pointer; ldy: ptr PRIMME_INT; blockSize: ptr cint; transpose: ptr cint;
        primme_svds: ptr primme_svds_params; ierr: ptr cint) {.noconv.}
    applyPreconditioner* {.importc: "applyPreconditioner".}: proc (x: pointer;
        ldx: ptr PRIMME_INT; y: pointer; ldy: ptr PRIMME_INT; blockSize: ptr cint;
        transpose: ptr cint; primme_svds: ptr primme_svds_params; ierr: ptr cint) {.noconv.} ##  Input for the following is only required for parallel programs
    numProcs* {.importc: "numProcs".}: cint
    procID* {.importc: "procID".}: cint
    mLocal* {.importc: "mLocal".}: PRIMME_INT
    nLocal* {.importc: "nLocal".}: PRIMME_INT
    commInfo* {.importc: "commInfo".}: pointer
    globalSumReal* {.importc: "globalSumReal".}: proc (sendBuf: pointer;
        recvBuf: pointer; count: ptr cint; primme_svds: ptr primme_svds_params;
        ierr: ptr cint) {.noconv.} ##  Though primme_svds_initialize will assign defaults, most users will set these
    numSvals* {.importc: "numSvals".}: cint
    target* {.importc: "target".}: primme_svds_target
    numTargetShifts* {.importc: "numTargetShifts".}: cint ##  For primme_svds_augmented method, user has to
    targetShifts* {.importc: "targetShifts".}: ptr cdouble ##  make sure  at least one shift must also be set
    `method`* {.importc: "method".}: primme_svds_operator ##  one of primme_svds_AtA, primme_svds_AAt or primme_svds_augmented
    methodStage2* {.importc: "methodStage2".}: primme_svds_operator ##  hybrid second stage method; accepts the same values as method
                                                                ##  These pointers are not for users but for d/zprimme_svds function
    intWorkSize* {.importc: "intWorkSize".}: cint
    realWorkSize* {.importc: "realWorkSize".}: csize
    intWork* {.importc: "intWork".}: ptr cint
    realWork* {.importc: "realWork".}: pointer ##  These pointers may be used for users to provide matrix/preconditioner
    matrix* {.importc: "matrix".}: pointer
    preconditioner* {.importc: "preconditioner".}: pointer ##  The following will be given default values depending on the method
    locking* {.importc: "locking".}: cint
    numOrthoConst* {.importc: "numOrthoConst".}: cint
    aNorm* {.importc: "aNorm".}: cdouble
    eps* {.importc: "eps".}: cdouble
    precondition* {.importc: "precondition".}: cint
    initSize* {.importc: "initSize".}: cint
    maxBasisSize* {.importc: "maxBasisSize".}: cint
    maxBlockSize* {.importc: "maxBlockSize".}: cint
    maxMatvecs* {.importc: "maxMatvecs".}: PRIMME_INT
    iseed* {.importc: "iseed".}: array[4, PRIMME_INT]
    printLevel* {.importc: "printLevel".}: cint
    outputFile* {.importc: "outputFile".}: ptr FILE
    stats* {.importc: "stats".}: primme_svds_stats
    convTestFun* {.importc: "convTestFun".}: proc (sval: ptr cdouble;
        leftsvec: pointer; rightsvec: pointer; rNorm: ptr cdouble; isconv: ptr cint;
        primme: ptr primme_svds_params; ierr: ptr cint) {.noconv.}
    convtest* {.importc: "convtest".}: pointer
    monitorFun* {.importc: "monitorFun".}: proc (basisSvals: pointer;
        basisSize: ptr cint; basisFlags: ptr cint; iblock: ptr cint; blockSize: ptr cint;
        basisNorms: pointer; numConverged: ptr cint; lockedSvals: pointer;
        numLocked: ptr cint; lockedFlags: ptr cint; lockedNorms: pointer;
        inner_its: ptr cint; LSRes: pointer; event: ptr primme_event; stage: ptr cint;
        primme_svds: ptr primme_svds_params; err: ptr cint) {.noconv.}
    monitor* {.importc: "monitor".}: pointer

  primme_svds_params_label* {.size: sizeof(cint).} = enum
    PRIMME_SVDS_primme = 0, PRIMME_SVDS_primmeStage2 = 1, PRIMME_SVDS_m = 2,
    PRIMME_SVDS_n = 3, PRIMME_SVDS_matrixMatvec = 4,
    PRIMME_SVDS_applyPreconditioner = 5, PRIMME_SVDS_numProcs = 6,
    PRIMME_SVDS_procID = 7, PRIMME_SVDS_mLocal = 8, PRIMME_SVDS_nLocal = 9,
    PRIMME_SVDS_commInfo = 10, PRIMME_SVDS_globalSumReal = 11,
    PRIMME_SVDS_numSvals = 12, PRIMME_SVDS_target = 13,
    PRIMME_SVDS_numTargetShifts = 14, PRIMME_SVDS_targetShifts = 15,
    PRIMME_SVDS_method = 16, PRIMME_SVDS_methodStage2 = 17,
    PRIMME_SVDS_intWorkSize = 18, PRIMME_SVDS_realWorkSize = 19,
    PRIMME_SVDS_intWork = 20, PRIMME_SVDS_realWork = 21, PRIMME_SVDS_matrix = 22,
    PRIMME_SVDS_preconditioner = 23, PRIMME_SVDS_locking = 24,
    PRIMME_SVDS_numOrthoConst = 25, PRIMME_SVDS_aNorm = 26, PRIMME_SVDS_eps = 27,
    PRIMME_SVDS_precondition = 28, PRIMME_SVDS_initSize = 29,
    PRIMME_SVDS_maxBasisSize = 30, PRIMME_SVDS_maxBlockSize = 31,
    PRIMME_SVDS_maxMatvecs = 32, PRIMME_SVDS_iseed = 33, PRIMME_SVDS_printLevel = 34,
    PRIMME_SVDS_outputFile = 35, PRIMME_SVDS_stats_numOuterIterations = 36,
    PRIMME_SVDS_stats_numRestarts = 37, PRIMME_SVDS_stats_numMatvecs = 38,
    PRIMME_SVDS_stats_numPreconds = 39, PRIMME_SVDS_stats_elapsedTime = 40,
    PRIMME_SVDS_monitorFun = 41, PRIMME_SVDS_monitor = 42,
    PRIMME_SVDS_stats_numGlobalSum = 391, PRIMME_SVDS_stats_volumeGlobalSum = 392,
    PRIMME_SVDS_stats_numOrthoInnerProds = 393, PRIMME_SVDS_stats_timeMatvec = 401,
    PRIMME_SVDS_stats_timePrecond = 402, PRIMME_SVDS_stats_timeOrtho = 403,
    PRIMME_SVDS_stats_timeGlobalSum = 404, PRIMME_SVDS_convTestFun = 405,
    PRIMME_SVDS_convtest = 406





proc sprimme_svds*(svals: ptr cfloat; svecs: ptr cfloat; resNorms: ptr cfloat;
                  primme_svds: ptr primme_svds_params): cint {.
    importc: "sprimme_svds", header: "primme.h".}
proc cprimme_svds*(svals: ptr cfloat; svecs: ptr PRIMME_COMPLEX_FLOAT;
                  resNorms: ptr cfloat; primme_svds: ptr primme_svds_params): cint {.
    importc: "cprimme_svds", header: "primme.h".}
proc dprimme_svds*(svals: ptr cdouble; svecs: ptr cdouble; resNorms: ptr cdouble;
                  primme_svds: ptr primme_svds_params): cint {.
    importc: "dprimme_svds", header: "primme.h".}
proc zprimme_svds*(svals: ptr cdouble; svecs: ptr PRIMME_COMPLEX_DOUBLE;
                  resNorms: ptr cdouble; primme_svds: ptr primme_svds_params): cint {.
    importc: "zprimme_svds", header: "primme.h".}
proc primme_svds_initialize*(primme_svds: ptr primme_svds_params) {.
    importc: "primme_svds_initialize", header: "primme.h".}
proc primme_svds_set_method*(`method`: primme_svds_preset_method;
                            methodStage1: primme_preset_method;
                            methodStage2: primme_preset_method;
                            primme_svds: ptr primme_svds_params): cint {.
    importc: "primme_svds_set_method", header: "primme.h".}
proc primme_svds_display_params*(primme_svds: primme_svds_params) {.
    importc: "primme_svds_display_params", header: "primme.h".}
proc primme_svds_free*(primme_svds: ptr primme_svds_params) {.
    importc: "primme_svds_free", header: "primme.h".}
proc primme_svds_get_member*(primme_svds: ptr primme_svds_params;
                            label: primme_svds_params_label; value: pointer): cint {.
    importc: "primme_svds_get_member", header: "primme.h".}
proc primme_svds_set_member*(primme_svds: ptr primme_svds_params;
                            label: primme_svds_params_label; value: pointer): cint {.
    importc: "primme_svds_set_member", header: "primme.h".}
proc primme_svds_member_info*(label: ptr primme_svds_params_label;
                             label_name: cstringArray; `type`: ptr primme_type;
                             arity: ptr cint): cint {.
    importc: "primme_svds_member_info", header: "primme.h".}
proc primme_svds_constant_info*(label_name: cstring; value: ptr cint): cint {.
    importc: "primme_svds_constant_info", header: "primme.h".}
