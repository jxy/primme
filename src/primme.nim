##  Nim interface for PRIMME.
##
##  The C library can be found at https://github.com/primme/primme

const
  primmeDir {.strdefine.} = "/usr"
  primmeLib {.strdefine.} = primmeDir&"/lib/libprimme.a"
  lapackLib {.strdefine.} = "/usr/lib/libopenblas.a -lm -lgfortran"
{.passC: "-I"&primmeDir&"/include".}
{.passL: primmeLib&" "&lapackLib.}

import primme/complex

type PRIMME_INT* = int64
type PRIMME_COMPLEX_FLOAT*{.importc: "float complex", header: "complex.h".} = complex[cfloat]
type PRIMME_COMPLEX_DOUBLE*{.importc: "double complex", header: "complex.h".} = complex[cdouble]

import primme/primme_eigs, primme/primme_svds

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

export complex
export primme_eigs, primme_svds
