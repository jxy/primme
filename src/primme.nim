##  Nim interface for PRIMME.
##
##  The C library can be found at https://github.com/primme/primme

import macros

const PRIMME_INT_SIZE* {.strdefine.} = "64"

macro gen_primme_int:untyped =
  template primme_int_type(t:typed):untyped {.dirty.} =
    type PRIMME_INT* = t
  if PRIMME_INT_SIZE == "64":
    result = getast primme_int_type(int64)
  elif PRIMME_INT_SIZE == "0":
    result = getast primme_int_type(cint)
  elif PRIMME_INT_SIZE == "32":
    result = getast primme_int_type(int32)
  else:
    result = getast primme_int_type(ident(PRIMME_INT_SIZE))

gen_primme_int()

const
  primmeDir {.strdefine.} = "/usr"
  primmeLib {.strdefine.} = primmeDir&"/lib/libprimme.a"
  lapackLib {.strdefine.} = "/usr/lib/libopenblas.a -lm -lgfortran"
{.passC: "-I"&primmeDir&"/include" & " -DPRIMME_INT_SIZE="&PRIMME_INT_SIZE.}
{.passL: primmeLib&" "&lapackLib.}

type PRIMME_COMPLEX_FLOAT*{.importc: "PRIMME_COMPLEX_FLOAT", header: "primme.h".} = object
type PRIMME_COMPLEX_DOUBLE*{.importc: "PRIMME_COMPLEX_DOUBLE", header: "primme.h".} = object

import primme/ccomplex

type
  complex_float = PRIMME_COMPLEX_FLOAT | ccomplex[cfloat] | ccomplex[float32]
  complex_double = PRIMME_COMPLEX_DOUBLE | ccomplex[cdouble] | ccomplex[float]

import primme/primme_eigs, primme/primme_svds

proc run*(primme: var primme_params;
          evals: var openarray[cfloat]; evecs: var openarray[cfloat];
          resNorms: var openarray[cfloat]): int =
  sprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cfloat]; evecs: var seq[complex_float];
          resNorms: var openarray[cfloat]): int =
  cprimme(evals[0].addr, cast[ptr PRIMME_COMPLEX_FLOAT](evecs[0].addr), resNorms[0].addr, primme.addr)
proc run*[I](primme: var primme_params;
          evals: var openarray[cfloat]; evecs: var array[I,complex_float];
          resNorms: var openarray[cfloat]): int =
  cprimme(evals[0].addr, cast[ptr PRIMME_COMPLEX_FLOAT](evecs[0].addr), resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cdouble]; evecs: var openarray[cdouble];
          resNorms: var openarray[cdouble]): int =
  dprimme(evals[0].addr, evecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_params;
          evals: var openarray[cdouble]; evecs: var seq[complex_double];
          resNorms: var openarray[cdouble]): int =
  zprimme(evals[0].addr, cast[ptr PRIMME_COMPLEX_DOUBLE](evecs[0].addr), resNorms[0].addr, primme.addr)
proc run*[I](primme: var primme_params;
          evals: var openarray[cdouble]; evecs: var array[I,complex_double];
          resNorms: var openarray[cdouble]): int =
  zprimme(evals[0].addr, cast[ptr PRIMME_COMPLEX_DOUBLE](evecs[0].addr), resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cfloat]; svecs: var openarray[cfloat];
          resNorms: var openarray[cfloat]): int =
  sprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cfloat]; svecs: var seq[complex_float];
          resNorms: var openarray[cfloat]): int =
  cprimme_svds(svals[0].addr, cast[ptr PRIMME_COMPLEX_FLOAT](svecs[0].addr), resNorms[0].addr, primme.addr)
proc run*[I](primme: var primme_svds_params;
          svals: var openarray[cfloat]; svecs: var array[I,complex_float];
          resNorms: var openarray[cfloat]): int =
  cprimme_svds(svals[0].addr, cast[ptr PRIMME_COMPLEX_FLOAT](svecs[0].addr), resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cdouble]; svecs: var openarray[cdouble];
          resNorms: var openarray[cdouble]): int =
  dprimme_svds(svals[0].addr, svecs[0].addr, resNorms[0].addr, primme.addr)
proc run*(primme: var primme_svds_params;
          svals: var openarray[cdouble]; svecs: var seq[complex_double];
          resNorms: var openarray[cdouble]): int =
  zprimme_svds(svals[0].addr, cast[ptr PRIMME_COMPLEX_DOUBLE](svecs[0].addr), resNorms[0].addr, primme.addr)
proc run*[I](primme: var primme_svds_params;
          svals: var openarray[cdouble]; svecs: var array[I,complex_double];
          resNorms: var openarray[cdouble]): int =
  zprimme_svds(svals[0].addr, cast[ptr PRIMME_COMPLEX_DOUBLE](svecs[0].addr), resNorms[0].addr, primme.addr)
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

template asarray*[T](p:pointer):untyped =
  type A = Uncheckedarray[T]
  cast[ptr A](p)

export primme_eigs, primme_svds
