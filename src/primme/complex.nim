{.passl:"-lm".}
type
  ccomplex*[T] = object
    a: array[2,T]
  ccomplex_float{.importc:"complex float", header:"complex.h", nodecl.} = object
  ccomplex_double{.importc:"complex double", header:"complex.h", nodecl.} = object
  ccomplex_types = ccomplex_float or ccomplex_double

func re*[T](x:ccomplex[T]):T = x.a[0]
func im*[T](x:ccomplex[T]):T = x.a[1]

proc `re=`*[T](x:var ccomplex[T], v:T) = x.a[0] = v
proc `im=`*[T](x:var ccomplex[T], v:T) = x.a[1] = v

func creal(x:ccomplex_double):cdouble {.importc.}
func cimag(x:ccomplex_double):cdouble {.importc.}
func crealf(x:ccomplex_float):cfloat {.importc.}
func cimagf(x:ccomplex_float):cfloat {.importc.}

func CMPLX(re, im:float):ccomplex_double {.importc, nodecl.}
func CMPLXF(re, im:float32):ccomplex_float {.importc, nodecl.}

func cmplx*[T](re, im:T):ccomplex[T] {.inline, noinit.} =
  result.a[0] = re
  result.a[1] = im

func toccomplex_type(x:ccomplex[float]):ccomplex_double {.inline, noinit.} =
  CMPLX(x.re, x.im)
func toccomplex_type(x:ccomplex[float32]):ccomplex_float {.inline, noinit.} =
  CMPLXF(x.re, x.im)

func toccomplex(x:ccomplex_double):ccomplex[float] {.inline, noinit.} =
  cmplx(creal(x), cimag(x))
func toccomplex(x:ccomplex_float):ccomplex[float32] {.inline, noinit.} =
  cmplx(crealf(x), cimagf(x))

func `+`(x:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = `x`;".}
func `-`(x:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = -`x`;".}

func `+`(x, y:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = `x` + `y`;".}
func `-`(x, y:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = `x` - `y`;".}
func `*`(x, y:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = `x` * `y`;".}
func `/`(x, y:ccomplex_types):ccomplex_types {.inline, noinit.} =
  {.emit: "`result` = `x` / `y`;".}

proc `+=`(x:var ccomplex_types, y:ccomplex_types) {.inline.} =
  {.emit: "*`x` += `y`;".}
proc `-=`(x:var ccomplex_types, y:ccomplex_types) {.inline.} =
  {.emit: "*`x` -= `y`;".}
proc `*=`(x:var ccomplex_types, y:ccomplex_types) {.inline.} =
  {.emit: "*`x` *= `y`;".}
proc `/=`(x:var ccomplex_types, y:ccomplex_types) {.inline.} =
  {.emit: "*`x` /= `y`;".}

func `+`*[T](x:ccomplex[T]):ccomplex[T] {.inline, noinit.} = x
func `-`*[T](x:ccomplex[T]):ccomplex[T] {.inline, noinit.} =
  result = toccomplex(-toccomplex_type(x))

func `+`*[T](x, y:ccomplex[T]):ccomplex[T] {.inline, noinit.} =
  result = toccomplex(toccomplex_type(x) + toccomplex_type(y))
func `-`*[T](x, y:ccomplex[T]):ccomplex[T] {.inline, noinit.} =
  result = toccomplex(toccomplex_type(x) - toccomplex_type(y))
func `*`*[T](x, y:ccomplex[T]):ccomplex[T] {.inline, noinit.} =
  result = toccomplex(toccomplex_type(x) * toccomplex_type(y))
func `/`*[T](x, y:ccomplex[T]):ccomplex[T] {.inline, noinit.} =
  result = toccomplex(toccomplex_type(x) / toccomplex_type(y))

proc `+=`*[T](x:var ccomplex[T], y:ccomplex[T]) {.inline.} =
  var z = toccomplex_type(x)
  z += toccomplex_type(y)
  x = toccomplex(z)
proc `-=`*[T](x:var ccomplex[T], y:ccomplex[T]) {.inline.} =
  var z = toccomplex_type(x)
  z -= toccomplex_type(y)
  x = toccomplex(z)
proc `*=`*[T](x:var ccomplex[T], y:ccomplex[T]) {.inline.} =
  var z = toccomplex_type(x)
  z *= toccomplex_type(y)
  x = toccomplex(z)
proc `/=`*[T](x:var ccomplex[T], y:ccomplex[T]) {.inline.} =
  var z = toccomplex_type(x)
  z /= toccomplex_type(y)
  x = toccomplex(z)

converter conv_toccomplex*[T](x:T):ccomplex[T] {.inline, noinit.} = cmplx(x, 0)

func `$`[T](x:ccomplex[T]):string =
  let s = if x.im > 0: "+" else: ""
  result = $x.re & s & $x.im & "i"

when isMainModule:
  var x: ccomplex[float]
  x.re = -3
  x.im = 4
  var y = cmplx(-0.4, 0.3)
  # ((x+y)/y - y)*y = x+y - y*y
  x += y
  x /= y
  x -= y
  x *= y
  echo "x: ",x
  echo "y: ",y
  # ((x+y-y*y)/y+y)*y - y = x
  let z = (x / y + y) * y - y
  echo "z: ",z
