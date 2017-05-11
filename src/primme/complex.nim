type complex*[T] = tuple[re,im:T]
converter toComplex*[T](x:T):complex[T] {.inline.} =
  result.re = x
  result.im = 0
proc `-`*[T](x:complex[T]):complex[T] {.inline.} =
  result.re = -x.re
  result.im = -x.im
proc `+`*[T](x:T, y:complex[T]):complex[T] {.inline.} =
  result.re = x+y.re
  result.im = y.im
proc `+`*[T](x:complex[T], y:T):complex[T] {.inline.} =
  result.re = x.re+y
  result.im = x.im
proc `+`*[T](x:complex[T], y:complex[T]):complex[T] {.inline.} =
  result.re = x.re+y.re
  result.im = x.im+y.im
proc `-`*[T](x:T, y:complex[T]):complex[T] {.inline.} =
  result.re = x-y.re
  result.im = -y.im
proc `-`*[T](x:complex[T], y:T):complex[T] {.inline.} =
  result.re = x.re-y
  result.im = x.im
proc `-`*[T](x:complex[T], y:complex[T]):complex[T] {.inline.} =
  result.re = x.re-y.re
  result.im = x.im-y.im
proc `*`*[T](x:T, y:complex[T]):complex[T] {.inline.} =
  result.re = x*y.re
  result.im = x*y.im
proc `*`*[T](x:complex[T], y:T):complex[T] {.inline.} =
  result.re = x.re*y
  result.im = x.im*y
proc `*`*[T](x:complex[T], y:complex[T]):complex[T] {.inline.} =
  result.re = x.re*y.re - x.im*y.im
  result.im = x.im*y.re + x.re*y.im
proc `/`*[T](x:T, y:complex[T]):complex[T] {.inline.} =
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
proc `/`*[T](x:complex[T], y:T):complex[T] {.inline.} =
  result.re = x.re/y
  result.im = x.im/y
proc `/`*[T](x:complex[T], y:complex[T]):complex[T] {.inline.} =
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
proc `+=`*[T](x:var complex[T], y:T) {.inline.} =
  x.re += y
proc `+=`*[T](x:var complex[T], y:complex[T]) {.inline.} =
  x.re += y.re
  x.im += y.im
proc `-=`*[T](x:var complex[T], y:T) {.inline.} =
  x.re -= y
proc `-=`*[T](x:var complex[T], y:complex[T]) {.inline.} =
  x.re -= y.re
  x.im -= y.im
proc `*=`*[T](x:var complex[T], y:T) {.inline.} =
  x.re *= y
  x.im *= y
proc `*=`*[T](x:var complex[T], y:complex[T]) {.inline.} =
  let im = x.re*y.im + x.im*y.re
  x.re = x.re*y.re - x.im*y.im
  x.im = im
proc `/=`*[T](x:var complex[T], y:T) {.inline.} =
  x.re /= y
  x.im /= y
proc `/=`*[T](x:var complex[T], y:complex[T]) {.inline.} =
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
