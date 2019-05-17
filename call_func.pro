pro call_func, p, kargs=kargs, logL=logL, data=data
  x = p
  ivar = kargs
  logL = -0.5 * total(ivar * x^2)
end
