pro call_func, p, args=args, logL=logL
  x = p
  ivar = args
  logL = -0.5 * total(ivar * x^2)
end
