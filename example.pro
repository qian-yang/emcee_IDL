pro example
  ndim = 10
  nwalkers = 100
  ivar = 1./randomu(seed, ndim, nwalkers)
  p0 = randomu(seed, ndim, nwalkers)
  emcee, p0, kargs=ivar, p=p, lnprob=lnprob
  print, p, lnprob
end
