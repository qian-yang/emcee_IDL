;+
; NAME:
;   emcee
;
; PURPOSE:
;   run MCMC (emcee method, learn from the emcee Python package)
;
; REFERENCE:
;   Foreman-Mackey et al. 2013
;
; INPUT:
;   p0, data, iterations, lp0
;   p0: initial positions of the walkers (ndim, nwalkers)
;
; OPTIONAL INPUTS:
;   lp0: the log posterior probabilities for the walkers at positions given by p0 (optional, (ndim, nwalkers))
;   data: data or data used to calculate the probabilities
;   iterations: optional, default to 1
;   kargs: kargs needs in call_func.pro
;
; KEYWORD: save_chains (output chains, allp)
;   Keywords for the photo-RM function
;   log: the parameters are in log space (alog10)
;   cont: if cont is set, there are only two parameters (sigma and tau)
;   scale: if scale is set, sigma = sigma_up * sqrt(tau/2.0)
;   lprior: if lprior is set, using lprior probability
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;  p: final positions of the walkers (ndim, nwalkers)
;  lnprob: the log posterior probabilities for the walkers at positions given by p
;  naccepted: number accepted for walkers (ndim, nwalkers)
;  chains: all chains when save_chains is set (ndim, nwalkers * iterations)
;  allp: the log posterior probabilities for the walkers at positions given by chains
;
; COMMENTS:
;   Need a file called call_func.pro to calculate the probabilities
;
; REVISION HISTORY:
;   28-Jul-2017  Written by Qian Yang, qianyang.astro@gmail.com
;
; Tudo:
;   get_autocorr_time
;-

function propose_stretch, p0, p1, lp0, data=data, a=a, newlnprob=newlnprob, accept=accept, nacc=nacc, kargs=kargs
  ; Propose a new position for one sub-ensemble given the positions of another. (see _propose_stretch in Python package.)
  ; input: p0, p1, lp0
  ; p0: The positions from which to jump. (ndim, Ns)
  ; p1: The positions of the other ensemble. (ndim, Nc)
  ; lp0: The log-probabilities at ``p0``. (Ns)
  ; return: q, newlnprob, accept
  ; q: The new proposed positions for the walkers in ``ensemble``. (ndim, Nc)
  ; newlnprob: The vector of log-probabilities at the positions given by ``q``. (Nc)
  ; `accept`` - A vector of 1 or 0 indicating whether or not the proposed position for each walker should be accepted. (Nc)
  ; The probability callculation function is in call_func.pro
  if ~keyword_set(a) then a=2.0
  Ns = n_elements(p0[0, *])
  Nc = n_elements(p1[0, *])
  ndim = n_elements(p0[*, 0])
  ; Generate the vectors of random numbers that will produce the proposal.
  zz = ((a - 1.0) * randomu(seed, Ns) + 1)^2/a ;between 1/a and a
  rint = floor(randomu(seed, Nc) * Ns)
  zzt = transpose(rebin(zz, Ns, ndim))
  q = p1[*, rint] + zzt * (p0 - p1[*, rint]) ; Equation 7 in FM2013
  call_func, q, logL = newlnprob, data=data, kargs=kargs
  ; decide whether or not accept the proposals.
  lnpdiff = (ndim - 1.0) * alog(zz) + newlnprob - lp0
  r01 = alog(randomu(seed, Ns))
  accept = where(lnpdiff gt r01, nacc)
  return, q
end

function update_step, p, lnprob, naccepted, S0, S1, p1=p1, lp1=lp1, nc1=nc1, data=data, kargs=kargs
  q = propose_stretch(p[*, S0], p[*, S0], lnprob[S0], data=data, newlnprob=newlnprob, accept=accept, nacc=nacc, kargs=kargs)
  if (nacc gt 0) then begin
    ; update the accepted walkers
    lnprob[S0[accept]] = newlnprob[accept]
    p[*, S0[accept]] = q[*, accept]
    naccepted[accept] = naccepted[accept] + 1
  endif
  p1 = p
  lp1 = lnprob
  nc1 = naccepted
  return, p1
end

pro emcee, p0, iterations=iterations, data=data, lp0=lp0, p=p, lnprob=lnprob, save_chains=save_chains, chains=chains, allp=allp, kargs=kargs
  if ~keyword_set(iterations) then iterations = 1
  nwalkers = n_elements(p0[0, *])
  ndim = n_elements(p0[*, 0])
  naccepted = lonarr(ndim, nwalkers)
  p = p0
  if keyword_set(save_chains) then begin
    nall = nwalkers * iterations
    chains = fltarr(ndim, nall)
    allp = fltarr(nall)
  endif
  halfk = long(nwalkers/2)
  first = lindgen(halfk)
  second = lindgen(halfk) + halfk
  ; calculate lp0 if it's not provided
  if ~keyword_set(lp0) then begin
    call_func, p, logL = lp0, data=data, kargs=kargs
  endif
  lnprob = lp0
  for i = 0, iterations - 1 do begin
    S0 = first
    S1 = second
    p = update_step(p, lnprob, naccepted, S0, S1, data=data, p1=p1, lp1=lp1, nc1=nc1, kargs=kargs)
    p = p1
    lnprob = lp1
    naccepted = nc1
    ; inverse
    S0 = second
    S1 = first
    p = update_step(p, lnprob, naccepted, S0, S1, data=data, p1=p1, lp1=lp1, nc1=nc1, kargs=kargs)
    p = p1
    lnprob = lp1
    naccepted = nc1
    ; output
    if keyword_set(save_chains) then begin
      chains[*, i*nwalkers : (i+1)*nwalkers-1] = p
      allp[i*nwalkers : (i+1)*nwalkers-1] = lnprob
    endif
  endfor
end
