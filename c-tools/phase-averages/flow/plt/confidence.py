#!/usr/bin/env python2
import numpy as np
from scipy import stats

def fit_param(x, y):
  # From "Simple linear regression -- numerical example" wiki page
  n = np.size(x)
  Sx = np.sum(x)
  Sy = np.sum(y)
  Sxx = np.sum(x**2)
  Sxy = np.sum(x*y)
  Syy = np.sum(x**2)

  # y = alpha + beta x
  beta = (n*Sxy - Sx*Sy)/(n*Sxx - Sx*Sx)
  alpha = Sy/n - beta*Sx/n

  s_eps2 = np.sum((y - alpha - beta*x)**2)
  s_beta = np.sqrt((s_eps2/(n-2))/(Sxx - Sx*Sx/n))
  s_alpha = s_beta*np.sqrt(Sxx/n)

  # assume confidence is 0.05
  confidence = 0.05
  t = stats.t.ppf(1.-0.5*confidence, n)

  alpha_err = t*s_alpha
  beta_err = t*s_beta
  
  alphaExp = np.exp(alpha)
  alphaExp1 = np.exp(alpha)*np.exp(alpha_err)
  alphaExp2 = np.exp(alpha)*np.exp(-alpha_err)
  meanAlphaError = 0.5*np.abs(alphaExp1 - alphaExp2)

  # return slope, intercept, errorSlope, errorIntercept
  return (beta, alphaExp, beta_err, meanAlphaError)
