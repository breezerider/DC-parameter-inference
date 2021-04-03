#!/bin/env python3

import copy
import numpy as np
import scipy as sp

import scipy.optimize as spopt
import scipy.special as spspec

import sys
sys.path.append('/home/sanyi/GIT/LpOptimization/python/fig1/')

import pyPSSAlib


def moments(A,order):
  return pyPSSAlib.moments(A,order)

def bootstrap(result, num_repetitions, bootstrap_size):
  mu_array = np.zeros((num_repetitions, result.shape[1], 2))
  for idx in range(num_repetitions):
    res_idx = np.random.random_integers(0, result.shape[0]-1, size=bootstrap_size)
    mu_array[idx, :, 0] = np.mean(result[res_idx, :], axis=0)
    mu_array[idx, :, 1] = np.std(result[res_idx, :], axis=0)
  return np.mean(mu_array,axis=0)

def min_eps(num_repetitions, num_species, alpha=0.05):
  return spspec.stdtrit(num_repetitions-1, 1.0 - alpha/(2.0*num_species))


class BaseModel:
  x0 = None   # Initial molecule numbers
  k0 = None   # Initial parameter values

  def _get_ki(self, par):
    return copy.deepcopy(self.k0)

  def __init__(self):
    self.k0 = None
    self.engine = pyPSSAlib.pSSAlib()
    #self.engine.test_case = '...'
    self.engine.time_start = 3000.0
    self.engine.time_step = 1.0
    self.engine.time_end = 4000.0

    self.num_repetitions = 100
    self.bootstrap_size = 1000

  def __call__(self, par=None):
    ki = self._get_ki(par)

    self.simulate(ki)

    result = bootstrap(self.engine.result, self.num_repetitions, self.bootstrap_size)
    result = result[:, 0]

    return { 'mean': result, 'ki' : ki }

  def simulate(self,ki):
    #dk = (np.asarray(ki[1:7]) - np.asarray(self.k0[1:7])) / np.asarray(self.k0[1:7]) * 100
    #print(f'simulate: ' + '; '.join(['{:.3e}'.format(x) for x in ki]) + ' -> ' + '; '.join(['{:.3e}%'.format(x) for x in dk]))
    print(f'simulate: ' + '; '.join(['{:.3e}'.format(x) for x in ki]))
    self.engine.params = ki
    if(not self.engine.run()):
      raise(RuntimeError("simulation failed"))

  def get_FanoVar(self, x):
    self.engine.params = x['ki']
    result = spopt.least_squares(lambda a: self.engine.odes(list(a)), x['mean'], bounds=(0,np.inf))

    if(any(result.x < 1.0)):
      print(f'FanoVar = inf; ki = ' + str(x['ki']))
      return np.inf

    Q = self.engine.lyapunovQ(list(result.x))
    J = self.engine.jacobian(list(result.x))
    C = sp.linalg.solve_lyapunov(J, -Q)
    ss= np.divide(np.asarray(x['mean']), np.asarray(result.x))
    return np.sqrt(ss * np.diag(C))

  def distance(self, x, y):
    Fano_Var = get_FanoVar(y)

    t = (np.asarray(x['mean']) - np.asarray(y['mean'])) / Fano_Var
    print(f'distance: {t} -> {np.max(np.abs(t.flatten()))}')
    return np.max(np.abs(t.flatten()))

  def result(self):
    return self.engine.result


class CAModel(BaseModel):
  __name__ = "Colloidal Aggregation Model"

  def _get_ki(self, par):
    ki = copy.deepcopy(self.k0)

    if(not (par is None)):
      ki[1] = par['k_1on']
      ki[2] = par['k_11']
      ki[3] = par['k_11m']
      ki[4] = par['k_1off']
      ki[5] = par['k_2off']
      ki[6] = par['Omega']

    return ki

  def __init__(self):
    BaseModel.__init__(self)
    self.k0 = [2, 2.1, 0.1, 1.0, 0.01, 0.1, 15, 0, 0]
    self.engine.test_case = 'ca'

  def min_eps(self):
    return min_eps(self.num_repetitions, self.k0[0])

class SingleBirthDeathModel(BaseModel):
  __name__ = "Single Birth Death Model"

  def _get_ki(self, par):
    ki = copy.deepcopy(self.k0)

    if(not (par is None)):
      ki[0] = par['k_g']
      ki[1] = par['k_d']

    return ki

  def __init__(self):
    BaseModel.__init__(self)
    self.k0 = [1, 1, 10, 0] # k_g, k_d, Omega
    self.engine.test_case = 'sbd'

  def min_eps(self):
    return min_eps(self.num_repetitions, 1)

class HomoReactionModel(BaseModel):
  __name__ = "Homo Reaction Model"

  def _get_ki(self, par):
    ki = copy.deepcopy(self.k0)

    if(not (par is None)):
      ki[0] = par['k_1']
      ki[1] = par['k_2']

    return ki

  def __init__(self):
    BaseModel.__init__(self)
    self.k0 = [1, 1, 10, 0] # k_1, k_2, Omega
    self.engine.test_case = 'homo'

  def min_eps(self):
    return min_eps(self.num_repetitions, 1)
