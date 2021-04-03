#!/usr/bin/env python3

import os
import sys

import copy
#import matplotlib.pyplot as plt

from models import moments
from models import bootstrap
from models import min_eps

from models import CAModel
from models import HomoReactionModel
from models import SingleBirthDeathModel

import numpy as np
import numpy.linalg as npalg

import multiprocessing

def costfcn_Mueller2012c(x, y):
  x_den = copy.deepcopy(x)
  x_den[x_den == 0] = 1.0

  return npalg.norm( (x - y)/x_den, ord=2, axis=1 )

def costfcn_Mueller2012c(x, y):
  x_den = copy.deepcopy(x)
  x_den[x_den == 0] = 1.0

  return npalg.norm( (x - y)/x_den, ord=2, axis=1 )

def start_process(data,MU,BS,RES,n,data_loaded):
  var_dict['data'] = data
  var_dict['MU']   = MU
  var_dict['BS']   = BS
  var_dict['RES']  = RES
  var_dict['n']    = n
  var_dict['data_loaded'] = data_loaded

  print('Starting ', multiprocessing.current_process().name)

# simulation worker process
def worker_simulate(args):
  global num_test_repetitions

  model = args[0]
  ki    = args[1]
  idx   = args[2]

  data_raw, data_shape = var_dict['data']
  #MU_raw, MU_shape     = var_dict['MU']
  #BS_raw, BS_shape     = var_dict['BS']
  n                    = var_dict['n']
  data_loaded          = var_dict['data_loaded']

  data_np = np.frombuffer(data_raw, dtype=np.float).reshape(data_shape)

  if(data_loaded):
    return
  else:
    for i in range(num_test_repetitions):
      print(f'simulation iteration {idx}, subiteration {i}...')
      model.simulate(ki)

      data_np[:, (idx + i*n + 1):(idx + i*n + n + 1)] = model.result()[:, :]

# compute worker process
def worker_compute(args):
  global bootstrap_repetitions, bootstrap_size, num_moments

  idx   = args[0]

  data_raw, data_shape = var_dict['data']
  MU_raw, MU_shape     = var_dict['MU']
  BS_raw, BS_shape     = var_dict['BS']
  n                    = var_dict['n']
  data_loaded          = var_dict['data_loaded']

  data_np = np.frombuffer(data_raw, dtype=np.float).reshape(data_shape)
  MU_np = np.frombuffer(MU_raw, dtype=np.float).reshape(MU_shape)
  BS_np = np.frombuffer(BS_raw, dtype=np.float).reshape(BS_shape)

  if(not data_loaded):
    raise RuntimeError('data must be pre-loaded')

  dt = np.zeros((data_shape[0],n),dtype=np.float)

  mns = list(range(1,num_moments+1))
  for i in range(num_test_repetitions):
    print(f'compute iteration {idx}, subiteration {i}...')

    dt[:] = data_np[:, (idx + i*n + 1):(idx + i*n + n + 1)]

    MU_np[idx, i, :, :] = moments(dt, mns).T
    BS_np[idx, i, :, :] = bootstrap(dt, bootstrap_repetitions, bootstrap_size)

def worker_compare(args):
  global bootstrap_repetitions, num_test_repetitions, num_param_points, alpha
  
  model = args[0]
  ki    = args[1]
  idx   = args[2]

  #data_raw, data_shape = var_dict['data']
  MU_raw, MU_shape     = var_dict['MU']
  BS_raw, BS_shape     = var_dict['BS']
  RES_raw, RES_shape   = var_dict['RES']
  n                    = var_dict['n']
  data_loaded          = var_dict['data_loaded']
  
  if(not data_loaded):
    raise RuntimeError('data must be pre-loaded')
  
  MU_np = np.frombuffer(MU_raw, dtype=np.float).reshape(MU_shape)
  BS_np = np.frombuffer(BS_raw, dtype=np.float).reshape(BS_shape)
  RES_np = np.frombuffer(RES_raw, dtype=np.uint).reshape(RES_shape)
  
  eps = min_eps(bootstrap_repetitions, n, alpha)
  bFano = True
    
  # true-positive comparison loop
  for i in range(num_test_repetitions):

    print(f'true-positive iteration {idx}, subiteration {i}...')

    var = model.get_FanoVar({'mean': BS_np[idx, i, :, 0], 'ki' : ki})
    if(np.isinf(var)):
      bFano = False
      #var = BS_np[idx, i, :, 1]
      #print(f'std var = {var}')

    for j in range(num_test_repetitions):
      
      if(i == j):
        continue

      c = costfcn_Mueller2012c(MU_np[idx, i, :, :], MU_np[idx, j, :, :])

      if(c <= alpha):
        RES_np[idx, 0] += 1
      
      if(not bFano):
        var = np.sqrt((BS_np[idx, i, :, 1]**2 + BS_np[idx, j, :, 1]**2)/num_test_repetitions)
        print(f'std var = {var}')

      t = np.max(np.abs((BS_np[idx, i, :, 0] - BS_np[idx, j, :, 0]) / var))

      if(t <= eps):
        RES_np[idx, 1] += 1

      #print(f'{i} vs {j}: cost = {c} >= {alpha}; oracle = {t} >= {eps};')
        
  # true-negative comparison loop
  i = np.random.randint(0, num_test_repetitions)
  
  print(f'true-negative iteration {idx}, subiteration {i}...')

  var = model.get_FanoVar({'mean': BS_np[idx, i, :, 0], 'ki' : ki})
  if(np.isinf(var)):
    bFano = False
    #var = BS_np[idx, i, :, 1]
    #print(f'std var = {var}')

  for idx2 in range(num_param_points*num_param_points):

    c = costfcn_Mueller2012c(MU_np[idx, i, :, :], MU_np[idx2, i, :, :])

    if(c > alpha):
      RES_np[idx, 2] += 1
      
    if(not bFano):
      var = np.sqrt((BS_np[idx, i, :, 1]**2 + BS_np[idx2, i, :, 1]**2)/num_test_repetitions)
      print(f'std var = {var}')

    t = np.max(np.abs((BS_np[idx, i, :, 0] - BS_np[idx2, i, :, 0]) / var))

    if(t > eps):
      RES_np[idx, 3] += 1

    #print(f'{i} vs {j}: cost = {c} >= {alpha}; oracle = {t} >= {eps};')


# global definitions
num_param_points = 10
num_test_repetitions = 10

bootstrap_repetitions = 1000
bootstrap_size = 1000

num_moments = 4
alpha = 0.15

# global variable dictionary
var_dict = {}

if __name__ == '__main__':
  pool_size       = multiprocessing.cpu_count() * 2
  
  #if((num_param_points % 2 == 0)):
  K1 = np.concatenate((np.logspace(-2, 0, num=num_param_points//2 - 1, endpoint=False), np.logspace(0, 2, num=num_param_points//2 + 1)))
  K2 = copy.deepcopy(K1)
  #else:
  #  K1 = np.concatenate((np.logspace(-2, 0, num=num_param_points//2, endpoint=False), np.logspace(0, 2, num=num_param_points//2 + 1)))
  #  K2 = copy.deepcopy(K1)

  KX, KY = np.meshgrid(K1, K2, sparse=False, indexing='ij')
  K = np.column_stack((np.ravel(KX), np.ravel(KY)))
  
  #print('SingleBirthDeathModel')
  #print((10.0*K[:,0]/K[:,1]) < 1.0)
  
  #print('HomoReactionModel')
  #print((10.0*np.sqrt(K[:,0]/(2.0*K[:,1]))) < 1.0)
  
  #quit()

  fname = f'K1_K2.csv'.replace(' ', '-')
  np.savetxt(fname, K, delimiter=',')

  # model loop
  for model, n in [(HomoReactionModel(), 1), (SingleBirthDeathModel(), 1)]:
    
    model.k0 = [1e-1, (1/5)**2 * 1e-1/2.0, 10.0, 0]
    
    data_shape = (model.engine.num_timepoints, n*num_test_repetitions*num_param_points*num_param_points + 1)
    data_raw = multiprocessing.RawArray('d', data_shape[0] * data_shape[1])
    data_np = np.frombuffer(data_raw).reshape(data_shape)

    data_loaded = None
    fname = f'data_{model.__name__}.csv'.replace(' ', '-')
    if(os.path.isfile(fname)):
      try:
        data_loaded = np.loadtxt(fname, delimiter=',')
        data_np[:] = data_loaded[:]
        data_loaded = True
        print(f'loaded data from {fname}...')
      except:
        data_loaded = False

    BS_shape = (num_param_points*num_param_points, num_test_repetitions, n, 2)
    BS_raw = multiprocessing.RawArray('d', int(np.prod(BS_shape)))
    BS_np = np.frombuffer(BS_raw, dtype=np.float).reshape(BS_shape)

    MU_shape = (num_param_points*num_param_points, num_test_repetitions, n, num_moments)
    MU_raw = multiprocessing.RawArray('d', int(np.prod(MU_shape)))
    MU_np = np.frombuffer(MU_raw, dtype=np.float).reshape(MU_shape)

    # results by column: true-positives (costfcn, oracle); true-negatives (costfcn, oracle)
    RES_shape = (num_param_points*num_param_points, 4)
    RES_raw = multiprocessing.RawArray('L', int(np.prod(RES_shape)))
    RES_np = np.frombuffer(RES_raw, dtype=np.uint).reshape(RES_shape)

    # simulate data
    if(not data_loaded):
      data_np[:, 0] = np.linspace(model.engine.time_start, model.engine.time_end, num=data_shape[0], endpoint=True)
      with multiprocessing.Pool(processes=pool_size,initializer=start_process, initargs=( (data_raw, data_shape), (MU_raw, MU_shape), (BS_raw, BS_shape), (RES_raw, RES_shape), n, data_loaded) ) as pool:
        # simulation arguments loop
        args_sim = list()
        for idx in range(K.shape[0]):
          #print(f'simulation iteration {idx}...')

          args_sim.append( [copy.deepcopy(model), [K[idx,0]*model.k0[0],K[idx,1]*model.k0[1],10,0], idx] ) 

        pool.map(worker_simulate, args_sim)

      fname = f'data_{model.__name__}.csv'.replace(' ', '-')
      np.savetxt(fname, data_np, delimiter=',')
      data_loaded = True

    # precompute comparison
    with multiprocessing.Pool(processes=pool_size,initializer=start_process, initargs=( (data_raw, data_shape), (MU_raw, MU_shape), (BS_raw, BS_shape), (RES_raw, RES_shape), n, data_loaded) ) as pool:
      # compute arguments loop
      args_comp = list()
      for idx in range(K.shape[0]):
        #print(f'compute iteration {idx}...')

        args_comp.append( [idx] ) 

      pool.map(worker_compute, args_comp)

    # compare
    with multiprocessing.Pool(processes=pool_size,initializer=start_process, initargs=( (data_raw, data_shape), (MU_raw, MU_shape), (BS_raw, BS_shape), (RES_raw, RES_shape), n, data_loaded) ) as pool:
      # compare arguments loop
      args_comp = list()
      for idx in range(K.shape[0]):
        #print(f'compute iteration {idx}...')

        args_comp.append( [copy.deepcopy(model), [K[idx,0]*model.k0[0],K[idx,1]*model.k0[1],10,0], idx] ) 

      pool.map(worker_compare, args_comp)
      
    # save results
    fname = f'comparison_{model.__name__}.csv'.replace(' ', '-')
    np.savetxt(fname, RES_np, delimiter=',')
