import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf

import numpy as np
import time
from .dictClass import dictClass
import pickle


def init_population_latinhypercube(num_population, len_of_vec, bounds, rng=None):
    """
    Initializes the population with Latin Hypercube Sampling.
    Latin Hypercube Sampling ensures that each parameter is uniformly
    sampled over its range.
    """
    rng = np.random
    bounds_np = np.array(bounds)
        
    segsize = 1.0 / num_population
    samples = (segsize * rng.uniform(size=(num_population,len_of_vec))
               + np.linspace(0., 1., num_population,
                             endpoint=False)[:, np.newaxis])

    population = np.zeros_like(samples)
    
    # Initialize population of candidate solutions by permutation
    for j in range(len_of_vec):
        order = rng.permutation(range(num_population))
        population[:, j] = samples[order, j]

    return population*(bounds_np[:,1]-bounds_np[:,0])[None,:]+bounds_np[:,0][None,:]



def minimize_GP_LCB(func, x0, bounds, args=(), initial_boundary_ratio = 0.8,
                    n_init_pop = None, maxnfev=10000, patience=1000,
                    LCB_nsig=1, trainGP_everyEpoch=False,
                    acquisition_optimize_options = {"num_restarts":2, "raw_samples":20},
                    scipy_minimize_options = {"ftol":2e-3},
                    returnGP=False, 
                    verbose=1, 
                    set_func_to_best = True,
                    prev_result=None ):
    '''
    func: function to obtimize
    x0: (numpy vector) initial function argument
    bounds: limit of function arguments so that [ bounds[i][0] <= x[i] <= bounds[i][1] for i in range(len(x)) ]
    maxnfev: (int) maximum number of function evaluation
    LCB_nsig: lower confidence bound in unit of stantard deviation of GP prediction 
    patience: quit optimization if the function value is not improving for patience number of function evaluation
    '''
    if not prev_result is None:
        return _minimize_GP_LCB_continue(prev_result=prev_result,
                                         func=func, bounds=bounds, args=args, 
                                         maxnfev=maxnfev, patience=patience, 
                                         LCB_nsig=LCB_nsig,
                                         acquisition_optimize_options = acquisition_optimize_options,
                                         scipy_minimize_options=scipy_minimize_options,
                                         returnGP = returnGP,   
                                         verbose = verbose, 
                                         set_func_to_best=set_func_to_best )
    ndim = len(x0)
    if n_init_pop == None:
        n_init_pop = min(2*ndim,maxnfev)
    
    bounds_np = np.array(bounds)
    bounds_torch = torch.tensor(np.transpose(bounds_np))

    x = np.zeros([n_init_pop,ndim])
    x[0 ,:] = x0[:]
    if n_init_pop>1:
        x[1:,:] = init_population_latinhypercube(n_init_pop-1, ndim, bounds_np*initial_boundary_ratio)
    y = np.array([-func(x_,*args) for x_ in x])
    x = torch.tensor(x)
    y = torch.tensor(y[:,None])    

    result = dictClass()
    result.history = dictClass()    
    result.history.x = x.detach().numpy()
    result.history.y = -y.detach().numpy().flatten()
    
    imin = result.history.y.argmin()    
    result.x = result.history.x[imin]
    y_min = result.history.y[imin]
    result.history.y_min = [y_min]
    result.surrogate = dictClass()
    result.surrogate.training_time = []
    result.surrogate.acqusition_optimization_time = []

    gp = SingleTaskGP(x,y)
    i_patience = 0 
    data_size = 0
    for i_eval in range(n_init_pop,maxnfev):
        if verbose >= 1:
            print("epoch: {0:3d}/{1:3d}, min_fun :{2:.3e}".format(i_eval, maxnfev, y_min))
            print("decision x: ", result.x)
        if i_patience >= patience and patience>0:       
            result.message = "maximum patience reached. optimization halt..."
            break
        
        t0 = time.time()
        # train GP hyper-parameters only when number of training data is increased by more than 10%
        if trainGP_everyEpoch:
            gp = SingleTaskGP(x,y)
        else:
            if data_size*1.1 < len(x):
                data_size = len(x)
                mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
                fit_gpytorch_model(mll,options=scipy_minimize_options)
        t1 = time.time()

        # increase bounds over optimization epochs starting from reduced bounds by factor of "initial_boundary_ratio"
        reduction = initial_boundary_ratio + (1-initial_boundary_ratio)*(i_eval+1)/maxnfev 
        bounds_reduced = bounds_torch*reduction
        UCB = UpperConfidenceBound(gp, beta=LCB_nsig**0.5)
        x1, y1acqu = optimize_acqf(UCB, 
                                   bounds=bounds_reduced, 
                                   q=1, num_restarts=2, raw_samples=20,
                                   options=scipy_minimize_options)
        t2 = time.time()
        if verbose==2:
            print("surrogate model training time (sec):",t1-t0,", acqusition optimization time (sec):",t2-t1)

        y1 = torch.tensor([-func(x1.detach().numpy()[0,:],*args)])[:,None]
        x = torch.concat((x,x1),axis=0)
        y = torch.concat((y,y1),axis=0)
        # add GP training data
        if not trainGP_everyEpoch:
            gp = gp.get_fantasy_model(x1,y1)

        result.surrogate.training_time.append(t1-t0)
        result.surrogate.acqusition_optimization_time.append(t2-t1)
        result.history.x = x.detach().numpy()
        result.history.y = -y.detach().numpy().flatten()
        imin = result.history.y.argmin()    
        result.x = result.history.x[imin]
        new_y_min = result.history.y[imin]
        result.history.y_min.append(new_y_min)
        if new_y_min <= y_min:
            i_patience = 0
            y_min = new_y_min
        else:
            i_patience += 1
        
    imin = result.history.y.argmin()    
    result.x = result.history.x[imin]
    result.fun = result.history.y[imin]
    result.nfev = i_eval+1

    if result.nfev == maxnfev:
        result.message = "maximum number of function evaluation reached"
    if returnGP:
        result.GP = gp
    if set_func_to_best:
        func(result.x)
    
    return result




def _minimize_GP_LCB_continue(prev_result, func, bounds, args=(), 
                              maxnfev=10000, patience=1000,
                              LCB_nsig=1,
                              acquisition_optimize_options = {"num_restarts":2, "raw_samples":20},
                              scipy_minimize_options = {"ftol":2e-3},
                              returnGP = True, verbose=1, 
                              set_func_to_best = True):
    '''
    func: function to obtimize
    x0: (numpy vector) initial function argument
    bounds: limit of function arguments so that [ bounds[i][0] <= x[i] <= bounds[i][1] for i in range(len(x)) ]
    maxnfev: (int) maximum number of function evaluation
    LCB_nsig: lower confidence bound in unit of stantard deviation of GP prediction 
    patience: quit optimization if the function value is not improving for patience number of function evaluation
    '''
    result = prev_result
    x = torch.tensor(result.history.x)
    y = torch.tensor(-result.history.y)[:,None]
    data_size, ndim = x.shape

    bounds_np = np.array(bounds)
    bounds_torch = torch.tensor(np.transpose(bounds_np))

    if hasattr(result,"GP"):
        gp = result.GP
    else:
        gp = SingleTaskGP(x,y)
    i_patience = 0 
    data_size = 0
    y_min = result.history.y_min[-1]
    for i_eval in range(maxnfev):
        if verbose >= 1:
            print("epoch: {0:3d}/{1:3d}, min_fun :{2:.3e}".format(i_eval, maxnfev, y_min))
            print("decision x: ", result.x)
        if i_patience >= patience and patience>0:       
            result.message = "maximum patience reached. optimization halt..."
            break
        
        t0 = time.time()
        # train GP hyper-parameters only when number of training data is increased by more than 10%
        if data_size*1.1 < len(x):
            data_size = len(x)
            mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
            fit_gpytorch_model(mll,options=scipy_minimize_options)
        t1 = time.time()

        # increase bounds over optimization epochs starting from reduced bounds by factor of "initial_boundary_ratio"
        UCB = UpperConfidenceBound(gp, beta=LCB_nsig**0.5)
        x1, y1acqu = optimize_acqf(UCB, 
                                   bounds=bounds_torch, 
                                   q=1, **acquisition_optimize_options,
                                   options=scipy_minimize_options)
        t2 = time.time()
        if verbose==2:
            print("surrogate model training time (sec):",t1-t0,", acqusition optimization time (sec):",t2-t1)

        y1 = torch.tensor([-func(x1.detach().numpy()[0,:],*args)])[:,None]
        x = torch.concat((x,x1),axis=0)
        y = torch.concat((y,y1),axis=0)
        # add GP training data
        gp = gp.get_fantasy_model(x1,y1)

        result.surrogate.training_time.append(t1-t0)
        result.surrogate.acqusition_optimization_time.append(t2-t1)
        result.history.x = x.detach().numpy()
        result.history.y = -y.detach().numpy().flatten()
        imin = result.history.y.argmin()    
        result.x = result.history.x[imin]
        new_y_min = result.history.y[imin]
        result.history.y_min.append(new_y_min)
        if new_y_min <= y_min:
            i_patience = 0
            y_min = new_y_min
        else:
            i_patience += 1
        
    imin = result.history.y.argmin()    
    result.x = result.history.x[imin]
    result.fun = result.history.y[imin]
    result.nfev = i_eval+1

    if result.nfev == maxnfev:
        result.message = "maximum number of function evaluation reached"
    if returnGP:
        result.GP = gp
    if set_func_to_best:
        func(result.x)
    
    return result

