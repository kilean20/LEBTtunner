import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf
import numpy as np


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



def minimize_GP_LCB(func, x0, bounds,
                    maxnfev=None, initial_boundary_ratio = 1.0,
                    n_init_pop = None, args=(), LCB_nsig=1, patience=0, save_optimization_history=True, returnGP=False, verbose=True, set_func_to_best = True ):
    '''
    func: function to obtimize
    x0: (numpy vector) initial function argument
    bounds: limit of function arguments so that [ bounds[i][0] <= x[i] <= bounds[i][1] for i in range(len(x)) ]
    maxnfev: (int) maximum number of function evaluation
    LCB_nsig: lower confidence bound in unit of stantard deviation of GP prediction 
    patience: quit optimization if the function value is not improving for patience number of function evaluation
    '''
    ndim = len(x0)
    if maxnfev == None:
        maxnfev = 10000
    if n_init_pop == None:
        n_init_pop = min(ndim,maxnfev)
    
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
    
    i_patience = 0 
    for i_eval in range(n_init_pop,maxnfev):
        if verbose:
            print("epoch: {0:3d}/{1:3d}, min_fun :{2:.3e}".format(i_eval, maxnfev, y_min))
            print("decision x: ", result.x)
        if i_patience >= patience and patience>0:       
            result.message = "maximum patience reached. optimization halt..."
            break
            
        gp = SingleTaskGP(x,y)
        mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
        fit_gpytorch_model(mll)
        UCB = UpperConfidenceBound(gp, beta=LCB_nsig**0.5)
        reduction = initial_boundary_ratio + (1-initial_boundary_ratio)*(i_eval+1-n_init_pop)/(maxnfev-n_init_pop) 
        bounds_tmp = bounds_torch*reduction
        
        x1, y1acqu = optimize_acqf(UCB, 
                                   bounds=bounds_tmp, 
                                   q=1, num_restarts=20, raw_samples=20)
 
        y1 = torch.tensor([-func(x1.detach().numpy()[0,:],*args)])[:,None]
        
        x = torch.concat((x,x1),axis=0)
        y = torch.concat((y,y1),axis=0)
        
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
