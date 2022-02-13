import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf
import numpy as np


from .dictClass import dictClass
import pickle


def init_population_latinhypercube(len_of_vec,num_population,bounds,rng=None):
    """
    Initializes the population with Latin Hypercube Sampling.
    Latin Hypercube Sampling ensures that each parameter is uniformly
    sampled over its range.
    """
    if rng is None:
        rng = np.random
        
    segsize = 1.0 / num_population
    samples = (segsize * rng.uniform(size=(num_population,len_of_vec))
               + np.linspace(0., 1., num_population,
                             endpoint=False)[:, np.newaxis])

    population = np.zeros_like(samples)
    
    # Initialize population of candidate solutions by permutation
    for j in range(len_of_vec):
        order = rng.permutation(range(num_population))
        population[:, j] = samples[order, j]

    return population*(bounds[:,1]-bounds[:,0])+bounds[:,0]



def minimize_GP_LCB(func, bounds,
                    x0=None, y0=None, 
                    maxnfev=None, 
                    n_init_pop=None, args=(), LCB_nsig=1, patience=0, save_optimization_history=True, returnGP=False, verbose=True, set_func_to_best = True ):
    '''
    func: function to obtimize
    x0: (numpy vector) initial function argument
    bounds: limit of function arguments so that [ bounds[i][0] <= x[i] <= bounds[i][1] for i in range(len(x)) ]
    maxnfev: (int) maximum number of function evaluation
    LCB_nsig: lower confidence bound in unit of stantard deviation of GP prediction 
    patience: quit optimization if the function value is not improving for patience number of function evaluation
    '''
    if n_init_pop is None:
        x0_shape = np.array(x0).shape
        if len(x0_shape)>1:
            n_init_pop = x0_shape[0]
        n_init_pop = len(bounds)
    
    x = init_population_latinhypercube(len(bounds),n_init_pop,bounds)
    if x0 is not None:
        x0_ = np.array(x0)
        if len(x0_.shape)==1:
            x[0,:] = x0_
        else:
            
            
        
    m,n = x0.shape
    ndim = len(x0)
    bounds_torch = torch.tensor(bounds)
    result = dictClass
    result.history = dictClass
    
    x = torch.tensor(x0)[None,:]
    y = torch.tensor([-func(x0,*args)])[None,:]
    print("GPBO got y")
    
    result.history.x = x.detach().numpy()
    result.history.y = -y.detach().numpy()
    y_min = result.history.y[0]
    result.history.y_min = [y_min]   
    result.x = result.history.x[0]
    
    i_patience = 0 
    for i_eval in range(1,maxnfev):
        if verbose:
            print("epoch :{0:3d}/{0:3d}, min_fun :{1:.3e}".format(i_eval, maxnfev, y_min))
            print("decision x: ", result.x)
        print(i_patience, patience)
        if i_patience >= patience and patience>0:       
            result.message = "maximum patience reached"
            break
            
        print('x.shape,y.shape',x.shape,y.shape)
        gp = SingleTaskGP(x,y)
        mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
        fit_gpytorch_model(mll)
        UCB = UpperConfidenceBound(gp, beta=LCB_nsig**0.5)
        x1, y1acqu = optimize_acqf(UCB, bounds=bounds_torch, q=1, num_restarts=20, raw_samples=20)
        y1 = torch.unsqueeze(torch.tensor(-func(x1.detach().numpy(),*args)),0)
        y1 = torch.tensor(y1)
        
        x = torch.concat((x,x1),axis=0)
        y = torch.concat((y,y1),axis=0)
        
        result.history.x = x.detach().numpy()
        result.history.y = -y.detach().numpy()
   
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
    result.nfev = i_eval
    if result.nfev == maxnfev:
        result.message = "maximum function evaluation reached"

    if returnGP:
        result.GP = gp
        
    if set_func_to_best:
        func(result.x)
    
    return result
