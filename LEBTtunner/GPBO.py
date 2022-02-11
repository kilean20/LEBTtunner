import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf


from .dictClass import dictClass
import pickle


def minimize(func,x0,bounds,maxnfev, args=(), LCB_nsig=1, patience=None, save_optimization_history=True, returnGP=False, verbose=True, set_func_to_best = True ):
    '''
    func: function to obtimize
    x0: (numpy vector) initial function argument
    bounds: limit of function arguments so that [ bounds[i][0] <= x[i] <= bounds[i][1] for i in range(len(x)) ]
    maxnfev: (int) maximum number of function evaluation
    LCB_nsig: lower confidence bound in unit of stantard deviation of GP prediction 
    patience: quit optimization if the function value is not improving for patience number of function evaluation
    '''
    ndim = len(x0)
    bounds_torch = torch.tensor(bounds)
    result = dictClass
    result.history = dictClass
    
    x = torch.unsqueeze(torch.tensor(x0),0)
    y = torch.unsqueeze(torch.tensor(-func(x0,*args)),0)
    
    result.history.x = x.detach().numpy()
    result.history.y = -y.detach().numpy()
    y_min = result.history.y[0]
    result.history.y_min = [y_min]   
    result.x = result.history.x[0]
    
    i_patience = 0 
    for i_eval in range(1,maxnfev):
        
        if verbose:
            print("epoch :{0:3d}/{0:3d}, min_fun :{1:.3e}".format(i_eval, maxfev, y_min))
            print("decision x: ", result.x)
            
        if i_patience >= patience:            
            result.message = "maximum patience reached"
            break
            
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