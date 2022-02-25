from typing import Any, Dict, List, Optional, Union
import torch
from botorch import settings
from botorch.models.gpytorch import BatchedMultiOutputGPyTorchModel
from botorch.models.transforms.input import InputTransform
from botorch.models.transforms.outcome import OutcomeTransform
from botorch.models.utils import validate_input_scaling
from botorch.utils.containers import TrainingData
from gpytorch.constraints.constraints import GreaterThan
from gpytorch.distributions.multivariate_normal import MultivariateNormal
from gpytorch.kernels.matern_kernel import MaternKernel
from gpytorch.kernels.scale_kernel import ScaleKernel
from gpytorch.likelihoods.gaussian_likelihood import GaussianLikelihood
from gpytorch.likelihoods.likelihood import Likelihood
from gpytorch.means.constant_mean import ConstantMean
from gpytorch.models.exact_gp import ExactGP
from gpytorch.module import Module
from gpytorch.priors.torch_priors import GammaPrior
from torch import Tensor


#from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from botorch.acquisition import UpperConfidenceBound
from botorch.optim import optimize_acqf
from gpytorch.mlls import ExactMarginalLogLikelihood

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



MIN_INFERRED_NOISE_LEVEL = 1e-4
class SingleTaskGP(BatchedMultiOutputGPyTorchModel, ExactGP):
    r"""A single-task exact GP model.

    A single-task exact GP using relatively strong priors on the Kernel
    hyperparameters, which work best when covariates are normalized to the unit
    cube and outcomes are standardized (zero mean, unit variance).

    This model works in batch mode (each batch having its own hyperparameters).
    When the training observations include multiple outputs, this model will use
    batching to model outputs independently.

    Use this model when you have independent output(s) and all outputs use the
    same training data. If outputs are independent and outputs have different
    training data, use the ModelListGP. When modeling correlations between
    outputs, use the MultiTaskGP.
    """

    def __init__(
        self,
        train_X: Tensor,
        train_Y: Tensor,
        likelihood: Optional[Likelihood] = None,
        covar_module: Optional[Module] = None,
        outcome_transform: Optional[OutcomeTransform] = None,
        input_transform: Optional[InputTransform] = None,
    ) -> None:
        r"""A single-task exact GP model.

        Args:
            train_X: A `batch_shape x n x d` tensor of training features.
            train_Y: A `batch_shape x n x m` tensor of training observations.
            likelihood: A likelihood. If omitted, use a standard
                GaussianLikelihood with inferred noise level.
            covar_module: The module computing the covariance (Kernel) matrix.
                If omitted, use a `MaternKernel`.
            outcome_transform: An outcome transform that is applied to the
                training data during instantiation and to the posterior during
                inference (that is, the `Posterior` obtained by calling
                `.posterior` on the model will be on the original scale).
            input_transform: An input transform that is applied in the model's
                forward pass.

        Example:
            >>> train_X = torch.rand(20, 2)
            >>> train_Y = torch.sin(train_X).sum(dim=1, keepdim=True)
            >>> model = SingleTaskGP(train_X, train_Y)
        """
        with torch.no_grad():
            transformed_X = self.transform_inputs(
                X=train_X, input_transform=input_transform
            )
        if outcome_transform is not None:
            train_Y, _ = outcome_transform(train_Y)
        self._validate_tensor_args(X=transformed_X, Y=train_Y)
        ignore_X_dims = getattr(self, "_ignore_X_dims_scaling_check", None)
        validate_input_scaling(
            train_X=transformed_X, train_Y=train_Y, ignore_X_dims=ignore_X_dims
        )
        self._set_dimensions(train_X=train_X, train_Y=train_Y)
        train_X, train_Y, _ = self._transform_tensor_args(X=train_X, Y=train_Y)
        if likelihood is None:
            noise_prior = GammaPrior(1.1, 0.05)
            noise_prior_mode = (noise_prior.concentration - 1) / noise_prior.rate
            likelihood = GaussianLikelihood(
                noise_prior=noise_prior,
                batch_shape=self._aug_batch_shape,
                noise_constraint=GreaterThan(
                    MIN_INFERRED_NOISE_LEVEL,
                    transform=None,
                    initial_value=noise_prior_mode,
                ),
            )
        else:
            self._is_custom_likelihood = True
        ExactGP.__init__(self, train_X, train_Y, likelihood)
        self.mean_module = ConstantMean(batch_shape=self._aug_batch_shape)
        if covar_module is None:
            self.covar_module = ScaleKernel(
                MaternKernel(
                    nu=2.5,
                    ard_num_dims=transformed_X.shape[-1],
                    batch_shape=self._aug_batch_shape,
                    lengthscale_prior=GammaPrior(3.0, 6.0),
                ),
                batch_shape=self._aug_batch_shape,
                outputscale_prior=GammaPrior(2.0, 0.15),
            )
            self._subset_batch_dict = {
                "likelihood.noise_covar.raw_noise": -2,
                "mean_module.constant": -2,
                "covar_module.raw_outputscale": -1,
                "covar_module.base_kernel.raw_lengthscale": -3,
            }
        else:
            self.covar_module = covar_module
        # TODO: Allow subsetting of other covar modules
        if outcome_transform is not None:
            self.outcome_transform = outcome_transform
        if input_transform is not None:
            self.input_transform = input_transform
        self.to(train_X)

    def forward(self, x: Tensor) -> MultivariateNormal:
        # == Kilean == 
        #if self.training:
        #    x = self.transform_inputs(x)
        x = self.transform_inputs(x)
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return MultivariateNormal(mean_x, covar_x)


    @classmethod
    def construct_inputs(
        cls, training_data: TrainingData, **kwargs: Any
    ) -> Dict[str, Any]:
        r"""Construct kwargs for the `Model` from `TrainingData` and other options.

        Args:
            training_data: `TrainingData` container with data for single outcome
                or for multiple outcomes for batched multi-output case.
            **kwargs: None expected for this class.
        """
        return {"train_X": training_data.X, "train_Y": training_data.Y}




def minimize_GP_LCB(func, x0, bounds, args=(),
                    n_init_pop = None, maxnfev=10000, patience=1000,
                    LCB_nsig = 1.0, 
                    trainGP_everyEpoch=True,
                    input_transform = None,
                    acquisition_optimize_options = {"num_restarts":20, "raw_samples":20},
                    scipy_minimize_options={},# = {"ftol":2e-3},
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
                                         LCB_nsig=LCB_nsig, trainGP_everyEpoch=trainGP_everyEpoch,
                                         acquisition_optimize_options = acquisition_optimize_options,
                                         scipy_minimize_options=scipy_minimize_options,
                                         returnGP = returnGP,   
                                         verbose = verbose, 
                                         set_func_to_best=set_func_to_best )
    ndim = len(x0)
    if n_init_pop == None:
        n_init_pop = min(2*ndim,maxnfev)
    
    bounds_np = np.array(bounds).astype(np.float64)
    bounds_torch = torch.tensor(np.transpose(bounds_np))

    x = np.zeros([n_init_pop,ndim])
    x[0 ,:] = x0[:]
    if n_init_pop>1:
        x[1:,:] = init_population_latinhypercube(n_init_pop-1, ndim, bounds_np)
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

    gp = SingleTaskGP(x,y,input_transform=input_transform)
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
            gp = SingleTaskGP(x,y,input_transform=input_transform)
            mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
            fit_gpytorch_model(mll,options=scipy_minimize_options)
        elif data_size*1.1 < len(x):
            data_size = len(x)
            mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
            fit_gpytorch_model(mll,options=scipy_minimize_options)
        t1 = time.time()

        if LCB_nsig=="NoRegret":
            beta = (1*(2*np.log((i_eval**(ndim/2.+2))*(np.pi**2)/0.3)))**0.5
        elif LCB_nsig=="Greedy":
            if i_eval < 0.8*maxnfev:
                beta = 2*(1-i_eval/(0.8*maxnfev))
            else:
                beta = 2e-3
        else:
            beta = LCB_nsig**0.5
        UCB = UpperConfidenceBound(gp, beta=beta)
        x1, y1acqu = optimize_acqf(UCB, 
                                   bounds=bounds_torch, 
                                   q=1, **acquisition_optimize_options,
                                   options = scipy_minimize_options)
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
                              LCB_nsig=1.0,trainGP_everyEpoch=True,
                              acquisition_optimize_options = {"num_restarts":20, "raw_samples":20},
                              scipy_minimize_options={},# = {"ftol":2e-3},
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

    bounds_np = np.array(bounds).astype(np.float64)
    bounds_torch = torch.tensor(np.transpose(bounds_np))

    if hasattr(result,"GP"):
        gp = result.GP
    else:
        gp = SingleTaskGP(x,y,input_transform=input_transform)
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
        if trainGP_everyEpoch:
            gp = SingleTaskGP(x,y,input_transform=input_transform)
            mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
            fit_gpytorch_model(mll,options=scipy_minimize_options)
        elif data_size*1.1 < len(x):
            data_size = len(x)
            mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
            fit_gpytorch_model(mll,options=scipy_minimize_options)
        t1 = time.time()

        if LCB_nsig=="NoRegret":
            beta = (1*(2*np.log((i_eval**(ndim/2.+2))*(np.pi**2)/0.3)))**0.5
        elif LCB_nsig=="Greedy":
            if i_eval < 0.8*maxnfev:
                beta = 2*(1-i_eval/(0.8*maxnfev))
            else:
                beta = 2e-3
        else:
            beta = LCB_nsig**0.5
        UCB = UpperConfidenceBound(gp, beta=beta)
        x1, y1acqu = optimize_acqf(UCB, 
                                   bounds=bounds_torch, 
                                   q=1, **acquisition_optimize_options,
                                   options = scipy_minimize_options)
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

