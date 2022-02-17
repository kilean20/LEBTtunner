# dictClass
from numpy import isinf

def is_notebook():
    try:
        from IPython import get_ipython
        if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def loss_regularizer(loss_func,loss_ifinf,args=(),lower_limit=None,upper_limit=None,flip_sign=False, loss_mean=0, loss_std=1):
    def regularized_lossfunc(x):
        loss = loss_func(x,*args)
        if flip_sign:
            loss = -loss
        if lower_limit is not None:
            loss = max([lower_limit,loss])
        if upper_limit is not None:
            loss = min([upper_limit,loss])
        if isinf(loss):
            loss = loss_ifinf

        return (loss-loss_mean)/loss_std

    return regularized_lossfunc
