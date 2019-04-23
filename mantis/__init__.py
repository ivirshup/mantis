import scanpy as sc
from anndata import AnnData

import pandas as pd
import numpy as np
import xarray as xr

from functools import partial, wraps
from inspect import isfunction, getmembers
from typing import Iterable, Optional

############################################
# Dot chaining
############################################

# Make scanpy api accessible through dot chaining

class Mapper(object):
    def __init__(self, parent, mod):
        self._parent = parent
        self._mod = mod
        self._modfuncs = getmembers(mod, isfunction)
        for _, f in self._modfuncs:
            add_method(self, self._parent)(f)


def add_method(obj, firstarg=None):
    if firstarg is None:
        firstarg = obj

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            ret = partial(func, firstarg)(*args, **kwargs)
            if ret is None:  # Hacky workaround for chaining inplace operations
                ret = firstarg
            return ret
        setattr(obj, func.__name__, wrapper)
        # Not binding func, but a wrapper which acts like it
        return func  # returning func means func can still be used normally, only meaningful for decorator usecase
    return decorator

AnnData.pl = property(partial(Mapper, mod=sc.pl))
AnnData.tl = property(partial(Mapper, mod=sc.tl))
AnnData.pp = property(partial(Mapper, mod=sc.pp))

# Groupby

class Grouped(object):
    """
    Generic grouped object for split-apply-combine style workflows.
    
    Designed to be used with chainable calls.
    """

    def __init__(self, groups: dict):
        self.groups = groups

    def __iter__(self):
        for k, v in self.groups.items():
            yield (k, v)

    def apply(self, func):
        """
        Apply func to each value, returning a Grouped result.
        """
        return Grouped({k: func(v) for k, v in self})

    def mutate(self, func):
        """
        Apply mutating func to each value.
        
        In contrast with apply, this method assumes passed `func` modifies 
        the values in place, and return the current Grouped.
        """
        for k, v in self:
            func(v)
        return self

    def combine(self, func):
        return func(self.groups)


def groupby(adata: AnnData, *, obs=None, var=None) -> Grouped:
    """
    Group anndata by obs or var dataframes.

    Example
    -------
    >>> pbmc = sc.datasets.pbmc68k_reduced()
    >>> means = (
            pmbc.groupby(obs="louvain")
            .apply(lambda x: pd.Series(np.ravel(x.raw.X.mean(axis=0)), index=x.var_names))
            .combine(pd.DataFrame)
        )
    """
    is_obs = obs is not None
    is_var = var is not None
    if is_obs == is_var:
        raise ValueError("Can only groupby obs or var")
    if is_obs:
        dim = 0
        grouped = adata.obs.groupby(obs)
        def make_idx(x): return (x, slice(None))
    elif is_var:
        dim = 1
        grouped = adata.var.groupby(var)
        def make_idx(x): return (slice(None), x)
    return Grouped({k: adata[make_idx(v.index)] for k, v in grouped})


def select(adata, obs=None, var=None):
    """
    Subset anndata, but in a chainable way.

    Example
    -------
    >>> (adata
            .select(obs=lambda x: x["batch"] == 1)
            .pl.umap()
        )
    """
    is_var = var is None
    is_obs = obs is None
    if is_var and is_obs:
        raise ValueError()
    if obs:
        idx = adata.obs.loc[obs].index
        return adata[idx, :]
    if var:
        idx = adata.var.loc[var].index
        return adata[:, idx]


setattr(AnnData, "groupby", groupby)
setattr(AnnData, "select", select)
