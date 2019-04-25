# mantis

Experimental APIs to give [Scanpy](https://github.com/theislab/scanpy) 5 times the colors and rocket hands.

## Current features

### Method chaining

```python
import scanpy as sc
import mantis

import pandas as pd
import numpy as np

pbmc = sc.datasets.pbmc3k()

pbmc_norm = (pbmc.copy()
    .pp.calculate_qc_metrics(inplace=True)
    # .pipe(lambda x: x.layers.__setitem__("counts", x.X))  # Not implemented yet
    .pp.normalize_per_cell(counts_per_cell_after=1000)
    .pp.log1p()
    .pp.pca()
)
pbmc_norm.pl.pca()
```

### Group by

Very early stages, likely to change. There's a good chance this example is out of date.

```python
pbmc = sc.datasets.pbmc68k()
# Get mean expression for a gene by cluster:
(pbmc
    .groupby(obs="louvain")
    .apply(lambda x: x.apply(lambda x: pd.Series(np.ravel(x.raw.X.mean(axis=0)), index=x.var_names)))
    .combine(pd.DataFrame)
)
```

## To-do 

* Merging anndata objects, the `combine` part of `split-apply-combine`
* Figure out how to allow mutations with callables
    * `.mutate_{attr}({attr_key}={callable})`?

## Ideas

### Mostly copy what xarray does

```python
import scanpy as sc
import xarray as xr
import pandas as pd
import numpy as np

adata = sc.datasets.pbmc68k_reduced()
adata.obs_names.name = "obs"
adata.var_names.name = "var"

ds = xr.Dataset(
    {
        "X": (["obs", "var"], adata.X),
        "X_raw": (["obs", "var"], adata.raw.X.toarray()),
        "louvain": (["obs"], adata.obs["louvain"]),
        "pca": (["obs", "pc"], adata.obsm["X_pca"]),
        "pcs": (["var", "pc"], adata.varm["PCs"]),
        "umap": (["obs", "umap_comp"], adata.obsm["X_umap"]),
    },
    coords={
        "var": adata.var_names,
        "obs": adata.obs_names
    }
)

groupedX = ds[["X", "louvain"]].groupby("louvain")
mean = groupedX.mean(dim="obs")
var = groupedX.var(dim="obs")
```