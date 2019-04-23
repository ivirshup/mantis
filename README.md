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
