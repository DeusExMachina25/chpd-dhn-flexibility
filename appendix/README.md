# The authors' online appendix

The case study in `data/` is built from the online appendix of

> L. Mitridati and J. A. Taylor, *Power systems flexibility from district
> heating networks*, PSCC 2018.

published at doi:[10.5281/zenodo.1195508](https://doi.org/10.5281/zenodo.1195508).

Only the downloaded archives are committed:

| File | Contents |
|---|---|
| `online-appendix.zip` | the Zenodo record (MD5 `a7388358c72b622a3397f99a620f8403`, verified) |
| `repo.zip` | the authors' reference Python implementation |

Unzip either one here to inspect them:

```bash
unzip online-appendix.zip
```

The extracted directories are deliberately **not** tracked. Their names run to
178 characters, which exceeds the 260-character `MAX_PATH` limit on Windows once
a clone prefix is added, and makes the repository fail to check out there unless
`git config --global core.longpaths true` is set. Keeping only the archives lets
the repository clone anywhere while preserving the exact bytes the data came
from.

## What was taken from where

`data/network.yaml` and `data/timeseries.csv` come from
`step 1.a initialization MILP MINI 3.py` in the reference implementation, which
is the only complete and self-consistent source: the appendix PDF's tables are
ambiguous once the PDF is de-laid-out, and its load and wind profiles appear
there only as figures.

Section 6 of `../NOTES.md` records the provenance of every parameter, and the
places where the appendix PDF and the authors' own code disagree: pressure
bounds, minimum heat exchanger flow, water density, and which generating units
are present. Where they conflict, the code was followed, because it is what
produces the published numbers.
