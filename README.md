# Online conformal inference for multi-step time series forecasting

Source for the paper and the code that produces its simulations, applications, and figures. The proposed methods are implemented in the [`conformalForecast`](https://github.com/xqnwang/conformalForecast/) R package; this repository uses version **[0.2.0.9000](https://github.com/xqnwang/conformalForecast/tree/e834ab8bbfd1efe6577b183af11a4b3e2a6ffc97)**.

The R scripts are in [`code/`](code), saved results and plots in [`result/`](result), and manuscript source in [`paper/cpts.qmd`](paper/cpts.qmd).

## Reproducing the results

Install [uvr](https://github.com/nbafrank/uvr#installation) for R and [uv](https://docs.astral.sh/uv/getting-started/installation/) for Python. The project uses R 4.4.3 and Python 3.9.25. If R 4.4.3 is not already available, install it with `uvr r install 4.4.3`. Then, from the repository root:

```sh
uv python install
uv sync --locked
uvr sync --frozen --no-companion
```

These commands restore the packages pinned in `uv.lock` and `uvr.lock` into `.venv/` and `.uvr/library/`.

Run all experiments from the repository root in the following order. The main scripts save forecasts to `result/*_fc.rds`, which the parameter sensitivity scripts reuse.

```sh
# Main simulations and applications
uvr run code/Sim_AR2.R
uvr run code/Sim_NL.R
uvr run code/Real_elec.R
uvr run code/Real_cafe.R

# Parameter sensitivity analyses
uvr run code/Sim_AR2_param.R
uvr run code/Sim_NL_param.R
uvr run code/Sim_elec_param.R
uvr run code/Sim_cafe_param.R

# Known versus forecast future predictors
uvr run code/Sim_predictor.R
```

The helper files `clip_acp.R` and `spci.R` are sourced by the experiment scripts and do not need to be run separately.

[`code/spci.R`](code/spci.R) adapts SPCI to the shared point forecasts and calibration windows, using out-of-sample errors and fixed tail probabilities of `alpha / 2`. It calls Python's `sklearn_quantile` through `reticulate`. Use `random_state` to control the forest's randomness; `set.seed()` alone is insufficient.

## Building the paper

Install [Quarto](https://quarto.org/) and a TeX distribution. From the repository root, in Bash or Zsh:

```sh
source .uvr/activate
quarto render paper/cpts.qmd
```
