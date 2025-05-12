# Miniproject Feedback and Assessment

## Report

**"Guidelines" below refers to the MQB report [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) [here](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html) (which were provided to the students in advance).**

**Title:** “Systematic Comparison of Six Microbial Growth Models: Evaluation of Convergence, Parameter Estimation, and Predictive Performance”

- **Introduction (15%)**  
  - **Score:** 12/15  
  - Provides context on predictive microbiology, referencing multiple models. Could further highlight the direct question or gap.

- **Methods (15%)**  
  - **Score:** 13/15  
  - Thorough approach: from 4,387 → 210 subsets, multi-start method. Additional param-bound detail could be included for completeness.

- **Results (20%)**  
  - **Score:** 15/20  
  - Gompertz tends to be top, Baranyi second, logistic failing more. Additional numeric or table detail recommended.

- **Tables/Figures (10%)**  
  - **Score:** 6/10  
  - The snippet references a “big table” but not deeply integrated. The [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) suggest explicit referencing in text.

- **Discussion (20%)**  
  - **Score:** 15/20  
  - Describes polynomial interpretability issues, Gompertz's advantage, logistic's issues. Could expand on data-limitation or next steps more.

- **Style/Structure (20%)**  
  - **Score:** 16/20  
  - Organized with minimal style errors, though references are sometimes broad.

**Summary:** A strong multi-model comparison, systematically done. More detailed numeric breakdown and thorough figure usage would be beneficial.

**Report Score:** 74  

---

## Computing

### Project Structure & Workflow

**Strengths**

* **Modular design**: Six focused R scripts each address a distinct stage of the analysisShell drivers (`run.Miniproject.sh` and `run_pipeline.sh`) enable the entire pipeline and LaTeX report to run with a single command (why two?).
* Dedicated `code/`, `data/`, and `results/` folders.

**Suggestions**

1. Consolidate shell drivers into one robust script (e.g. `run_pipeline.sh`), with:

   * Portable shebang: `#!/usr/bin/env bash`
   * Strict mode: `set -euo pipefail`
   * Directory guard: `cd "$(dirname "$0")"`
   * Timestamped logging: `|& tee "../results/pipeline_$(date +'%Y%m%d_%H%M').log"`
   * Clear exit codes on failure.

2. **Lock environments**:

   * **R**: Initialize `renv` in `code/`, commit `renv.lock`, and run `renv::restore()` in the driver.
   * **System tools**: Document required LaTeX and other OS dependencies in the README.

---

### README File

**Strengths**

* Good overview of goals, directory layout, prerequisites, and usage.
* Dependency listing for R and Python packages.

**Suggestions**

1. Verify script names and paths in the README match the repository (e.g. `run_pipeline.sh`).
2. Add environment setup:

   ```bash
   git clone <repo_url> && cd MiniProject/code
   Rscript -e "renv::restore()"
   bash run_pipeline.sh
   ```
3. **Embed directory tree** under **Directory Structure**.
4. **Summarize scripts** under **code/** with inputs and outputs.
5. **Include license and data citation** in a **Data** section.

---

## R Script Structure & Syntax

####  `01_data_preprocessing.R`

* **Eliminate global side effects**: Remove `rm(list = ls())`; use clean sessions or wrap in `main()`.
* **Portable paths**: Replace relative paths with `here::here('data', 'LogisticGrowthData.csv')`.
* **Vectorized IDs**: Use `mutate(ID = as.integer(factor(paste(Species, Temp, Medium, Citation))))`.
* **Functional design**: Encapsulate logic into functions (e.g. `clean_data()`) and apply with `purrr::map_df()`.
* **Progress logging**: Use `cli::cli_alert_info()` to report rows read, filtered, and saved.

####  `02_define_and_fit_models.R`

* Move model formulas into `models.R` and source once.
* Adopt `nls.multstart::nls_multstart()` to handle random starts and bounds; e.g.:

  ```r
  nls_multstart(
    logPopBio ~ logistic_model(Time, r_max, K, N_0),
    data, iter = 100,
    start_lower = c(r_max=0, K=min, N_0=min),
    start_upper = c(r_max=5, K=max, N_0=max)
  )
  ```
* Replace `for` loops with `split(cleaned_data, .$ID) %>% purrr::map()`.
* Call `set.seed()` once at the top.

####  `03_compare_models_and_akaike_weights.R`

* Consolidate shared functions in `utils.R` and `source()` it.
* Bind model fits into a long tibble, then use `slice_min(AICc)` / `slice_min(BIC)` and `MuMIn::Weights()`.
* Create a `theme_mproj <- theme_minimal()` object for consistent styling.

#### `04_visualize_parameters.R`

* Source multi‐start and model functions from `utils.R`.
* Use `broom::augment()` to combine fits with data.
* Consiude Looping with purrr:

  ```r
  param_df %>%
    group_split(term) %>%
    purrr::iwalk(~ ggsave(paste0(...), make_param_boxplot(.x)))
  ```

*Compute `coord_cartesian()` limits programmatically or use free scales in facets.

### `05_plot_each_id.R`

* **Function composition**: Define `plot_id(id, df)` returning a `ggplot` and then `purrr::walk(ids, ~ ggsave(..., plot = plot_id(.x)))`.
* **ggplot2**: Replace base‑device loops with `ggsave()` and layered geoms.
* **Panel assembly**: Use `patchwork::wrap_plots()` to collect multiple IDs and shared legends.

---

## Shell Driver Enhancements

**Consolidate** the two drivers into one `run_pipeline.sh` that:

* Sets `#!/usr/bin/env bash` and `set -euo pipefail`.
* Changes to the `code/` directory.
* Executes scripts 01–06 in order, logging output.
* Compiles the LaTeX report, directing the PDF to `../results/`.
* Accepts flags for `--data-dir`, `--results-dir`, and `--n-starts`.

---

## NLLS Fitting Approach

**Strengths**

* Comprehensive suite: Logistic, Gompertz, Baranyi, Three‑Phase, plus polynomial benchmarks.
* Multi‑start strategy improves convergence across heterogeneous data.
* In‑sample metrics (AICc, BIC, R²) provide a balanced selection framework.

**Opportunities**

1. **nls.multstart**: Streamline multi‑start with built‑in bounds, diagnostics, and CIs.
2. **Parallel execution**: Use `furrr` or `future.apply` to distribute fits across cores.
3. **Biological bounds**: Enforce `r_max ≥ 0`, `K ≥ N0`, and flag boundary hits.
4. **Cross‑validation**: Implement leave‑one‑timepoint‑out CV for predictive assessment.
5. **Akaike weights plots**: Visualize weight distributions with density or ridge plots.
6. **Diagnostics log**: Export convergence codes, errors, and runtimes per model/ID.

---

## Summary

Your pipeline is well‑structured; To further enhance:

* Lock environments with `renv` and lockfiles.
* Refactor repetitive code into shared utility functions.
* Leverage specialized packages (`nls.multstart`, `purrr`, `ggplot2`, `patchwork`).
* Implement robust logging and parallelization.

### Score: 66

---

## Overall Score: (74+66)/2 = 70