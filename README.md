# Fair Logistic Regression via Penalized Mean-Field Variational Inference

Python implementation for the experiments in the accompanying paper.

- `main.py` runs cross-validation for German Credit and Student Performance.
- `optimizer.py`, `vb_logit_fit.py` and `vb_logit_pred.py` handle hyperparameter selection, CAVI fitting and prediction. The data loaders and `cross_validation.py` prepare the inputs.
- `fair_mfvi_revised_pipeline.ipynb` also runs ACSIncome experiments. In Jupyter, set `RUN_FRESH` to `True` and run the setup and ACSIncome cells. This runs the current optimizer. The paper's historical ACS fits are represented by the saved results.
- `graphs.py` draws figures from saved results. `figures/` holds the exact published panels.

## Run

Use Python 3.14 in a virtual environment. Run these commands from the repository root.

```sh
python -m pip install -r requirements.txt
python main.py --dataset german --results_dir runs/german
python main.py --dataset student --results_dir runs/student
python graphs.py
```

Fitting can take substantial time. `python graphs.py` redraws the supported figures under `runs/figures/` without fitting. Use `python graphs.py --dataset german --results-dir runs/german` to plot a fresh run.

## Citation and data

Please cite Emilio Carrizosa, Rafael Jiménez-Llamas and Pepa Ramírez-Cobo, *Fair Logistic Regression via Penalized Mean-Field Variational Inference*. The variational Bayesian logistic regression routines build on Jan Drugowitsch's [VBLinLogit](https://github.com/DrugowitschLab/VBLinLogit).

The UCI datasets are [German Credit](https://archive.ics.uci.edu/dataset/144/statlog-german-credit-data) (Hans Hofmann, 1994, processed for this study) and [Student Performance](https://archive.ics.uci.edu/dataset/320/student%2Bperformance) (Paulo Cortez, 2014), both under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The Nevada person file is from the U.S. Census Bureau's [2018 ACS 1-Year PUMS](https://www.census.gov/programs-surveys/acs/microdata/access.2018.html).
