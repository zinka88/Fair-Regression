# Fair-Regression

Code and data for "Fair Regression for Health Care Spending" by Anna Zink and Sherri Rose (2019), *Biometrics*, [doi:10.1111/biom.13206](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13206). 

All code is written in R and was run on R version 3.5.1. Required packages are CVXR (the analysis was run on version 0.99) and dplyr (the analysis was run on version 0.7.6).

1. **Analysis**

   **Code:** *fairness_for_risk_adjustment.R* runs the main analysis and depends on functions defined in *fairness_functions.R*.

   **Data:** *simulated_analysis_data.csv* (variable descriptions in *simulated_analysis_data_dictionary.xlsx*)

2. **Simulation**

  	**Code:** *run_simulation.R* runs the simulation and depends on functions defined in *sim_function.R* and         *fairness_functions_sim.R*.

  	**Data:** *simulation_data.csv* (variable descriptions in *simulation_data_dictionary.xlsx*)

