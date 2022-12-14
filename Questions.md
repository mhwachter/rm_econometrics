# RM Econometrics
## Topics
- ML for selection of OLS control variables, sensitivity analysis
- ML for choice of instruments for IV estimation
- ML for selection of control variables in IV models
---
## Structure
1. Introduction
2. Theory: Theory of Lasso. Evaluation of properties of lasso.
3. Simulations: 
4. Applications
5. Conclusion
---
## Questions
- All three topics?
- R or Stata for simulation and application?
- Replicate random forest results?
- Find the real data?
- Code from CompStat?
- Data generating process for simulation study? How many different should we consider?
---
## Meeting
Focus on Lasso. Group that works on Lasso estimation. Brief introduction in Lasso estimation. Properties of Lasso. How to choose the tuning parameter. 

Go through steps of covariance selection. Lasso for IV selection not doing good jobs and then pin down issues.

Connection between theory and application part. Assumptions for Lasso estimation. Sparsity assumption and low dimension of coefficients. 
---
## Allocation of Topics
- Cristian: Theory
    - Should I try a simulation with normally distributed regressors or uniformly distributed is enough?
    - What happens under different penalties?
        - According to literature prediction does not improve. This simulation still needs to be done. Did Marcel figured out how to model the penalties from Springler codes ?
    - How should we defined the tunning parameter (lambda)?
    - How do glmnet package manage deliver results for OLS under          high-dimensionality ?
- Marcel: Simulation
## To Do:
- Run Code for Lasso
- Start with the presentation (do we have a template?)
    - Motivation: Why Lasso in Labor Economics ? -> Use Angrist 2022 as example (Cristian)
    - Theory (Explanation + Formulas + Graphs)
    - Simulation (Marcel)

## TO DO UNTIL (08.01)
- Cristian: Finish simulations for penalty selection
- Group: At least a draft of the presentation to discuss with the Professor
- Marcel: Upload codes and visualize the results

