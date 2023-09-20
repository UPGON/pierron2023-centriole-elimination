# Centriole-elimination (Pierron et al., 2023, under review)

## ImageJ macro for dual colour quantifications
The script is under src/ij_scripts/

## Estimation of the parameters for the signal intensity decay of two markers in centrioles
### Data
The data consist of signal intensities of centriolar proteins (Tubulin, SAS-4) for different genotypes (WT, sas-1ts) at different phases (EP, LP, Diplo), with EP < LP < Diplo. 
The intensities are obtained by computing the mean of the ROI, from which the background signal has been subtracted.
Furthermore, all intensities of each phase for each genotype-protein has been scaled by the mean of the EP, assumed to be full initial intensity.

### Model
We modelled the intensities for each comparison (phase_{ref} vs phase_{next}) as a Poisson regression, with a log link function. 
We encoded the phase indicator as a contrast binary vector that is set to zero for the reference phase and set to one for the next phase. 
Thus, the intensity is $exp^{\alpha}$ for the reference phase and $exp^{\alpha+\beta}$ for the next phase. 
We imposed two normal priors N(0, 10) over $\alpha$ and $\beta$.

$$ Intensity \sim Poisson(\lambda) $$

$$ log(\lambda)=\alpha+\beta X $$

$$ \alpha \sim Normal(0,10) $$

$$ \beta \sim Normal(0,10) $$

### Sampling and Estimation
We sampled the posterior and estimated the parameters alpha and beta using the Stan program (Stan Development Team, 2023) via its Python binding PyStan (Riddell, Hartikainen and Carter, 2021). 
Effectively, we ran four chains in parallel for 8000 iterations each. The code for the model and the sampling is available at in the GitHub repo (https://github.com/UPGON/pierron2023-centriole-elimination) under samples/.

### Interpretation of the coefficients and their uncertainty
The coefficient of decay $\beta_1$, that is the decrease between the two phases, is summarised with the 95 % credible interval, hereafter referred to as 95 %-C.I.
Tubulin decays earlier in the WT than in the sas-1ts (WT: $\beta_1$ = -0.33, 95 %-C.I. = \[-0.422, -0.241\]; sas-1ts: $\beta_1$ = -1.08, 95 %-C.I. = \[-1.154, -1.011\]).
SAS-4 decays likewise, yet with slightly less confidence (WT: $\beta_1 = -0.29$, 95 %-C.I. = \[-0.380, -0.199\]; sas-1ts: $\beta_1$ = -0.42, 95 %-C.I. = \[-0.480, -0.365\]).


### References
- Riddell, A., Hartikainen, A., & Carter, M. (2021). PyStan (3.0.0). https://pypi.org/project/pystan
- Stan Development Team. 2023. Stan Modeling Language Users Guide and Reference Manual, Version 2.33. https://mc-stan.org
