# Centriole-elimination (Pierron et al., 2023, under review)

## ImageJ macro for dual colour quantifications
The script is under src/ij_scripts/

## Estimation of the parameters for the signal intensity decay of two markers in centrioles
### Data
The data consist of signal intensities of centriole markers (Tubulin, SAS-4) for different genotypes (WT, sas-1ts) at different phases (EP, LP, Diplo), with EP < LP < Diplo. 
The intensities are obtained by computing the mean of the ROI, from which the background signal has been subtracted.
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
Effectively, we ran four chains in parallel for 8000 iterations each. The code for the model and the sampling is available at <url>.

### References
- Riddell, A., Hartikainen, A., & Carter, M. (2021). PyStan (3.0.0). https://pypi.org/project/pystan
- Stan Development Team. 2023. Stan Modeling Language Users Guide and Reference Manual, Version 2.33. https://mc-stan.org
