# MultiSKAT
MultiSKAT is an R-package focused at rare-variant analysis of continuous multiple phenotype data. 
This project contains the R-codes/functions (including an example dataset) to carry out the MultiSKAT tests.

Main functions:

MultiSKAT_NULL: Builds a MultiSKAT null model extracting the sufficient statistics from the phenotype and the covariate data under the 
                assumption of no association for independent samples. Options are available (is.fast) for a faster approximation to the null model and is set to 
                TRUE by default. 
            
MultiSKAT_NULL_Kins: Builds a MultiSKAT null model extracting the sufficient statistics from the phenotype and the covariate data under the 
                     assumption of no association for samples with relatedness. For details on the function arguments please refer to the 
                     package vignette.
                     
MultiSKAT: Performs the MultiSKAT test for a given phenotype kernel (Sigma_p) for related or unrelated samples. For details on the function
           arguments please refer to the package vignette.                     
           
minP: Performs the minP omnibus test for a given genotype kernel (SKAT or Burden) over a given list of phenotype kernels.

minPcom: Performs the minPcom omnibus test with the genotype kernels being SKAT and Burden and a given list of phenotype kernels. It is not
         advisable to use more than 5 phenotype kernels since the tail approximation can be unstable for copula-based methods.
          
