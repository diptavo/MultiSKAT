# MultiSKAT `v(1.0)`
MultiSKAT is an R-package focused at rare-variant analysis of continuous multiple phenotype data. 
This project contains the R-codes/functions (including an example dataset) to carry out the MultiSKAT tests. Please note that this package is still under some developement. Any questions or bug-reports should be mailed to **diptavo@umich.edu**

# Citation

If you use MultiSKAT R package for data analysis, please consider citing our manuscript on biorXiv.
 https://t.co/pxA9W8CZ0S
 
[Diptavo Dutta, Laura Scott, Michael Boehnke, Seunggeun Lee. *Multi-SKAT: General framework to test multiple phenotype associations of rare variants*]( https://t.co/pxA9W8CZ0S)


# Installation

MultiSKAT R-package can be installed from the downloaded gzipped tarball as
```R CMD INSTALL MultiSKAT_1.0.tar.gz ```
with the following packages pre-installed: SKAT, copula, nlme

# Main functions:

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
          

MultiSKAT provides a generalized framework for performing rare-variant based tests of associations. Many published methods (e.g. GAMuT, MSKAT, MAAUSS, MK-FM) can be viewed as a special case of MultiSKAT with certain specific choice of kernels. In addition, MultiSKAT includes a set of omnibus tests that can aggregate results over a multiplicity of kernels, producing robust results.
