context("test-missing-genotype-data")



test_that("multiplication works", {

  data(MultiSKAT.example.data)
  attach(MultiSKAT.example.data)
  
  ## Introduce NAs in Genotypes
  Genotypes[1,1] <- NA
  
  expect_warning(MultiSKAT(MultiSKAT_NULL(Phenotypes, Cov), Genotypes, Sigma_p = cov(Phenotypes), verbose = FALSE), 
                 "The missing genotype rate is 0.000004. Imputation is applied.")

  detach(MultiSKAT.example.data)
  
})
