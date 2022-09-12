if(requireNamespace("testthat", quietly = TRUE)) {
  library("predmod", quietly = TRUE)
  testthat::test_check("predmod")
}
