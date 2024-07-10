test_that("the mrca() function works", {
  ex <- example_recPhyloXML()
  ex_clade <- ex$spTree$clade
  expect_equal(find_mrca_clade(ex_clade, c("A"))$name, "A")
  expect_equal(find_mrca_clade(ex_clade, c("D", "C"))$name, "CD")
  expect_equal(find_mrca_clade(ex_clade, c("D", "B"))$name, "BCD")
  expect_equal(find_mrca_clade(ex_clade, c("D", "A"))$name, "ABCD")
  expect_true(is.na(find_mrca_clade(ex$spList, c("D", "E"))))
})
