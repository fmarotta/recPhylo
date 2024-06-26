test_that("the mrca() function works", {
  ex <- RecPhylo$new(recphylo_example("example_1.recphyloxml"), use_branch_length = 10)
  expect_equal(mrca(ex$spList, c("A")), "A")
  expect_equal(mrca(ex$spList, c("D", "C")), "CD")
  expect_equal(mrca(ex$spList, c("D", "B")), "BCD")
  expect_equal(mrca(ex$spList, c("D", "A")), "ABCD")
  expect_true(is.na(mrca(ex$spList, c("D", "E"))))
})
