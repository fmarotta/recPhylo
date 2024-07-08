# NOTE: use testthat::snapshot_review() to inspect failures

skip_if_not_installed("vdiffr")
test_that("the plots with `branch_length` look fine", {
  ex <- example_recPhyloXML()
  vdiffr::expect_doppelganger(
    "scale = 3, padding = 3, y_shift = T, branch_length = T",
    RecPhyloLayout$new(ex,
      branch_length_scale = 3,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = T
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = T, branch_length = T",
    RecPhyloLayout$new(ex,
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = T
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 5, padding = 3, y_shift = T, branch_length = T",
    RecPhyloLayout$new(ex,
      branch_length_scale = 5,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = T
    )$plot()
  )
})

skip_if_not_installed("vdiffr")
test_that("the plots without `branch_length` look fine", {
  ex <- example_recPhyloXML()
  vdiffr::expect_doppelganger(
    "scale = 3, padding = 3, y_shift = T, branch_length = F",
    RecPhyloLayout$new(ex,
      branch_length_scale = 3,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = T, branch_length = F",
    RecPhyloLayout$new(ex,
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 5, padding = 3, y_shift = T, branch_length = F",
    RecPhyloLayout$new(ex,
      branch_length_scale = 5,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
})

skip_if_not_installed("vdiffr")
test_that("the plots without `y_shift` look fine", {
  ex <- example_recPhyloXML()
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = F, branch_length = T",
    RecPhyloLayout$new(ex,
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = F,
      use_branch_length = T
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = F, branch_length = F",
    RecPhyloLayout$new(ex,
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = F,
      use_branch_length = F
    )$plot()
  )
})
