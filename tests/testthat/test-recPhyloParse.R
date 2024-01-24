# NOTE: use testthat::snapshot_review() to inspect failures

skip_if_not_installed("vdiffr")
test_that("the plots with `branch_length` look fine", {
  ex <- RecPhylo$new(recphylo_example("example_1.recphyloxml"))
  vdiffr::expect_doppelganger(
    "scale = 3, padding = 3, y_shift = T, branch_length = T",
    ex$redraw(
      branch_length_scale = 3,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = T, branch_length = T",
    ex$redraw(
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 5, padding = 3, y_shift = T, branch_length = T",
    ex$redraw(
      branch_length_scale = 5,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
})

skip_if_not_installed("vdiffr")
test_that("the plots without `branch_length` look fine", {
  ex <- RecPhylo$new(recphylo_example("example_1.recphyloxml"))
  vdiffr::expect_doppelganger(
    "scale = 3, padding = 3, y_shift = T, branch_length = F",
    ex$redraw(
      branch_length_scale = 3,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 4, padding = 3, y_shift = T, branch_length = F",
    ex$redraw(
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "scale = 5, padding = 3, y_shift = T, branch_length = F",
    ex$redraw(
      branch_length_scale = 5,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = F
    )$plot()
  )
})

#
#   # y_shift = F
#   ex$redraw(branch_length_scale = 3, x_padding = 3, use_y_shift = F, use_branch_length = "branch_length")
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 4, x_padding = 3, use_y_shift = F, use_branch_length = "branch_length")
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 5, x_padding = 3, use_y_shift = F, use_branch_length = "branch_length")
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 6, x_padding = 3, use_y_shift = F, use_branch_length = "branch_length")
#   plot(ex)
#
#
#
#   # withOUT branch_length
#   ex$redraw(branch_length_scale = 3, x_padding = 3, use_y_shift = T, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 4, x_padding = 3, use_y_shift = T, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 5, x_padding = 3, use_y_shift = T, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 6, x_padding = 3, use_y_shift = T, use_branch_length = F)
#   plot(ex)
#
#   # y_shift = F
#   ex$redraw(branch_length_scale = 3, x_padding = 3, use_y_shift = F, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 4, x_padding = 3, use_y_shift = F, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 5, x_padding = 3, use_y_shift = F, use_branch_length = F)
#   plot(ex)
#
#   ex$redraw(branch_length_scale = 6, x_padding = 3, use_y_shift = F, use_branch_length = F)
#   plot(ex)
