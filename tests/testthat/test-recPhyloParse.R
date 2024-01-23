# NOTE: use testthat::snapshot_review() to inspect failures

skip_if_not_installed("vdiffr")
test_that("the plots with `branch_length` look fine", {
  ex <- RecPhylo$new(recphylo_example("example_1.recphyloxml"))
  vdiffr::expect_doppelganger(
    "Gibberish",
    ex$redraw(
      branch_length_scale = 3,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "A bit better",
    ex$redraw(
      branch_length_scale = 4,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
  vdiffr::expect_doppelganger(
    "OK",
    ex$redraw(
      branch_length_scale = 5,
      x_padding = 3,
      use_y_shift = T,
      use_branch_length = "branch_length"
    )$plot()
  )
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
})
