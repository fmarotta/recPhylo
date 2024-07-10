test_that("Phylogenies have idx field", {
  # One of our business logics is that lists of phylogeny elements (i.e.
  # phyloXML and recPhyloXML$recGeneTrees) can be identified by a $idx
  # field in each phylogeny, representing the index.
  ex_phylo <- example_phyloXML()
  expect_equal(ex_phylo[[1]]$idx, 1)

  ex_recPhylo <- example_recPhyloXML()
  expect_null(ex_recPhylo$spTree$idx)
  expect_equal(ex_recPhylo$recGeneTrees[[1]]$idx, 1)
})

test_that("Phylogenies have unique names", {
  # Another business logic is that the <name> elements' content must be unique
  # across phyloXML and recGeneTree.
  ex_phylo <- example_phyloXML()
  all_names <- traverse_clades(ex_phylo[[1]]$clade, function(x) x$name)
  expect_equal(all_names, unique(all_names))

  ex_recPhylo <- example_recPhyloXML()
  all_names <- traverse_clades(ex_recPhylo$recGeneTrees[[1]]$clade, function(x) x$name)
  expect_equal(all_names, unique(all_names))
})
