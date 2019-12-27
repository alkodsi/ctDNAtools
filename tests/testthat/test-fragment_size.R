library(purrr)
library(assertthat)

data("mutations", package = "ctDNAtools")
data("targets", package = "ctDNAtools")
bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
bamT2 <- system.file("extdata", "T2.bam", package = "ctDNAtools")
bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
bamN2 <- system.file("extdata", "N2.bam", package = "ctDNAtools")
bamN3 <- system.file("extdata", "N3.bam", package = "ctDNAtools")


fs <- list(
  fs1 = get_fragment_size(bamN3),
  fs2 = get_fragment_size(bamN3,
    mutations = mutations, isProperPair = T,
    ignore_trimmed = F, different_strands = F
  ),
  fs3 = get_fragment_size(bamN3, max_size = 200, min_size = 100),
  fs4 = get_fragment_size(bamN3, tag = "ID1", targets = targets)
)

test_that("get_fragment_size works", {
  map(fs, ~ expect_is(.x, "data.frame"))

  map(fs, ~ expect_true(assertthat::has_name(.x, c("Sample", "ID", "start", "end", "size"))))

  map(fs, ~ expect_true(assertthat::noNA(.x$size)))

  map(fs, ~ expect_true(is.numeric(.x$size)))

  expect_true(assertthat::has_name(fs$fs2, "category"))
})

bfs <- list(
  bfs1 = bin_fragment_size(bamN3),
  bfs2 = bin_fragment_size(bamT1, mutations = mutations),
  bfs3 = bin_fragment_size(bamN3, normalized = TRUE, targets = targets),
  bfs4 = bin_fragment_size(bamN3, normalized = TRUE, custom_bins = c(50, 100, 150))
)

test_that("bin_fragment_size works", {
  map(bfs, ~ expect_is(.x, "data.frame"))

  map(bfs, ~ expect_true(has_name(.x, c("Breaks"))))

  map(bfs, ~ expect_true(noNA(.x[[2]])))

  map(bfs, ~ expect_true(is.numeric(.x[[2]])))

  expect_equal(sum(bfs$bfs4[[2]]), 1)

  expect_equal(ncol(bfs$bfs2), 4)
})


sfs <- list(
  sfs1 = summarize_fragment_size(bam = bamN3, regions = targets),
  sfs2 = summarize_fragment_size(bam = bamN3, regions = targets, summary_functions = list(sd = sd)),
  sfs3 = summarize_fragment_size(bam = bamN3, regions = targets, max_size = 200)
)

test_that("summarize_fragment_size works", {
  map(sfs, ~ expect_is(.x, "data.frame"))

  map(sfs, ~ expect_true(ncol(.x) > 1))

  map(sfs, ~ expect_true(nrow(.x) == nrow(targets)))

  expect_false(all(sfs$sfs1[[2]] == sfs$sfs3[[2]]))
})

af <- list(
  af1 = analyze_fragmentation(bamN3, targets, step_size = 5),
  af2 = analyze_fragmentation(bamN3, targets, step_size = 10)
)

test_that("analyze_fragmentation works", {
  map(af, ~ expect_is(.x, "data.frame"))

  map(af, ~ expect_true(has_name(
    .x,
    c("WPS", "WPS_adjusted", "n_fragment_ends", "n_fragment_ends_adjusted", "n_reads")
  )))

  expect_gt(nrow(af$af1), nrow(af$af2))
})

test_that("get_hist_bins works", {
  ghb <- get_hist_bins(c(0:10), from = 0, to = 10, by = 2, normalized = T)

  expect_is(ghb, "list")

  expect_true(has_name(ghb, c("counts", "breaks")))

  expect_is(ghb$counts, "numeric")

  expect_is(ghb$breaks, "character")

  expect_equal(sum(ghb$counts), 1)
})

test_that("get_mutations_fragment_size works", {
  mutfs <- get_mutations_fragment_size(bam = bamT1, mutations = mutations[1:3, ])

  expect_is(mutfs, "list")

  expect_length(mutfs, 3)

  map(mutfs, ~ expect_length(.x, 2))

  expect_is(unlist(mutfs), "integer")
})
