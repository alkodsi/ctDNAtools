library(purrr)
library(assertthat)
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

data("mutations", package = "ctDNAtools")
data("targets", package = "ctDNAtools")
bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
bamT2 <- system.file("extdata", "T2.bam", package = "ctDNAtools")
bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
bamN2 <- system.file("extdata", "N2.bam", package = "ctDNAtools")
bamN3 <- system.file("extdata", "N3.bam", package = "ctDNAtools")

tests1 <- list(
  test1 = test_ctDNA(
    mutations = mutations, bam = bamT1,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, informative_reads_threshold = 100
  ),
  test2 = test_ctDNA(
    mutations = mutations, bam = bamT1, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, ID_column = "PHASING"
  ),
  test3 = test_ctDNA(
    mutations = mutations, bam = bamT1, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, ID_column = "PHASING", black_list = "chr14_106327474_C_G"
  ),
  test4 = test_ctDNA(
    mutations = mutations, bam = bamT1, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19, substitution_specific = FALSE,
    n_simulation = 100, ID_column = "PHASING", black_list = "chr14_106327474"
  )
)

tests2 <- list(
  test1 = test_ctDNA(
    mutations = mutations, bam = bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, informative_reads_threshold = 100
  ),
  test2 = test_ctDNA(
    mutations = mutations, bam = bamT2, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, ID_column = "PHASING"
  ),
  test3 = test_ctDNA(
    mutations = mutations, bam = bamT2, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    n_simulation = 100, ID_column = "PHASING", black_list = "chr14_106327474_C_G"
  ),
  test4 = test_ctDNA(
    mutations = mutations, bam = bamT2, informative_reads_threshold = 100,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19, substitution_specific = FALSE,
    n_simulation = 100, ID_column = "PHASING", black_list = "chr14_106327474"
  ),
  test5 = test_ctDNA(
    mutations = mutations, bam = bamT2, tag = "ID1", targets = targets,
    n_simulation = 100, ID_column = "PHASING", bam_list = c(bamN1, bamN2, bamN3),
    bam_list_tags = c("ID1", "", ""),
    informative_reads_threshold = 100, reference = BSgenome.Hsapiens.UCSC.hg19
  )
)

test_that("test_ctDNA works", {
  map(tests1, ~ expect_is(.x, "data.frame"))

  map(tests2, ~ expect_is(.x, "data.frame"))

  map(tests1, ~ expect_true(has_name(.x, c("sample", "pvalue", "decision", "informative_reads", "background_rate"))))

  map(tests1, ~ expect_equal(.x$decision, factor("positive", levels = c("positive", "negative", "undetermined"))))

  map(tests2, ~ expect_equal(.x$decision, factor("negative", levels = c("positive", "negative", "undetermined"))))

  expect_warning(test_ctDNA(
    mutations = mutations[mutations$CHROM != "chr14", ], bam = bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19
  ))

  expect_warning(test_ctDNA(
    mutations = mutations[1, ], bam = bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    black_list = "chr14_106327474", substitution_specific = FALSE
  ))

  expect_error(test_ctDNA(
    mutations = mutations, bam = bamT2,
    targets = targets[targets$chr != "chr14", ], reference = BSgenome.Hsapiens.UCSC.hg19
  ))
})

bgr <- list(
  bgr1 = get_background_rate(bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    min_base_quality = 20
  ),
  bgr2 = get_background_rate(bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    min_base_quality = 30
  ),
  bgr3 = get_background_rate(bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    black_list = "chr14_106327474_C_G"
  ),
  bgr4 = get_background_rate(bamT2,
    targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
    black_list = "chr14_106327474", substitution_specific = F
  )
)

test_that("get_background_rate works", {
  map(bgr, ~ expect_is(.x, "list"))

  map(bgr, ~ expect_true(has_name(.x, c("rate"))))

  expect_gt(bgr$bgr1$rate, bgr$bgr2$rate)
})

bams <- c(bamN1, bamN2, bamN3)

bgp <- list(
  bgp1 = create_background_panel(bams, bam_list_tags = c("", "", "ID1"), targets, BSgenome.Hsapiens.UCSC.hg19),
  bgp2 = create_background_panel(bams, targets, BSgenome.Hsapiens.UCSC.hg19,
    min_base_quality = 20, substitution_specific = FALSE
  )
)


test_that("create_background_panel works", {
  map(bgp, ~ expect_is(.x, "list"))

  map(bgp, ~ expect_true(has_name(.x, c("depth", "alt", "vaf"))))

  map(bgp, ~ expect_true(nrow(.x$depth) == nrow(.x$vaf)))

  map(bgp, ~ expect_true(nrow(.x$depth) == nrow(.x$alt)))

  expect_true(all(map_dbl(strsplit(bgp$bgp1$depth$Locus, "_"), length) == 4))

  expect_true(all(map_dbl(strsplit(bgp$bgp2$depth$Locus, "_"), length) == 2))

  expect_is(create_black_list(bgp$bgp1), "character")

  expect_is(create_black_list(bgp$bgp2), "character")

  expect_is(create_black_list(bgp$bgp1, mean_vaf_quantile = 0.5, min_samples_n_reads = 1, n_reads = 3), "character")
})

phase <- merge_mutations_in_phase(mutations, bamT1, ID_column = "PHASING")

test_that("merge_mutations_in_phase works", {
  expect_is(phase, "list")

  expect_true(has_name(phase, c("out", "purification_prob", "multi_support", "informative_reads")))

  expect_is(phase$out, "data.frame")

  expect_true(has_name(phase$out, c("Phasing_id", "ref", "alt", "n_reads_multi_mutation", "all_reads", "multi_support")))

  expect_equal(nrow(phase$out), 7)
})


fm <- list(
  fm1 = filter_mutations(mutations, black_list = "chr14_106327474_C_G"),
  fm2 = filter_mutations(mutations, black_list = "chr14_106327474", substitution_specific = FALSE),
  fm3 = filter_mutations(mutations, bams = bams)
)

test_that("filter_mutations works", {
  map(fm, ~ expect_is(.x, "data.frame"))
})


test_that("extract_trinucleotide_context works", {
  tri <- list(
    tri1 = extract_trinucleotide_context(mutations, BSgenome.Hsapiens.UCSC.hg19),
    tri2 = extract_trinucleotide_context(mutations, BSgenome.Hsapiens.UCSC.hg19, destrand = FALSE)
  )

  map(tri, ~ expect_is(.x, "data.frame"))

  map(tri, ~ expect_true(ncol(.x) == 2))

  map(tri, ~ expect_true(has_name(.x, c("substitution", "context"))))
})

test_that("get_mutations_read_counts works", {
  rc <- get_mutations_read_counts(mutations, bamT2)

  expect_is(rc, "list")

  expect_length(rc, 2)

  expect_true(has_name(rc, c("ref", "alt")))

  expect_true(noNA(rc$ref))

  expect_true(noNA(rc$alt))

  expect_true(is.numeric(rc$ref))

  expect_true(is.numeric(rc$alt))
})

test_that("get_mutations_read_names works", {
  rn <- get_mutations_read_names(mutations = mutations, bam = bamT2)

  expect_is(rn, "list")

  expect_length(rn, nrow(mutations))

  expect_true(all(map_dbl(rn, length) == 2))

  expect_is(unlist(rn), "character")
})


test_that("positivity_test works", {
  expect_equal(simulator(depths = 5, rate = bgr$bgr1, alt_reads = 1, seed = 123), 0)

  expect_equal(simulator(depths = 5, rate = bgr$bgr1, alt_reads = 0, seed = 123), 1)

  expect_equal(positivity_test(depths = 5, rate = bgr$bgr1, alt_reads = 0, n_simulations = 100), 1)

  expect_equal(positivity_test(depths = 5, rate = bgr$bgr1, alt_reads = 1, n_simulations = 100), 1 / 101)
})

bamNoSM <- system.file("extdata", "ex1.bam", package = "Rsamtools")

test_that("bamutils works", {
  expect_true(verify_tag(bamT1, "ID1"))

  expect_false(verify_tag(bamT1, "ID11"))

  expect_equal(get_bam_SM(bamT1), "T1")

  expect_equal(get_bam_SM(bamNoSM), "ex1.bam")

  expect_true(all(paste0("chr", c(1:22)) %in% get_bam_chr(bamT1)))
})

test_that("vcf_to_mutations_df works", {
  
  vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  df <- vcf_to_mutations_df(vcf, sample_name = "HG00096")
  
  expect_is(df, "data.frame")
})
