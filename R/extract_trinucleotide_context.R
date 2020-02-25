#' Extracts the trinucleotide context for a set of mutations

#' @param mutations A data frame having the mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param reference the reference genome in BSgenome format
#' @param destrand logical, whether to destrand mutations
#' @return A data frame with two columns having the substitutions and the trinucleotide context
#' @export
#' @examples
#' \donttest{
#' data("mutations", package = "ctDNAtools")
#' ## Use human reference genome from BSgenome.Hsapiens.UCSC.hg19 library
#' suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#'
#' ## with destranding
#' extract_trinucleotide_context(mutations, BSgenome.Hsapiens.UCSC.hg19)
#'
#' ## without destranding
#' extract_trinucleotide_context(mutations, BSgenome.Hsapiens.UCSC.hg19,
#'   destrand = FALSE
#' )
#' }
extract_trinucleotide_context <- function(mutations, reference, destrand = TRUE) {
  assertthat::assert_that(
    is.data.frame(mutations),
    assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT"))
  )

  assertthat::assert_that(all(mutations$CHROM %in% GenomeInfoDb::seqnames(reference)),
    msg = "Chromosomes in mutations and/or targets don't match the specified reference"
  )
  
  assertthat::assert_that(class(reference) == "BSgenome")

  assertthat::assert_that(all(nchar(mutations$REF) == 1), all(nchar(mutations$ALT) == 1),
    msg = "Only SNVs are supported"
  )

  assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
    all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
    msg = "REF and ALT in mutations should be characters having basepairs"
  )

  assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

  assertthat::assert_that(is.logical(destrand), length(destrand) == 1)

  context <- BSgenome::getSeq(
    reference,
    GenomicRanges::GRanges(
      mutations$CHROM,
      IRanges::IRanges(mutations$POS - 1, mutations$POS + 1)
    )
  )


  if (destrand) {
    context <- ifelse(mutations$REF %in% c("C", "T"), context, Biostrings::reverseComplement(context))

    ref <- Biostrings::DNAStringSet(mutations$REF)

    ref <- ifelse(mutations$REF %in% c("C", "T"), ref, Biostrings::reverseComplement(ref))

    alt <- Biostrings::DNAStringSet(mutations$ALT)

    alt <- ifelse(mutations$REF %in% c("C", "T"), alt, Biostrings::reverseComplement(alt))

    substitutions <- factor(paste0(ref, alt), levels = c("CA", "CG", "CT", "TA", "TC", "TG"))
  } else {
    context <- as.character(context)

    substitutions <- factor(paste0(mutations$REF, mutations$ALT),
      levels = c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG")
    )
  }

  combinations <- expand.grid(c("C", "G", "A", "T"), c("C", "G", "A", "T"))

  context <- factor(paste0(substring(context, 1, 1), "_", substring(context, 3, 3)),
    levels = sort(paste0(combinations[, 1], "_", combinations[, 2]))
  )

  out <- data.frame(substitution = substitutions, context = context)

  return(out)
}
