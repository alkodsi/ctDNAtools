#' Collapses mutations in phase into one event
#'
#' Given a mutations data frame and a bam file, this function collapses mutations in phase identified by the ID_column into one event.
#' While doing that, it ignores the reads that support both the reference and alternative alleles for different mutations in phase.
#' When use_unique_molecules is TRUE, it counts the unique molecules only once, i.e. if two mutations are close to each others, the
#' counts for the second mutation include only the reads that do not map to the first one.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param ID_column The name of the column in mutations data.frame that has the IDs for mutations in phase.
#' NA values will be filled automatically by unique mutation identifiers.
#' @param min_base_quality minimum base quality for a read to be counted
#' @return A data frame that has the ref and alt counts for the mutations/events, as well as, 
#' the count of reads that support multiple mutations in phase, and the total number of reads.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%

merge_mutations_in_phase <- function(mutations, bam, tag = "", ID_column = "phasingID", min_base_quality = 20, use_unique_molecules = T) {
	
	assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
       assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT", ID_column)))
    
    mutation_id <- paste0(mutations$CHROM, ":",mutations$POS, "_", mutations$REF, "_", mutations$ALT)
    mutations[,ID_column] <- as.character(ifelse(is.na(mutations[,ID_column]), mutation_id, mutations[,ID_column]))

    IDs <- data.frame(mutations_id = mutation_id, phasing_id = mutations[,ID_column], stringsAsFactors = F)
    IDs_list <- purrr::map(unique(IDs$phasing_id), ~ IDs[IDs$phasing_id == .x, "mutations_id"])

    read_names <- get_mutations_read_names(mutations = mutations, bam = bam,
        tag = tag, min_base_quality = min_base_quality)

    out <- purrr::map(IDs_list, function(.x) {
    	## extract mutations belonging to the phase ID
    	read_names_to_merge <- read_names[.x]  

    	## get the ref and alt reads
    	alt <- unlist(purrr::map(read_names_to_merge, "alt"))	
    	ref <- unlist(purrr::map(read_names_to_merge, "ref"))
    	
        ## number of reads that map to more than one mutation
        all_reads <- c(ref, alt)
        n_reads_multi_mutation <- length(unique(all_reads[duplicated(all_reads)]))

        ## purify mismatches mapping only to one of the mutations in phase
    	alt <- alt[! alt %in% ref ]

    	## count how many reads support more than one mutation
    	multi_support <- length(unique(alt[duplicated(alt)]))

        ## put the counts together
    	df <- data.frame(ref = length(unique(ref)),
    		alt = length(unique(alt)),
            n_reads_multi_mutation = n_reads_multi_mutation,
            all_reads = length(unique(all_reads)),
    		multi_support = multi_support)
    	ref <- unique(ref)
    	alt <- unique(alt)
        list(df = df, ref = ref, alt = alt)    
    })
   
    df <- purrr::map_dfr(out,"df")

    purification_prob <- sum(df$n_reads_multi_mutation)/sum(df$all_reads) 
    
    out <- data.frame(Phasing_id = unique(IDs$phasing_id),
  	    df, stringsAsFactors = F)

    return(list(out = out, purification_prob = purification_prob))

}