
#' @export
merge_mutations_in_phase <- function(mutations, bam, tag = "", ID_column = "phasingID", min_base_quality = 20) {
	
	assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
       assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT", ID_column)))
    
    mutation_id <- paste(mutations$CHROM, mutations$POS, mutations$REF, mutations$ALT, sep = "_")
    mutations[,ID_column] <- ifelse(is.na(mutations[,ID_column]), mutation_id, mutations[,ID_column])

    IDs <- data.frame(mutations_id = mutation_id, phasing_id = mutations[,ID_column])
    IDs_list <- purrr::map(unique(IDs$phasing_id), ~ IDs[IDs$phasing_id == .x, "mutations_id"])

    read_names <- get_mutations_read_names(mutations = mutations, bam = bam,
        tag = tag, min_base_quality = min_base_quality)

    out <- purrr::map_dfr(IDs_list, function(.x){
    	## extract mutations belonging to the phase ID
    	read_names_to_merge <- read_names[.x]  

    	## get the ref and alt reads
    	alt <- unlist(purrr::map(read_names_to_merge, "alt"))	
    	ref <- unlist(purrr::map(read_names_to_merge, "ref"))
    	
        ## purify mismatches mapping only to one of the mutations in phase
    	alt <- alt[! alt %in% ref ]

    	## count how many reads support more than one mutation
    	multi_support <- length(unique(alt[duplicated(alt)]))

        ## put the counts together
    	data.frame(ref = length(unique(ref)),
    		alt = length(unique(alt)),
    		multi_support = multi_support)
    	})

    out <- data.frame(Phasing_id = unique(IDs$phasing_id),
     	out, stringsAsFactors = F)

    return(out)

}