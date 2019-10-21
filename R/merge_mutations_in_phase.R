#' @export

merge_mutations_in_phase <- function(mutations, bam, tag = "", ID_column = "phasingID", min_base_quality = 20, use_unique_molecules = T) {
	
	assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
       assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT", ID_column)))
    
    mutation_id <- paste0(mutations$CHROM, ":",mutations$POS, "_", mutations$REF, "_", mutations$ALT)
    mutations[,ID_column] <- as.character(ifelse(is.na(mutations[,ID_column]), mutation_id, mutations[,ID_column]))

    IDs <- data.frame(mutations_id = mutation_id, phasing_id = mutations[,ID_column], stringsAsFactors = F)
    IDs_list <- purrr::map(unique(IDs$phasing_id), ~ IDs[IDs$phasing_id == .x, "mutations_id"])

    read_names <- get_mutations_read_names(mutations = mutations, bam = bam,
        tag = tag, min_base_quality = min_base_quality)

    out <- purrr::map(IDs_list, function(.x){
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

    if(use_unique_molecules) {
        ref <- purrr::map(out, "ref")
    
        alt <- purrr::map(out, "alt")
       
        ref_df <- data.frame(ID = factor(unlist(map2(unique(IDs$phasing_id),ref,
    	    ~rep(.x, each = length(.y)))), levels = unique(IDs$phasing_id)) , ref = unlist(ref)) %>%
            dplyr::filter(!duplicated(.data$ref)) %>%
            dplyr::count(.data$ID, .drop = F)

        alt_df <- data.frame(ID = factor(unlist(map2(unique(IDs$phasing_id),alt,
    	    ~rep(.x, each = length(.y)))), levels = unique(IDs$phasing_id)) , alt = unlist(alt)) %>%
            dplyr::filter(!duplicated(.data$alt)) %>%
            dplyr::count(.data$ID, .drop = F)

        df <- df %>%
            dplyr::mutate(ref = ref_df$n, alt = alt_df$n)
    }
    
    out <- data.frame(Phasing_id = unique(IDs$phasing_id),
  	    df, stringsAsFactors = F)

    return(out)

}