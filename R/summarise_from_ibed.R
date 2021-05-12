#' Summarize ibed
#'
#' summary statistics of interactions in the context of bait, otherEnd and interaction
#'
#' @param ibed ibed file from Chicago or data.frame with ibed header
#' ibed header requirement: bait_chr, bait_start, bait_end, bait_name, otherEnd_chr, otherEnd_start,
#' otherEnd_end, otherEnd_name, int_id
#' @param baitmap bait design file corresponding to interaction being called in Chicago or data.frame
#' containing column "chr", "start", "end" and "name"|"anno"
#'
#'
#' @return data.frame of summary stats
#' @export
#'
#' @import dplyr tidyr
#' @importFrom stats median
#'
#' @examples
#' data(example)
#' \donttest{summarise_from_ibed(ibed, baitmap)}
#'
#'
summarise_from_ibed <- function(ibed, baitmap){
        # depends on baitmap to determine gene/transcript annotation
        if ("data.frame" %in% class(baitmap)){
                baitmap_df = baitmap
                colnames(baitmap_df)[grepl("chr",tolower(colnames(baitmap_df)))] = "bait_chr"
                colnames(baitmap_df)[grepl("start",tolower(colnames(baitmap_df)))] = "bait_start"
                colnames(baitmap_df)[grepl("end",tolower(colnames(baitmap_df)))] = "bait_end"
                colnames(baitmap_df)[grepl("name|anno",tolower(colnames(baitmap_df)))] = "bait_name"
                baitmap_df = baitmap_df %>% select(bait_chr, bait_start, bait_end, bait_name)
        }else{
                baitmap_df = readr::read_delim(baitmap, delim="\t",
                                               col_names=c("bait_chr","bait_start","bait_end","frag_id","bait_name"),
                                               col_types = "ciiic"
                ) %>%
                select(bait_chr, bait_start, bait_end, bait_name)
        }


        if("data.frame" %in% class(ibed)){
                file="ibed"
        }else if(file.exists(ibed)){
                file = ibed
        }

        if("data.frame" %in% class(ibed)){
                ibed=ibed
        }else if(file.exists(ibed)){
                ibed = read_ibed_with_int_id(ibed)
        }else{
                stop("ibed can be either a ibed file or ibed data.frame")
        }

        if (!"int_id" %in% colnames(ibed)){
                ibed = read_ibed_with_int_id(ibed)
        }

        # all bait fragment number in design
        bait_desig_n = baitmap_df %>% count() %>% pull(n)

        # all transcript number in design
        transcript_design_n = baitmap_df %>% select(bait_name) %>%
                tidyr::separate_rows(bait_name, sep="\\|") %>%
                distinct() %>% filter(!is.na(bait_name)) %>%
                count() %>% pull(n)

        # total bait number in interaction
        bait_n = ibed %>% distinct(bait_chr, bait_start, bait_end) %>% count() %>% pull()

        # total bait ratio in interaction
        bait_ratio = bait_n/bait_desig_n

        # total interaction number
        int_n = ibed %>% distinct(int_id) %>% count() %>% pull(n)

        # trans interaction number
        int_trans_n = ibed %>% filter(bait_chr!=otherEnd_chr) %>% distinct(int_id) %>% count() %>% pull(n)

        # trans interaction ratio
        int_trans_ratio = int_trans_n/int_n

        # bait2bait interaction number
        int_b2b_n = ibed %>% filter(otherEnd_name!=".") %>% distinct(int_id) %>% count() %>% pull(n)

        # bait2bait interaction ratio
        int_b2b_ratio = int_b2b_n/int_n

        # bait2bait cis interaction number
        int_b2b_cis_n = ibed %>%
                filter(otherEnd_name!=".", bait_chr==otherEnd_chr) %>%
                distinct(int_id) %>% count() %>% pull(n)

        # median interaction distance (including bait2bait)
        median_dist = ibed %>% filter(bait_chr==otherEnd_chr) %>%
                mutate(dist=abs((bait_start+bait_end)/2 - (otherEnd_start+otherEnd_end)/2)) %>%
                distinct(int_id, dist) %>%
                summarise(dist=median(dist)) %>% pull(dist)

        # median interaction number per bait
        int_n_per_bait = ibed %>% group_by(bait_chr, bait_start, bait_end) %>%
                summarise(int_n = n_distinct(int_id)) %>%
                ungroup() %>%
                summarise(int_n=median(int_n)) %>% pull(int_n)

        # median cis interaction number per bait
        int_cis_n_per_bait = ibed %>%
                filter(bait_chr==otherEnd_chr) %>%
                group_by(bait_chr, bait_start, bait_end) %>%
                summarise(int_n = n_distinct(int_id)) %>%
                ungroup() %>%
                summarise(int_n=median(int_n)) %>% pull(int_n)

        # total other end number (exclude bait)
        oe_n = ibed %>% filter(otherEnd_name==".") %>%
                distinct(otherEnd_chr, otherEnd_start, otherEnd_end) %>%
                count() %>% pull(n)

        # median bait number per other end
        bait_n_per_oe = ibed %>% filter(otherEnd_name==".") %>%
                mutate(bait=glue::glue("{bait_chr}:{bait_start}-{bait_end}")) %>%
                group_by(otherEnd_chr, otherEnd_start, otherEnd_end) %>%
                summarise(bait_n = n_distinct(bait)) %>%
                ungroup() %>%
                summarise(bait_n=median(bait_n)) %>% pull(bait_n)

        # other end number with only 1 bait
        oe_n_1bait_per_oe = ibed %>% filter(otherEnd_name==".") %>%
                mutate(bait=glue::glue("{bait_chr}:{bait_start}-{bait_end}")) %>%
                group_by(otherEnd_chr, otherEnd_start, otherEnd_end) %>%
                summarise(bait_n = n_distinct(bait)) %>%
                ungroup() %>% filter(bait_n==1) %>%
                count() %>% pull(n)

        # other end ratio with only 1 bait
        oe_ratio_1bait_per_oe = oe_n_1bait_per_oe/oe_n

        # other end number with more than 4 baits
        oe_n_more4bait_per_oe = ibed %>% filter(otherEnd_name==".") %>%
                mutate(bait=glue::glue("{bait_chr}:{bait_start}-{bait_end}")) %>%
                group_by(otherEnd_chr, otherEnd_start, otherEnd_end) %>%
                summarise(bait_n = n_distinct(bait)) %>%
                ungroup() %>% filter(bait_n >= 4) %>%
                count() %>% pull(n)

        # other end ratio with more than 4 baits
        oe_ratio_more4bait_per_oe = oe_n_more4bait_per_oe/oe_n


        tibble(
                file=file,
                bait_desig_n=bait_desig_n,
                transcript_design_n=transcript_design_n,
                bait_n=bait_n,
                bait_ratio=bait_ratio,
                int_n=int_n,
                int_trans_n=int_trans_n,
                int_trans_ratio=int_trans_ratio,
                int_b2b_n=int_b2b_n,
                int_b2b_ratio=int_b2b_ratio,
                int_b2b_cis_n=int_b2b_cis_n,
                median_dist=median_dist,
                int_n_per_bait=int_n_per_bait,
                int_cis_n_per_bait=int_cis_n_per_bait,
                oe_n=oe_n,
                bait_n_per_oe=bait_n_per_oe,
                oe_n_1bait_per_oe=oe_n_1bait_per_oe,
                oe_ratio_1bait_per_oe=oe_ratio_1bait_per_oe,
                oe_n_more4bait_per_oe=oe_n_more4bait_per_oe,
                oe_ratio_more4bait_per_oe=oe_ratio_more4bait_per_oe
        )
}
