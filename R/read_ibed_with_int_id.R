#' Read ibed file to data.frame and add interaction id as new column
#'
#' A ibed file contains duplicates interaction because bait-to-bait interaction may be
#' doubled for both baits. Some bait-to-bait show significant at only one direction.
#' This function allows to detect bi-directional or uni-directional significant interactions,
#' and assign bi-directional with same interaction id.
#'
#' @param ibed ibed file from Chicago or data.frame with ibed header
#' ibed header requirement: bait_chr, bait_start, bait_end, bait_name, otherEnd_chr, otherEnd_start,
#' otherEnd_end, otherEnd_name
#' @param parse_b2b convert all bait-to-bait interactions to bi-directional ("double all"), or
#' or only keep bi-directional bait-to-bait interactions ("only bi-direction"). default = NULL
#' @param max_score whether convert bait-to-bait interactions to max score of both directions.
#' default = F. In case of parse_b2b="double all" and max_score=T, NA will fill score and N_reads in bait-to-bait
#' interactions in which only one direction is significant.
#'
#' @return data.frame with extra column called "int_id"
#' @export
#' @import dplyr tidyr
#'
#' @examples
#'
#' \donttest{
#' read_ibed_with_int_id(
#' system.file("extdata", "chicagoResults.ibed", package = "parseIbed"),
#' parse_b2b="only bi-direction"
#' )
#' }
#'
#' data(example)
#' \donttest{read_ibed_with_int_id(ibed, parse_b2b="double all")}
#'
read_ibed_with_int_id <- function(ibed, parse_b2b=NULL, max_score=F){

        if("data.frame" %in% class(ibed)){
                ibed=ibed
        }else if(file.exists(ibed)){
                ibed <- readr::read_delim(ibed, delim="\t")
        }else{
                stop("ibed can be either a ibed file or ibed data.frame")
        }

        if(!'score' %in% colnames(ibed)){
                ibed = ibed %>% mutate(score=5)
        }
        if(!'N_reads' %in% colnames(ibed)){
                ibed = ibed %>% mutate(N_reads=1)
        }
        nob2b=ibed %>% filter(otherEnd_name==".")
        b2b=ibed %>% filter(otherEnd_name!=".")

        b2b_Alldouble = bind_rows(
                b2b %>% select(-N_reads, -score),
                b2b %>% select(-N_reads, -score) %>%
                        dplyr::rename(
                                bait_chr=otherEnd_chr, bait_start=otherEnd_start, bait_end=otherEnd_end, bait_name=otherEnd_name,
                                otherEnd_chr=bait_chr, otherEnd_start=bait_start, otherEnd_end=bait_end, otherEnd_name=bait_name
                        )
        ) %>% distinct()

        b2b_Alldouble = bind_rows(
                b2b_Alldouble %>%
                        filter(bait_start > otherEnd_start) %>%
                        mutate(int_id = row_number()),
                b2b_Alldouble %>%
                        filter(bait_start > otherEnd_start) %>%
                        mutate(int_id = row_number()) %>%
                        dplyr::rename(
                                bait_chr=otherEnd_chr, bait_start=otherEnd_start, bait_end=otherEnd_end, bait_name=otherEnd_name,
                                otherEnd_chr=bait_chr, otherEnd_start=bait_start, otherEnd_end=bait_end, otherEnd_name=bait_name
                        )
        )

        b2b = left_join(
                b2b,
                b2b_Alldouble
        )

        nob2b = nob2b %>%
                mutate(int_id=row_number()+nrow(b2b_Alldouble)/2)


        if (is.null(parse_b2b)){
                if(max_score==T){
                        ibed = bind_rows(
                                left_join(
                                        b2b,
                                        b2b %>%
                                                group_by(int_id) %>%
                                                dplyr::slice(which.max(score)) %>%
                                                ungroup() %>%
                                                select(int_id, N_reads, score)
                                ),
                                nob2b
                        )
                }else{
                        ibed = bind_rows(b2b, nob2b)
                }
        }else if (parse_b2b == "double all"){
                if(max_score==T){
                        ibed = bind_rows(
                                left_join(
                                        b2b_Alldouble,
                                        b2b %>%
                                                group_by(int_id) %>%
                                                dplyr::slice(which.max(score)) %>%
                                                ungroup() %>%
                                                select(int_id, N_reads, score)
                                ),
                                nob2b
                        )
                }else{
                        ibed = bind_rows(
                                left_join(
                                        b2b_Alldouble,
                                        b2b
                                ),
                                nob2b
                        )
                }
        }else if (parse_b2b == "only bi-direction"){
                if(max_score==T){
                        ibed = bind_rows(
                                left_join(
                                        semi_join(
                                                b2b_Alldouble,
                                                b2b %>% select(int_id) %>%
                                                        count(int_id) %>% filter(n>1)
                                        ),
                                        b2b %>% group_by(int_id) %>%
                                                dplyr::slice(which.max(score)) %>%
                                                ungroup() %>%
                                                select(int_id, N_reads, score)
                                ),
                                nob2b
                        )
                }else{
                        ibed = bind_rows(
                                semi_join(
                                        b2b,
                                        b2b %>% select(int_id) %>%
                                                count(int_id) %>% filter(n>1)
                                ),
                                nob2b
                        )
                }
        }
        ibed
}
