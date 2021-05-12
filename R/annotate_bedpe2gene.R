#' annotate_bedpe2gene
#'
#' annotate bedpe loop with one end as gene promoter
#'
#' @param bedpe bedpe data.frame
#' @param prom_bed prom_bed data.frame
#' @param gene_side NULL or bedpe column name match (eg."bait");
#' if NULL, search genes on both ends of bedpe.
#' @param bedpe_chr_a character
#' @param bedpe_start_a character
#' @param bedpe_end_a character
#' @param bedpe_chr_b character
#' @param bedpe_start_b character
#' @param bedpe_end_b character
#' @param prom_chr character
#' @param prom_start character
#' @param prom_end character
#' @param prom_anno character
#'
#' @return data.frame with a and b end annotated to pro_anno
#' @export
#'
#' @import dplyr tidyr DFbedtools
#' @examples
#' # annotate bedpe on both sides
#'
#' \donttest{data(bedpe)
#' data(prom_bed)
#' annotate_bedpe2gene(bedpe, prom_bed)
#' }
#'
#' # annotate ibed file on one side
#'
#' \donttest{data(ibed)
#' data(prom_bed)
#' annotate_bedpe2gene(
#' ibed %>% mutate_at(vars(ends_with("_start")), function(x){x=x-1}),
#' prom_bed, gene_side="bait",
#' bedpe_chr_a="bait_chr", bedpe_start_a="bait_start", bedpe_end_a="bait_end",
#' bedpe_chr_b="otherEnd_chr", bedpe_start_b="otherEnd_start", bedpe_end_b="otherEnd_end"
#' )
#' }
#'
#'
annotate_bedpe2gene <- function(bedpe, prom_bed, gene_side=NULL,
                                   bedpe_chr_a="chr_a", bedpe_start_a="start_a", bedpe_end_a="end_a",
                                   bedpe_chr_b="chr_b", bedpe_start_b="start_b", bedpe_end_b="end_b",
                                   prom_chr="pro_chr", prom_start="pro_start", prom_end="pro_end", prom_anno="pro_anno"
){
        ## convert name
        bedpe = bedpe %>%
                dplyr::rename(
                        chr_a=!!rlang::sym(bedpe_chr_a),
                        start_a=!!rlang::sym(bedpe_start_a),
                        end_a=!!rlang::sym(bedpe_end_a),
                        chr_b=!!rlang::sym(bedpe_chr_b),
                        start_b=!!rlang::sym(bedpe_start_b),
                        end_b=!!rlang::sym(bedpe_end_b)
                )

        prom_bed = prom_bed %>%
                dplyr::rename(
                        pro_chr=!!rlang::sym(prom_chr),
                        pro_start=!!rlang::sym(prom_start),
                        pro_end=!!rlang::sym(prom_end),
                        pro_anno=!!rlang::sym(prom_anno)
                )

        bedpe_cols = c(
                "chr_a"=bedpe_chr_a,
                "start_a"=bedpe_start_a,
                "end_a"=bedpe_end_a,
                "chr_b"=bedpe_chr_b,
                "start_b"=bedpe_start_b,
                "start_b"=bedpe_end_b
        )


        ## gene
        # a side
        a2prom=overlap_df(
                bedpe %>% distinct(chr_a,start_a,end_a),
                prom_bed,
                df1_chr_col="chr_a", df1_start_col="start_a", df1_end_col="end_a",
                df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end",
                df1_0base=T, df2_0base=T, minoverlap=1L
        )
        a2prom=bind_cols(a2prom$overlap_df1, a2prom$overlap_df2)

        # b site
        b2prom=overlap_df(
                bedpe %>% distinct(chr_b,start_b,end_b),
                prom_bed,
                df1_chr_col="chr_b", df1_start_col="start_b", df1_end_col="end_b",
                df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end",
                df1_0base=T, df2_0base=T, minoverlap=1L
        )
        b2prom=bind_cols(b2prom$overlap_df1, b2prom$overlap_df2)

        # add annotation
        if (is.null(gene_side)){
                bedpe_gene = bedpe %>%
                        left_join(
                                a2prom %>%
                                        select(contains("_a"), pro_anno) %>%
                                        dplyr::rename(anno_a=pro_anno)
                        ) %>%
                        left_join(
                                b2prom %>%
                                        select(contains("_b"), pro_anno) %>%
                                        dplyr::rename(anno_b=pro_anno)
                        ) %>%
                        filter(!is.na(anno_a)|!is.na(anno_b))
        }else{
                gene_side = unique(gsub("chr|start|end","",names(bedpe_cols[grepl(gene_side, bedpe_cols)])))
                gene_anno_name = paste0("anno",gene_side)
                if(length(gene_side)==1){
                        if (gene_side=="_a"){
                                half2prom=a2prom
                        }else if(gene_side=="_b"){
                                half2prom=b2prom
                        }
                        bedpe_gene = bedpe %>%
                                left_join(
                                        half2prom %>%
                                                select(contains(gene_side), pro_anno) %>%
                                                dplyr::rename_at(c("pro_anno"), function(x){gene_anno_name})
                                ) %>%
                                filter_at(gene_anno_name, function(x){!is.na(x)})
                }else{
                        stop("cannot find gene_side match")
                }
        }
        bedpe_gene
}
