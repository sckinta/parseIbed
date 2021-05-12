#' annotate_bedpe2geneOCR
#'
#' annotate bedpe loop with one end as gene promoter and the other end as OCR
#'
#' @param bedpe bedpe data.frame
#' @param ocr_bed ocr_bed data.frame
#' @param prom_bed prom_bed data.frame
#' @param gene_side NULL or bedpe column name match (eg."bait");
#' if NULL, search genes on both ends of bedpe.
#' if gene_side is not NULL, ocr_side must not be NULL
#' @param ocr_side NULL or bedpe column name match (eg."otherEnd")
#' if NULL, search ocrs on both ends of bedpe
#' @param bedpe_chr_a character
#' @param bedpe_start_a character
#' @param bedpe_end_a character
#' @param bedpe_chr_b character
#' @param bedpe_start_b character
#' @param bedpe_end_b character
#' @param ocr_chr character
#' @param ocr_start character
#' @param ocr_end character
#' @param ocr_id character
#' @param prom_chr character
#' @param prom_start character
#' @param prom_end character
#' @param prom_anno character
#'
#' @return data.frame with a and b end annotated to pro_anno and ocr_coord:ocr_id
#' @export
#'
#' @import dplyr tidyr DFbedtools
#' @examples
#' # annotate bedpe on both sides
#'
#' \donttest{data(bedpe)
#' data(prom_bed)
#' annotate_bedpe2geneOCR(bedpe, ocr_bed, prom_bed)
#' }
#'
#' # annotate ibed file on one side
#'
#' \donttest{data(ibed)
#' annotate_bedpe2geneOCR(
#' ibed %>% mutate_at(vars(ends_with("_start")), function(x){x=x-1}),
#' ocr_bed,
#' prom_bed, gene_side="bait", ocr_side="otherEnd",
#' bedpe_chr_a="bait_chr", bedpe_start_a="bait_start", bedpe_end_a="bait_end",
#' bedpe_chr_b="otherEnd_chr", bedpe_start_b="otherEnd_start", bedpe_end_b="otherEnd_end"
#' )
#' }
annotate_bedpe2geneOCR <- function(bedpe, ocr_bed, prom_bed, gene_side=NULL, ocr_side=NULL,
                                   bedpe_chr_a="chr_a", bedpe_start_a="start_a", bedpe_end_a="end_a",
                                   bedpe_chr_b="chr_b", bedpe_start_b="start_b", bedpe_end_b="end_b",
                                   ocr_chr="chr", ocr_start="start", ocr_end="end", ocr_id="id",
                                   prom_chr="pro_chr", prom_start="pro_start", prom_end="pro_end", prom_anno="pro_anno"
                                   ){
        ## check gene_side, ocr_side
        if(as.integer(is.null(gene_side)) != as.integer(is.null(ocr_side))){
                stop("gene_side and ocr_side need to be set either both null or neither null")
        }

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

        ocr_bed = ocr_bed %>%
                dplyr::rename(
                        chr=!!rlang::sym(ocr_chr),
                        start=!!rlang::sym(ocr_start),
                        end=!!rlang::sym(ocr_end),
                        id=!!rlang::sym(ocr_id)
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


        ## OCR
        # a side
        a2ocr=overlap_df(
                bedpe_gene %>% distinct(chr_a,start_a,end_a),
                ocr_bed,
                df1_chr_col="chr_a", df1_start_col="start_a", df1_end_col="end_a",
                df2_chr_col="chr", df2_start_col="start", df2_end_col="end",
                df1_0base=T, df2_0base=T, minoverlap=1L
        )
        a2ocr=bind_cols(a2ocr$overlap_df1, a2ocr$overlap_df2)

        # b site
        b2ocr=overlap_df(
                bedpe_gene %>% distinct(chr_b,start_b,end_b),
                ocr_bed,
                df1_chr_col="chr_b", df1_start_col="start_b", df1_end_col="end_b",
                df2_chr_col="chr", df2_start_col="start", df2_end_col="end",
                df1_0base=T, df2_0base=T, minoverlap=1L
        )
        b2ocr=bind_cols(b2ocr$overlap_df1, b2ocr$overlap_df2)

        # add annotation
        if(is.null(ocr_side)){
                bedpe_gene2ocr = bedpe_gene %>%
                        left_join(
                                a2ocr %>%
                                        mutate_at(c("start","end"), function(x){format(x,scientific=F)}) %>%
                                        mutate(ocr_a=glue::glue("{chr}:{start}:{end}:{id}")) %>%
                                        mutate(ocr_a=gsub(" ","",ocr_a)) %>%
                                        select(-chr, -start, -end, -id)
                        ) %>%
                        left_join(
                                b2ocr %>%
                                        mutate_at(c("start","end"), function(x){format(x,scientific=F)}) %>%
                                        mutate(ocr_b=glue::glue("{chr}:{start}:{end}:{id}")) %>%
                                        mutate(ocr_b=gsub(" ","",ocr_b)) %>%
                                        select(-chr, -start, -end, -id)
                        ) %>%
                        filter((!is.na(anno_a) & !is.na(ocr_b)) | (!is.na(anno_b) & !is.na(ocr_a)))
        }else{
                ocr_side = unique(gsub("chr|start|end","",names(bedpe_cols[grepl(ocr_side, bedpe_cols)])))
                ocr_anno_name = paste0("ocr",ocr_side)
                if(length(ocr_side)==1){
                        if (ocr_side=="_a"){
                                half2ocr=a2ocr
                        }else if(ocr_side=="_b"){
                                half2ocr=b2ocr
                        }
                        bedpe_gene2ocr = bedpe_gene %>%
                                left_join(
                                        half2ocr %>%
                                                mutate_at(c("start","end"), function(x){format(x,scientific=F)}) %>%
                                                mutate(ocr=glue::glue("{chr}:{start}:{end}:{id}")) %>%
                                                mutate(ocr=gsub(" ","",ocr)) %>%
                                                select(-chr, -start, -end, -id) %>%
                                                dplyr::rename_at(c("ocr"), function(x){ocr_anno_name})
                                ) %>%
                                filter_at(c(gene_anno_name, ocr_anno_name), function(x){!is.na(x)})
                }else{
                        stop("cannot find ocr_side match")
                }
        }


        bedpe_gene2ocr
}
