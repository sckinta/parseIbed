#' Example promoter bed data
#'
#' prom_bed example data list the coordinates of non-pseudogene promoter region
#' (-1500 ~ 500bp of TSS in gencodeV19 annotation), and corresponding annotation (gene_name+transcript_id)
#' and gene_name
#'
#' @docType data
#'
#' @usage data(prom_bed)
#'
#' @format a data.frame with column names
#' pro_chr, pro_start, pro_end, pro_anno, strand and gene_id
#'
#' @author Chun Su
#'
#'
#'
#' @examples
#' data(prom_bed)
#' \donttest{
#' head(prom_bed)}
"prom_bed"
