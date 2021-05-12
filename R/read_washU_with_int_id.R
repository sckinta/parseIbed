#' Read chicago interaction in washU format
#'
#' @param file chicago output in washU format
#'
#' @export
#'
#' @import dplyr tidyr
#'
read_washU_with_int_id <- function(file){
        df = readr::read_delim(file, delim="\t", col_names=c("bait","oe","score"))
        df = df %>%
                separate(bait, c("bait_chr","bait_start","bait_end")) %>%
                separate(oe,c("otherEnd_chr","otherEnd_start","otherEnd_end")) %>%
                mutate_at(c("bait_start","bait_end","otherEnd_start","otherEnd_end"), as.integer)
        df %>% mutate(int_id = row_number())
}
