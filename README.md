# parseIbed
parse Ibed bedpe R package

## installation
devtools::install_github("sckinta/parseIbed", ref = "main")

## functions
read_ibed_with_int_id(ibed, parse_b2b=NULL, max_score=F) # ibed can be df or file

read_washU_with_int_id(file)

summarise_from_ibed(ibed, baitmap) # ibed can be df or file

annotate_bedpe2gene(bedpe, prom_bed) # also work for ibed. ?annotate_bedpe2gene

annotate_bedpe2geneOCR(bedpe, ocr_bed, prom_bed) # also work for ibed. ?annotate_bedpe2geneOCR
