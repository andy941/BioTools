library(dplyr)
library(biomaRt)
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-H", "--header"), type="logical", default=TRUE, 
              help="header TRUE or FALSE [default= %default]", metavar="character"),
  make_option(c("-s", "--sep"), type="character", default="\t", 
              help="separator [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Annotated_", 
              help="output name prefix (Use empty string to modify 'in place') [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for -f", call.=FALSE)
}

# SET variables for analysis
file = opt$file
sep = opt$sep 
out = opt$out
header = opt$header

# add gene names from biomaRt
ensembl <- useMart(host="https://plants.ensembl.org", biomart = "plants_mart",
                   dataset = "athaliana_eg_gene")

t2g <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "transcript_biotype",
                            "description"), 
             mart = ensembl)

t2g <- dplyr::rename(t2g,
                     ens_gene = ensembl_gene_id, 
                     ext_gene = external_gene_name)
gene_info <- t2g %>%
  dplyr::select(ens_gene, ext_gene, transcript_biotype, description) %>%
  unique.data.frame()

dat <- read.table(file = file,
           header = header,
           sep = sep)
gene_col <- colnames(dat)[1]
by = c("ens_gene")
names(by) = gene_col

# Save the table with info about the genes --------------------------
dat <- left_join(dat, gene_info, by = by)
if (!header) colnames(dat)[1] = "Gene_id"
write.table(dat, 
            file = paste0(out, file),
            sep = sep,
            row.names = F)
