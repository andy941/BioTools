library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(viridis)
library(optparse)

option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="reference gene", metavar="character"),
  make_option(c("-A", "--ConditionA"), type="character", default=NULL, 
              help="first condition to compare", metavar="character"),
  make_option(c("-B", "--ConditionB"), type="character", default=NULL, 
              help="second condition to compare", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--sep"), type="character", default="\t", 
              help="separator [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="qPCR_", 
              help="output file base name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) || is.null(opt$ref)){
  print_help(opt_parser)
    stop("At least one argument must be supplied for -f and -r", call.=FALSE)
}

# SET variables for analysis
file = opt$file
reference = opt$ref
CondA = opt$ConditionA 
CondB = opt$ConditionB
sep = opt$sep 
out = opt$out

# ggplot2 theme
theme <- theme(axis.text.x = element_text(angle=0, vjust=0.5, colour="black", size = 12),
               axis.text.y = element_text(colour="black", size = 10),
               axis.title.x = element_text(size=14),
               axis.title.y = element_text(size=14),
               strip.text.x = element_text(size=14),
               strip.text.y = element_text(size=14),
               axis.ticks=element_line(colour="black"),
               panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               plot.title = element_text(size = 16, colour="black", hjust=0.5),
               plot.subtitle = element_text(size = 12, colour="black", hjust=0.5),
               plot.margin=unit(c(1,1,1,1),"line"))


# ------------------------------------------------------------------------------------------


dat <- read.csv(file, 
                 sep = sep, 
                 header = T) %>%
  data.frame()

if (mean(c("SampleName", "CrossingPoint", "Replica", "Primers") %in% colnames(dat)) != 1) {
  print_help(opt_parser)
	stop("Dataset HAS to have this columns:\n'SampleName' 'CrossingPoint' 'Replica' 'Primers'\n")
}

# ------------------------------------------------------------------------------------------


ggplot(dat, aes(Primers, CrossingPoint)) +
  facet_wrap(~SampleName) +
  geom_boxplot(aes(fill=Primers),
               outlier.shape = NA, 
               show.legend=F) +
  geom_jitter(aes(shape=Replica), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
              size = 2,
              alpha = 0.7) +
  labs(title = paste0("Overview of Cp Values"),
       subtitle = paste0("Raw values")) +
  theme +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T)
ggsave(filename = paste0(out,"Overview.pdf"), device = "pdf", width = 12, height = 8)

#__________________________________________________________________________________________________________________________

dat <- dat %>%
  group_by(SampleName, Primers, Replica) %>%
  summarise(CrossingPoint=mean(CrossingPoint)) %>%
  ungroup()

ggplot(dat, aes(Primers, CrossingPoint, fill=Primers)) +
  facet_wrap(~SampleName) +
  geom_boxplot(aes(fill=Primers),
               outlier.shape = NA, 
               show.legend=F) +
  geom_jitter(aes(shape=Replica), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
              size = 2,
              alpha = 0.7) +
  labs(title = paste0("Overview of Cp Values"),
       subtitle = paste0("Means of Technical Replicas")) +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  theme
ggsave(filename = paste0(out,"Overview_TechMeans.pdf"), device = "pdf", width = 12, height = 8)

#__________________________________________________________________________________________________________________________

dat <- dat %>%
  group_by(SampleName) %>%
  mutate(CrossingPoint=CrossingPoint - mean(CrossingPoint[Primers==reference])) %>%
  filter(Primers!=reference)

ggplot(dat, aes(Primers, 2**-CrossingPoint, fill=Primers)) +
  facet_wrap(~SampleName) +
  geom_boxplot(aes(fill=Primers),
               outlier.shape = NA, 
               show.legend=F) +
  geom_jitter(aes(shape=Replica), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
              size = 2,
              alpha = 0.7) +
  labs(title = paste0("Normalized Expression Values"),
       subtitle = paste0("Reference Gene: ", reference)) +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  theme
ggsave(filename = paste0(out,"Normalized.pdf"), device = "pdf", width = 12, height = 8)

#__________________________________________________________________________________________________________________________

if (!(is.null(CondA) && is.null(CondB))) {

	dat <- dat %>%
		group_by(SampleName) %>%
		summarise(CrossingPoint=mean(CrossingPoint[Primers==CondA]) - mean(CrossingPoint[Primers==CondB]))


	ggplot(dat, aes(SampleName, 2**-CrossingPoint, fill=SampleName)) +
		geom_bar(stat="identity", position = "dodge", show.legend=F) +
		geom_text(aes(label= format(2**-CrossingPoint, nsmall = 0, digits=3, scientific = FALSE)), 
				  position=position_dodge(.9), hjust=.5, vjust=-.4) +
				   labs(title = paste0("Compare Expression Values"),
						subtitle = paste0("Comparison: ", CondA, " / ", CondB)) +
				   scale_fill_viridis(discrete = T) +
				   scale_color_viridis(discrete = T) +
	   theme
ggsave(filename = paste0(out,"Compare.pdf"), device = "pdf", width = 12, height = 8)
}
