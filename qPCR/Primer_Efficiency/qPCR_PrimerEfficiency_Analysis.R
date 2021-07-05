library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--sep"), type="character", default="\t", 
              help="separator [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="qPCR_", 
              help="output file base name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for -f", call.=FALSE)
}

# Dataset HAS to have this columns: 
# "SampleName" "CrossingPoint" "Replica" "Primers" 

# SET variables for analysis
file = opt$file
sep = opt$sep
out = opt$out
Dilution_Factor = 10  # It doesn't really matter fro the resulting efficiency but can be tweaked

# ------------------------------------------------------------------------------------------

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
  dplyr::select(SampleName, CrossingPoint, Primers, Template) %>%
  data.frame()
dat <- read.csv(file, 
                sep = sep, 
                header = T) %>%
  data.frame()

if (mean(c("SampleName", "CrossingPoint", "Primers", "Template") %in% colnames(dat)) != 1) {
  print_help(opt_parser)
  stop("Dataset HAS to have this columns:\n'SampleName' 'CrossingPoint' 'Primers' 'Template'\n")
}

# ------------------------------------------------------------------------------------------


ggplot(dat, aes(SampleName, CrossingPoint)) +
  facet_wrap(~Primers, dir = "v") +
  # geom_boxplot(aes(fill=SampleName),
  #              outlier.shape = NA, 
  #              show.legend=F) +
  geom_jitter(aes(color=SampleName), 
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
  group_by(SampleName, Template, Primers) %>%
  summarise(CrossingPoint=mean(CrossingPoint)) %>%
  mutate(Template = log(Template, Dilution_Factor)) %>%
  ungroup()

ggplot(dat, aes(Template, CrossingPoint, color=Primers)) +
  geom_point(size = 2,
             alpha = 0.7) +
  geom_smooth(method = 'lm', se=F) +
  stat_regline_equation(label.x.npc = "middle", label.y.npc = "top") +
  labs(title = paste0("Overview of Cp Values"),
       subtitle = paste0("Means of Technical Replicas")) +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  theme
ggsave(filename = paste0(out,"Overview_TechMeans.pdf"), device = "pdf", width = 12, height = 8)

#__________________________________________________________________________________________________________________________

dat <- dat %>%
  group_by(Primers) %>%
  summarise(Efficiency=Dilution_Factor^(-1/summary(lm(formula = CrossingPoint ~ Template))$coefficients[2,1]) -1)


ggplot(dat, aes(Primers, Efficiency, fill=Primers)) +
  geom_bar(stat="identity", position = "dodge", show.legend=F) +
  geom_text(aes(label= format(Efficiency, nsmall = 0, digits=3, scientific = FALSE)), 
           position=position_dodge(.9), hjust=.5, vjust=-.4) +
  labs(title = "Primer Efficiencies",
       subtitle = "Formula: DilutionFactor^( - 1/slope ) - 1") +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  theme
ggsave(filename = paste0(out,"Primer_Efficiencies.pdf"), device = "pdf", width = 12, height = 8)

