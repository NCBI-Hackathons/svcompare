library(reshape)
library(ggplot2)

#difference from median start and end
bpvar <- read.table("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/giab.breakpoint_variation.tsv", header = TRUE, sep = "\t", na.strings = c("NA", NA, "NaN"))
x<-levels(bpvar$CALLSET)
levels(bpvar$CALLSET) <- c("IllFB","IllGATK","PBMultibreak","IllCNVnator","PBAssmblyt1","Illcommonlaw","Illfermikit","IllMetaSV","PBParlAss","PBParlPB","PBHoneytails","PBsmrtsvdip","Bionano","PBAssemblyt2","CGvcfbeta")

ggplot(bpvar, aes(x=START_DIST)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=5,colour="white") + facet_grid(CALLSET ~ TYPE) + xlim(0,500) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/giab.breakpoint_variation_histstartbycallset.jpg")
ggplot(bpvar, aes(x=END_DIST)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=5,colour="white") + facet_grid(CALLSET ~ TYPE) + xlim(0,500) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/giab.breakpoint_variation_histendbycallset.jpg")
#ggplot(bpvar[bpvar$CALLSET=="Bionano",], aes(x=START_DIST)) +  geom_histogram(binwidth=10,colour="white") + facet_grid(CALLSET ~ .)
sum(bpvar$CALLSET=="Bionano")


#count per callset by type
countbytech <- read.table("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/giab.total_callset_count_by_tech.tsv", header = FALSE, sep = "\t", na.strings = c("NA", NA, "NaN"))
colnames(countbytech) <- c("CALLSET","TYPE","COUNT")
levels(countbytech$CALLSET) <- c("IllFB","IllGATK","PBMultibreak","IllCNVnator","PBAssmblyt1","Illcommonlaw","Illfermikit","IllMEI","IllMetaSV","PBParlAss","PBParlPB","PBHoneytails","PBsmrtsvdip","Bionano","PBAssemblyt2","CGvcfbeta")
ggplot(countbytech, aes(x=CALLSET,y=COUNT)) +  geom_bar(stat="identity") + facet_grid(. ~ TYPE) + theme(text = element_text(size=12),axis.text.x = element_text(angle=90, vjust=1)) 
ggsave("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/giab.total_callset_count_by_tech_barplot.jpg")


#callsets per event by type
callsetsperevent <- read.csv("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/callsets_per_event_type_aware.txt", header = TRUE, na.strings = c("NA", NA, "NaN"))
callsetspereventmelt <- melt(callsetsperevent)
ggplot(callsetspereventmelt, aes(x=value)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=1,colour="white") + facet_grid(. ~ variable) + xlim(0,20) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/callsets_per_event_type_aware_histbytype.jpg")


#Max difference between breakpoints between callsets
maxdiffbreakpoint <- read.csv("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/coorcordance_break_points_3b_i.txt", header = TRUE, na.strings = c("NA", NA, "NaN"))
callsetspereventmelt <- melt(callsetsperevent)
ggplot(maxdiffbreakpoint, aes(x=StartDifference)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=10,colour="black")+ xlim(0,2000) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/AJTrio/SVs/hackathon_Aug2016/coorcordance_break_points_hist.jpg")
