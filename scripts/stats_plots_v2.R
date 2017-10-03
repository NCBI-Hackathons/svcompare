library(reshape)
library(ggplot2)

#difference from median start and end
bpvar <- read.table("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation.tsv", header = TRUE, sep = "\t", na.strings = c("NA", NA, "NaN"))
x<-levels(bpvar$CALLER)
levels(bpvar$CALLER) <- c("IllFB","IllGATK","PBMultibreak","IllCNVnator","PBAssmblyt1","Illcommonlaw","Illfermikit","IllMetaSV","PBParlAss","PBParlPB","PBHoneytails","PBsmrtsvdip","Bionano","PBAssemblyt2","CGvcfbeta")

ggplot(bpvar, aes(x=START_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=2,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,500) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histstartbycallset.jpg")
ggplot(bpvar, aes(x=STOP_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=2,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,500) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histendbycallset.jpg")
#ggplot(bpvar[bpvar$CALLER=="Bionano",], aes(x=START_DIST)) +  geom_histogram(binwidth=10,colour="white") + facet_grid(CALLER ~ .)
sum(bpvar$CALLER=="Bionano")

sizes <- c(0,50,100,300,400,1000,3000,10000)
sizenames<- c("salt50","sb50to100","sc100to300","sd300to400","se400to1000","sf1000to3000","sg3000to10000","shgt10000")
sizeno <- length(sizes)
for (k in 1:(sizeno-1)) {
  ggplot(bpvar[bpvar$LENGTH>=sizes[k] & bpvar$LENGTH<sizes[k+1],], aes(x=START_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,1) + theme(strip.text.y = element_text(size = 6, colour = "black", angle = 270))
  ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histstartbycallset_",sizes[k],"-",sizes[k+1],"bp.jpg"))
  ggplot(bpvar[bpvar$LENGTH>=sizes[k] & bpvar$LENGTH<sizes[k+1],], aes(x=STOP_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.5) + theme(strip.text.y = element_text(size = 6, colour = "black", angle = 270))
  ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histendbycallset_",sizes[k],"-",sizes[k+1],"bp.jpg"))
  
}
ggplot(bpvar[bpvar$LENGTH>=sizes[sizeno],], aes(x=START_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.3) + theme(strip.text.y = element_text(size = 6, colour = "black", angle = 270))
ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histstartbycallset_gt",sizes[sizeno],"bp.jpg"))
ggplot(bpvar[bpvar$LENGTH>=sizes[sizeno],], aes(x=STOP_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.3) + theme(strip.text.y = element_text(size = 6, colour = "black", angle = 270))
ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histendbycallset_gt",sizes[sizeno],"bp.jpg"))

#only assemblytics and MetaSV
for (k in 1:(sizeno-1)) {
  ggplot(bpvar[bpvar$LENGTH>=sizes[k] & bpvar$LENGTH<sizes[k+1] & (bpvar$CALLER=="PBAssmblyt1" | bpvar$CALLER=="IllMetaSV" ),], aes(x=START_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,1) + theme(strip.text.y = element_text(size = 10, colour = "black", angle = 270)) + labs(x="Distance from median start",y="Proportion of calls")
  ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histstartbycallset_PBAssmblyt1andIllMetaSV_",sizes[k],"-",sizes[k+1],"bp.jpg"))
  ggplot(bpvar[bpvar$LENGTH>=sizes[k] & bpvar$LENGTH<sizes[k+1] & (bpvar$CALLER=="PBAssmblyt1" | bpvar$CALLER=="IllMetaSV" ),], aes(x=STOP_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.5) + theme(strip.text.y = element_text(size = 10, colour = "black", angle = 270)) + labs(x="Distance from median end",y="Proportion of calls")
  ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histendbycallset_PBAssmblyt1andIllMetaSV_",sizes[k],"-",sizes[k+1],"bp.jpg"))
  
}
ggplot(bpvar[bpvar$LENGTH>=sizes[sizeno] & (bpvar$CALLER=="PBAssmblyt1" | bpvar$CALLER=="IllMetaSV" ),], aes(x=START_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.3) + theme(strip.text.y = element_text(size = 10, colour = "black", angle = 270)) + labs(x="Distance from median start",y="Proportion of calls")
ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histstartbycallset_PBAssmblyt1andIllMetaSV_gt",sizes[sizeno],"bp.jpg"))
ggplot(bpvar[bpvar$LENGTH>=sizes[sizeno] & (bpvar$CALLER=="PBAssmblyt1" | bpvar$CALLER=="IllMetaSV" ),], aes(x=STOP_DIST_TO_MEDIAN)) +  geom_histogram(aes(y = (..density..)),binwidth=1,colour="black") + facet_grid(CALLER ~ TYPE) + xlim(0,100) + ylim(0,0.3) + theme(strip.text.y = element_text(size = 10, colour = "black", angle = 270)) + labs(x="Distance from median end",y="Proportion of calls")
ggsave(paste0("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.breakpoint_variation_histendbycallset_PBAssmblyt1andIllMetaSV_gt",sizes[sizeno],"bp.jpg"))

#count per callset by type
countbytech <- read.table("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.callset_total_count_by_type.tsv", header = TRUE, sep = "\t", na.strings = c("NA", NA, "NaN"))
colnames(countbytech) <- c("CALLER","TYPE","COUNT")
levels(countbytech$CALLER) <- c("IllFB","IllGATK","PBMultibreak","IllCNVnator","PBAssmblyt1","Illcommonlaw","Illfermikit","IllMEI","IllMetaSV","PBParlAss","PBParlPB","PBHoneytails","PBsmrtsvdip","Bionano","PBAssemblyt2","CGvcfbeta")
ggplot(countbytech, aes(x=CALLER,y=COUNT)) +  geom_bar(stat="identity") + facet_grid(. ~ TYPE) + theme(text = element_text(size=12),axis.text.x = element_text(angle=90, vjust=1)) 
ggsave("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.total_callset_count_by_tech_barplot.jpg")


#callsets per event by type
callsetsperevent <- read.table("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.callset_total_count_by_type.tsv", header = TRUE, sep = "\t", na.strings = c("NA", NA, "NaN"))
callsetspereventmelt <- melt(callsetsperevent)
ggplot(callsetspereventmelt, aes(x=value)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=1,colour="white") + facet_grid(. ~ variable) + xlim(0,20) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/callsets_per_event_type_aware_histbytype.jpg")


#Max difference between breakpoints between callsets
maxdiffbreakpoint <- read.csv("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/testgiab2.coorcordance_break_points_3b_i.txt", header = TRUE, na.strings = c("NA", NA, "NaN"))
callsetspereventmelt <- melt(callsetsperevent)
ggplot(maxdiffbreakpoint, aes(x=StartDifference)) +  geom_histogram(aes(y = (..ncount..)/sum(..ncount..)),binwidth=10,colour="black")+ xlim(0,2000) + theme(strip.text.y = element_text(size = 4, colour = "black", angle = 270))
ggsave("/Users/jzook/Documents/OneDrive - National Institute of Standards and Technology (NIST)/AJTrio/SVs/hackathon_Aug2016/test2/coorcordance_break_points_hist.jpg")
