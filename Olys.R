## OLYS TEMPORAL SCRIPT

### WORKING DIRECTORY AND PACKAGES
#################################################################################
setwd('/Users/jmartine/Desktop/Olys/')

library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(RColorBrewer) #display.brewer.all()
library(zoo)
library(corrplot) 
library(lme4) 
library(lubridate)
library(car)
library(wesanderson)
library(viridis)
library(ggeffects)
library(plotrix)

### LOADING DATA
data <- read.csv('Oly_full.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') 
str(data)
data$y <- as.numeric(data$y)
data$worms <- as.numeric(data$worms)

### CALCULATING PREVALENCE 
data_pop <- data %>%
        group_by(population, valve) %>% #
        summarize_at(vars(worms), list(mean = mean), na.rm=TRUE)

data_summ <- data_pop %>%
        group_by(population) %>% #
        summarize_at(vars(mean), list(mean = mean), na.rm=TRUE)

data_summ <- as.data.frame(data_summ)


## range of sizes fossil pop and std error
min(data$y, na.rm=T) # 0.8
max(data$y, na.rm=T) # 6.3
mean(data$y, na.rm=T) # 3.21
std.error(data$y, na.rm=T) # 0.0378

### PLOTTING PREVALENCE
plot <- ggplot(data_summ, aes(x= population, y= mean, fill= population, show.legend = FALSE))  +
        scale_fill_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
        geom_bar(stat = "identity", col = "black", show.legend = FALSE) + 
        ylim(0,1) +
        labs(x = 'Population', y = 'Prevalence', size=16) 
plot + theme_classic(base_size = 18) 

data_y <- data %>%
        group_by(population) %>% #
        summarize_at(vars(y), list(mean = mean), na.rm=TRUE)

data_summ <- as.data.frame(data_summ)

### PLOTTING SHELL HEIGHT
ploty <- ggplot(data, aes(x= population, y= y, fill= population, show.legend = FALSE))  +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(size=4, shape = 21, show.legend = FALSE) +
        scale_fill_manual(values=wes_palette("GrandBudapest1", n = 4)) + 
        ylim(0,7) +
        labs(x = 'Population', y = 'Shell height (cm)', size=16) 
ploty + theme_classic(base_size = 18) 


level_order <- c('Modern','Archeo','Recent F', 'Fossil') 

histdens <- ggplot(data, aes(x= y, fill= factor(population, level=level_order), show.legend = FALSE))  +
        geom_histogram(aes(y=..density..), alpha=0.9, show.legend = FALSE) +
        geom_density(alpha=.9) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 4)) +
        labs(x = 'Shell height (cm)', y = 'Frequency', size=16) 
histdens + theme_classic(base_size = 18) 

hist <- ggplot(data, aes(x= y, fill= population))  +
        geom_histogram(aes(y=..density..), alpha=0.8) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 4)) +
        labs(x = 'Shell height (cm)', y = 'Frequency', size=16) 
hist + theme_classic(base_size = 18) 

dens <- ggplot(data, aes(x= y, fill= population, show.legend = FALSE))  +
        geom_density(alpha=.8) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 4)) +
        labs(x = 'Shell height (mm)', y = 'Frequency', size=16) 
dens + theme_classic(base_size = 18) 

### LOADING TRACES DATA
traces <- read.csv('Oly_traces.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') 
str(traces)
traces$TracesPresent <- as.numeric(traces$TracesPresent)
traces$TracesA <- as.numeric(traces$TracesA)
traces$TracesB <- as.numeric(traces$TracesB)

### PLOTTING TRACE WIDTH

tracewidth <- ggplot(traces, aes(x= factor(Population, level=level_order), y= mmWidth, fill= Population, show.legend = FALSE))  +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(size=4, shape = 21, show.legend = FALSE) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 3)) + 
        labs(x = 'Population', y = 'Burrow width (mm)', size=16) 
tracewidth + theme_classic(base_size = 18) 


denswidth <- ggplot(traces, aes(x= mmWidth, fill= Population, show.legend = FALSE))  +
        geom_density(alpha=.8) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 3)) +
        labs(x = 'Burrow width (mm)', y = 'Frequency', size=16) 
denswidth + theme_classic(base_size = 18) 

### PLOTTING TRACE LENGTH
tracelength <- ggplot(traces, aes(x= factor(Population, level=level_order), y= mmLength, fill= Population, show.legend = FALSE))  +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(size=4, shape = 21, show.legend = FALSE) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 3)) + 
        labs(x = 'Population', y = 'Burrow length (mm)', size=16) +
        ylim(0,30)
tracelength + theme_classic(base_size = 18) 

denslength <- ggplot(traces, aes(x= mmLength, fill= Population, show.legend = FALSE))  +
        geom_density(alpha=.8) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 3)) +
        labs(x = 'Burrow length (mm)', y = 'Frequency', size=16) 
denslength + theme_classic(base_size = 18) 


## PLOTTING TRACE INTENSITY, AND CHI-SQ TEST
traceA <- traces %>%
        group_by(Population) %>% #
        summarize_at(vars(TracesA), list(sum = sum), na.rm=TRUE)

traceB <- traces %>%
        group_by(Population) %>% #
        summarize_at(vars(TracesB), list(sum = sum), na.rm=TRUE)

traceC <- traces %>%
        group_by(Population) %>% #
        summarize_at(vars(TracesC), list(sum = sum), na.rm=TRUE)

level_order <- c('Recent F', 'Archeo', 'Modern') 

traceabundA <- ggplot(traceA, aes(x= factor(Population, level=level_order), y=sum, fill= Population, show.legend = FALSE))  +
        scale_fill_manual(values=wes_palette("Darjeeling2")) + 
        geom_bar(stat = "identity", col = "black", show.legend = FALSE) +
        ylim(0,30) +
        labs(x = 'Population', y = 'No. of oysters with 1-4 burrows', size=16) 
traceabundA + theme_classic(base_size = 18) 

traceabundB <- ggplot(traceB, aes(x= factor(Population, level=level_order), y=sum, fill= Population, show.legend = FALSE))  +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 4)) + 
        geom_bar(stat = "identity", col = "black", show.legend = FALSE) +
        ylim(0,30) +
        labs(x = 'Population', y = 'No. of oysters with 5-10 burrows', size=16) 
traceabundB + theme_classic(base_size = 18) 

traceabundC <- ggplot(traceC, aes(x= factor(Population, level=level_order), y=sum, fill= Population, show.legend = FALSE))  +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 4)) + 
        geom_bar(stat = "identity", col = "black", show.legend = FALSE) +      
        ylim(0,30) +
        labs(x = 'Population', y = 'No. of oysters with >11 burrows', size=16) 
traceabundC + theme_classic(base_size = 18) 

traceA <- as.data.frame(traceA)
traceB <- as.data.frame(traceB)
traceC <- as.data.frame(traceC)

chisq.test(x=traceA$sum) # X-squared = 8.3158, df = 2, p-value = 0.01564
chisq.test(x=traceB$sum) # X-squared = 0.89655, df = 2, p-value = 0.6387
chisq.test(x=traceC$sum) # X-squared = 9.4146, df = 2, p-value = 0.009029

# Modern vs Recent Fossil
ModT <- subset(traces, traces$Population == "Modern")
ArchT <- subset(traces, traces$Population == "Archeo")
FosT <- subset(traces, traces$Population == "Recent F")

hist(ModT$TracesC)

t.test(ModT$TracesC, y=FosT$TracesC) # t = 4.3695, df = 59.616, p-value = 5.054e-05

# Modern vs Archeo
t.test(ModT$TracesC, y=ArchT$TracesC) # t = 2.0483, df = 58.954, p-value = 0.04499

# R Fossil vs Archeo
t.test(FosT$TracesC, y=ArchT$TracesC) # t = -1.9731, df = 55.529, p-value = 0.05347

### NON-PARAMETRIC WILCOX TEST, SAME RESULTS
wilcox.test(ModT$TracesC, FosT$TracesC) # W = 714, p-value = 0.000136
wilcox.test(ModT$TracesC, ArchT$TracesC) # W = 602, p-value = 0.04631
wilcox.test(FosT$TracesC, y=ArchT$TracesC) # W = 345, p-value = 0.0552


### LOADING RADIOCARBON DATA
radiocarbon <- read.csv('radiocarbon.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') 
str(radiocarbon)

### PLOTTING 
boxages <- ggplot(radiocarbon, aes(x= population, y= C14, fill= population, show.legend = FALSE))  +
        geom_boxplot(alpha=0.7, lwd=1, outlier.shape = NA, show.legend = FALSE) + 
        geom_point(size=4, shape = 21, show.legend = FALSE) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 2)) + 
        labs(x = 'Population', y = 'C14 age (BP)', size=16) 
boxages + theme_classic(base_size = 18) 

densage <- ggplot(radiocarbon, aes(x= C14, fill= population, show.legend = FALSE))  +
        geom_density(alpha=.8) +
        scale_fill_manual(values=wes_palette("Darjeeling2", n = 2)) +
        labs(x = '14 age (BP)', y = 'Frequency', size=16) 
densage + theme_classic(base_size = 18) 

### MODEL
model <- glmer(worms ~ y + population + (1|site), family="binomial", data = data)
summary(model)
anova(model)
car::Anova(model, type=3) # getting p-values 

plotmod <- ggpredict(model,c("y","population"))
plotmod2 <- ggpredict(model,c("population"))



allpops_plot <- ggplot(plotmod,aes(x,predicted,color=group), color=group) +
        scale_color_manual(values=wes_palette("Darjeeling2", n = 4)) + 
        geom_point(size=4) +
        geom_errorbar(data=plotmod, mapping=aes(x=x, ymin=conf.low, ymax=conf.high), width=0.18) +
        geom_line(aes(group=group)) +
        xlab("Shell height (cm)") +
        ylab(expression(paste("Predicted infestation"))) +
        theme_classic() +
        guides(color=guide_legend("Population")) +
        theme(plot.title=element_text(size=14,hjust=0.5,face="plain"), axis.text.y=element_text(size=14), 
              axis.title.y=element_text(size=14), axis.text.x=element_text(size=14), axis.title.x=element_text(size=14),
              panel.grid.minor=element_line(color=NA))
allpops_plot

level_order <- c('Fossil','Recent F','Archeo', 'Modern') 
#aes(x= factor(season, level= level_order), y= predicted,

time_plot <- ggplot(plotmod2,aes(x,predicted,color=group), color=group) +
        scale_color_manual(values=wes_palette("Darjeeling2", n = 4)) + 
        geom_point(size=4) +
        geom_errorbar(data=plotmod2, mapping=aes(x=x, ymin=conf.low, ymax=conf.high), width=0.18) +
        #geom_line(aes(group=group)) +
        xlab("Shell height (cm)") +
        ylab(expression(paste("Predicted infestation"))) +
        theme_classic() +
        guides(color=guide_legend("Population")) +
        theme(plot.title=element_text(size=14,hjust=0.5,face="plain"), axis.text.y=element_text(size=14), 
              axis.title.y=element_text(size=14), axis.text.x=element_text(size=14), axis.title.x=element_text(size=14),
              panel.grid.minor=element_line(color=NA))
time_plot


##################
### More LOADING DATA
dataplot <- read.csv('Oly plot.csv', header=TRUE, row.names = NULL, stringsAsFactors=FALSE,  sep=',') 
heastr(dataplot)

#modern <- subset(dataplot, dataplot$Population=='Modern')
#historical <- subset(dataplot, dataplot$Population=='Historical')

level_order <- c('Fossil','Recent F','Archeological', 'Modern - MB', 'Modern - SB', 'Modern - SB2', 'Modern - DH') 
#"#ECCBAE" "#046C9A" "#D69C4E" "#ABDDDE"

plotprev <- ggplot(dataplot, aes(x= factor(Location, level= level_order), y= Prevalence, fill= Population, show.legend = FALSE))  +
        scale_fill_manual(values = c("#D69C4E", "#046C9A")) +
        geom_bar(stat = "identity", col = "black", show.legend = FALSE) + 
        ylim(0,100) +
        labs(x = 'Population', y = 'Prevalence', size=16)
plotprev + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


geom_hline(yintercept=40, linetype='dotted', col = 'red')+
annotate("text", x = "Feb", y = 40, label = "Previous Level", vjust = -0.5)