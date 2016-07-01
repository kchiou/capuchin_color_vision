#!/usr/bin/env Rscript

load.packages <- function(required.packages) {
	for (i in required.packages)  if (!i %in% installed.packages()[,'Package']) install.packages(i,repos='http://cran.wustl.edu')
	suppressPackageStartupMessages( for (i in required.packages) require(i,character.only=TRUE) )
}

# Load required packages (installing if necessary)
load.packages(required.packages=c('ggplot2','lme4','MASS','car','lsmeans','reshape2'))

# ========================================================================================
# === Import data
# ========================================================================================

# Read in dataset exported from MS Access
# "CapuchinForagingDataCorrected.txt" is a tab-delimited table of behavioral data with 23 columns
# Each row contains information about a new state behavior
# The following are descriptions of relevant columns

# StateBegin             : timestamp of the start time of the state behavior
# StateDuration          : duration of the state behavior
# Name                   : name of the animal observed performing the state behavior
# CurrentGroup           : group membership of the observed animal
# AgeClass               : age class of the observed animal
# Sex                    : sex of the observed animal
# ColorVisionType        : color vision status (di/trichromacy) of the observed animal
# Behavior               : observed state behavior
# Event                  : concatenated single-letter codes of foraging events that took
#                          place within the state
# PhenologyID            : phenological identifier for a foraging patch
# ScientificName         : binomial Latin name of the plant species
# BoutBegin              : timestamp of the start time of a feeding bout
# BoutDuration           : duration of a feeding bout
# CorrectedPhenologyID   : phenological identifier (corrected after scanning for errors)
# SubBoutDuration        : duration of subbouts (bouts that were necessarily subdivided
#                          because they traversed multiple foraging patches). These are
#                          later calculated for all bouts, but sometimes were done manually, in which
#                          case they are entered here

raw.data <- read.delim('data/CapuchinForagingDataCorrected.txt')

# "FruitColoration.txt" is a tab-delimited table containing information about taxonomy
# and visual characteristics of fruits

fruit.coloration <- read.delim('data/FruitColoration.txt')

# The following are descriptions of relevant columns

# Code                   : project code for the fruit taxon
# Genus                  : generic name of the fruit taxon (=dwc:genus)
# Species                : specific epithet of the fruit taxon (=dwc:specificEpithet)
# FruitColoration        : color category of the fruit

# ========================================================================================
# === Define functions
# ========================================================================================

# * * * Credit to Stack Overflow user caracal * * *
# http://stackoverflow.com/questions/13847936/in-r-plotting-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot

# Visualize random effects (re = object of class ranef.mer)
ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
    require(ggplot2)
    f <- function(x) {
        pv   <- attr(x, "postVar")
        cols <- 1:(dim(pv)[1])
        se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
        pDf  <- data.frame(y=unlist(x)[ord],
                           ci=1.96*se[ord],
                           nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                           ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                           ind=gl(ncol(x), nrow(x), labels=names(x)))

        if(QQ) {  ## normal QQ-plot
            p <- ggplot(pDf, aes(nQQ, y))
            p <- p + facet_wrap(~ ind, scales="free")
            p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
        } else {  ## caterpillar dotplot
            p <- ggplot(pDf, aes(ID, y)) + coord_flip()
            if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                p <- p + facet_wrap(~ ind)
            } else {           ## different scales for random effects
                p <- p + facet_grid(ind ~ ., scales="free_y")
            }
            p <- p + xlab("Levels") + ylab("Random effects")
        }

        p <- p + theme(legend.position="none")
        p <- p + geom_hline(yintercept=0)
        p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
        p <- p + geom_point(aes(size=1.2), colour="blue")
        return(p)
    }

    lapply(re, f)
}

# ========================================================================================
# === Data cleanup
# ========================================================================================

cleaned.data <- raw.data

# Strip whitespaces from the Events
levels(cleaned.data$Event) <- gsub(' ','',levels(cleaned.data$Event))

# Correct dates and times
cleaned.data$Date <- as.Date(cleaned.data$Date,format='%m/%d/%Y')
cleaned.data$StateBegin <- as.POSIXlt(cleaned.data$StateBegin,tz='America/Costa_Rica',format='%m/%d/%Y %I:%M:%S %p')
cleaned.data$BoutBegin <- as.POSIXlt(cleaned.data$BoutBegin,tz='America/Costa_Rica',format='%m/%d/%Y %I:%M:%S %p')

# Retain only events that involve eating (E) or rejecting (R)
cleaned.data <- cleaned.data[grep('[ER]',cleaned.data$Event),]

if (length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event)))) {
	## If matches both E and R, ignore that row and flip a warning

	warning('Found ',length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event))),' ambiguous (E or R) event',if (length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event))) > 1) 's','.',sep='')
	cleaned.data <- cleaned.data[-intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event)),]
}


# Correct phenology IDs and disambiguate bout durations for bouts that had multiple phenology IDs
cleaned.data$CorrectedPhenologyID[cleaned.data$CorrectedPhenologyID == ''] <- NA
cleaned.data$CorrectedPhenologyID <- cleaned.data$CorrectedPhenologyID[,drop=TRUE]

# Clean up factor levels of phenology IDs
cleaned.data$PhenologyID <- factor(cleaned.data$PhenologyID,levels=union(levels(cleaned.data$PhenologyID),levels(cleaned.data$CorrectedPhenologyID)))

# Rewrite phenology IDs with corrected phenology IDs
cleaned.data$PhenologyID[!is.na(cleaned.data$CorrectedPhenologyID)] <- as.character(cleaned.data$CorrectedPhenologyID[!is.na(cleaned.data$CorrectedPhenologyID)])

# Reorder the cleaned dataset
cleaned.data <- cleaned.data[order(cleaned.data$BoutBegin,cleaned.data$StateBegin),]

# Some remaining entries remain confusing/ambiguous

# DECISION D1: Amanda made the executive decision 9/24/2015 to remove the Jacquinia nervosa record from Nutella/W07AM63/13:13:02 (ambiguous entry)
cleaned.data <- cleaned.data[!(cleaned.data$Name %in% 'Nutella' & cleaned.data$PhenologyID %in% 'W07AM63' & cleaned.data$BoutBegin == '2007-01-26 13:13:02' & cleaned.data$ScientificName %in% 'Jacquinia nervosa'),]

# DECISION D2: Kenny made the decision to change phenology ID of Marmite/W07AM1131/11:11:08 to W07AM1132 where the species is Zuelania guidonia (apparent user error)
cleaned.data[cleaned.data$Name %in% 'Marmite' & cleaned.data$BoutBegin %in% '2007-04-13 11:11:08' & cleaned.data$ScientificName %in% 'Zuelania guidonia',]$PhenologyID <- 'W07AM1132'

# Add a bout ID to cleaned data (no such ID exists so must be inferred)
# The combination of an animal name and the bout timestamp will uniquely identify the bout
multi.bout.id <- c('Name','BoutBegin')

# Extract unique combinations of Name-BoutBegin
unique.bouts <- unique(cleaned.data[,multi.bout.id])
unique.bouts <- unique.bouts[order(unique.bouts$BoutBegin,unique.bouts$Name),]

# Assign integers to serve as unique identifiers for each bout
unique.bouts$BoutID <- 1:nrow(unique.bouts)

# Merge new bout IDs with the working data
cleaned.data <- merge(cleaned.data,unique.bouts,by=multi.bout.id,all.x=TRUE,sort=FALSE)

# DECISION D3: Bouts as defined thus far are insufficient because it is assumed that each
# bout can only encompass one set of phenological conditions (=one unique PhenologyID)
# In practice this is not the case so we must identify subbouts and determine their
# durations (which do not yet exist in the data)

# Bouts are subdivided into subbouts if the bout traversed multiple patches
# The combination of animal name, phenology conditions, and the bout timestamp will uniquely identify the subbout
multi.subbout.id <- c('Name','PhenologyID','BoutBegin','BoutDuration')

# Extract unique combinations of Name-PhenologyID-BoutBegin
unique.subbouts <- unique(cleaned.data[,multi.subbout.id])
unique.subbouts <- unique.subbouts[order(unique.subbouts$BoutBegin,unique.subbouts$Name,unique.subbouts$PhenologyID),]

# Assign subbout IDs
unique.subbouts$SubBoutID <- 1:nrow(unique.subbouts)

# Merge new subbout IDs with the working data
cleaned.data <- merge(cleaned.data,unique.subbouts,by=multi.subbout.id,all.x=TRUE,sort=FALSE)

# - - - - - - - - - - - - Calculate bout durations - - - - - - - - - - - - #

bouts <- split(cleaned.data,cleaned.data$BoutID)

# Calculate CorrectedBoutDuration (note that this is different than the SubBoutDuration earlier, which was intended to deal with a different case)
cleaned.data$CorrectedBoutDuration <- NA

# For each bout, calculate the number of unique phenology IDs per bout
uIDs <- do.call(c,lapply(bouts,function(x) length(unique(x$PhenologyID))))

# Get indices of "singular" bouts (iS) that already have only one phenology ID
iS <- cleaned.data$BoutID %in% as.integer(names(uIDs[which(uIDs == 1)]))

# Use the existing bout duration for cases where the bout already had only one phenology ID
cleaned.data$CorrectedBoutDuration[iS] <- cleaned.data$BoutDuration[iS]

# Extract bouts that need to be split (uIDs > 1)
split.bouts <- bouts[as.integer(names(uIDs[which(uIDs > 1)]))]

# For bouts that need splitting, go through and calculate subbout durations
# Sub-bout durations can be calculated on the basis of timed state behaviors that occurred within each bout
bouts.split.by.phenologies <- do.call(rbind,lapply(split.bouts,function(x) {
	# Get subject, duration, and begin time of the bout (only one possible for each)
	name <- unique(x$Name)
	bout.duration <- unique(x$BoutDuration)
	bout.begin <- unique(x$BoutBegin)

	# Duration should be calculated on the basis of states. The first state begin is not always the same as the bout begin time.
	first.state.begin <- min(x$StateBegin)
	
	# Split each bout into subbouts defined by unique phenology IDs
	phenologies <- do.call(rbind,lapply(split(x,x$PhenologyID[,drop=TRUE]),function(y) {
		phenology.id <- unique(y$PhenologyID)
		# Calculate metadata for the subbouts (most importantly, the state begin time)
		data.frame(PhenologyID=phenology.id,Name=name,BoutDuration=bout.duration,BoutBegin=bout.begin,StateBegin=min(y$StateBegin))
	}))
	
	# Order subbouts by the state begin time
	phenologies <- phenologies[order(phenologies$StateBegin),]

	# For each row (except the last), calculate the end time as the begin time of the subsequent row
	# The end time for the last row should be the bout begin time + the bout duration, taking into account the difference between the first state and the beginning of the bout
	phenologies$NextStateBegin <- c(phenologies$StateBegin[2:nrow(phenologies)],bout.begin + bout.duration + as.integer(difftime(first.state.begin,bout.begin,units='secs')))
	
	# The duration of the subbout is the difference between the state begin time and the state end time calculated above
	phenologies$PhenologyDuration <- as.integer(difftime(phenologies$NextStateBegin,phenologies$StateBegin,units='secs'))

	# Return subbouts
	phenologies
}))

row.names(bouts.split.by.phenologies) <- NULL

# Rename StateBegin in the split bouts to SubBoutBegin
names(bouts.split.by.phenologies)[names(bouts.split.by.phenologies) %in% 'StateBegin'] <- 'SubBoutBegin'

# Grab the relevant columns
phenology.durations <- bouts.split.by.phenologies[,c('Name','PhenologyID','BoutBegin','PhenologyDuration','SubBoutBegin')]

# Merge the newly calculated SubBout information back into cleaned data
cleaned.data <- merge(cleaned.data,phenology.durations,by=c('Name','PhenologyID','BoutBegin'),all.x=TRUE,sort=FALSE)

# Use the calculated PhenologyDuration (= sub-bout duration) as the CorrectedBoutDuration (which were left NA for bouts that needed splitting)
cleaned.data$CorrectedBoutDuration[is.na(cleaned.data$CorrectedBoutDuration)] <- cleaned.data$PhenologyDuration[is.na(cleaned.data$CorrectedBoutDuration)]

# The BoutDuration is now the CorrectedBoutDuration
cleaned.data$BoutDuration <- cleaned.data$CorrectedBoutDuration

# DECISION D4: We subsequently decided to recalculate durations across the board based on state times rather than bout times
#
# A bout is now defined thusly:
# (1) a bout starts at the time of the StateBegin that contains the first foraging event
# (2) a bout ends at the time of the StateEnd that contains the last foraging event
# (3) a new bout duration is calculated between the time of the first state begin and last state end containing foraging events within a bout.

# Calculate the state end time
cleaned.data$StateEnd <- cleaned.data$StateBegin + cleaned.data$StateDuration

# Calculate for each subbout the duration (temporarily named "ForagingDuration")
foraging.states <- do.call(rbind,lapply(split(cleaned.data,cleaned.data$SubBoutID),function(x) {
	data.frame(SubBoutID=unique(x$SubBoutID),ForagingBegin=min(x$StateBegin),ForagingEnd=max(x$StateEnd),ForagingDuration=as.numeric(difftime(max(x$StateEnd),min(x$StateBegin),units='secs')),stringsAsFactors=FALSE)
}))

# Merge this information back into cleaned.data
cleaned.data <- merge(cleaned.data,foraging.states,by='SubBoutID',all.x=TRUE,sort=FALSE)

# - - - - - - - - - - - - Tabulate foraging outcomes - - - - - - - - - - - - #

# Set up the column as a factor
cleaned.data$ForagingOutcome <- factor(NA,levels=c('Eat','Reject'))

# If the Event contains the code "E", the fruit was eaten
cleaned.data$ForagingOutcome[grep('E',cleaned.data$Event)] <- 'Eat'

# If the Event contains the code "R", the fruit was rejected
cleaned.data$ForagingOutcome[grep('R',cleaned.data$Event)] <- 'Reject'

# - - - - - - - - - - - - Assemble full dataset - - - - - - - - - - - - #

# The full dataset for analysis is named "analysis.table"

# These are the important metadata fields for each bout/subbout
multi.bout.subbout.id <- c('BoutID','SubBoutID','BoutDuration','ForagingDuration','Name','PhenologyID','ScientificName','BoutBegin','SubBoutBegin')

# Each row of the dataset is equivalent to a bout or subbout
analysis.table <- unique(cleaned.data[,multi.bout.subbout.id])

# SubBoutBegin does not exist for proper bouts that did not require splitting
# Set it to BoutBegin for these cases
analysis.table$SubBoutBegin[is.na(analysis.table$SubBoutBegin)] <- analysis.table$BoutBegin[is.na(analysis.table$SubBoutBegin)]

# Sort this by time and then animal
analysis.table <- analysis.table[order(analysis.table$SubBoutBegin,analysis.table$Name),]

# Convert the animal data to strings to avoid unexpected behaviors
analysis.table$Name <- as.character(analysis.table$Name)


# These are the important metadata fields for each individual
multi.animal.id <- c('Name','CurrentGroup','AgeClass','Sex','ColorVisionType')

unique.animals <- unique(cleaned.data[,multi.animal.id])

unique.animals$Name <- as.character(unique.animals$Name)

# Maturity is a redefined AgeClass for analysis
unique.animals$Maturity <- as.character(unique.animals$AgeClass)

# Adults and subadults are "Mature"
unique.animals$Maturity[unique.animals$AgeClass %in% c('Adult','Subadult')] <- 'Mature'

# Infants are excluded from analysis
unique.animals$Maturity[unique.animals$AgeClass %in% c('Infant')] <- NA

# Finalize this factor
unique.animals$Maturity <- factor(unique.animals$Maturity,levels=c('Mature','LargeImmature','SmallImmature'))

# Merge animal information into the dataset
analysis.table <- merge(analysis.table,unique.animals,by='Name',all.x=TRUE,sort=FALSE)

# Replace BoutDuration with ForagingDuration and backup BoutDuration (see Decision D4)
analysis.table$BoutDurationOriginal <- analysis.table$BoutDuration
analysis.table$BoutDuration <- analysis.table$ForagingDuration

# Drop ForagingDuration
analysis.table$ForagingDuration <- NULL


# Foraging outcomes must be merged back into the dataset

# Calculate number of eat/reject per subbout
foraging.outcomes <- table(cleaned.data$ForagingOutcome,cleaned.data$SubBoutID)
analysis.table$NumberEaten <- as.numeric(foraging.outcomes['Eat',][analysis.table$SubBoutID])
analysis.table$NumberRejected <- as.numeric(foraging.outcomes['Reject',][analysis.table$SubBoutID])

# Eat/Reject rates are calculated as counts divided by the bout duration (see decision D4)
analysis.table$RateEaten <- analysis.table$NumberEaten / analysis.table$BoutDuration
analysis.table$RateRejected <- analysis.table$NumberRejected / analysis.table$BoutDuration

# The investigation rate is the combined eat/reject rate
analysis.table$InvestigationRate <- (analysis.table$NumberEaten + analysis.table$NumberRejected) / analysis.table$BoutDuration

# The acceptance index is the fraction of inspected fruits that are eaten
analysis.table$AcceptanceIndex <- analysis.table$NumberEaten / (analysis.table$NumberEaten + analysis.table$NumberRejected)


# Keep only bouts where at least one thing was eaten and it lasted more than zero seconds
analysis.table <- analysis.table[as.logical(analysis.table$NumberEaten) & as.logical(analysis.table$BoutDuration) & as.logical(analysis.table$BoutDurationOriginal),]


# Incorporate metadata on fruits

# Concatenate genus + species
fruit.coloration$ScientificName <- factor(paste(fruit.coloration$Genus,fruit.coloration$Species),levels=levels(analysis.table$ScientificName))
fruit.coloration <- fruit.coloration[,c('ScientificName','Code','FruitColoration')]

# Rename column to SpeciesCode for clarity
names(fruit.coloration)[names(fruit.coloration) == 'Code'] <- 'SpeciesCode'

# Merge fruit metadata back into the dataset
analysis.table <- merge(analysis.table,fruit.coloration,by='ScientificName',all.x=TRUE,sort=FALSE)


# Preferred sorting of columns in the dataset
column.sort.order <- c('BoutID','SubBoutID','BoutBegin','SubBoutBegin','BoutDuration','BoutDurationOriginal','Name','CurrentGroup','AgeClass','Maturity','Sex','ColorVisionType','PhenologyID','ScientificName','SpeciesCode','FruitColoration','NumberEaten','NumberRejected','RateEaten','RateRejected','InvestigationRate','AcceptanceIndex')

# Preferring sorting of rows in the dataset
row.sort.order <- order(analysis.table$BoutBegin,analysis.table$SubBoutBegin,analysis.table$Name,analysis.table$PhenologyID)

analysis.table <- analysis.table[row.sort.order,column.sort.order]
rownames(analysis.table) <- NULL

# Use factors when appropriate and prettify their codes
analysis.table$Name <- factor(analysis.table$Name)
levels(analysis.table$AgeClass)[levels(analysis.table$AgeClass)=='LargeImmature'] <- 'Large Immature'
levels(analysis.table$AgeClass)[levels(analysis.table$AgeClass)=='SmallImmature'] <- 'Small Immature'
analysis.table$AgeClass <- factor(analysis.table$AgeClass,levels=c('Adult','Subadult','Large Immature','Small Immature','Infant'))
levels(analysis.table$Maturity)[levels(analysis.table$Maturity)=='LargeImmature'] <- 'Large Immature'
levels(analysis.table$Maturity)[levels(analysis.table$Maturity)=='SmallImmature'] <- 'Small Immature'
analysis.table$Maturity <- factor(analysis.table$Maturity,levels=c('Mature','Large Immature','Small Immature'))
analysis.table <- droplevels(analysis.table)

# Save the dataset if necessary
# write.table(analysis.table,file="analysis_table.txt",sep='\t',row.names=FALSE,quote=FALSE)

# ========================================================================================
# === Exploratory data visualization
# ========================================================================================

# Histogram of the feeding rate
ggplot(analysis.table,aes(RateEaten)) +
	geom_histogram(binwidth=0.05,position='dodge') +
	facet_wrap(~AgeClass,ncol=2) +
	xlab('Feeding Rate') + ylab('Count') + ggtitle('Feeding rate distribution')
ggsave(filename='output/DataExploration_feeding_rate_distribution.pdf')

# Boxplot + jitter plot of the feeding rate by color vision phenotype
ggplot(analysis.table,aes(ColorVisionType,RateEaten,color=Sex)) +
	geom_boxplot(outlier.shape=21) +
	geom_jitter(alpha=0.2) +
	facet_wrap(~AgeClass,nrow=1) +
	theme(axis.text.x=element_text(angle=-30,hjust=0),panel.grid=element_blank()) +
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ColorVisionType)-0.5,1),color='white',size=0.5) +
	xlab('Color Vision Type') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type')
ggsave(filename='output/DataExploration_feeding_rate_by_color_vision_type.pdf')

# The same plot by sex (dichromats only)
ggplot(analysis.table[analysis.table$ColorVisionType %in% 'Dichromat',],aes(Sex,RateEaten,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(legend.title = element_blank(),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$Sex)-0.5,1),color='white',size=0.5) + 
	xlab('Sex') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by sex')
ggsave(filename='output/DataExploration_feeding_rate_by_sex_dichromats_only.pdf')

# Boxplot of feeding rate by color vision phenotype and fruit species
ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(ScientificName,RateEaten,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ScientificName)-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Plant Species') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type and plant species')
ggsave(filename='output/DataExploration_feeding_rate_by_color_vision_type_and_species.pdf')

# Boxplot of feeding rate by color vision phenotype and fruit coloration
ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(FruitColoration,RateEaten,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(panel.grid=element_blank()) + geom_vline(xintercept=seq(1.5,nlevels(analysis.table$FruitColoration)-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Fruit Coloration') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type and fruit coloration')
ggsave(filename='output/DataExploration_feeding_rate_by_color_vision_type_and_fruit_coloration.pdf')

# Scatterplot of feeding rate by the acceptance index
ggplot(analysis.table,aes(round(AcceptanceIndex,2),RateEaten)) + 
	geom_point() + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~Name,ncol=floor(sqrt(nlevels(analysis.table$Name)))) + 
	theme(axis.text.x=element_text(angle=-45,size=6),axis.text.y=element_text(size=6)) + 
	xlab('Acceptance Index') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by acceptance index (by animal)')
ggsave(filename='output/DataExploration_feeding_rate_by_acceptance_index.pdf')

# Histogram of the acceptance index distribution by age class
ggplot(analysis.table,aes(AcceptanceIndex)) + 
	geom_histogram(binwidth=0.05,position='dodge') + 
	facet_wrap(~AgeClass,ncol=2) + 
	xlab('Acceptance Index') + ylab('Count') + ggtitle('Acceptance index distribution')
ggsave(filename='output/DataExploration_acceptance_index_distribution.pdf')

# Boxplot + jitter plot of the acceptance index by color vision phenotype
ggplot(analysis.table,aes(ColorVisionType,AcceptanceIndex,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(axis.text.x=element_text(angle=-30,hjust=0),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ColorVisionType)-0.5,1),color='white',size=0.5) + 
	xlab('Color Vision Type') + ylab('Acceptance Index') + ggtitle('Acceptance index by color vision type')
ggsave(filename='output/DataExploration_acceptance_index_by_color_vision_type.pdf')

# Boxplot + jitter plot of the acceptance index by sex (dichromats only)
ggplot(analysis.table[analysis.table$ColorVisionType %in% 'Dichromat',],aes(Sex,AcceptanceIndex,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(legend.title = element_blank(),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$Sex)-0.5,1),color='white',size=0.5) + 
	xlab('Sex') + ylab('Acceptance Index') + ggtitle('Acceptance index by sex')
ggsave(filename='output/DataExploration_acceptance_index_by_sex_dichromats_only.pdf')

# Boxplot of the acceptance index by color vision phenotype and fruit species
ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(ScientificName,AcceptanceIndex,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ScientificName)-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Plant Species') + ylab('Acceptance Index') + ggtitle('Acceptance index by color vision type and plant species')
ggsave(filename='output/DataExploration_acceptance_index_by_color_vision_type_and_species.pdf')

# Boxplot + jitter plot of bout duration by age class
ggplot(analysis.table,aes(AgeClass,BoutDuration,color=AcceptanceIndex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~ColorVisionType+Sex) + 
	scale_color_continuous(name='Acceptance Index') + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	xlab('Age Class') + ylab('Length of Bout') + ggtitle('Lengths of bouts by age class')
ggsave(filename='output/DataExploration_bout_duration_by_age_class.pdf')

# Scatterplot of the bout duration by the acceptance index
ggplot(analysis.table,aes(AcceptanceIndex,BoutDuration,color=ColorVisionType)) + 
	geom_point(alpha=0.5) + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Acceptance Index') + ylab('Length of Bout') + ggtitle('Lengths of bouts by acceptance index')
ggsave(filename='output/DataExploration_bout_duration_by_acceptance_index.pdf')

# Scatterplot of the bout duration by the rejection rate ("pickiness")
ggplot(analysis.table,aes(RateRejected,BoutDuration,color=ColorVisionType)) + 
	geom_point(alpha=0.5) + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~AgeClass+Sex,nrow=1) + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Rate Rejected') + ylab('Length of Bout') + ggtitle('Lengths of bouts by pickiness')
ggsave(filename='output/DataExploration_bout_duration_by_pickiness.pdf')

# ========================================================================================
# === Finalize the dataset
# ========================================================================================

# Drop infants from the analysis
analysis.table <- analysis.table[!is.na(analysis.table$Maturity),]


# ========================================================================================
# === Explore data distributions
# ========================================================================================

RateEaten <- analysis.table$RateEaten

# Add one to allow analysis of nonzero nonnegative distributions
RateEatenPlusOne <- RateEaten + 1

# Calculate parameters for Poisson and gamma distributions
po <- fitdistr(RateEatenPlusOne,'Poisson')
ga <- fitdistr(RateEatenPlusOne,'gamma')

# Test gaussian, log-normal, exponential, and gamma distributions
pdf(file='output/DistributionTests.pdf')
	layout(t(matrix(1:4,nrow=2)))
	qqp(RateEaten,distribution='norm',ylab='Feeding rate',xlab='Quantiles',main='Gaussian distribution')
	qqp(RateEaten,distribution='lnorm',ylab='Feeding rate',xlab='Quantiles',main='Log-normal distribution')
	qqp(RateEaten,distribution='exp',ylab='Feeding rate',xlab='Quantiles',main='Exponential distribution')
	qqp(RateEatenPlusOne,distribution='gamma',shape=ga$estimate[[1]],rate=ga$estimate[[2]],ylab='Feeding rate + 1',xlab='Quantiles',main='Gamma distribution')
dev.off()


# ========================================================================================
# === Model fitting (main model)
# ========================================================================================

# Linear mixed model with:
# Feeding rate as dependent variable
# Maturity, color vision phenotype, and fruit species as fixed effects
# Individual animal and phenological ID as random effects

mixed.model <- lmer(RateEaten ~ Maturity + ColorVisionType + ScientificName + (1|Name) + (1|PhenologyID),data=analysis.table)

## Type 2 ANOVA table from package 'car'
mixed.model.anova <- Anova(mixed.model)

# Calculate least-square means (color vision phenotype)
cv.lsmeans <- summary(lsmeans(mixed.model,'ColorVisionType'))

# Calculate least-square means (maturity)
mat.lsmeans <- summary(lsmeans(mixed.model,'Maturity'))

# Write results to file
sink(file='output/model_fitting.log')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- Main model\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
#	cat('# - - - - - - - - - - - - - - - - - - - - Summary - - - - - - - - - - - - - - - - - - - - \n')
#	cat('\n')
#	print(summary(mixed.model))
#	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(mixed.model.anova)
	cat('\n')
	cat('# - - - - - - - - - - - Least-square means (color vision phenotype) - - - - - - - - - - - \n')
	cat('\n')
	print(cv.lsmeans)
	cat('\n')
	cat('# - - - - - - - - - - - - - -  Least-square means (maturity)  - - - - - - - - - - - - - - \n')
	cat('\n')
	print(mat.lsmeans)
	cat('\n')
sink()

# Visualize random effects

model.random.effects <- ranef(mixed.model,condVar=TRUE)

# Plot using the ggCaterpillar function
random.effects.plot <- ggCaterpillar(model.random.effects,QQ=FALSE)

# Plot random effects (phenological ID)
random.effects.plot[[1]] + xlab('Phenology ID') + ggtitle('Random effects: phenology ID')
ggsave(filename='output/GLMM_random_effects_phenology_id.pdf')

# Plot random effects (individual ID)
random.effects.plot[[2]] + xlab('Animal Name') + ggtitle('Random effects: animal name')
ggsave(filename='output/GLMM_random_effects_animal_id.pdf')


# ========================================================================================
# === Model fitting (dichromats only)
# ========================================================================================

# Re-run the linear mixed model only males and females

# Limit to dichromats only because trichromats cannot be male
dichromat.table <- analysis.table[analysis.table$ColorVisionType == 'Dichromat',]

# Linear mixed model with:
# Feeding rate as dependent variable
# Sex, maturity, and fruit species as fixed effects
# Individual animal and phenological ID as random effects

sex.model <- lmer(RateEaten ~ Sex + Maturity + ScientificName + (1|Name) + (1|PhenologyID),data=dichromat.table)

# Calculate least square means (sex)
sex.lsmeans <- summary(lsmeans(sex.model,'Sex'))

# Write results to file
sink(file='output/model_fitting.log',append=TRUE)
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- Sex model (dichromats only)\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
#	cat('# - - - - - - - - - - - - - - - - - - - - Summary - - - - - - - - - - - - - - - - - - - - \n')
#	cat('\n')
#	print(summary(sex.model))
#	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(sex.model))
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - Least-square means (sex)  - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(sex.lsmeans)
	cat('\n')
sink()


# ========================================================================================
# === Model fitting (investigation rate and acceptance index)
# ========================================================================================

# Linear mixed model with:
# Investigation rate as dependent variable
# Maturity, color vision phenotype, and fruit species as fixed effects
# Individual animal and phenological ID as random effects

ir.model <- lmer(InvestigationRate ~ Maturity + ColorVisionType + ScientificName + (1|Name) + (1|PhenologyID),data=analysis.table)

# Linear mixed model with:
# Acceptance index as dependent variable
# Maturity, color vision phenotype, and fruit species as fixed effects
# Individual animal and phenological ID as random effects

ai.model <- lmer(AcceptanceIndex ~ Maturity + ColorVisionType + ScientificName + (1|Name) + (1|PhenologyID),data=analysis.table)

# Write results to file
sink(file='output/model_fitting.log',append=TRUE)
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- Investigation rate model\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
#	cat('# - - - - - - - - - - - - - - - - - - - - Summary - - - - - - - - - - - - - - - - - - - - \n')
#	cat('\n')
#	print(summary(ir.model))
#	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(ir.model))
	cat('\n')
	cat('# - - - - - - - - - - - Least-square means (color vision phenotype) - - - - - - - - - - - \n')
	cat('\n')
	print(summary(lsmeans(ir.model,'ColorVisionType')))
	cat('\n')
	cat('# - - - - - - - - - - - - - -  Least-square means (maturity)  - - - - - - - - - - - - - - \n')
	cat('\n')
	print(summary(lsmeans(ir.model,'Maturity')))
	cat('\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- Acceptance index model\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
#	cat('# - - - - - - - - - - - - - - - - - - - - Summary - - - - - - - - - - - - - - - - - - - - \n')
#	cat('\n')
#	print(summary(ai.model))
#	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(ai.model))
	cat('\n')
	cat('# - - - - - - - - - - - Least-square means (color vision phenotype) - - - - - - - - - - - \n')
	cat('\n')
	print(summary(lsmeans(ai.model,'ColorVisionType')))
	cat('\n')
	cat('# - - - - - - - - - - - - - -  Least-square means (maturity)  - - - - - - - - - - - - - - \n')
	cat('\n')
	print(summary(lsmeans(ai.model,'Maturity')))
	cat('\n')
sink()

# ========================================================================================
# === Post-hoc tests
# ========================================================================================

# Run version of model without ColorVisionType or ScientificName to examine effect of color vision on different species

# Linear mixed model with:
# Feeding rate as dependent variable
# Maturity as fixed effect
# Individual animal and phenological ID as random effects
less.model <- lmer(RateEaten ~ Maturity + (1|Name) + (1|PhenologyID),data=analysis.table)

# Analyze residuals from this model
resid.analysis <- analysis.table
resid.analysis$Residuals = resid(less.model)

# Visualize these residuals
ggplot(resid.analysis,aes(ScientificName,Residuals,fill=ColorVisionType)) + 
	geom_boxplot() + 
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,face='italic'),panel.grid=element_blank()) + 
	geom_hline(yintercept=0,color='black',size=0.5)
ggsave(filename='output/LessModel_residuals.pdf')

# Assemble p-value, test statistic, and degrees of freedom of a linear model of the residuals by color vision phenotype for each fruit species separately
test <- as.data.frame.matrix(do.call(rbind,lapply(split(resid.analysis,resid.analysis$ScientificName),function(x) {
	if (all(table(x$ColorVisionType) > 0)) {
		# Only do this if both color vision phenotypes are represented
		c(
			p=anova(lm(Residuals ~ ColorVisionType,x))$`Pr(>F)`[1],
			F=anova(lm(Residuals ~ ColorVisionType,x))$`F value`[1],
			df=anova(lm(Residuals ~ ColorVisionType,x))$`Df`[1]
		)
	} else {
		# If either color vision phenotype is not represented, return NA
		rep(NA,3)
	}
})))
names(test) <- c('cv.res.p.value','cv.res.F','cv.res.df')

# Rerun the main model on each fruit species separately
# Then calculate p-values, test statistics, and degrees of freedom for maturity and color vision phenotype and assemble the data
test <- cbind(test,data.frame(do.call(rbind,lapply(split(analysis.table,analysis.table$ScientificName),function(x) {
	test <- try(Anova(lmer(RateEaten ~ Maturity + ColorVisionType + (1|Name) + (1|PhenologyID),data=x)))
	if ('try-error' %in% class(test)) {
		c('mt.p.value'=NA,'mt.Chisq'=NA,'mt.df'=NA,'cv.p.value'=NA,'cv.Chisq'=NA,'cv.df'=NA)
	} else {
		p.values <- test$`Pr(>Chisq)`[which(attr(test,'row.names') %in% c('Maturity','ColorVisionType'))]
		chisq.values <- test$`Chisq`[which(attr(test,'row.names') %in% c('Maturity','ColorVisionType'))]
		df.values <- test$`Df`[which(attr(test,'row.names') %in% c('Maturity','ColorVisionType'))]
		values <- c(p.values[1],chisq.values[1],df.values[1],p.values[2],chisq.values[2],df.values[2])
		names(values) <- c('mt.p.value','mt.Chisq','mt.df','cv.p.value','cv.Chisq','cv.df')
		values
	}
}))))

# Calculate number of trichromats represented in bouts for each fruits species
test$trichromat.n <-  do.call(c,lapply(split(resid.analysis,resid.analysis$ScientificName),function(x) {
	table(x$ColorVisionType)['Trichromat']
}))

# Calculate number of dichromats represented in bouts for each fruits species
test$dichromat.n <-  do.call(c,lapply(split(resid.analysis,resid.analysis$ScientificName),function(x) {
	table(x$ColorVisionType)['Dichromat']
}))

# Adjust for multiple comparisons using the Holm adjustment method
test$cv.res.holm <- p.adjust(test$cv.res.p.value,'holm')
test$mt.holm <- p.adjust(test$mt.p.value,'holm')
test$cv.holm <- p.adjust(test$cv.p.value,'holm')

# Write posthoc test results to file
write.csv(test,file='output/posthoc.csv')

# ========================================================================================
# === Calculate table information for publication
# ========================================================================================

# Summarize bout counts, sum durations, feeding counts (events), and reject counts for each fruit species
bouts.by.species <- melt(table(analysis.table$ScientificName))
duration.by.species <- melt(tapply(analysis.table$BoutDuration,analysis.table$ScientificName,sum))
events.by.species <- melt(tapply(analysis.table$NumberEaten,analysis.table$ScientificName,sum))
rejected.by.species <- melt(tapply(analysis.table$NumberRejected,analysis.table$ScientificName,sum))

# Tally up summary info for each bout (by plant species)

# Start with number of bouts
bout.summary <- bouts.by.species
names(bout.summary) <- c('species','bouts')

# Add duration, feeding events, and reject events
bout.summary <- data.frame(
	bout.summary,
	duration = duration.by.species$value ,
	feeding_events = events.by.species$value ,
	reject_events = rejected.by.species$value
)

# Order of plants given in the table
plant.order <- c("Sciadodendron excelsum", "Cordia guanacastensis", "Cordia panamensis",
"Bursera simaruba", "Diospyros salicifolia", "Sloanea terniflora",
"Erythroxylum havanense", "Vachellia collinsii", "Casearia arguta",
"Casearia sylvestris", "Muntingia calabura", "Zuelania guidonia",
"Ficus cotinifolia", "Ficus hondurensis", "Ficus morazaniana",
"Ficus obtusifolia", "Ficus ovalis", "Maclura tinctoria", "Trophis racemosa",
"Karwinskia calderoni", "Krugiodendron ferreum", "Genipa americana",
"Randia monantha", "Randia thurberi", "Allophylus occidentalis",
"Dipterodendron costaricense", "Manilkara chicle", "Simarouba glauca",
"Jacquinia nervosa")

# Order summary accordingly
bout.summary <- bout.summary[match(plant.order,bout.summary$species),]

# Write info to file
write.csv(bout.summary,file='output/bout_summary.csv')

# ========================================================================================
# === Prepare figures for publication
# ========================================================================================

# vis.table is the main analysis table modified for visualization
vis.table <- analysis.table

# Order the fruit species by increasing mean feeding rate
sorted.scientific.name <- names(sort(tapply(vis.table$RateEaten,vis.table$ScientificName,mean)))
vis.table$ScientificName <- factor(vis.table$ScientificName,levels=sorted.scientific.name)

# Bring in information on conspicuous fruits
# Derived from Excel sheet "Species color categorizations_June12 2016.xlsx"
conspicuous <- structure(c(1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1), .Names = c("acol", "aocc",
"bsim", "carg", "cgua", "cpan", "csyl", "dcos", "dsal", "ehav",
"fcot", "fhon", "fmor", "fobt", "fova", "game", "jpun", "kcal",
"kfer", "mcal", "mchi", "mtin", "rmon", "rthu", "sexc", "sgla",
"ster", "trac", "zgui"))

# Conspicuity is coded using species codes, so need to translate this information to 
# binomial scientific names used in the dataset
plant.species <- unique(vis.table[,c('ScientificName','SpeciesCode')])
plant.species <- plant.species[order(plant.species$SpeciesCode),]
rownames(plant.species) <- NULL

# Check that plant.species codes and conspicuous codes are the same
if (!identical(as.character(plant.species$SpeciesCode),names(conspicuous))) stop('Species codes do not match!')

# Merge plant.species and conspicuous
conspicuous <- data.frame(plant.species,Conspicuous = conspicuous)

# Sort this information in the same order as the species order in the plot
conspicuous <- conspicuous[match(levels(conspicuous$ScientificName),conspicuous$ScientificName),]

# Check that the order of conspicuous codes and their appearance in the plot are the same
if (!identical(as.character(conspicuous$ScientificName),levels(conspicuous$ScientificName))) stop('Species orders do not match!')

# Calculate max and min x-axis coordinates for rectangles
conspicuous$xmin <- seq(0.5,nrow(conspicuous) - 0.5,1)
conspicuous$xmax <- seq(1.5,nrow(conspicuous) + 0.5,1)

# Recode conspicuity (0 == 'not conspicuous' ; 1 == 'conspicuous')
conspicuous$Conspicuous <- c('not conspicuous','conspicuous')[as.numeric(conspicuous$Conspicuous)+1]

# Reorder the factor levels
conspicuous$Conspicuous <- factor(conspicuous$Conspicuous,levels=c('conspicuous','not conspicuous'))

theme_set(theme_bw(base_size = 6))

# Figure: Feeding rate by color vision phenotype
ggplot(vis.table,aes(ColorVisionType,RateEaten,color=Sex)) + 
	geom_boxplot(outlier.shape=21,outlier.size=0.75,size=0.75) + 
	geom_jitter(alpha=0.2,size=0.75,pch=16) + 
	facet_wrap(~Maturity,nrow=1) + 
	theme(panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(vis.table$ColorVisionType)-0.5,1),color='white',size=0.5) + 
	scale_color_manual(name='Sex',values=c('#db8aa1','#00bfc4')) + 
	scale_x_discrete(labels=c('Dichr.','Trichr.')) +
	xlab('Color Vision Type') + ylab('Feeding Rate')
ggsave(filename='output/fig_feeding.rate.by.color.vision.pdf',width=3.2,height=2.4)

# Figure: Least-square means for color vision phenotype
ggplot(cv.lsmeans,aes(ColorVisionType,lsmean,color=ColorVisionType)) + 
	geom_point(size=2) + 
	geom_errorbar(aes(ymin=lsmean-SE,ymax=lsmean+SE),width=0.5,size=0.5) + 
	xlab('Color Vision Type') + ylab('Feeding Rate (Least-Square Means)') + 
	theme(legend.position='none') + 
	scale_color_manual(values=c('#000000','#000000'))
ggsave(filename='output/fig_lsmeans.color.vision.pdf',width=1.6,height=2.4,useDingbats=FALSE)

# Figure: Least-square means for maturity
ggplot(mat.lsmeans,aes(Maturity,lsmean,color=Maturity)) + 
	geom_point(size=2) + 
	geom_errorbar(aes(ymin=lsmean-SE,ymax=lsmean+SE),width=0.5,size=0.5) + 
	xlab('Age Class') + ylab('Feeding Rate (Least-Square Means)') + 
	theme(legend.position='none') + 
	scale_color_manual(values=c('#000000','#000000','#000000'))
ggsave(filename='output/fig_lsmeans.maturity.pdf',width=1.6,height=2.4,useDingbats=FALSE)

# Figure: Feeding rate by fruit species
theme_set(theme_grey(base_size = 8))
ggplot() +
	geom_rect(data=conspicuous,aes(NULL,xmin=xmin,xmax=xmax,ymin=0,ymax=1,fill=Conspicuous)) +
	scale_fill_manual(values=c('#cccccc','#ebebeb'),name='Fruit Conspicuity') +
	geom_boxplot(data=vis.table,aes(ScientificName,RateEaten,color=ColorVisionType),outlier.shape=21) +
	geom_vline(xintercept=seq(1.5,nlevels(vis.table$ScientificName)-0.5,1),color='white',size=0.5) +
	scale_color_manual(name='Color Vision Type',values=c('#fdc086','#af8dc3'),labels=c('Dichromat','Trichromat')) +
	scale_y_continuous(limits=c(0,1),breaks=c(0,1)) +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,face='italic'),panel.grid=element_blank()) +
	xlab('Plant Species') +
	ylab('Feeding Rate')
ggsave(filename='output/fig_feeding.rate.by.fruit.species.pdf',width=7,height=4)

message('Analysis complete!')