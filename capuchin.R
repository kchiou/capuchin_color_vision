#!/usr/bin/env Rscript

load.packages = function(required.packages) {
	for (i in required.packages)  if (!i %in% installed.packages()[,'Package']) install.packages(i,repos='http://cran.wustl.edu')
	suppressPackageStartupMessages( for (i in required.packages) require(i,character.only=TRUE) )
}

# Load required packages (installing if necessary)
load.packages(required.packages=c('ggplot2','lme4','MASS','car','lsmeans','lmerTest','reshape2','gdata'))

# ========================================================================================
# === Import data
# ========================================================================================

# Read in dataset exported from MS Access
# "CapuchinForagingData.txt" is a tab-delimited table of behavioral data with 23 columns
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

raw.data = read.delim('data/CapuchinForagingData.txt')

# "FruitColoration.txt" is a tab-delimited table containing information about taxonomy
# and visual characteristics of fruits

fruit.coloration = read.delim('data/FruitColoration.txt')

# The following are descriptions of relevant columns

# Code                   : project code for the fruit taxon
# Genus                  : generic name of the fruit taxon (=dwc:genus)
# Species                : specific epithet of the fruit taxon (=dwc:specificEpithet)
# FruitColoration        : color category of the fruit

# Separate information on fruit conspicuity

conspicuous = structure(c(1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1), .Names = c("acol", "aocc",
"bsim", "carg", "cgua", "cpan", "csyl", "dcos", "dsal", "ehav",
"fcot", "fhon", "fmor", "fobt", "fova", "game", "jpun", "kcal",
"kfer", "mcal", "mchi", "mtin", "rmon", "rthu", "sexc", "sgla",
"ster", "trac", "zgui")) # cpan changed to conspicuous

# Import rank information

dom.rank = read.xls('data/2007_2008Dominance.xlsx',stringsAsFactors=FALSE)

# ========================================================================================
# === Define functions
# ========================================================================================

# * * * Credit to Stack Overflow user caracal * * *
# http://stackoverflow.com/questions/13847936/in-r-plotting-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot

# Visualize random effects (re = object of class ranef.mer)
ggCaterpillar = function(re, QQ=TRUE, likeDotplot=TRUE) {
    require(ggplot2)
    f = function(x) {
        pv   = attr(x, "postVar")
        cols = 1:(dim(pv)[1])
        se   = unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ord  = unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
        pDf  = data.frame(y=unlist(x)[ord],
                           ci=1.96*se[ord],
                           nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                           ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                           ind=gl(ncol(x), nrow(x), labels=names(x)))

        if(QQ) {  ## normal QQ-plot
            p = ggplot(pDf, aes(nQQ, y))
            p = p + facet_wrap(~ ind, scales="free")
            p = p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
        } else {  ## caterpillar dotplot
            p = ggplot(pDf, aes(ID, y)) + coord_flip()
            if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                p = p + facet_wrap(~ ind)
            } else {           ## different scales for random effects
                p = p + facet_grid(ind ~ ., scales="free_y")
            }
            p = p + xlab("Levels") + ylab("Random effects")
        }

        p = p + theme(legend.position="none")
        p = p + geom_hline(yintercept=0)
        p = p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
        p = p + geom_point(aes(size=1.2), colour="blue")
        return(p)
    }

    lapply(re, f)
}

Sidak = function(vecP)
#
# This function corrects a vector of probabilities for multiple testing
# using the Bonferroni (1935) and Sidak (1967) corrections.
#
# References: Bonferroni (1935), Sidak (1967), Wright (1992).
#
# Bonferroni, C. E. 1935. Il calcolo delle assicurazioni su gruppi di teste. 
# Pp. 13-60 in: Studi in onore del Professore Salvatore Ortu Carboni. Roma.
#
# Sidak, Z. 1967. Rectangular confidence regions for the means of multivariate 
# normal distributions. Journal of the American Statistical Association 62:626-633.
#
# Wright, S. P. 1992. Adjusted P-values for simultaneous inference. 
# Biometrics 48: 1005-1013. 
#
#                  Pierre Legendre, May 2007
{
k = length(vecP)

vecPB = 0
vecPS = 0

for(i in 1:k) {
   bonf = vecP[i]*k
   if(bonf > 1) bonf=1
   vecPB = c(vecPB, bonf)
   vecPS = c(vecPS, (1-(1-vecP[i])^k))
   }
#
return(list(OriginalP=vecP, BonfP=vecPB[-1], SidakP=vecPS[-1]))
}

# ========================================================================================
# === Data cleanup
# ========================================================================================

# Clean dominance rank table
names(dom.rank)[1] = 'Date'

convert.date = function(date) {
	# 01-Jan-07 to 2007-01-01
	dates = strsplit(date,'-')
	dates = lapply(dates,function(dmy) {
		d = dmy[1]
		m = dmy[2]
		y = dmy[3]
		y = if (as.numeric(y) < 50) paste0('20',y) else paste0(19,y)
		m = formatC(match(m,month.abb),width=2,flag='0')
		paste(y,m,d,sep='-')
	})
	as.Date(do.call(c,dates))
}

dom.rank$Date = convert.date(dom.rank$Date)

# Divide into before Jan 1, 2008 and after

dom.rank$PeriodID = as.numeric(dom.rank$Date >= '2008-01-01') + 1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

cleaned.data = raw.data

# Strip whitespaces from the Events
levels(cleaned.data$Event) = gsub(' ','',levels(cleaned.data$Event))

# Correct dates and times
cleaned.data$Date = as.Date(cleaned.data$Date,format='%m/%d/%Y')
cleaned.data$StateBegin = as.POSIXlt(cleaned.data$StateBegin,tz='America/Costa_Rica',format='%m/%d/%Y %I:%M:%S %p')
cleaned.data$BoutBegin = as.POSIXlt(cleaned.data$BoutBegin,tz='America/Costa_Rica',format='%m/%d/%Y %I:%M:%S %p')

# Retain only events that involve eating (E) or rejecting (R)
cleaned.data = cleaned.data[grep('[ER]',cleaned.data$Event),]

if (length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event)))) {
	## If matches both E and R, ignore that row and flip a warning

	warning('Found ',length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event))),' ambiguous (E or R) event',if (length(intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event))) > 1) 's','.',sep='')
	cleaned.data = cleaned.data[-intersect(grep('E',cleaned.data$Event),grep('R',cleaned.data$Event)),]
}


# Correct phenology IDs and disambiguate bout durations for bouts that had multiple phenology IDs
cleaned.data$CorrectedPhenologyID[cleaned.data$CorrectedPhenologyID == ''] = NA
cleaned.data$CorrectedPhenologyID = cleaned.data$CorrectedPhenologyID[,drop=TRUE]

# Clean up factor levels of phenology IDs
cleaned.data$PhenologyID = factor(cleaned.data$PhenologyID,levels=union(levels(cleaned.data$PhenologyID),levels(cleaned.data$CorrectedPhenologyID)))

# Rewrite phenology IDs with corrected phenology IDs
cleaned.data$PhenologyID[!is.na(cleaned.data$CorrectedPhenologyID)] = as.character(cleaned.data$CorrectedPhenologyID[!is.na(cleaned.data$CorrectedPhenologyID)])

# Reorder the cleaned dataset
cleaned.data = cleaned.data[order(cleaned.data$BoutBegin,cleaned.data$StateBegin),]

# Some remaining entries remain confusing/ambiguous

# DECISION D1: Amanda made the executive decision 9/24/2015 to remove the Jacquinia nervosa record from Nutella/W07AM63/13:13:02 (ambiguous entry)
cleaned.data = cleaned.data[!(cleaned.data$Name %in% 'Nutella' & cleaned.data$PhenologyID %in% 'W07AM63' & cleaned.data$BoutBegin == '2007-01-26 13:13:02' & cleaned.data$ScientificName %in% 'Jacquinia nervosa'),]

# DECISION D2: Kenny made the decision to change phenology ID of Marmite/W07AM1131/11:11:08 to W07AM1132 where the species is Zuelania guidonia (apparent user error)
cleaned.data[cleaned.data$Name %in% 'Marmite' & cleaned.data$BoutBegin %in% '2007-04-13 11:11:08' & cleaned.data$ScientificName %in% 'Zuelania guidonia',]$PhenologyID = 'W07AM1132'

# Add a bout ID to cleaned data (no such ID exists so must be inferred)
# The combination of an animal name and the bout timestamp will uniquely identify the bout
multi.bout.id = c('Name','BoutBegin')

# Extract unique combinations of Name-BoutBegin
unique.bouts = unique(cleaned.data[,multi.bout.id])
unique.bouts = unique.bouts[order(unique.bouts$BoutBegin,unique.bouts$Name),]

# Assign integers to serve as unique identifiers for each bout
unique.bouts$BoutID = 1:nrow(unique.bouts)

# Merge new bout IDs with the working data
cleaned.data = merge(cleaned.data,unique.bouts,by=multi.bout.id,all.x=TRUE,sort=FALSE)

# DECISION D3: Bouts as defined thus far are insufficient because it is assumed that each
# bout can only encompass one set of phenological conditions (=one unique PhenologyID)
# In practice this is not the case so we must identify subbouts and determine their
# durations (which do not yet exist in the data)

# Bouts are subdivided into subbouts if the bout traversed multiple patches
# The combination of animal name, phenology conditions, and the bout timestamp will uniquely identify the subbout
multi.subbout.id = c('Name','PhenologyID','BoutBegin','BoutDuration')

# Extract unique combinations of Name-PhenologyID-BoutBegin
unique.subbouts = unique(cleaned.data[,multi.subbout.id])
unique.subbouts = unique.subbouts[order(unique.subbouts$BoutBegin,unique.subbouts$Name,unique.subbouts$PhenologyID),]

# Assign subbout IDs
unique.subbouts$SubBoutID = 1:nrow(unique.subbouts)

# Merge new subbout IDs with the working data
cleaned.data = merge(cleaned.data,unique.subbouts,by=multi.subbout.id,all.x=TRUE,sort=FALSE)

# - - - - - - - - - - - - Calculate bout durations - - - - - - - - - - - - #

bouts = split(cleaned.data,cleaned.data$BoutID)

# Calculate CorrectedBoutDuration (note that this is different than the SubBoutDuration earlier, which was intended to deal with a different case)
cleaned.data$CorrectedBoutDuration = NA

# For each bout, calculate the number of unique phenology IDs per bout
uIDs = do.call(c,lapply(bouts,function(x) length(unique(x$PhenologyID))))

# Get indices of "singular" bouts (iS) that already have only one phenology ID
iS = cleaned.data$BoutID %in% as.integer(names(uIDs[which(uIDs == 1)]))

# Use the existing bout duration for cases where the bout already had only one phenology ID
cleaned.data$CorrectedBoutDuration[iS] = cleaned.data$BoutDuration[iS]

# Extract bouts that need to be split (uIDs > 1)
split.bouts = bouts[as.integer(names(uIDs[which(uIDs > 1)]))]

# For bouts that need splitting, go through and calculate subbout durations
# Sub-bout durations can be calculated on the basis of timed state behaviors that occurred within each bout
bouts.split.by.phenologies = do.call(rbind,lapply(split.bouts,function(x) {
	# Get subject, duration, and begin time of the bout (only one possible for each)
	name = unique(x$Name)
	bout.duration = unique(x$BoutDuration)
	bout.begin = unique(x$BoutBegin)

	# Duration should be calculated on the basis of states. The first state begin is not always the same as the bout begin time.
	first.state.begin = min(x$StateBegin)
	
	# Split each bout into subbouts defined by unique phenology IDs
	phenologies = do.call(rbind,lapply(split(x,x$PhenologyID[,drop=TRUE]),function(y) {
		phenology.id = unique(y$PhenologyID)
		# Calculate metadata for the subbouts (most importantly, the state begin time)
		data.frame(PhenologyID=phenology.id,Name=name,BoutDuration=bout.duration,BoutBegin=bout.begin,StateBegin=min(y$StateBegin))
	}))
	
	# Order subbouts by the state begin time
	phenologies = phenologies[order(phenologies$StateBegin),]

	# For each row (except the last), calculate the end time as the begin time of the subsequent row
	# The end time for the last row should be the bout begin time + the bout duration, taking into account the difference between the first state and the beginning of the bout
	phenologies$NextStateBegin = c(phenologies$StateBegin[2:nrow(phenologies)],bout.begin + bout.duration + as.integer(difftime(first.state.begin,bout.begin,units='secs')))
	
	# The duration of the subbout is the difference between the state begin time and the state end time calculated above
	phenologies$PhenologyDuration = as.integer(difftime(phenologies$NextStateBegin,phenologies$StateBegin,units='secs'))

	# Return subbouts
	phenologies
}))

row.names(bouts.split.by.phenologies) = NULL

# Rename StateBegin in the split bouts to SubBoutBegin
names(bouts.split.by.phenologies)[names(bouts.split.by.phenologies) %in% 'StateBegin'] = 'SubBoutBegin'

# Grab the relevant columns
phenology.durations = bouts.split.by.phenologies[,c('Name','PhenologyID','BoutBegin','PhenologyDuration','SubBoutBegin')]

# Merge the newly calculated SubBout information back into cleaned data
cleaned.data = merge(cleaned.data,phenology.durations,by=c('Name','PhenologyID','BoutBegin'),all.x=TRUE,sort=FALSE)

# Use the calculated PhenologyDuration (= sub-bout duration) as the CorrectedBoutDuration (which were left NA for bouts that needed splitting)
cleaned.data$CorrectedBoutDuration[is.na(cleaned.data$CorrectedBoutDuration)] = cleaned.data$PhenologyDuration[is.na(cleaned.data$CorrectedBoutDuration)]

# The BoutDuration is now the CorrectedBoutDuration
cleaned.data$BoutDuration = cleaned.data$CorrectedBoutDuration

# DECISION D4: We subsequently decided to recalculate durations across the board based on state times rather than bout times
#
# A bout is now defined as:
# (1) a bout starts at the time of the StateBegin that contains the first foraging event
# (2) a bout ends at the time of the StateEnd that contains the last foraging event
# (3) a new bout duration is calculated between the time of the first state-begin and last state-end containing foraging events within a bout.

# Calculate the state end time
cleaned.data$StateEnd = cleaned.data$StateBegin + cleaned.data$StateDuration

# Calculate for each subbout the duration (temporarily named "ForagingDuration")
foraging.states = do.call(rbind,lapply(split(cleaned.data,cleaned.data$SubBoutID),function(x) {
	data.frame(SubBoutID=unique(x$SubBoutID),ForagingBegin=min(x$StateBegin),ForagingEnd=max(x$StateEnd),ForagingDuration=as.numeric(difftime(max(x$StateEnd),min(x$StateBegin),units='secs')),stringsAsFactors=FALSE)
}))

# Merge this information back into cleaned.data
cleaned.data = merge(cleaned.data,foraging.states,by='SubBoutID',all.x=TRUE,sort=FALSE)

# - - - - - - - - - - - - Tabulate foraging outcomes - - - - - - - - - - - - #

# Set up the column as a factor
cleaned.data$ForagingOutcome = factor(NA,levels=c('Eat','Reject'))

# If the Event contains the code "E", the fruit was eaten
cleaned.data$ForagingOutcome[grep('E',cleaned.data$Event)] = 'Eat'

# If the Event contains the code "R", the fruit was rejected
cleaned.data$ForagingOutcome[grep('R',cleaned.data$Event)] = 'Reject'

# - - - - - - - - - - - - Assemble full dataset - - - - - - - - - - - - #

# The full dataset for analysis is named "analysis.table"

# These are the important metadata fields for each bout/subbout
multi.bout.subbout.id = c('BoutID','SubBoutID','BoutDuration','ForagingDuration','Name','PhenologyID','ScientificName','BoutBegin','SubBoutBegin')

# Each row of the dataset is equivalent to a bout or subbout
analysis.table = unique(cleaned.data[,multi.bout.subbout.id])

# SubBoutBegin does not exist for proper bouts that did not require splitting
# Set it to BoutBegin for these cases
analysis.table$SubBoutBegin[is.na(analysis.table$SubBoutBegin)] = analysis.table$BoutBegin[is.na(analysis.table$SubBoutBegin)]

# Sort this by time and then animal
analysis.table = analysis.table[order(analysis.table$SubBoutBegin,analysis.table$Name),]

# Convert the animal data to strings to avoid unexpected behaviors
analysis.table$Name = as.character(analysis.table$Name)


# These are the important metadata fields for each individual
multi.animal.id = c('Name','CurrentGroup','AgeClass','Sex','ColorVisionType')

unique.animals = unique(cleaned.data[,multi.animal.id])

unique.animals$Name = as.character(unique.animals$Name)

# Maturity is a redefined AgeClass for analysis
unique.animals$Maturity = as.character(unique.animals$AgeClass)

# Adults and subadults are "Mature"
unique.animals$Maturity[unique.animals$AgeClass %in% c('Adult','Subadult')] = 'Mature'

# Infants are excluded from analysis
unique.animals$Maturity[unique.animals$AgeClass %in% c('Infant')] = NA

# Finalize this factor
unique.animals$Maturity = factor(unique.animals$Maturity,levels=c('Mature','LargeImmature','SmallImmature'))

# Merge animal information into the dataset
analysis.table = merge(analysis.table,unique.animals,by='Name',all.x=TRUE,sort=FALSE)

# Replace BoutDuration with ForagingDuration and backup BoutDuration (see Decision D4)
analysis.table$BoutDurationOriginal = analysis.table$BoutDuration
analysis.table$BoutDuration = analysis.table$ForagingDuration

# Drop ForagingDuration
analysis.table$ForagingDuration = NULL


# Foraging outcomes must be merged back into the dataset

# Calculate number of eat/reject per subbout
foraging.outcomes = table(cleaned.data$ForagingOutcome,cleaned.data$SubBoutID)
analysis.table$NumberEaten = as.numeric(foraging.outcomes['Eat',][analysis.table$SubBoutID])
analysis.table$NumberRejected = as.numeric(foraging.outcomes['Reject',][analysis.table$SubBoutID])

# Eat/Reject rates are calculated as counts divided by the bout duration (see decision D4)
analysis.table$RateEaten = analysis.table$NumberEaten / analysis.table$BoutDuration
analysis.table$RateRejected = analysis.table$NumberRejected / analysis.table$BoutDuration

# The investigation rate is the combined eat/reject rate
analysis.table$InvestigationRate = (analysis.table$NumberEaten + analysis.table$NumberRejected) / analysis.table$BoutDuration

# The acceptance index is the fraction of inspected fruits that are eaten
analysis.table$AcceptanceIndex = analysis.table$NumberEaten / (analysis.table$NumberEaten + analysis.table$NumberRejected)


# Keep only bouts where at least one thing was eaten and it lasted more than zero seconds
analysis.table = analysis.table[as.logical(analysis.table$NumberEaten) & as.logical(analysis.table$BoutDuration) & as.logical(analysis.table$BoutDurationOriginal),]



# DECISION D5: Get rid of BH group altogether
# Get rid of large immature ranks
analysis.table = analysis.table[!analysis.table$CurrentGroup %in% 'BH',]

# DECISION D6: Get rid of a couple of plant species (too little data)
analysis.table = analysis.table[!analysis.table$ScientificName %in% c('Ficus obtusifolia','Trophis racemosa'),]
analysis.table$ScientificName = analysis.table$ScientificName[,drop=TRUE]

# Incorporate metadata on rank
analysis.table$PeriodID = as.numeric(as.Date(analysis.table$BoutBegin) >= '2008-01-01') + 1

# Names must match (get rid of spaces globally)
analysis.table$Name = gsub(' ','',analysis.table$Name)
dom.rank$IndividualID = gsub(' [FM]$','',dom.rank$IndividualID)
dom.rank$IndividualID = gsub(' ','',dom.rank$IndividualID)

# Some animals have slight inconsistencies between their names in the feeding and rank datasets
# Correct them here
analysis.table = within(analysis.table, {
	Name[Name %in% 'AlbusDumbledoor'] = 'AlbusDumbledore'
	Name[Name %in% 'Blankita'] = 'Blanquita'
	Name[Name %in% 'Artimis'] = 'Artemis'
	Name[Name %in% 'Barty'] = 'BartyCrouch'
})

# Redo this with rank and conspicuity (interaction), separate males and females

dom.rank = dom.rank[!dom.rank$IndividualID %in% sort(as.character(unique(analysis.table$Name[!analysis.table$Maturity %in% 'Mature']))),]

# Reclassify dominance ranks into categories

dom.hierarchies = split(dom.rank,list(dom.rank$GroupID,dom.rank$Sex,dom.rank$Date))

dom.hierarchies = lapply(dom.hierarchies,function(x) {
	x = x[order(x$Rank),]
	ranks = x$Rank
	if (length(ranks) == 1) {
		result = 'high'
	} else if (length(ranks) == 2) {
		result = c('high','low')
	} else if (length(ranks) == 3) {
		result = c('high','medium','low')
	} else if (length(ranks) == 4) {
		result = c('high','medium','medium','low')
	} else {
		result = character(length(ranks))
		result[1:2] = 'high'
		result[(length(ranks)-1):(length(ranks))] = 'low'
		result[result %in% ''] = 'medium'
	}
	x$RankClass = result
	x
})

dom.rank = do.call(rbind,dom.hierarchies)
rownames(dom.rank) = NULL

# Merge dominance rank info
dom.rank.merge = dom.rank[,c('IndividualID','PeriodID','Rank','RankClass')]
names(dom.rank.merge)[names(dom.rank.merge) %in% 'IndividualID'] = 'Name'

analysis.table = merge(analysis.table,dom.rank.merge,by=c('Name','PeriodID'),all.x=TRUE)

# There are some missing or inconsistent ranks following the merge. Sort out these few cases manually

# DECISION D6: Buzz's rank is always low, Lavendar and Mrs W should always be low, Albus always mid
analysis.table$RankClass[analysis.table$Name %in% c('Buzz','Lavender','MrsWeasley')] = 'low'
analysis.table$RankClass[analysis.table$Name %in% c('AlbusDumbledore')] = 'medium'

# Pull out the year
analysis.table$Year = as.numeric(substr(analysis.table$BoutBegin,1,4))

# Make a table of information about each monkey and save it to a file
monkeys = unique(analysis.table[,c('Name','CurrentGroup','Sex','AgeClass','Maturity','RankClass','Year','ColorVisionType')])
monkeys = monkeys[order(monkeys$Name,monkeys$Year),]
rownames(monkeys) = NULL
write.csv(monkeys,file='output/monkey_info.csv',row.names=FALSE)

# Incorporate metadata on fruits

# Concatenate genus + species
fruit.coloration$ScientificName = factor(paste(fruit.coloration$Genus,fruit.coloration$Species),levels=levels(analysis.table$ScientificName))
fruit.coloration = fruit.coloration[,c('ScientificName','Code','FruitColoration')]

# Rename column to SpeciesCode for clarity
names(fruit.coloration)[names(fruit.coloration) == 'Code'] = 'SpeciesCode'

# Add conspicuity information
conspicuous[as.logical(conspicuous)] = 'conspicuous'
conspicuous[conspicuous == '0'] = 'nonconspicuous'

fruit.coloration$Conspicuity = conspicuous[fruit.coloration$SpeciesCode]
fruit.coloration$Conspicuity[fruit.coloration$FruitColoration %in% 'dark'] = 'dark'

# Merge fruit metadata back into the dataset
analysis.table = merge(analysis.table,fruit.coloration,by='ScientificName',all.x=TRUE,sort=FALSE)

# Preferred sorting of columns in the dataset
column.sort.order = c('BoutID','SubBoutID','BoutBegin','SubBoutBegin','BoutDuration','BoutDurationOriginal','Name','CurrentGroup','AgeClass','Maturity','Sex','Rank','RankClass','ColorVisionType','PhenologyID','ScientificName','SpeciesCode','Conspicuity','NumberEaten','NumberRejected','RateEaten','RateRejected','InvestigationRate','AcceptanceIndex')

# Preferred sorting of rows in the dataset
row.sort.order = order(analysis.table$BoutBegin,analysis.table$SubBoutBegin,analysis.table$Name,analysis.table$PhenologyID)

analysis.table = analysis.table[row.sort.order,column.sort.order]
rownames(analysis.table) = NULL

# Use factors when appropriate and prettify their codes
analysis.table$Name = factor(analysis.table$Name)
levels(analysis.table$AgeClass)[levels(analysis.table$AgeClass)=='LargeImmature'] = 'Large Immature'
levels(analysis.table$AgeClass)[levels(analysis.table$AgeClass)=='SmallImmature'] = 'Small Immature'
analysis.table$AgeClass = factor(analysis.table$AgeClass,levels=c('Adult','Subadult','Large Immature','Small Immature','Infant'))
levels(analysis.table$Maturity)[levels(analysis.table$Maturity)=='LargeImmature'] = 'Large Immature'
levels(analysis.table$Maturity)[levels(analysis.table$Maturity)=='SmallImmature'] = 'Small Immature'
analysis.table$Maturity = factor(analysis.table$Maturity,levels=c('Mature','Large Immature','Small Immature'))
analysis.table = droplevels(analysis.table)

# Save the dataset if necessary
# write.table(analysis.table,file="analysis_table.txt",sep='\t',row.names=FALSE,quote=FALSE)

# ========================================================================================
# === Check energy dataset
# ========================================================================================

# Check the energy dataset to test for differences between conspicuous, cryptic, and dark fruits
energy = read.delim('data/energy.txt',na.strings='')
kruskal.test(Diameter~Conspicuity,data=energy)
kruskal.test(DryEnergy~Conspicuity,data=energy)
kruskal.test(FinalEnergyRate~Conspicuity,data=energy)

# ========================================================================================
# === Exploratory data visualization
# ========================================================================================

# Histogram of the feeding rate
p = ggplot(analysis.table,aes(RateEaten)) +
	geom_histogram(binwidth=0.05,position='dodge') +
	facet_wrap(~AgeClass,ncol=2) +
	xlab('Feeding Rate') + ylab('Count') + ggtitle('Feeding rate distribution')
ggsave(p,filename='output/DataExploration_feeding_rate_distribution.pdf')

# Histogram of feeding rate split by species
p = ggplot(analysis.table,aes(RateEaten)) +
	geom_histogram(bins=10,position='dodge') +
	facet_wrap(~ScientificName,scales='free') +
	xlab('Feeding Rate') + ylab('Count') + ggtitle('Feeding rate distribution (by plant)')
ggsave(p,filename='output/DataExploration_feeding_rate_distribution_by_species.pdf')

# Boxplot + jitter plot of the feeding rate by color vision phenotype
p = ggplot(analysis.table,aes(ColorVisionType,RateEaten,color=Sex)) +
	geom_boxplot(outlier.shape=21) +
	geom_jitter(alpha=0.2) +
	facet_wrap(~AgeClass,nrow=1) +
	theme(axis.text.x=element_text(angle=-30,hjust=0),panel.grid=element_blank()) +
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ColorVisionType)-0.5,1),color='white',size=0.5) +
	xlab('Color Vision Type') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type')
ggsave(p,filename='output/DataExploration_feeding_rate_by_color_vision_type.pdf')

# The same plot by sex (dichromats only)
p = ggplot(analysis.table[analysis.table$ColorVisionType %in% 'Dichromat',],aes(Sex,RateEaten,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(legend.title = element_blank(),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$Sex)-0.5,1),color='white',size=0.5) + 
	xlab('Sex') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by sex')
ggsave(p,filename='output/DataExploration_feeding_rate_by_sex_dichromats_only.pdf')

# Boxplot of feeding rate by color vision phenotype and fruit species
p = ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(ScientificName,RateEaten,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ScientificName)-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Plant Species') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type and plant species')
ggsave(p,filename='output/DataExploration_feeding_rate_by_color_vision_type_and_species.pdf')

# Boxplot of feeding rate by color vision phenotype and fruit coloration
p = ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(Conspicuity,RateEaten,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(panel.grid=element_blank()) + geom_vline(xintercept=seq(1.5,length(unique(analysis.table$Conspicuity))-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Fruit Coloration') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by color vision type and fruit coloration')
ggsave(p,filename='output/DataExploration_feeding_rate_by_color_vision_type_and_fruit_coloration.pdf')

# Scatterplot of feeding rate by the acceptance index
p = ggplot(analysis.table,aes(round(AcceptanceIndex,2),RateEaten)) + 
	geom_point() + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~Name,ncol=floor(sqrt(nlevels(analysis.table$Name)))) + 
	theme(axis.text.x=element_text(angle=-45,size=6),axis.text.y=element_text(size=6)) + 
	xlab('Acceptance Index') + ylab('Feeding Rate') + ggtitle('Feeding efficiency by acceptance index (by animal)')
ggsave(p,filename='output/DataExploration_feeding_rate_by_acceptance_index.pdf')

# Histogram of the acceptance index distribution by age class
p = ggplot(analysis.table,aes(AcceptanceIndex)) + 
	geom_histogram(binwidth=0.05,position='dodge') + 
	facet_wrap(~AgeClass,ncol=2) + 
	xlab('Acceptance Index') + ylab('Count') + ggtitle('Acceptance index distribution')
ggsave(p,filename='output/DataExploration_acceptance_index_distribution.pdf')

# Boxplot + jitter plot of the acceptance index by color vision phenotype
p = ggplot(analysis.table,aes(ColorVisionType,AcceptanceIndex,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(axis.text.x=element_text(angle=-30,hjust=0),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ColorVisionType)-0.5,1),color='white',size=0.5) + 
	xlab('Color Vision Type') + ylab('Acceptance Index') + ggtitle('Acceptance index by color vision type')
ggsave(p,filename='output/DataExploration_acceptance_index_by_color_vision_type.pdf')

# Boxplot + jitter plot of the acceptance index by sex (dichromats only)
p = ggplot(analysis.table[analysis.table$ColorVisionType %in% 'Dichromat',],aes(Sex,AcceptanceIndex,color=Sex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(legend.title = element_blank(),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$Sex)-0.5,1),color='white',size=0.5) + 
	xlab('Sex') + ylab('Acceptance Index') + ggtitle('Acceptance index by sex')
ggsave(p,filename='output/DataExploration_acceptance_index_by_sex_dichromats_only.pdf')

# Boxplot of the acceptance index by color vision phenotype and fruit species
p = ggplot(analysis.table[!is.na(analysis.table$Maturity),],aes(ScientificName,AcceptanceIndex,color=ColorVisionType)) + 
	geom_boxplot(outlier.shape=21) + 
	facet_wrap(~Maturity,nrow=nlevels(analysis.table$Maturity)) + 
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5),panel.grid=element_blank()) + 
	geom_vline(xintercept=seq(1.5,nlevels(analysis.table$ScientificName)-0.5,1),color='white',size=0.5) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Plant Species') + ylab('Acceptance Index') + ggtitle('Acceptance index by color vision type and plant species')
ggsave(p,filename='output/DataExploration_acceptance_index_by_color_vision_type_and_species.pdf')

# Boxplot + jitter plot of bout duration by age class
p = ggplot(analysis.table,aes(AgeClass,BoutDuration,color=AcceptanceIndex)) + 
	geom_boxplot(outlier.shape=21) + 
	geom_jitter(alpha=0.2) + 
	facet_wrap(~ColorVisionType+Sex) + 
	scale_color_continuous(name='Acceptance Index') + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	xlab('Age Class') + ylab('Length of Bout') + ggtitle('Lengths of bouts by age class')
ggsave(p,filename='output/DataExploration_bout_duration_by_age_class.pdf')

# Scatterplot of the bout duration by the acceptance index
p = ggplot(analysis.table,aes(AcceptanceIndex,BoutDuration,color=ColorVisionType)) + 
	geom_point(alpha=0.5) + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~AgeClass,nrow=1) + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Acceptance Index') + ylab('Length of Bout') + ggtitle('Lengths of bouts by acceptance index')
ggsave(p,filename='output/DataExploration_bout_duration_by_acceptance_index.pdf')

# Scatterplot of the bout duration by the rejection rate ("pickiness")
p = ggplot(analysis.table,aes(RateRejected,BoutDuration,color=ColorVisionType)) + 
	geom_point(alpha=0.5) + 
	geom_smooth(method=lm,se=FALSE) + 
	facet_wrap(~AgeClass+Sex,nrow=1) + 
	theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
	scale_color_discrete(name='Color Vision Type') + 
	xlab('Rate Rejected') + ylab('Length of Bout') + ggtitle('Lengths of bouts by pickiness')
ggsave(p,filename='output/DataExploration_bout_duration_by_pickiness.pdf')

# ========================================================================================
# === Finalize the dataset
# ========================================================================================

# Drop infants from the analysis
analysis.table = analysis.table[!is.na(analysis.table$Maturity),]

# ========================================================================================
# === Explore data distributions
# ========================================================================================

# NOTE: SPLIT BY SCIENTIFIC NAME
shapiro.results = data.frame(
	ScientificName = unique(analysis.table$ScientificName),
	NoTransform = round(do.call(c,lapply(1:nlevels(analysis.table$ScientificName),function(x) shapiro.test(analysis.table$RateEaten[analysis.table$ScientificName %in% unique(analysis.table$ScientificName)[x]])$p.value)),5),
	Transform = round(do.call(c,lapply(1:nlevels(analysis.table$ScientificName),function(x) shapiro.test(log(analysis.table$RateEaten)[analysis.table$ScientificName %in% unique(analysis.table$ScientificName)[x]])$p.value)),5)
)

shapiro.results$reg.passed = shapiro.results$NoTransform >= 0.05
shapiro.results$log.passed = shapiro.results$Transform >= 0.05

print(shapiro.results)

RateEaten = analysis.table$RateEaten

# Add one to allow analysis of nonzero nonnegative distributions
RateEatenPlusOne = RateEaten + 1

# Calculate parameters for Poisson and gamma distributions
po = fitdistr(RateEatenPlusOne,'Poisson')
ga = fitdistr(RateEatenPlusOne,'gamma')

# Test gaussian, log-normal, exponential, and gamma distributions
pdf(file='output/DistributionTests.pdf')
	layout(t(matrix(1:4,nrow=2)))
	qqp(RateEaten,distribution='norm',ylab='Feeding rate',xlab='Quantiles',main='Gaussian distribution')
	qqp(RateEaten,distribution='lnorm',ylab='Feeding rate',xlab='Quantiles',main='Log-normal distribution')
	qqp(RateEaten,distribution='exp',ylab='Feeding rate',xlab='Quantiles',main='Exponential distribution')
	qqp(RateEatenPlusOne,distribution='gamma',shape=ga$estimate[[1]],rate=ga$estimate[[2]],ylab='Feeding rate + 1',xlab='Quantiles',main='Gamma distribution')
dev.off()

AcceptanceIndex = analysis.table$AcceptanceIndex

# Add one to allow analysis of nonzero nonnegative distributions
AcceptanceIndexPlusOne = AcceptanceIndex + 1

# Calculate parameters for Poisson and gamma distributions
po = fitdistr(AcceptanceIndexPlusOne,'Poisson')
ga = fitdistr(AcceptanceIndexPlusOne,'gamma')

# Test gaussian, log-normal, exponential, and gamma distributions
pdf(file='output/DistributionTestsAI.pdf')
	layout(t(matrix(1:4,nrow=2)))
	qqp(AcceptanceIndex,distribution='norm',ylab='Acceptance index',xlab='Quantiles',main='Gaussian distribution')
	qqp(AcceptanceIndex,distribution='lnorm',ylab='Acceptance index',xlab='Quantiles',main='Log-normal distribution')
	qqp(AcceptanceIndex,distribution='exp',ylab='Acceptance index',xlab='Quantiles',main='Exponential distribution')
	qqp(AcceptanceIndexPlusOne,distribution='gamma',shape=ga$estimate[[1]],rate=ga$estimate[[2]],ylab='Acceptance index + 1',xlab='Quantiles',main='Gamma distribution')
dev.off()

# ========================================================================================
# === Model fitting
# ========================================================================================

# Drop extraneous factor levels prior to model-building
analysis.table = droplevels(analysis.table)

# - - - - - Main model - - - - - #
full.model.no.rank = lmer(log(RateEaten) ~ Maturity * ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=analysis.table)

# - - - - - Mature-only model - - - - - #
mature.model.no.rank = lmer(log(RateEaten) ~ ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Maturity %in% 'Mature')))

# - - - - - Rank model - - - - - #
mature.model.with.rank = lmer(log(RateEaten) ~ ColorVisionType * Conspicuity * RankClass + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Maturity %in% 'Mature')))

# - - - - - Maturity model (immatures only) - - - - - #
immature.model.with.maturity  = lmer(log(RateEaten) ~ Maturity * ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (1|PhenologyID),data=droplevels(subset(analysis.table,!Maturity %in% 'Mature')))
sex.model  = lmer(log(RateEaten) ~ Maturity * Sex + Conspicuity + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=subset(analysis.table,ColorVisionType %in% 'Dichromat'))
mature.female.model.with.rank = lmer(log(RateEaten) ~ ColorVisionType * Conspicuity * RankClass + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Maturity %in% 'Mature' & Sex %in% 'Female')))
mature.male.model.with.rank = lmer(log(RateEaten) ~ Conspicuity * RankClass + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=droplevels(subset(analysis.table,Maturity %in% 'Mature' & Sex %in% 'Male')))
female.model.no.rank = lmer(log(RateEaten) ~ Maturity * ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=subset(analysis.table,Sex %in% 'Female'))

# Because these names are getting long and out of hand, use the following shorthand codes for major models
# fm : full.model.no.rank
# rm : mature.model.with.rank
# im : immature.model.with.maturity
# sm : sex.model

# Give aliases using these shorthand names
f.m = full.model.no.rank
r.m = mature.model.with.rank
i.m = immature.model.with.maturity
s.m = sex.model

# Calculate least-square means (variable names use shorthand prefixes and "ls" suffix)
fm.ls = lsmeansLT(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous'))),c('Maturity:ColorVisionType'))
#rm.ls = lsmeansLT(lmer(log(RateEaten) ~ ColorVisionType * RankClass + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous' & Maturity %in% 'Mature'))),c('ColorVisionType:RankClass'))
im.ls = lsmeansLT(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,!Maturity %in% 'Mature' & Conspicuity %in% 'conspicuous'))),'Maturity:ColorVisionType')
sm.ls = lsmeansLT(lmer(log(RateEaten) ~ Maturity + Sex + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=subset(analysis.table,ColorVisionType %in% 'Dichromat')),'Sex')

# Calculate differences of least-square means (variable names use shorthand prefixes and "dm" suffix)
fm.dm = difflsmeans(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous'))),c('Maturity:ColorVisionType'))
#rm.dm = difflsmeans(lmer(log(RateEaten) ~ ColorVisionType * RankClass + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous' & Maturity %in% 'Mature'))),c('ColorVisionType:RankClass'))
im.dm = difflsmeans(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName)  + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,!Maturity %in% 'Mature' & Conspicuity %in% 'conspicuous'))),'Maturity:ColorVisionType')
sm.dm = difflsmeans(lmer(log(RateEaten) ~ Maturity + Sex  + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=subset(analysis.table,ColorVisionType %in% 'Dichromat')),'Sex')

# Redo models above with Conspicuity as interaction term when allowed (some models will throw an error)

fm.ls = lsmeansLT(lmer(log(RateEaten) ~ Maturity * ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,TRUE))),c('Maturity:ColorVisionType:Conspicuity'))
#rm.ls = lsmeansLT(lmer(log(RateEaten) ~ ColorVisionType * RankClass * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Maturity %in% 'Mature'))),c('ColorVisionType:RankClass:Conspicuity'))
sm.ls = lsmeansLT(lmer(log(RateEaten) ~ Maturity + Sex + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=subset(analysis.table,ColorVisionType %in% 'Dichromat')),'Sex')

fm.dm = difflsmeans(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous'))),c('Maturity:ColorVisionType'))
# rm.dm = difflsmeans(lmer(log(RateEaten) ~ ColorVisionType * RankClass + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous' & Maturity %in% 'Mature'))),c('ColorVisionType:RankClass'))
sm.dm = difflsmeans(lmer(log(RateEaten) ~ Maturity + Sex  + (1|Name) + (1|ScientificName) + (1|PhenologyID),data=subset(analysis.table,ColorVisionType %in% 'Dichromat')),'Sex')

# Do a quality check on the full model (e.g., check for homoscedasticity)
pdf(file='output/qualitycheck.pdf',width=6,height=6,useDingbats=FALSE)
	layout(matrix(1:4,ncol=2))

	# QQ plot
	qqnorm(resid(full.model.no.rank),main='QQ plot')
	qqline(resid(full.model.no.rank),col='red')

	# Residual plot
	plot(resid(full.model.no.rank),main='Residual plot',ylab='Residuals')

	# Autocorrelation
	acf(resid(full.model.no.rank),main='Autocorrelation plot')

	# Fitted vs. residuals
	d = data.frame(f=fitted(full.model.no.rank),r=resid(full.model.no.rank))
	plot(r~f,d,main='Fitted vs. residuals',xlab='Fitted',ylab='Residuals')
	j = order(d$f) # indices marking the order of points on x-axis
	lines(d$f[j],predict(loess(r~f,d))[j],col='red',lwd=2)
	abline(h=0,lty=2)
dev.off()

# Post-hoc analysis

# Main model
main.diffmeans = difflsmeans(lmer(log(RateEaten) ~ Maturity + ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,TRUE))),c('ColorVisionType:Conspicuity'))$diffs.lsmeans.table

# Planned post-hoc comparisons (dichromats vs. trichromats for each conspicuity level)
i = c(
	"ColorVisionType:Conspicuity  Dichromat conspicuous -  Trichromat conspicuous",
	"ColorVisionType:Conspicuity  Dichromat nonconspicuous -  Trichromat nonconspicuous",
	"ColorVisionType:Conspicuity  Dichromat dark -  Trichromat dark" 
)

main.diffmeans.subset = main.diffmeans[i,]

# Adjust p-values for multiple comparisons
main.diffmeans.subset$holm = p.adjust(main.diffmeans.subset$`p-value`,'holm')

# ========================================================================================
# === Model reporting
# ========================================================================================

# Write results to file
sink(file='output/model_fitting.log')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- "Main model"\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(f.m))
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - -  Least-square means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(fm.ls)
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - Difference of means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(fm.dm)
	cat('\n')

	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- "Dominance rank model" (adults only)\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(r.m))
	cat('\n')

	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- "Maturity model" (immatures only)\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(i.m))
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - -  Least-square means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(im.ls)
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - Difference of means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(im.dm)
	cat('\n')

	cat('# ----------------------------------------------------------------------------------------\n')
	cat('# --- "Sex model" (dichromats only)\n')
	cat('# ----------------------------------------------------------------------------------------\n')
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - -  ANOVA results  - - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(Anova(s.m))
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - -  Least-square means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(sm.ls)
	cat('\n')
	cat('# - - - - - - - - - - - - - - - - - Difference of means - - - - - - - - - - - - - - - - - \n')
	cat('\n')
	print(sm.dm)
	cat('\n')
sink()

# ========================================================================================
# === Prepare figures for publication
# ========================================================================================

# Run the following models for figure plotting

# Least-square means for the three maturity classes (conspicuous fruits only)
maturity.lsmeans = lsmeansLT(lmer(log(RateEaten) ~ Maturity * ColorVisionType + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,Conspicuity %in% 'conspicuous'))),c('Maturity'))$lsmeans.table

# Least-square means for the two color vision types (full dataset)
di.tri.lsmeans = lsmeansLT(lmer(log(RateEaten) ~ Maturity + ColorVisionType * Conspicuity + (ColorVisionType|Name) + (ColorVisionType|ScientificName) + (ColorVisionType|PhenologyID),data=droplevels(subset(analysis.table,TRUE))),c('ColorVisionType:Conspicuity'))$lsmeans.table


# Order plant species by feeding rate
sorted.scientific.name = names(sort(tapply(analysis.table$RateEaten,analysis.table$ScientificName,mean)))

conspicuity = unique(analysis.table[,c('ScientificName','SpeciesCode','Conspicuity')])

conspicuity = conspicuity[match(sorted.scientific.name,conspicuity$ScientificName),]
conspicuity$ScientificName = factor(conspicuity$ScientificName,levels=sorted.scientific.name)

rownames(conspicuity) = NULL

# Calculate max and min x-axis coordinates for rectangles
conspicuity$xmin = seq(0.5,nrow(conspicuity) - 0.5,1)
conspicuity$xmax = seq(1.5,nrow(conspicuity) + 0.5,1)

# Sort factor levels
conspicuity$Conspicuity = factor(conspicuity$Conspicuity,levels=c('conspicuous','nonconspicuous','dark'))
di.tri.lsmeans$Conspicuity = factor(di.tri.lsmeans$Conspicuity,levels=c('conspicuous','nonconspicuous','dark'))
maturity.lsmeans$Maturity = factor(maturity.lsmeans$Maturity,levels=c('Mature','Large Immature','Small Immature'))

# Rename "nonconspicuous" to "cryptic"
levels(conspicuity$Conspicuity)[levels(conspicuity$Conspicuity) == 'nonconspicuous'] = 'cryptic'
levels(di.tri.lsmeans$Conspicuity)[levels(di.tri.lsmeans$Conspicuity) == 'nonconspicuous'] = 'cryptic'

theme_set(theme_classic(base_size = 8))

# Recode some variables and calculate SE intervals
maturity.lsmeans = within(maturity.lsmeans,{
	LogEstimate = Estimate
	UpperSE = Estimate + `Standard Error`
	LowerSE = Estimate - `Standard Error`
	Estimate = exp(Estimate)
	LowerCI = `Lower CI`
	UpperCI = `Upper CI`
	`Lower CI` = exp(`Lower CI`)
	`Upper CI` = exp(`Upper CI`)
})

# Recode some variables and calculate SE intervals
di.tri.lsmeans = within(di.tri.lsmeans,{
	LogEstimate = Estimate
	UpperSE = Estimate + `Standard Error`
	LowerSE = Estimate - `Standard Error`
	Estimate = exp(Estimate)
	LowerCI = `Lower CI`
	UpperCI = `Upper CI`
	`Lower CI` = exp(`Lower CI`)
	`Upper CI` = exp(`Upper CI`)
})

p = ggplot(maturity.lsmeans,aes(Maturity,LogEstimate,color=Maturity)) + 
	geom_point(size=2) + 
#	geom_errorbar(aes(ymin=lsmean-SE,ymax=lsmean+SE),width=0.5,size=0.5) + 
	geom_errorbar(aes(ymin=LowerSE,ymax=UpperSE),width=0.5,size=0.5) + 
	geom_errorbar(aes(ymin=LowerCI,ymax=UpperCI),width=0.5,size=0.2,linetype=3) + 
	xlab('Age Class') + ylab('log(Feeding Rate)') + 
	theme(legend.position='none') + 
	scale_color_manual(values=c('#000000','#000000','#000000'))
ggsave(p,filename='output/fig_lsmeans_maturity.pdf',width=2.4,height=2.4,useDingbats=FALSE)

p = ggplot(di.tri.lsmeans,aes(ColorVisionType,LogEstimate,color=ColorVisionType)) + 
	geom_point(size=2) + 
	geom_errorbar(aes(ymin=LowerSE,ymax=UpperSE),width=0.5,size=0.5) + 
	geom_errorbar(aes(ymin=LowerCI,ymax=UpperCI),width=0.5,size=0.2,linetype=3) + 
	ylim(c(min(di.tri.lsmeans$LowerCI),max(di.tri.lsmeans$UpperCI)+0.02)) +
	facet_wrap(~Conspicuity,nrow=1) +
	ylab('log(Feeding Rate)') + 
	theme(legend.position='bottom',
		axis.title.x=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank()) +
	scale_color_manual(name='Color vision type',values=c('#7fbf7b','#af8dc3'))
ggsave(p,filename='output/fig_lsmeans_colorvision_conspicuity.pdf',width=3.0,height=2.4,useDingbats=FALSE)

# Create clone of analysis.table with reordered plant species
vis.table = analysis.table
vis.table$ScientificName = factor(vis.table$ScientificName,levels=sorted.scientific.name)

p = ggplot() +
	geom_rect(data=conspicuity,aes(NULL,xmin=xmin,xmax=xmax,ymin=-0.5,ymax=1.5,fill=Conspicuity)) +
	scale_fill_manual(values=c('#cccccc','#eeeeee','#aaaaaa'),name='Fruit Conspicuity') +
	geom_boxplot(data=vis.table,aes(ScientificName,RateEaten,color=ColorVisionType),outlier.shape=21) +
	geom_vline(xintercept=seq(1.5,nlevels(vis.table$ScientificName)-0.5,1),color='white',size=0.5) +
	scale_color_manual(name='Color Vision Type',values=c('#fdc086','#af8dc3'),labels=c('Dichromat','Trichromat')) +
	scale_y_continuous(breaks=c(0,1)) +
	coord_cartesian(ylim=c(0,1)) +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,face='italic'),panel.grid=element_blank()) +
	xlab('Plant Species') +
	ylab('Feeding Rate')
ggsave(p,filename='output/fig_feeding.rate.by.fruit.species.pdf',width=7,height=4)

# ========================================================================================
# === Prepare tables for publication
# ========================================================================================

# Summarize bout counts, sum durations, feeding counts (events), and reject counts for each fruit species
bouts.by.species = melt(table(analysis.table$ScientificName))
duration.by.species = melt(tapply(analysis.table$BoutDuration,analysis.table$ScientificName,sum))
events.by.species = melt(tapply(analysis.table$NumberEaten,analysis.table$ScientificName,sum))
rejected.by.species = melt(tapply(analysis.table$NumberRejected,analysis.table$ScientificName,sum))

# Tally up summary info for each bout (by plant species)

# Start with number of bouts
bout.summary = bouts.by.species
names(bout.summary) = c('species','bouts')

# Add duration, feeding events, and reject events
bout.summary = data.frame(
	bout.summary,
	duration = duration.by.species$value ,
	feeding_events = events.by.species$value ,
	reject_events = rejected.by.species$value
)

# Order of plants given in the table (see Decision D6 for plants removed)
plant.order = c("Sciadodendron excelsum", "Cordia guanacastensis", "Cordia panamensis", 
"Bursera simaruba", "Diospyros salicifolia", "Sloanea terniflora", 
"Erythroxylum havanense", "Vachellia collinsii", "Casearia arguta", 
"Casearia sylvestris", "Muntingia calabura", "Zuelania guidonia", 
"Ficus cotinifolia", "Ficus hondurensis", "Ficus morazaniana", 
"Ficus ovalis", "Maclura tinctoria", "Karwinskia calderoni", 
"Krugiodendron ferreum", "Genipa americana", "Randia monantha", 
"Randia thurberi", "Allophylus occidentalis", "Dipterodendron costaricense", 
"Manilkara chicle", "Simarouba glauca", "Jacquinia nervosa")

# Order summary accordingly
bout.summary = bout.summary[match(plant.order,bout.summary$species),]

# Write info to file
write.csv(bout.summary,file='output/bout_summary.csv')

# ========================================================================================
# === Coda
# ========================================================================================

message('Analysis complete!')
