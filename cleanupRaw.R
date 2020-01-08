# THIS FILE CLEANS AND SUMMARIZES WETLAND AND INFIELD DATA FROM ECOLOGICS

library(tidyverse)

#Load raw data
load('rawData.Rdata')

#NOTES ON COLUMNS
# BBID = Unique ID for speciment
# BTID = Unique ID for trap
# BLID = Unique site (can have multiple traps)
# Notes on Replicates
#A,B,C = Largely in S Calgary area, used in 2015 and 2016 year, all ditch, B & C occasionally used
#WXX = wetland exclusively, largely in S Calgary 2016 (W = wetland, number = m from wetland) 
#DBV, DPF = ditch sampling from 2017, wider geographic sampling (D = ditch, BV = blue vane, PF = pitfall) - ask Hailey about processing
#WCCXX, WBVXX, WPFXX = anything from into a field, NOT just wetland, 2017 (CC = coloured cup, BV = blue vane, PF = pitfall) 

#Timing of samples:
# 2015 all ditch,
# 2016 ditch + infield
# 2017 ditch  + infield
# 2018 restoration chronosequence + ditch, all WBV,WPF,WCC are all 25 m from wetland

#Notes on locationType 
#locationType  = ditch, infield_wetland, infield_pivot, infield_control
#   within field, positioned at distances from features (control has no distance)

#Other notes:
#Check whether zeros are real! Later material may not have been entered. Dialictus missing from 2018.
#Some species haven't been IDed (e.g. Lasioglossum_spp. 6)
#Arthropod table needs cleaning (prior to joining) due to typing errors
#GDD is a more relevant measure of time than Jday (this is what the bees would "feel" rather than DOY)

# Clean up site data ------------------------------------------------------

site <- site %>% select(-aliasBLID1:-locality,-country,-created_at:-expt_yield_Lakeland) %>% 
  mutate_if(is.factor,as.character)

# Clean up trap data ------------------------------------------------------

# #PROBLEM: SOME BLIDS FROM TRAP DF ARE NOT PRESENT IN SITE DF (either in BLID or aliases)
# badBLID <- trap %>% select(BTID:trapType) %>% filter(!(BLID %in% site$BLID),!is.na(BLID)) %>%
#   pull(BLID) %>% unique()

#Bad BLID values:
# 20427 - duplicate of 24027, nothing in site or arth df: REMOVE FROM TRAP DF
# 23460 - should be 13460: CHANGE IN TRAP AND ARTH DF
# 13395 - should be 11395: CHANGE IN TRAP AND ARTH DF
# 13048 - should be 10348: CHANGE IN TRAP AND ARTH DF
# 15280 - duplicate of 15208, nothing in arth df: REMOVE FROM TRAP DF
# 13752 - should be 13751: CHANGE IN TRAP DF
# 25095 - 25091 or 25035, but were recorded on different days (June 22, 23), nothing in site or arth df: REMOVE FROM TRAP DF
# 13903 - duplicate of 13902: REMOVE FROM TRAP DF
# 27029 - present in trap table, but not arth table: REMOVE FROM TRAP DF
#allBad <- c(20427,23460,13395,13048,15280,13752,25095,13903,27029)

removeBLID <- c(20427,15280,25095,13903,27029) #BLID values to remove from trap df
#BLID values to change
changeBLID <- data.frame(BLID=c(23460,13395,13048,13752),to=c(13460,11395,10348,13751)) 

# Notes: 
# "Coloured cups" consist of blue, white & yellow cup, but were sometimes all amalgamated into single category, so 3x effort makes sense.
trap <- trap %>% mutate_if(is.factor,as.character) %>% #Convert all factors to characters
  filter(startMonth!=0) %>% #Removes pass 0 (where traps are deployed)
  filter(!(BLID %in% removeBLID)) %>% #Removes bad BLIDs
  mutate(change=(BLID %in% changeBLID$BLID)) %>% #Marks BLIDs for changes
  left_join(changeBLID,by='BLID') %>% rowwise() %>% 
  mutate(BTID=ifelse(change,gsub(as.character(BLID),to,BTID),BTID),BLID=ifelse(change,to,BLID)) %>% 
  ungroup() %>% select(-change,-to,-startNESWPhoto_DSC,-endNESWPhoto_DSC) %>% #Cleanup
  filter(!is.na(trapType),BLID!=0) %>% #Get rid of NA and 0 rows
  #Fixes typos
  mutate(locationType=gsub('Ditch','ditch',locationType)) %>% #Fix typo in locationtype
  #Fixes reversal in year labels
  mutate(BTID=gsub('-2007','-2017',BTID)) %>% mutate(BTID=gsub('-20107','-2017',BTID)) %>% 
  mutate(BTID=gsub('-2107','-2017',BTID)) %>% mutate(BTID=gsub('-2071','-2017',BTID)) %>% 
  mutate(replicate=gsub('WCCC','WCC',replicate),BTID=gsub('WCCC','WCC',BTID)) %>%
  mutate(replicate=gsub('WCC2R1','WCC25R1',replicate),BTID=gsub('WCC2R1','WCC25R1',BTID)) %>% 
  mutate(collector=case_when(collector=='J. VIckruck' ~ 'J. Vickruck',
                             collector=='M.Gavin' ~ 'M. Gavin',
                             TRUE ~ collector)) %>% 
  # "Coloured cups" consist of blue, white & yellow cup, but were sometimes all amalgamated into single category, so 3x effort makes sense.
  #3x deployedhours for "coloured cups" category (3 colour types)
  mutate(deployedhours=case_when(trapType=='Coloured Cups' ~ deployedhours*3,TRUE ~ deployedhours)) %>% 
  #Adds deployedhours from separate cups to each other
  group_by(BTID) %>% mutate(deployedhours=sum(deployedhours)) %>% ungroup() %>% 
  filter(trapType!='Yellow Cup'& trapType!='White Cup') %>%  #Removes yellow and white cup rows
  mutate(trapType=gsub('Blue Cup','Coloured Cups',trapType)) %>%  #Changes blue cup rows to "coloured cup"
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>% mutate(replicate=gsub('W[BWY]C','WCC',replicate)) %>% 
  #Convert -99 in mowed data to NA
  mutate(adjMowed=ifelse(adjMowed==-99,NA,adjMowed),oppMowed=ifelse(oppMowed==-99,NA,oppMowed)) %>% 
  #Beginning to reorganize trap dataframe. Goal:
  #dist = distance from some "source" feature (usually wetland/SNL if in field)
  #distFrom = type of "source" feature
  rename('trapLoc'='locationType') %>% #trapLoc = location of trap (what feature is the trap sitting in?)
  mutate(trapLoc=case_when( #Deal with empty trap locations from 2016
    trapLoc=='' & grepl('W[027]',replicate) ~ 'wetland',
    trapLoc=='' & grepl('A',replicate) ~ 'ditch',
    TRUE ~ trapLoc
  )) %>% 
  mutate(dist=as.numeric(gsub('[A-Z]','0',replicate)),distFrom=trapLoc) %>% #Create distance and distFrom column
  mutate(trapLoc=case_when( #Sets trapLoc to adjCrop if trap was placed at a distance into some feature
    (distFrom=='wetland' | distFrom=='pivot') & (startYear==2016 | startYear==2017) ~ adjCrop, 
    trapLoc=='control' ~ adjCrop, #I think "control" just means a canola field without any wetland/pivot nearby (i.e. dist = Inf)
    TRUE ~ trapLoc
  )) %>% 
  mutate(trapType=case_when(
    trapType=='infield' & replicate=='WBV' ~ 'Blue Vane',
    trapType=='infield' & replicate=='WPF' ~ 'Pit Fall',
    trapType=='infield' & replicate=='WCC' ~ 'Coloured Cups',
    TRUE ~ trapType
  )) 

# #Looks OK
# trap %>% select(startYear,trapLoc,replicate,dist,distFrom,trapType) %>% distinct() %>% as.data.frame()
# 
# trap %>% group_by(trapType,endYear) %>% summarize(nPasses=n(),daysDeployed=round(sum(deployedhours)/24),
#         nLocs=length(unique(BLID))) %>% ungroup() %>%
#   mutate(daysPerLoc=daysDeployed/nLocs)

# Clean up arthropod data -------------------------------------------------

#Questions for Paul:
# What does ascension_date mean? Date that it was uploaded to database? Some dates are in the future (2022); does this matter?
# More info on missing dates for spiders: 12490-DPF-2017, 12490-X-DPF-2017, 17774-DPF-2017, 25070-DPF-2017
# Are all spiders within these from the same pass? If not, how should these be dealt with?
# Potential Excel-related problem: dates on many BTIDs have been "dragged down" resulting in increased years. Easy to fix when year is >2019, but could be a lot of material with erronious 2019 labels on them (should be 2018)

temp <- arth %>% mutate_if(is.factor,as.character) %>% #Convert all factors to characters
  mutate(change=(BLID %in% changeBLID$BLID)) %>%  #Marks BLIDs for changes
  left_join(changeBLID,by='BLID') %>% rowwise() %>% 
  mutate(BTID=ifelse(change,gsub(as.character(BLID),to,BTID),BTID),BLID=ifelse(change,to,BLID)) %>% 
  ungroup() %>% select(-change,-to) %>% #Cleanup
  select(-specimenID) %>% #Get rid of specimenID (all zeros) and BBID (other specimen ID column?)
  filter(BLID!=0,!is.na(BLID),nchar(arthOrder)>0) %>%  #Get rid of NA and 0 rows, and what look like "test" rows
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>% 
  #Fixes typos:
  mutate(BTID=gsub('-DBV-DBV','-DBV',BTID)) %>% mutate(BTID=gsub('-a-','-A-',BTID)) %>% 
  mutate(BTID=gsub('-WPF26-','-WPF25-',BTID)) %>% 
  mutate(BTID=gsub('-2007','-2017',BTID)) %>% mutate(BTID=gsub('-20107','-2017',BTID)) %>% 
  mutate(BTID=gsub('-2107','-2017',BTID)) %>% mutate(BTID=gsub('-2071','-2017',BTID)) %>% 
  mutate(BTID=gsub('WCCC','WCC',BTID)) %>% mutate(BTID=gsub('WCC2R1','WCC25R1',BTID)) %>% 
  mutate(BTID=gsub('14455-8-DPF2017','14455-8-DPF-2017',BTID)) %>%
  mutate(BTID=gsub('22576-6WPF25-2017','22576-6-WPF25-2017',BTID)) %>%
  mutate(BTID=gsub('31217-40DBV-2018','31217-4-DBV-2018',BTID)) %>%
  mutate(BTID=gsub('13902-7-DBV2017','13902-7-DBV-2017',BTID)) %>%
  mutate(BTID=gsub('_','-',BTID)) %>% mutate(BTID=gsub('115208','11520',BTID)) %>%
  mutate(BTID=gsub('15208-1-PF-2017','15208-1-DPF-2017',BTID)) %>%
  #Years that are greater than sampling period - Danielle's bees from 2019 haven't been added yet
  mutate(BTID=gsub('-20(19|20|21|22|23|24|24|25)','-2018',BTID)) %>% #Change all dates >2018 to 2018
  mutate(BTID=case_when( #BTIDs that weren't fixed by the date swap
    BTID=='11520-8-DBV-2018' ~ '11520-8-DBV-2017',
    BTID=='15208-7-DBV-2018' ~ '15208-7-DBV-2017',
    grepl('25238-7-DBV-2018',BTID) ~ '25238-7-DBV-2017',
    grepl('25224-6-DBV-2018',BTID) ~ '25224-6-DBV-2017',
    grepl('25091-5-DBV-2018',BTID) ~ '25091-5-DBV-2017',
    grepl('31128-1-WBV-2018',BTID) ~ '31028-1-WBV-2018',
    TRUE ~ BTID
  )) %>% 
  #Deals with missing pass labels (essentially guesses)
  mutate(BTID=gsub('11520-DBV-2017','11520-8-DBV-2017',BTID)) %>% #Should be 8, based on other B. insularis
  mutate(BTID=gsub('13825-DBV-2017','13825-7-DBV-2017',BTID)) %>% #Should be 7, based on other B. nevadensis
  mutate(BTID=gsub('17823-DBV-2017','17823-7-DBV-2017',BTID)) %>% #Should be 7, based on other B. nevadensis
  #I need more info about these spiders before assigning them to a discrete pass. For now, all missing passes marked as NA
  mutate(BTID=gsub('12490-DPF-2017','12490-NA-DPF-2017',BTID)) %>%
  mutate(BTID=gsub('12490-X-DPF-2017','12490-NA-DPF-2017',BTID)) %>%
  mutate(BTID=gsub('17774-DPF-2017','17774-NA-DPF-2017',BTID)) %>% 
  mutate(BTID=gsub('25070-DPF-2017','25070-NA-DPF-2017',BTID)) %>%
  mutate(BTID=gsub('15208-DPF-2017','15208-NA-DPF-2017',BTID)) %>% 
  #Create separate columns from BTID
  separate(BTID,c('BLID','pass','rep','year'),remove=F,convert=T,sep='-') %>% 
  #Tests merging from 2 different dfs:
  #Join in site type from site df
  left_join(select(site,BLID,siteType),by='BLID') %>%
  mutate(mismatchBLID=is.na(siteType)) %>% #Records mismatch in BLID (site label)
  #Join in location of specific trap from trap df
  left_join(select(trap,BTID,trapLoc),by='BTID') %>% 
  mutate(mismatchBTID=is.na(trapLoc)) #Records mismatch in BLID (trapping label)
  
#Problem: large number of NAs in traploc/siteType. First check trap df to makes sure that trapLoc is correctly specified, then decide what to do with remaining mismatches.
#This would take a large amount of time (I think) to go through.

#Number of mismatches per year - roughly 1-2% each year
temp %>% select(year,contains('mismatch')) %>% 
  summarize(total=n(),mismatchBLID=sum(mismatchBLID),mismatchBTID=sum(mismatchBTID)) %>% 
  mutate(propMismatchBTID=mismatchBTID/total)

#Function for comparing character-by-character similarity of a string to a vector of other strings
#Must be same length of strings
closestMatch <- function(target,possible,returnType=1){
  #Check inputs
  lt <- nchar(target)
  lp <- nchar(possible)
  np <- length(possible)
  if(length(lt)>1) stop('Multiple targets. Use sapply.')
  if(any(lp!=max(lp))) stop('Possible strings have different lengths.')
  
  out <- rep(NA,np) #Output vector
  splPoss <- strsplit(possible,'') #Splits into characters
  splTarg <- strsplit(target,'')[[1]]
  
  for(i in 1:length(possible)){
    out[i] <- sum(splTarg == splPoss[[i]])
  }
  if(returnType==1) return(out) #Return vector of character similarity
  if(returnType==2) {
    ind <- which(out==max(out)) #index of values to return
    if(length(ind)>1) warning('Multiple matches')
    return(possible[ind])
  }
}

testTarg <- 'abcde'
testPoss <- c('abcee','abcee','aeeee','aaade')
closestMatch(testTarg,testPoss,2)


# Save to file ------------------------------------------------------------

save(site,trap,arth,file='cleanData.Rdata') #Saves to cleaned data file


