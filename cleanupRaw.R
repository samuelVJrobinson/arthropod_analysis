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

# NOTE: aliases are not present in either trap or arthropod database (see commented code below to verify)
site <- site %>% select(-aliasBLID1:-locality,-country,-created_at:-expt_yield_Lakeland) %>% 
  mutate_if(is.factor,as.character)

# #Checking to see if any of the aliases are represented in the rest of the database
# aliases <- site %>% select(aliasBLID1,aliasBLID2) %>% 
#   pivot_longer(cols=c(aliasBLID1,aliasBLID2),names_to='type',values_to='alias') %>% 
#   filter(alias!=0,!is.na(alias)) %>% select(alias) %>% 
#   mutate(inTrapBLID=alias %in% trap$BLID,inArthBLID=alias %in% arth$BLID) %>% 
#   rowwise() %>% 
#   mutate(inTrapBTID=any(grepl(alias,trap$BTID))) %>%
#   mutate(inArthBTID=any(grepl(alias,arth$BTID)))


# Clean up trap data ------------------------------------------------------

#PROBLEM: SOME BLIDS FROM TRAP DF ARE NOT PRESENT IN SITE DF (either in BLID or aliases)
# badBLID <- trap %>% select(BTID:trapType) %>% filter(!(BLID %in% site$BLID),!is.na(BLID)) %>% 
#   pull(BLID) %>% unique()
#Bad BLID values:
# 23460 - should be 13460: CHANGE IN TRAP AND ARTH DF
# 13395 - should be 11395: CHANGE IN TRAP AND ARTH DF
# 13048 - should be 10348: CHANGE IN TRAP AND ARTH DF
# 13752 - should be 13751: CHANGE IN TRAP DF
# 20427 - duplicate of 24027, nothing in site or arth df: REMOVE FROM TRAP DF
# 15280 - duplicate of 15208, nothing in arth df: REMOVE FROM TRAP DF
# 25095 - 25091 or 25035, but were recorded on different days (June 22, 23), nothing in site or arth df: REMOVE FROM TRAP DF
# 13903 - duplicate of 13902: REMOVE FROM TRAP DF
# 27029 - present in site table, but not arth table: REMOVE FROM TRAP DF

removeBLID <- c(20427,15280,25095,13903,27029) #BLID values to remove from trap df
#BLID values to change
changeBLID <- data.frame(BLID=c(23460,13395,13048,13752),to=c(13460,11395,10348,13751)) 

# Notes: 
# "Coloured cups" consist of blue, white & yellow cup, but were sometimes all amalgamated into single category, so 3x effort makes sense.

#trap2 <- 
temp <- trap %>% mutate_if(is.factor,as.character) %>% #Convert all factors to characters
  filter(!(BLID %in% removeBLID)) %>% #Removes bad BLIDs
  mutate(change=(BLID %in% changeBLID$BLID)) %>% #Marks BLIDs for changes
  left_join(changeBLID,by='BLID') %>% rowwise() %>% 
  mutate(BTID=ifelse(change,gsub(as.character(BLID),to,BTID),BTID),BLID=ifelse(change,to,BLID)) %>% 
  ungroup() %>% select(-change,-to) %>% #Cleanup
  #Re-merge site data
  select(-lonTrap:-locationType,-created_at:-updated_at,-startNESWPhoto_DSC:-endNESWPhoto_DSC,-floralAdjacentNotes) %>% 
  left_join(select(site,BLID:lon,elevation,siteType),by='BLID') %>% #Join in lat,lon,elev, siteType
  filter(!is.na(trapType)) %>% #Get rid of NA row
  #Fixes typos
  mutate(siteType=gsub('Ditch','ditch',siteType)) %>% #Fix typo in locationtype
  mutate(BTID=gsub('2107','2017',BTID)) %>% #Fixes reversal in year label
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
  #Beginning to reorganize trap dataframe:
  rename('trapLoc'='siteType') %>% #trapLoc = location of trap (what feature is the trap sitting in?)
  #dist = distance from some "source" feature (usually wetland/SNL if in field)
  #distFrom = type of "source" feature
  rowwise() %>% mutate(dist=as.numeric(gsub('[A-Z]','0',replicate)),distFrom='NA') %>% ungroup() %>% 
  mutate(useYear=grepl('-201[6-8]',BTID)) %>%   #Fix 2016-2018 data
  mutate(distFrom=ifelse(useYear,gsub('inField_','',trapLoc),distFrom)) %>% 
  mutate(trapLoc=case_when(
    useYear & dist>0  ~ adjCrop,
    useYear & dist==0  ~ distFrom,
    TRUE ~ trapLoc)) %>% 
  select(-useYear) %>% 
  filter(!grepl('-2019',BTID)) #Material from 2019 not prepped yet, so filter out for now

#Other things to fix:
#Replace 'infield' with appropriate trap type in 2018 data
# Fix labels in arth database using appropriate BLIDs, re-merge data from trap df
# Create exclosure/non-exclosure column in 2019 data, remove "Exclosure Pit Fall" from trapType
    

  
  
  
  
  
  
  #in-field trapType recorded as "infield" for 2018, with type of trap recorded in replicate code. Confusing...
  mutate(isInfield=trapType=='infield',
         trapType=case_when(isInfield & replicate=='WBV'~'Blue Vane',
                            isInfield & replicate=='WPF'~'Pit Fall',
                            isInfield & replicate=='WCC'~'Coloured Cups',
                            TRUE ~ trapType))
  
# Clean up arthropod data -------------------------------------------------
# head(arth)

arth <- arth %>% mutate_if(is.factor,as.character) %>% #Converts all factors to character
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>% 
  #Spelling mistake
  mutate(BTID=gsub('WCCC','WCC',BTID))


# Save to file ------------------------------------------------------------

save(site,trap,arth,file='cleanData.Rdata') #Saves to cleaned data file


