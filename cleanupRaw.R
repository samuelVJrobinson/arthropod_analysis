# THIS FILE CLEANS AND SUMMARIZES WETLAND AND INFIELD DATA FROM ECOLOGICS

library(dplyr)
library(tidyr)

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
site <- site %>% select(-aliasBLID1:-locality,-country,-created_at:-expt_yield_Lakeland)

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
badBLID <- trap %>% select(BTID:trapType) %>% filter(!(BLID %in% site$BLID),!is.na(BLID)) %>% 
  pull(BLID) %>% unique()
#Bad BLID values:
# 20427 - duplicate of 24027, nothing in arth df: REMOVE
# 23460 - should be 13460, change in trap and arth
# 13395 - should be 11395, change in trap and arth
# 13048 - should be 10348, change in trap and arth
# 15280 - duplicate of 15208, nothing in arth df: REMOVE
# 13752 - should be 13751, change in trap df
# 25095 - 25091 or 25035, but were recorded on different days (June 22, 23), nothing in arth df
# 13903 - duplicate of 13902, but has correct deployedhours and startHour
# 27029 -


trap %>% filter(BLID==27029)
trap %>% filter(BLID==27029,endYear==2017,trapType=='Blue Vane') %>% arrange(BLID)
trap %>% filter(BLID==13903|BLID==13902,endYear==2017,trapType=='Blue Vane',pass==4) %>% arrange(BLID)

arth %>% filter(BLID==13902)

arth %>% filter(BTID=='13902-3-DBV-2017')




# Notes: 
# "Coloured cups" consist of blue, white & yellow cup, but were sometimes all amalgamated into single category, so 3x effort makes sense.
trap2 <- trap %>% 
  mutate_if(is.factor,as.character) %>% #Converts all factors to character
  select(-lonTrap:-locationType)
  
  #Fixes typos
  mutate(locationType=gsub('Ditch','ditch',locationType)) %>% #Fix typo in locationtype
  mutate(replicate=gsub('WCCC','WCC',replicate),BTID=gsub('WCCC','WCC',BTID)) %>%
  mutate(replicate=gsub('WCC2R1','WCC25R1',replicate),BTID=gsub('WCC2R1','WCC25R1',BTID)) %>% 
  mutate(collector=case_when(collector=='J. VIckruck' ~ 'J. Vickruck',
                             collector=='M.Gavin' ~ 'M. Gavin',
                             TRUE ~ collector)) %>% 
  #3x deployedhours for "coloured cups" category (3 colour types)
  mutate(deployedhours=case_when(trapType=='Coloured Cups' ~ deployedhours*3,TRUE ~ deployedhours)) %>% 
  #Adds deployedhours from separate cups to each other
  group_by(BTID) %>% mutate(deployedhours=sum(deployedhours)) %>% ungroup() %>% 
  filter(trapType!='Yellow Cup'& trapType!='White Cup') %>%  #Removes yellow and white cup rows
  mutate(trapType=gsub('Blue Cup','Coloured Cups',trapType)) %>%  #Changes blue cup rows to "coloured cup"
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>%
  mutate(replicate=gsub('W[BWY]C','WCC',replicate)) %>%
  data.frame()

#Other things to fix:
# Exclosure/non-exclosure column, removing "Exclosure Pit Fall" from trapType

#Beginning to reorganize trap dataframe. 
#trapLoc = location of trap (what feature is the trap sitting in?)
#dist = distance from some other feature (usually wetland/SNL if in field)
#distFrom = feature that 
  trap2 %>% select(trapType,collector,endYear,replicate,locationType,adjCrop) %>% distinct() %>% 
    #Creates dist column
    rowwise() %>% mutate(dist=as.numeric(gsub('[A-Z]','0',replicate))) %>% ungroup() %>% 
    mutate(locationType=case_when( #Fix NA locations
      is.na(locationType)&endYear==2016&grepl("W",replicate) ~ 'wetland', #2016 W traps were next to wetland
      is.na(locationType)&endYear==2016&replicate=='A' ~ 'ditch',
      is.na(locationType)&collector=='D. Clake' ~ 'alpine',
      is.na(locationType)&collector=='G. Sekulic' ~ 'ditch',
      TRUE ~ locationType)) %>%
    data.frame()
    
    
    mutate(trapLoc=case_when(dist==0 ~ locationType )) %>% data.frame()
  
  
  
  
  
  
  
  #in-field trapType recorded as "infield" for 2018, with type of trap recorded in replicate code. Confusing...
  mutate(isInfield=trapType=='infield',
         trapType=case_when(isInfield & replicate=='WBV'~'Blue Vane',
                            isInfield & replicate=='WPF'~'Pit Fall',
                            isInfield & replicate=='WCC'~'Coloured Cups',
                            TRUE ~ trapType))

  
  


  
  
  
  
  as.data.frame()


# Clean up arthropod data -------------------------------------------------
# head(arth)

arth <- arth %>% mutate_if(is.factor,as.character) %>% #Converts all factors to character
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>% 
  #Spelling mistake
  mutate(BTID=gsub('WCCC','WCC',BTID))


# Save to file ------------------------------------------------------------

save(site,trap,arth,file='cleanData.Rdata') #Saves to cleaned data file


