# THIS FILE CLEANS AND SUMMARIZES WETLAND AND INFIELD DATA FROM ECOLOGICS

library(dplyr)
library(tidyr)

#Load raw data
load('rawData.Rdata')

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

# Notes: 
# "Coloured cups" consist of blue, white & yellow cup, but were sometimes all amalgamated into single category, so 3x effort makes sense.
trap <- trap %>% 
  mutate_if(is.factor,as.character) %>% #Converts all factors to character
  #Fixes typos
  mutate(replicate=gsub('WCCC','WCC',replicate),BTID=gsub('WCCC','WCC',BTID)) %>%
  mutate(replicate=gsub('WCC2R1','WCC25R1',replicate),BTID=gsub('WCC2R1','WCC25R1',BTID)) %>% 
  #3x deployedhours for "coloured cups" category (3 colour types)
  mutate(deployedhours=case_when(trapType=='Coloured Cups' ~ deployedhours*3,TRUE ~ deployedhours)) %>% 
  #Adds deployedhours from separate cups to each other
  group_by(BTID) %>% mutate(deployedhours=sum(deployedhours)) %>% ungroup() %>% 
  filter(trapType!='Yellow Cup'& trapType!='White Cup') %>%  #Removes yellow and white cup rows
  mutate(trapType=gsub('Blue Cup','Coloured Cups',trapType)) %>%  #Changes blue cup rows to "coloured cup"
  #Convert blue, white, yellow cups to "coloured cups"
  mutate(BTID=gsub('W[BWY]C','WCC',BTID)) %>%
  mutate(replicate=gsub('W[BWY]C','WCC',replicate)) %>% 
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


