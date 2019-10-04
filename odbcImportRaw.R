#THIS FILE DOWNLOADS DATA FROM ECOLOGICS

library(RODBC)
conn<-odbcConnect('ecologics',case='nochange') #ODBC connection must be set up beforehand

#Queries 
site <- sqlQuery(conn,'SELECT * FROM site') #Site table
#Key info: BLID
trap <- sqlQuery(conn,'SELECT * FROM trap') #Trap table
arth <- sqlQuery(conn,'SELECT * FROM arthropod') #Arthropod table - animal records
#NOTE: THIS TAKES ABOUT 1 MINS TO DOWNLOAD
save(site,trap,arth,file='rawData.Rdata') #Saves to raw data file

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

#Notes on siteType
#siteType = ditch, infield_wetland, infield_pivot, infield_control
#   within field, positioned at distances from features (control has no distance)

#Other notes:
#Check whether zeros are real! Later material may not have been entered. Dialictus missing from 2018.
#Some species haven't been IDed (e.g. Lasioglossum_spp. 6)
#Arthropod table needs cleaning (prior to joining) due to typing errors
#GDD is a more relevant measure of time than Jday (this is what the bees would "feel" rather than DOY)

