#NEW METHOD
# Download data from ecologics using mySQL workbench, then import csv files
site <- read.table('./data/site.csv',header=T,sep=';',stringsAsFactors = F)
trap <- read.table('./data/trap.csv',header=T,sep=';',stringsAsFactors = F)
arth <- read.table('./data/arthropod.csv',header=T,sep=',',stringsAsFactors = F)

save(site,trap,arth,file='./data/rawData.Rdata') #Saves to raw data file

#OLD METHOD
# library(RODBC)
# 
# conn <- odbcConnect(DSN='ecologics') #ODBC connection must be set up beforehand
# 
# #Queries - takes about 1min to download
# site <- sqlQuery(conn,'SELECT * FROM site') #Site table
# #Key info: BLID
# trap <- sqlQuery(conn,'SELECT * FROM trap') #Trap table
# arth <- sqlQuery(conn,'SELECT * FROM arthropod') #Arthropod table - animal records
# odbcClose(conn) #Closes connection
# 
# save(site,trap,arth,file='rawData.Rdata') #Saves to raw data file
