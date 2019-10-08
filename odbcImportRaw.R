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
