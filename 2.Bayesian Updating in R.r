library(e1071)
library(data.table)
library(dplyr)

wd <- 'D:/Dropbox/Cardille lab/Work/ascii maps/Bayesian maps'
setwd(wd)



rowsstart <- 3000
rowsend <- 4000
rows <- rowsend - rowsstart

colstart <- 3000
colend <- 4000
cols <- colend - colstart

initial<- tbl_df(data.table( Class1 = rep(1/3, rows*cols),Class2 = rep(1/3, rows*cols),Class3 = rep(1/3, rows*cols)))

map1<-na.omit(tbl_df(fread(dir()[1]))[rowsstart:rowsend,colstart:colend])

map2<-na.omit(tbl_df(fread((dir()[2])))[rowsstart:rowsend,colstart:colend])

system.time(list1<- tbl_df(as.data.frame(unlist(map1))))#2 minutes


list2<- tbl_df(as.data.frame(unlist(map2)))

#list2<-list2)

names(list1)[1]<- "Classes"
names(list2)[1]<- "Classes"
landscape.class<- na.omit(unique(list1))


		
truthtable1<- list1 %.%
	mutate (Class1 = ((Classes%/%-9999)==1)*1, Class2 = ((Classes%/%2)==1)*1, Class3 = ((Classes%/%3)==1)*1) %.%
	select (Class1,Class2,Class3)
	
truthtable2<- list2 %.%
	mutate (Class1 = ((Classes%/%-9999)==1)*1, Class2 = ((Classes%/%2)==1)*1, Class3 = ((Classes%/%3)==1)*1) %.%
	select (Class1,Class2,Class3)

summarytable<- summarise (truthtable1,Class1sum = sum(Class1), Class2sum= sum(Class2), Class3sum = sum(Class3))
	
	
	
#probability<- mutate((truthtable*0.925)+0.025)


#test <- truthtable[1:100,]


	


#model <- naiveBayes(zz ~., data = 'aa')


