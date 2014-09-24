library(e1071)
library(data.table)
library(dplyr)


wd <- 'D:/Dropbox/Cardille lab/Work/ascii maps/Bayesian maps'
setwd(wd)

initial<- tbl_df(data.table( Class1 = rep(1/3, 6668*6668),Class2 = rep(1/3, 6668*6668),Class3 = rep(1/3, 6668*6668)))

map1<-tbl_df(fread(dir()[1]))

map2<-tbl_df(fread((dir()[2])))

system.time(list1<- tbl_df(as.data.frame(unlist(map1))))#2 minutes


list2<- tbl_df(as.data.frame(unlist(map2)))

#list2<-list2)

names(list1)[1]<- "Classes"

landscape.class<- na.omit(unique(list1))


ptm<-proc.time()

for (i in 1:length(list1[[1]])){
	if (match(list2[i,],list1[i,],nomatch= 0) == TRUE) {
		change.class<- match(list2[i,],landscape.class[[1]])
		prob.A<-initial[i,change.class]*0.95
		probability <- prob.A / (prob.A + sum(initial[i,-change.class])*0.05)
		initial[i,change.class] <- probability
		}
	}
	
proc.time()-ptm		
		
		
		
		
		
		
#truthtable<- list1 %.%
#	mutate (Class1 = ((Classes%/%-9999)==1)*1, Class2 = ((Classes%/%2)==1)*1, Class3 = ((Classes%/%3)==1)*1) %.%
#	select (Class1,Class2,Class3)
	
#probability<- mutate((truthtable*0.925)+0.025)


#test <- truthtable[1:100,]


	


#model <- naiveBayes(zz ~., data = 'aa')


