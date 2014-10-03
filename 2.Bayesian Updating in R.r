library(data.table)
library(dplyr)

wd <- 'D:/Dropbox/Cardille lab/Work/ascii maps/Bayesian maps'
setwd(wd)

#bayes<- function(x) {
#	x*


rowsstart <- 3000
rowsend <- 4000
rows <- rowsend - rowsstart

colstart <- 3000
colend <- 4000
cols <- colend - colstart



map1<-na.omit(tbl_df(fread(dir()[1]))[rowsstart:rowsend,colstart:colend])

map2<-na.omit(tbl_df(fread((dir()[2])))[rowsstart:rowsend,colstart:colend])

list1<- tbl_df(as.data.frame(unlist(map1)))#2 minutes

list2<- tbl_df(as.data.frame(unlist(map2)))


initial<- tbl_df(data.table( Class1 = rep(1/3, nrow(list1)),Class2 = rep(1/3, nrow(list1)),Class3 = rep(1/3, nrow(list1))))


#list2<-list2)

#names(list1)[1]<- "Classes"
#names(list2)[1]<- "Classes"
#landscape.class<- na.omit(unique(list1))

confusiontable<- table(list2[[1]],list1[[1]])


probabilitytable<- apply(confusiontable,2,function(x) x/sum(x))


prob.mat<- matrix(0,3,3)	
colnames(prob.mat)<- c(1,2,3)
rownames(prob.mat)<- c(1,2,3)

prob.mat[,match(colnames(probabilitytable),colnames(prob.mat))]<- probabilitytable[,na.omit(match(colnames(prob.mat),colnames(probabilitytable)))]


prob.mat<- prob.mat[match(list1[[1]],rownames(prob.mat)),]

rownames(prob.mat)<-rownames(list1)

prob.mat<-tbl_df(as.data.frame(prob.mat))

new.mat<- prob.mat*initial

posterior<-t(apply(new.mat,1,function(x) x/sum(x)))

initial[list1[[1]] == list2[[1]],] <- posterior[list1[[1]] == list2[[1]],]
 
new.list<- max.col(initial)

new.map<- matrix(new.list,dim(map1)[1],dim(map1[2]))


