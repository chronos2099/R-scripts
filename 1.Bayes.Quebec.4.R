##################################################
##################################################
################ BAYES R SCRIPT ##################
##################################################
##################################################
# Modified by Julie, july 2014
##################################################

#install.packages("e1071")
library(e1071)
#install.packages("gmodels")
library(gmodels)

#####################
##### FUNCTIONS #####
#####################
#####################




afn.checkfidelity = function(ourfinalresult, truthtocheck)
{
  
  confusiontable = afn.makesquareconfusiontable(truthtocheck, ourfinalresult) 
  squaretable = confusiontable[[1]]
  print(squaretable)

  kappastuff = kappa(squaretable)
  print(kappastuff)
  
  fidelityfile = paste(kappadir, "finalguess.vs.check.txt", sep="")
  
  for (ii in 1:(length(kappastuff)))
  {
    towrite = paste(names(kappastuff)[ii], kappastuff[ii], sep=": ")
    write(towrite, fidelityfile, append = TRUE)
  }
  
}


kappa = function(CM) 
{
  
  # FROM DGRossiter 2004: TECHNICAL NOTE: STATISTICAL METHODS FOR ACCURACY ASSESSMENT OF CLASSIFIED THEMATIC MAPS
  
  #convert both data frames and vectors to matrices
  cmx<-as.matrix(CM)
  #try to convert a vector to a square matrix
  if (ncol(cmx) == 1)
    cmx<-matrix(cmx, byrow=TRUE, nrow=sqrt(nrow(cmx)))
  nr<-nrow(cmx); nc<-ncol(cmx)
  if (nr != nc)
  { print("Error: matrix is not square"); break }
  n<-sum(cmx)
  d<-diag(cmx); dsum<-sum(d); th1<-dsum/n #diag (extract diagonal of matrix) sum of diagonal divided by sum of matrix
  th1v<-((th1*(1-th1))/n)
  csum<-apply(cmx,2,sum); rsum<-apply(cmx,1,sum)
  ua<-d/rsum; pa<-d/csum
  th2 <- sum(rsum*csum) / n^2; kh <- (th1-th2)/(1-th2)
  th3 <- sum( (csum + rsum) * d ) / n^2;
  th4 <- 0; for (i in 1:nr) for (j in 1:nc)
    th4 <- th4 + (cmx[i,j] * ((csum[i] + rsum[j])^2));
  th4 <- th4 / n^3;
  th1c <- 1 - th1; th2c <- 1 - th2;
  khv <- 1/n *
    ( ( ( th1 * th1c ) / th2c^2 )
      + ( ( 2 * th1c * ((2*th1*th2) - th3) ) / th2c^3 )
      + ( ( th1c^2 * ( th4 - (4 * th2^2 ) ) ) / th2c^4 )
    )
  #per-class kappa, user's accuracy...
  p <- cmx/n; uap <- apply(p,1,sum); pap <- apply(p,2,sum); dp<-diag(p);
  kpu <- (dp/uap - pap)/(1 - pap);
  #...and its variance
  t1 <- uap-dp; t2 <- (pap*uap)-dp; t3 <- dp*(1 - uap - pap + dp);
  kpuv <- ( (t1/(uap^3 * (1-pap)^3)) * ((t1*t2) + t3) )/n;
  #per-class kappa, producer's reliability...
  kpp <- (dp/pap - uap)/(1 - uap);
  #...and its variance
  t1 <- (pap-dp);
  kppv <- ( (t1/(pap^3 * (1-uap)^3)) * ((t1*t2) + t3) )/n;
  #return all statistics as a list
  return(list(sum.n=n, sum.naive=th1, sum.var=th1v, sum.kappa=kh, sum.kvar=khv,
         user.naive=ua, prod.naive=pa,
         user.kappa=kpu, user.kvar=kpuv, prod.kappa=kpp, prod.kvar=kppv))
}


afn.daybreak = function(dayvalue)
{
  
  # this function is for visual purposes when running R (easier to distinguish when a new day iteration starts)
  
  pounds = "##############################################################"
  for (ii in c(1:3)) print(pounds)
  print(paste("##########################    ", dayvalue,"   ##########################"))  
  for (ii in c(1:3)) print(pounds)
  for (ii in c(1:3)) print("...")
  
}

###START HERE####
afn.dotheday = function (onerootID, onedayname, oneeventdayname, priors, event.comparison)
{
  
  # this function writes out directories for output files and calls necessary functions
  
  rootID     = onerootID
  theday     = onedayname
  feventname = paste(eventdir, oneeventdayname, sep="")
    
  ftruthtable      = paste(truthdir,   "1.truthtable.",    rootID, ".txt",      sep="")
  ftruthtable.mini = paste(truthdir,   "1.truthtable.",    rootID, ".mini.txt", sep="")
  fnameposterior   = paste(statedir,   "themachine.post.", theday, ".txt",      sep="")
  fname            = paste(daycounter, ".event.trimmed",                        sep="")
  
  event.current = afn.readandtrimevent(feventname, trimthearea, thereplacementvalue)[[1]]
  afn.writeforarc(matrix(event.current, ncol = colsforanalysis), fname, eventdir, xllcorner, yllcorner, F)
  
  results.0vs1 = afn.processoneevent(theday, ftruthtable, event.current, event.comparison, priors, ftruthtable.mini, fnameposterior, thereplacementvalue)
  
  posteriors               = results.0vs1[[1]]
  therunningclassification = results.0vs1[[2]]
  therunningmaxprob        = results.0vs1[[3]]
  
  return(list(posteriors, therunningclassification, therunningmaxprob, event.current))
  
}


afn.writerunningagreements = function(thecols, therows, kappafile.running, descrip, kappaf2, numtoround = 3)
{

  # this function writes out the kappa agreements to a txt file

  print(kappafile.running)
  print(kappaf2)
  
  confusiontable = afn.makesquareconfusiontable(thecols, therows)
  squaretable = confusiontable[[1]]
  print(squaretable)
  
  write(paste("Accuracy assessment statistical analysis for ", descrip, sep=""), kappafile.running, append = TRUE)
  sink(kappafile.running, append = TRUE)
  kappa.0vs1 = kappa(squaretable)
  print(kappa.0vs1)
  sink()

  write.table(squaretable, kappaf2, row.names = T, col.names = T, sep=" ")
  
}


afn.makesquareconfusiontable = function(thefirst, thesecond)
{
  
  # this function makes sure that the confusion table is square (adds rows/cols as necessary)
  # using http://www.statmethods.net/stats/frequencies.html
  
  print("making square truth table...")
  mytable = table(thesecond, thefirst)
  print(mytable)
  print(dim(mytable))
  
  print("making cross-tabulation truth table...")
  allanames   = rownames(mytable)
  allbnames   = colnames(mytable)
  lena        = length(allanames)
  lenb        = length(allbnames)
  allnames    = sort(unique(c(allanames, allbnames)))
  lenallnames = length(allnames)
  
  #numrowsintruthtable = max(as.numeric(dimnames(thetable.partial)[[1]]))
  thetable.big = matrix(0, lenallnames, lenallnames)
  rownames(thetable.big) = allnames
  colnames(thetable.big) = allnames
  print(paste("the big table dimensions: ", dim(thetable.big)[[1]], dim(thetable.big)[[2]]))
  print(thetable.big)
  
  thetable.partial = mytable
  for(ii in c(1:dim(thetable.partial)[[1]] ))
  {
    for(jj in c(1:dim(thetable.partial)[[2]] )) 
    {
      bigrownames = rownames(thetable.big) 
      bigcolnames = colnames(thetable.big)
      i.in.big.table = which(bigrownames == rownames(thetable.partial)[ii])
      j.in.big.table = which(bigcolnames == colnames(thetable.partial)[jj])
      print(paste(ii, jj, rownames(thetable.partial)[ii], colnames(thetable.partial)[jj], thetable.partial[ii,jj], i.in.big.table, j.in.big.table))
      valtowrite = thetable.partial[ii,jj]
      whichrowinbigtable = as.numeric(rownames(thetable.partial)[ii])
      thetable.big[i.in.big.table,j.in.big.table] = valtowrite
    }
  }
  print(paste("partial table dimensions: ", dim(thetable.partial)[[1]], dim(thetable.partial)[[2]]))
  print(thetable.partial)

  col.sums <- apply(thetable.big, 2, sum)
  thetable.proportions =  sweep(thetable.big, 2, (0.01+col.sums), `/`)
  round(thetable.proportions, 3)
#   thetable.proportions = thetable.big / col.sums
#   thetable.proportions[dim(thetable.proportions)[[1]],]= 1  # any constant will do
  
  print("done afn.makesquareconfusiontable")

  return(list(thetable.big, thetable.proportions, thetable.partial))# why 3 tables?
  
}


afn.load.tile = function(filepath, numrowstoread = -1, na.strings = "-9999")
{
  
  # this function loads in an ascii tile
  
  # note that na -9999 of arcgis will be read in as na
  # the is.na test later is faster than the x<1 test, so this is good and can maybe point the way to removing the x<1 test.
  
  tile = as.matrix(read.table(filepath, header=F, skip=6, sep=" ", nrows=numrowstoread))
#   tile = tile[,1:4918]  # this gets rid of an extra col if there is one
#   tile = tile[1:2049,]  # this gets rid of an extra row if there is one

  numcolsonreadin = dim(tile)[[2]]
  numrowsonreadin = dim(tile)[[1]]
  print(paste("Read in", numcolsonreadin, "columns and", numrowsonreadin, "rows from", filepath, sep=" "))
  
  return(tile)
  
}


afn.getheader.tile = function(filepath)
{
  
  # this function retrieves the header of the tile
  
  tile = read.table(filepath, header=F, nrows=6, row.names=1)
  
  return(list(tile[1,1], tile[2,1], tile[3,1], tile[4,1]))  
  
}


afn.writeposteriorclassandmaxprobabilityforarc = function(listofclassificationandmaxprob, thepath, theroot = "day", whichday, thecoltitles, colsforanalysis)
{

  # this function writes out the posterior information into an esrigrid
  
  # first, write the classification
  fname = paste(theroot, whichday, ".class",   sep="")
  afn.writeforarc(matrix(listofclassificationandmaxprob[[1]], ncol = colsforanalysis), fname, thepath, xllcorner, yllcorner, F)
  
  # next, write the max prob
  fname = paste(theroot, whichday, ".maxprob", sep="")
  afn.writeforarc(matrix(listofclassificationandmaxprob[[2]], ncol = colsforanalysis), fname, thepath, xllcorner, yllcorner, T, 3)
  
}


afn.writeforarc = function(themat, theband, arcpath, xllcorner, yllcorner, roundit=F, theround=2)
{
  
  # this function writes out the esrigrid for arc
  
  fname = paste(arcpath, theband, ".esrigrid.asc", sep="")
  print(fname)
  t1  = paste("ncols    ",     dim(themat)[2])
  t2  = paste("nrows    ",     dim(themat)[1])
  t3  = paste("xllcorner  ",   xllcorner)
  t4  = paste("yllcorner   ",  yllcorner)
  t5  = paste("cellsize      30.0")
  t6  = paste("NODATA_value  -9999")
  ttt = c(t1, t2, t3, t4, t5, t6)
  write.table(ttt, fname, col.names=F, row.names=F, quote=F, eol="\r\n")
  if(!roundit)
  {
    write.table(themat, fname, append=T, col.names=F, row.names=F, quote=F, eol="\r\n")  
  }
  else
  {
    write.table(round(themat, theround), fname, append=T, col.names=F, row.names=F, quote=F, eol="\r\n")  
  }
  
  return(fname)
  
}


afn.rebalancebyrow = function(themat, ndigits = 8)
{
  
  # this function takes away from the max and spreads it out over the rest of the rows
  
  thetotal = apply(themat, 1, sum)
  
  return(( round( apply(themat, 2, function(x)  x/thetotal), ndigits)))
  
}


afn.eachposteriorgetsonepercent = function(themat, multfactor, addfactor)
{
  
  # this function adds a small number to each posterior value
  
  return((themat * multfactor) + addfactor) 
  
}


afn.cleanarcvector = function(arcvector, replacevalue)
{
  
  # this function takes care of NA pixels
  
  arcvector[is.na(arcvector)]          = replacevalue
  arcvector[which(arcvector == -9999)] = replacevalue
  if(debugging)
    summary(arcvector)
  if(debugging)
    print(length(arcvector))
  
#   arcvector[arcvector<1] = replacevalue
#   if(debugging)
#     summary(arcvector)

#   numrows = dim(event.matrix)[[1]] 
#   numcols = dim(event.matrix)[[2]]
  
  return(arcvector)
  
}


afn.calculateposteriorsafterevent = function(priors, theevent, truthtable, theday, fnameposterior)
{
  
  # this function applies BULC to get posteriors
  
  print("entering afn.calculateposteriorsafterevent")
  write.table(date(), "time.1.start.txt")
  
  #numberofpixels= 3e6  #numberofpixels= 15  #numberofpixels= numrows * numcols  #numberofpixels = length(event.current)/100
  print(theevent)
  print(truthtable)
  
  print("applying the BULC algorithm")
  print(date())
  probmat = truthtable[theevent,]
  if(debugging)
  {
    print(probmat[1:numtoprint,])
    print(summary(probmat))
  }
  print("did apply")
  print(date())
  
  print("multiplying to get intermediate posteriors")
  posteriors.intermediate = priors * probmat
  if(debugging)
  {
    print(posteriors.intermediate[1:numtoprint,])
  }

  print("rebalancing intermediate posteriors")
  posteriors.rebalanced = afn.rebalancebyrow(posteriors.intermediate, numberofclasses)
  if(debugging)
  {  
    print(cbind((theevent[1:numtoprint]), posteriors.rebalanced [1:numtoprint,]))
  }
  
  print("finalizing posteriors with low values")
  posteriors = afn.eachposteriorgetsonepercent(posteriors.rebalanced, (1-(numberofclasses/200)), (1/200))
  if(debugging)
  {  
    print(cbind((theevent[1:numtoprint]), round(posteriors [1:numtoprint,],3)))
  }
  
  print("posteriors now finalized for this day.")
  
  dimnames(posteriors)[[2]] = dimnames(truthtable)[[2]]
  
  fnameposterior

  write.table(date(), "time.2.end.txt")

  if(debugging)
  {
    print(posteriors[1:numtoprint,])
  }
  
  print("leaving afn.calculateposteriorsafterevent")
  
  return(posteriors)
 
}


afn.readandtrimevent = function(feventname, trimthearea=F,  replacementvalue = numberofclasses)
{
  
  # this function reads in a matrix from a file and trims it if required
  
  print(paste("reading arc ascii file:", feventname))
  event.matrix = afn.load.tile(feventname, numrowstoread = rowlimited)
  print("finished reading.")

  if(trimthearea) #the indexes are set outside the function global to it
  {
    event.matrix.sub = event.matrix[thetop:thebottom, theleft:theright]
    rm(event.matrix)
  } 
  else
  {
    event.matrix.sub = event.matrix
  }
  
  rowsforanalysis = dim(event.matrix.sub)[[1]]
  colsforanalysis = dim(event.matrix.sub)[[2]]
  event.current = as.vector(event.matrix.sub)
    
  if(debugging)
  {
    print(dim(event.current))
    print(summary(as.factor(event.current)))
  }
  
  print("cleaning...")
  event.current.cleaned = afn.cleanarcvector(event.current, replacementvalue)
  
  if(debugging)
  {
    print(event.current.cleaned[1:numtoprint])
    print(dim(event.current.cleaned))
    print(summary(as.factor(event.current.cleaned)))
  }
  
  return(list(event.current.cleaned, rowsforanalysis, colsforanalysis))
  
}


afn.recordclassification = function(pixelsbyclasses)
{
  
  # this function records the class and prob from BULC
  
  print("recording the classification...")
  
  maxclass = max.col(pixelsbyclasses)
  for (i in c(length(thecolindices):1))
  {
    maxclass[maxclass == i] = thecolindices[i]
  }
  
  maxprob = apply(pixelsbyclasses, 1, 'max')
  
  if(debugging)
  {
    print(maxclass[1:numtoprint])
    print(maxprob[1:numtoprint])
  }
  
  print("finished recording the classification")
  theclassification = list(maxclass, maxprob)
  
  return(theclassification)
  
}


afn.makeconfusiontablepriorvscurrent = function(thetruth, thecurrentestimate, replacementvalue = 100)
{
  
  # this function makes a confusion table between either two consecutive days or one day and an NLCD
  
  print("making cross-tabulation truth table...")
  
  thetable.partial = table(thecurrentestimate, thetruth)#contigency table
  print("the small truth table")
  print(thetable.partial)
  numrowsintruthtable = max(as.numeric(dimnames(thetable.partial)[[1]]))
  
  thetable.big = matrix(0, replacementvalue, numberofclasses)#probability table of 4 classes 100 rows
  dimnames(thetable.big)[[2]] = as.character(thecolindices) # but we also need the last number 100
  
  for(ii in c(1:dim(thetable.partial)[[1]]))
  {
    for(jj in c(1:dim(thetable.partial)[[2]])) 
    {
      bigcolnames        = colnames(thetable.big)
      whichcolinbigtable = which(bigcolnames == colnames(thetable.partial)[jj])

      partialrownames    = rownames(thetable.partial)
      whichrowinbigtable = as.numeric(rownames(thetable.partial)[ii])
      
      valtowrite = thetable.partial[ii,jj]
      thetable.big[whichrowinbigtable,whichcolinbigtable] = valtowrite
    }
  }
  
  print("the truth table in terms of pixel counts, not proportions:")
  print(thetable.big)
  
  col.sums <- apply(thetable.big, 2, sum)
  thetable.proportions = thetable.big / col.sums
  thetable.proportions = sweep(thetable.big, 2, (0.01+col.sums), `/`)#summary stats 

  # add a row for processing nodata, which will be set to be the replacement value elsewhere
  lastrow = dim(thetable.proportions)[[1]]
  standardvalue  = 1
  standardvector = rep(standardvalue, numberofclasses)
  thetable.proportions[lastrow,] = standardvector
  print("the truth table in terms of proportions:")
  print(thetable.proportions)
  
  return(list(thetable.big, thetable.proportions, thetable.partial))
  
}


afn.processoneevent = function(name.today, filename.truthtable, rowsfortruthtable, colsfortruthtable, thepriorprobabilities, ftruthtablesmall, filename.posterior, replacementvalue)
{
  
  # this function does truth table stuff and calls the function that runs BULC
  print("summary of rows for truth table (ie the event.today) is below:")
  print(summary(as.factor(rowsfortruthtable)))
  
  truthtable.prior.current = afn.makeconfusiontablepriorvscurrent(colsfortruthtable, rowsfortruthtable)

  print("Truth table of the running estimate vs the day")
  bigtruthtable   = truthtable.prior.current[[2]]
  print("big truth table is below")
  print(bigtruthtable)
  smalltruthtable = truthtable.prior.current[[3]]
  write.table(round(bigtruthtable,standardroundingdigits), paste(filename.truthtable,".initial.txt",sep=""))
  write.table(smalltruthtable, ftruthtablesmall)
  
  print("Keeping pixels that were given the same class two consecutive days; clearing the rest")
  rowsfortruthtable[which(rowsfortruthtable != colsfortruthtable)] = replacementvalue
  print(paste("truth table dimensions are: " , dim(bigtruthtable)))
  write.table(round(bigtruthtable,standardroundingdigits), filename.truthtable)
  
  print("Making posteriors now...")
  posteriors = afn.calculateposteriorsafterevent(thepriorprobabilities, rowsfortruthtable, bigtruthtable, name.today, filename.posterior)
  theclassification = afn.recordclassification(posteriors)
  runningclass      = as.vector(theclassification[[1]])
  runningmaxprob    = as.vector(theclassification[[2]])
  cbind(runningclass[1:numtoprint], runningmaxprob[1:numtoprint])
  afn.writeposteriorclassandmaxprobabilityforarc(theclassification, statedir, rootarcname, name.today, thecoltitles, colsforanalysis)
  
  results = list(posteriors, runningclass, runningmaxprob)
  return(results)
  
}


afn.makeinitialpriors.equalpriors = function(numberofpixels, columntitles, outputfile)
{
  
  onepixelprob = rep(1/numberofclasses, numberofclasses, byrow = T)
  print(sum(onepixelprob))
  
  thepriorsmatrix = matrix(rep(onepixelprob, numberofpixels), numberofpixels, byrow = T)
  colnames(thepriorsmatrix) = columntitles
  print(thepriorsmatrix[1:numtoprint,])
  write.table(thepriorsmatrix, outputfile)
  
  return(thepriorsmatrix)
  
}


############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ ############ ############ ############


####################
# Define variables #
####################
####################

rowlimited             = -1
numberofclasses        = 3
thereplacementvalue    = 100
thecellsize            = 30
numtoprint             = 15
standardroundingdigits = 3
debugging              = T

theleft     = 3000
theright    = 3299
thetop      = 4000
thebottom   = 4299
trimthearea = T

thecolindices = c(1, 2, 3)
outputcolnames = c(     
                    "1.water"
                  , "2.forest"
                  , "3.human"
                   )
thecoltitles = c(outputcolnames, "Class", "MaxProb")


######################
# Define directories #
######################
######################

truthdir = "../1.TruthTables/"
eventdir = "../2.Events/"
statedir = "../3.SystemStates/"
picdir   = "../4.Pictures/"
kappadir = "../6.Kappas/"

kappafile.eventbyevent    = paste(kappadir, "0.kappas.eventbyevent.txt",    sep="")

rootpicname = "2.post."
rootarcname = "2.post"
fstatesname = paste(statedir, "*.txt", sep="")


################################################
# Read an event to get preliminary information #
################################################
################################################

feventname = paste(eventdir,"1.CART/130602cart.asc", sep="")
theevent   = afn.readandtrimevent(feventname, trimthearea, thereplacementvalue) 

rowsforanalysis = theevent[[2]]
colsforanalysis = theevent[[3]]

numberofpixels = rowsforanalysis*colsforanalysis
print(paste(numberofpixels, "pixels"))

if(trimthearea == F)
{
  tileheader = afn.getheader.tile(feventname)  
  xllcorner = tileheader[[3]]
  yllcorner = tileheader[[4]]
} else 
{
  xllcorner = 0
  yllcorner = 0
}


#################################
# Make equal probability priors #
#################################
#################################

fstate0 = "../3.SystemStates/0.equalpriors.txt"
priors = afn.makeinitialpriors.equalpriors(numberofpixels, outputcolnames, fstate0)  # first time only


##################
# Read in events #
##################
##################

event.comparison = theevent[[1]]

rootIDnames    = c(
                      "0vs1" 
                    , "1vs2"
                    , "2vs3"
                    , "3vs4"
                    , "4vs5"
                    , "5vs6"
                    , "6vs7"
                    , "7vs8"
                    , "8vs9"
                    , "9vs10"
                  )

eventnames     = c(
                      "1.18jun2013"
                    , "2.26jun2013"
                    , "3.04jul2013"
                    , "4.12jul2013"
                    , "5.29aug2013"
                    , "6.06sep2013"
                    , "7.14sep2013"
                    , "8.22sep2013"
                    , "9.08oct2013"
                    , "10.18jun2014"
                  )
                 
eventlocations = c(
                      "1.CART/130618cart.asc"
                    , "1.CART/130626cart.asc"
                    , "1.CART/130704cart.asc"
                    , "1.CART/130712cart.asc"
                    , "1.CART/130829cart.asc"
                    , "1.CART/130906cart.asc"
                    , "1.CART/130914cart.asc"
                    , "1.CART/130922cart.asc"
                    , "1.CART/131008cart.asc"
                    , "1.CART/140618cart.asc"
                  )  
  


####################
# Execute each day #
####################
####################

posteriors = priors # first time only

numberofdaystoexecute = 1

for (daycounter in c(1:numberofdaystoexecute)) 
  {
    
    afn.daybreak(daycounter)
    onedayresult = afn.dotheday(rootIDnames[daycounter], eventnames[daycounter], eventlocations[daycounter], posteriors, event.comparison)
    # param 1: used for naming files
    # param 2: used for names
    # param 3: used to open an event file and use it.
    # param 4: used inside the function as the probabilities *prior* to the event taking place.
    
    posteriors = onedayresult[[1]]
    running    = onedayresult[[2]]
    yesterday  = onedayresult[[4]]
    
    cellswithdata = which(yesterday != thereplacementvalue)
    for (ii in 1:length(cellswithdata))
      event.comparison[cellswithdata[ii]] = yesterday[cellswithdata[ii]]
    
    print(paste("done day", daycounter))
    
  }


##################################################


print("done.")
