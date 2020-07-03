#Neural Nets
print("X");
library(randomForest)
library(rpart)
library(rattle)
library(caret)
print("Y");
input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  print("READING INPUT FILES...");
  t1 <<- read.table(toString(parameters["training",2]), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  t2 <<- read.table(toString(parameters["clinical",2]), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  print("DONE");
  prefix <<- toString(parameters["prefix", 2]);
}

run <- function() {
   t2[t2$STUDYID == "DEE4X H1N1",]$STUDYID <<- "H1N1"
   t2[t2$STUDYID == "DEE3 H1N1",]$STUDYID <<- "H1N1"
   t2[t2$STUDYID == "DEE2 H3N2",]$STUDYID <<- "H3N2"
   t2[t2$STUDYID == "DEE5 H3N2",]$STUDYID <<- "H3N2"
   t2[t2$STUDYID == "Rhinovirus Duke",]$STUDYID <<- "Rhinovirus"
   t2[t2$STUDYID == "Rhinovirus UVA",]$STUDYID <<- "Rhinovirus"

   t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
   colnames(t1)[1] <<- "CEL"
   x <<- as.data.frame(merge(t1, t2, by ="CEL", stringsAsFactors=FALSE))
   studyID <<- unique(x$STUDYID)
}

output <- function(outputfile) {
 for(virus in "DEE1 RSV"){
  
  v1 <- x[x[,"STUDYID"]==virus,]
  v1 <- v1[order(v1$SUBJECTID),]
  
  rf.label = as.factor(v1[v1[,"TIMEHOURS"]==-24,22286])
  
  times = c(-24, 0, 5, 12, 21.5, 36, 45.5)
  maxAc = 0
  for(t in times){
    print("TIME")
    print(t)
    for(threshold in 0.0005)
    {
      t.x = v1[v1[,"TIMEHOURS"]==t,]
      #t3 = read.csv(paste("ChosenGenes_Cumulative_",virus,"_",t,".csv",sep=""),header = TRUE)
      t3 = read.csv(paste(prefix,virus,"_",t,".csv",sep=""),header = TRUE)
      # newV = t.x[,c(1,as.numeric(unlist(t3[1]))+1,22279:dim(t.x)[2])]
      
      resCon = sapply(t3[2], function(x) x > threshold)
      fin = t3[2][resCon]
      if(t == -24)
      {
        newV = t.x[,c(as.numeric(unlist(t3[1][[1]][1:length(fin)]))+1)]
      }
      else
      {
      newV = cbind(newV, t.x[,c(as.numeric(unlist(t3[1][[1]][1:length(fin)]))+1)])
      }
      set.seed(1283)
      if(length(fin)!=0)
      {
        datX = data.matrix(newV)
        
        cv.folds <- createMultiFolds(rf.label, k=10, times = 10)
        fit  = trainControl(method = "repeatedcv", number = 10, repeats = 10, index = cv.folds)
        res = train(x = datX, y = rf.label, method = "rf", tuneLength = 3, ntree = 1000, trControl = fit)
        
        df = c(max(res$results$Accuracy), t, virus, length(fin), threshold)
        write(df, file = outputfile, sep = ",", append=TRUE)
        print(max(res$results$Accuracy))
        if(max(res$results$Accuracy)>maxAc)
        {
          maxAc = max(res$results$Accuracy)
          bestFit = fit
          bestRes = res
          bestTime = t
          bestThreshold = threshold
          bestNumberOfFeatures = length(fin)
        }
      }
    }
  }
 }
}

