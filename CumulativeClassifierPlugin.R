dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

#Neural Nets
library(randomForest)
library(rpart)
library(rattle)
library(caret)
input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
      pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

  print("READING INPUT FILES...");
  t1 <<- read.table(paste(pfix, toString(parameters["training",2]), sep=""), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  t2 <<- read.table(paste(pfix, toString(parameters["clinical",2]), sep=""), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  print("DONE");
  prefix <<- paste(pfix, toString(parameters["prefix", 2]), sep="");
  t3 <<- read.table(paste(pfix, toString(parameters["times", 2]), sep=""));
  joinby <<- toString(parameters["joinby", 2])
  orderby <<- toString(parameters["orderby", 2])
  target <<- toString(parameters["target", 2])
  id <<- toString(parameters["id", 2])
  myX <<- toString(parameters["x", 2])
  traincontrolmethod <<- toString(parameters["trainControl", 2])
  trainmethod <<- toString(parameters["train", 2])
  thresh <<- as.numeric(toString(parameters["threshold", 2]))
  classcol <<- toString(parameters["classcol", 2])
}

run <- function() {
   t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
   colnames(t1)[1] <<- joinby
   x <<- as.data.frame(merge(t1, t2, by =joinby, stringsAsFactors=FALSE))
   studyID <<- unique(as.character(unlist(x[id])))
   train_set_size <<- ncol(t1)
   class_index <<- grep(classcol, colnames(t2))

   #studyID <<- unique(x$STUDYID)
}

output <- function(outputfile) {
 for(virus in target){
  
  v1 <- x[x[,id]==virus,]
v1 <<- v1[order(as.character(unlist(v1[orderby]))),]
  #v1 <- v1[order(v1$SUBJECTID),]
  
   times <<- as.numeric(unlist(t3))
  rf.label = as.factor(v1[v1[,myX]==times[1],(train_set_size+class_index)])
  
  #times = c(-24, 0, 5, 12, 21.5, 36, 45.5)
  maxAc = 0
  for(t in times){
    print("TIME")
    print(t)
    for(threshold in thresh)
    {
      t.x = v1[v1[,myX]==t,]
      #t3 = read.csv(paste("ChosenGenes_Cumulative_",virus,"_",t,".csv",sep=""),header = TRUE)
      t3 = read.csv(paste(prefix,virus,"_",t,".csv",sep=""),header = TRUE)
      # newV = t.x[,c(1,as.numeric(unlist(t3[1]))+1,22279:dim(t.x)[2])]
      
      resCon = sapply(t3[2], function(x) x > threshold)
      fin = t3[2][resCon]
      #if(t == -24)
      if(t == times[1])
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
        fit  = trainControl(method = trainControl, number = 10, repeats = 10, index = cv.folds)
        res = train(x = datX, y = rf.label, method = train, tuneLength = 3, ntree = 1000, trControl = fit)
        
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

