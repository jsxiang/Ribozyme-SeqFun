setwd("~/Documents/MFS/RS_MFS_all/")

install.packages("gbm")
require(gbm)
require(MASS)

#separating training and test data
train=sample(1:506,size=374)


Boston.boost=gbm(medv ~ . ,data = Boston[train,],distribution = "gaussian",n.trees = 10000,
                 shrinkage = 0.01, interaction.depth = 4)
Boston.boost
summary(Boston.boost)


xandata<-read.csv(file="xan_onehotwlabel.csv")
xtrain=sample(1:5066,size=round(0.8*5066))

xandata.boost=gbm(RNA.DNA ~ . ,data = xandata[xtrain,],distribution = "gaussian",n.trees = 10000,
                  shrinkage = 0.01, interaction.depth = 10)
summary(xandata.boost)

theodata<-read.csv(file="theo_onehotwlabel.csv")
ttrain=sample(1:658,size=round(0.8*658))

theodata.boost=gbm(RNA.DNA ~ . ,data = theodata[ttrain,],distribution = "gaussian",n.trees = 10000,
                  shrinkage = 0.01, interaction.depth = 10)
summary(theodata.boost)


cdGdata<-read.csv(file="cdG_onehotwlabel.csv")
ctrain=sample(1:4966,size=round(0.8*4966))

cdGdata.boost=gbm(RNA.DNA ~ . ,data = cdGdata[ctrain,],distribution = "gaussian",n.trees = 10000,
                  shrinkage = 0.01, interaction.depth = 10)
summary(cdGdata.boost)

FAdata<-read.csv(file="FA_onehotwlabel.csv")
ftrain=sample(1:5117,size=round(0.8*5117))

FAdata.boost=gbm(RNA.DNA ~ . ,data = FAdata[ftrain,],distribution = "gaussian",n.trees = 10000,
                  shrinkage = 0.01, interaction.depth = 10)
summary(FAdata.boost)


#Plot of Response variable with lstat variable
plot(FAdata.boost,i="U6") 
#Inverse relation with lstat variable

plot(Boston.boost,i="rm") 
#as the average number of rooms increases the the price increases

n.trees = seq(from=100 ,to=10000, by=100) #no of trees-a vector of 100 values 

#Generating a Prediction matrix for each Tree
predmatrix<-predict(xandata.boost,xandata[-xtrain,],n.trees = n.trees)
dim(predmatrix) #dimentions of the Prediction Matrix

#Calculating The Mean squared Test Error
test.error<-with(xandata[-xtrain,],apply( (predmatrix-RNA.DNA)^2,2,mean))


head(test.error) #contains the Mean squared test error for each of the 100 trees averaged

#Plotting the test error vs number of trees

plot(n.trees , test.error , pch=19,col="blue",xlab="Number of Trees",ylab="Test Error", main = "Perfomance of Boosting on Test Set")

#adding the RandomForests Minimum Error line trained on same data and similar parameters
abline(h = min(test.err),col="red") #test.err is the test error of a Random forest fitted on same data
legend("topright",c("Minimum Test error Line for Random Forests"),col="red",lty=1,lwd=1)




