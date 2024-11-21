#For BayesLiNGAM
dir.create("./Bayes_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

df <- read.csv("FCPC_Bayes.csv")

# Please install the commands for BayesLiNGAM based on the URL (https://www.cs.helsinki.fi/group/neuroinf/lingam/bayeslingam/)
setwd("/Users/_____/Bayeslingamtest")

library('fastICA')
source('bayeslingam/main/loud.R')

loud()

#Calculation
x <- df$Test
y <- df$AIB
z <- df$Butyrate
w <- df$Blautia
w1 <- df$Lactobacillus
#w3 <- df$Tyramine
d <- data.frame(x1=x, x2=y,x3=z,x4=w,x5=w1)
result <- greedybayeslingam(d,model='GL')

#Data Export
sink('./Bayes_result/Bayes_Raw_data.txt', append = TRUE)
print (result)
sink()


#Visualize the plot
par("mar"=c(1,1,1,1))

#reference https://qiita.com/tomiyou/items/ca7032b1e0f1bf2b437b
library(igraph)
par(mfrow=c(2,3))
prob <- round(result$prob,digits=4) * 100
n.node <- max(result$components$node)
case <- nrow(result$DAGs)
node <- result$DAGs[,1:n.node]
res <- as.matrix(result$DAGs[,-c(1:n.node)],nrow=case)
name <-paste("X",1:n.node,sep="")

for(i in order(prob,decreasing=T)[1:6]){
  amat <- matrix(0,n.node,n.node)
  index <- node[i,]
  amat[lower.tri(amat,diag=FALSE)] <- res[i,]
  amat <- t(amat[index,index])
  g <- graph.adjacency(amat,mode="directed")
  E(g)$label <- NA
  pos <- layout.circle(g)
  rot <- matrix(c(0,1,-1,0),nrow=2,byrow=T)
  la <- pos %*% rot
  if(n.node == 2)la <- matrix(c(-1,0.5,1,0.5),,nrow=2,byrow=T)
  plot(g,layout=la,edge.arrow.size = 0.5,vertex.size=30,vertex.color = "blue",vertex.label=name)     #Change the color dependent upon the conditions
  mtext(paste(prob[i],"%"),line=0)
}

