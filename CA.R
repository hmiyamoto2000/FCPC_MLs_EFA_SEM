#For CA (Fig.S10)
#Create the folder for Correlation network
dir.create("./CA_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Make the file with the components selected by ELA, AA, RF, and XGBoost
#File name: FCPC_MLs.txt
d <- read.delim("FCPC_MLs.txt", header=TRUE, row.name=1, sep="\t", fileEncoding="UTF-8",check.names=F)
d$name <- NULL
dat_t <- t(d)
no_zero<-subset(d,apply(d, 1, sum)>0)
dat_spearman<-cor(dat_t,dat_t,method="spearman")
write.csv(dat_spearman, "./CA_result/spearman_Pra.csv")

#r>0.4 (positive correlation)
dir.create("./CA_result/Positive_Cor", showWarnings = TRUE, recursive = FALSE, mode = "0777")

selected_data_spearman<-read.csv("FCPC_plus.csv", row.names=1, header=T)

selected_data_cor_mtrx<-ifelse(selected_data_spearman>0.4,1,0) #>0.7:1, <0.7:0
selected_dataz <- as.data.frame(selected_data_spearman) 
selected_datazz <- colnames(selected_data_spearman) 
selected_dataa <- sprintf("dl n=%s\nformat = fullmatrix\nlabels:\n",nrow(selected_data_cor_mtrx)) 
cat(selected_dataa, file="Spe.dl", sep=",") 
selected_datab <- sprintf("%s",selected_datazz)
cat(selected_datab, file="Spe.dl", append=T, sep=",") 
selected_datac <- sprintf("\nselected_data:\n") 
cat(selected_datac, file="Spe.dl", append=T)
write.table(selected_data_cor_mtrx,"./CA_result/Positive_Cor/Spe.dl",col.names=FALSE, row.names=FALSE ,quote=FALSE, sep=" ", append=T) 

#Visualize the data stored in the file "Spe.dl" using the Gephi (http://gephi.org)

#|r|>0.4 (negative correlation), r>0: 0 
#Change as absolute number
selected_data_spearman<-read.csv("FCPC_minus.csv", row.names=1, header=T)

dir.create("./CA_result/Negative_Cor", showWarnings = TRUE, recursive = FALSE, mode = "0777")

selected_data_cor_mtrx<-ifelse(selected_data_spearman>0.4,1,0)
selected_dataz <- as.data.frame(selected_data_spearman) 
selected_datazz <- colnames(selected_data_spearman) 
selected_dataa <- sprintf("dl n=%s\nformat = fullmatrix\nlabels:\n",nrow(selected_data_cor_mtrx)) 
cat(selected_dataa, file="Spe.dl", sep=",") 
selected_datab <- sprintf("%s",selected_datazz)
cat(selected_datab, file="Spe.dl", append=T, sep=",") 
selected_datac <- sprintf("\nselected_data:\n") 
cat(selected_datac, file="Spe.dl", append=T)
write.table(selected_data_cor_mtrx,"./CA_result/Negative_Cor/Spe.dl",col.names=FALSE, row.names=FALSE ,quote=FALSE, sep=" ", append=T) 

#Visualize the data stored in the file "Spe.dl" using the Gephi (http://gephi.org)