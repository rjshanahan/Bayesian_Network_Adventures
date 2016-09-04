#Richard Shanahan  
#https://github.com/rjshanahan  
#6 August 2016

###### INFS 5094: Project Breast Cancer Genes

# load required packages
library(dplyr)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(e1071)
library(plotly)

library(gRain)
library(bnlearn)
library(pcalg)

library(network)
library(igraph)
library(ggnet)
library(intergraph)



#seet ggplot2 theme
theme = theme_set(theme_minimal())
theme = theme_update(legend.position="top")

###### 0.1 READ & PRE-PROCESS BRCA-50 DATASET ###### 
breast_cancer <- read.csv('YOUR_CSV.csv',
                          header=T,
                          sep=",",
                          quote='"',
                          # colClasses=c(
                          # ),
                          strip.white=T,
                          stringsAsFactors=F,
                          fill=T)

#inspect
str(breast_cancer)
describe(breast_cancer)

#check for duplicate records
nrow(breast_cancer) - length(unique(breast_cancer))

#check if there are any missing values
colSums(is.na(breast_cancer)) 

#recode 'class' as character
breast_cancer$class_char <- ifelse(breast_cancer$class == "C",
                                   "cancer",
                                   "normal")

#add id column for visualisations
breast_cancer$id <- 1:nrow(breast_cancer)

str(breast_cancer)



###### 0.2 CREATE BINARISED VERSION OF THE DATASET ###### 

#create new binarised dataset 

#function to discretise based on average value of gene
binariser <- function(myCol) {
  
  # myMean <- mean(myCol,
  #                na.rm=T)
  
  myMean <- sum(colSums(breast_cancer[1:50]))/(nrow(breast_cancer)*ncol(breast_cancer[1:50]))
  
  myCol <- ifelse(myCol > myMean,
                  1,
                  0)
  return(myCol)
}

#apply function to each gene variable
breast_cancer_dc <- breast_cancer %>%
  mutate_each(funs(binariser),
              -class, -class_char, -id)



#binarise class variable
breast_cancer_dc$class_bin <- ifelse(breast_cancer_dc$class_char == "cancer",
                                     1,
                                     0)

str(breast_cancer_dc)


###### 0.3 VISUALISE DISTRIBUTIONS ###### 

#class distribution
ggplot(data = breast_cancer, 
       aes(x=class_char,
           fill=class_char)) + 
  geom_bar() +
  ggtitle("BRCA-50 Breast Cancer 'Class' Distribution")


#reshape dataset for boxplot representation
breast_cancer.m <- melt(select(breast_cancer, -class, -class_char),
                        id.var="id")

#distribution of other variables
ggplot(data = breast_cancer.m, 
       aes(x=value),
       fill=class) + 
  #geom_histogram(aes(fill=factor(variable)),
  #               binwidth = 1) +
  geom_bar(aes(fill=factor(variable))) +
  facet_wrap("variable") +
  guides(fill=FALSE) +
  ggtitle("BRCFA-50 Breast Cancer Dataset Distributions (bin=1)")


#create boxplots of each variable in same graphic
p <- ggplot(data = breast_cancer.m, 
            aes(x=variable, y=log(value))) + 
  geom_boxplot(aes(fill=variable)) + 
  xlab("Breast Cancer Gene") + 
  ylab("Value") + 
  guides(fill=FALSE) +
  ggtitle("BRCA-50 Breast Cancer Dataset Boxplots - log transformed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#guides(fill=guide_legend(title="Attribute Legend"))

ggplotly(p)



#reshape dataset for boxplot representation - binarised
breast_cancer_dc.m <- melt(select(breast_cancer_dc, -class, -class_char),
                           id.var="id")

#distribution of other variables
ggplot(data = breast_cancer_dc.m, 
       aes(x=value),
       fill=class) + 
  #geom_histogram(aes(fill=factor(variable)),
  #               binwidth = 1) +
  geom_bar(aes(fill=factor(variable))) +
  facet_wrap("variable") +
  guides(fill=FALSE) +
  ggtitle("BRCFA-50 Breast Cancer Binarised Dataset")



###### 0.4 CORRELATIONS ###### 

#correlations for standard dataset
breast_cancer.cor <- round(cor(select(breast_cancer, -class, -class_char), 
                               use = "complete.obs",
                               y=NULL,
                               method = "pearson"), 2)

breast_cancer.cor

#correlations for binarised dataset
breast_cancer_dc.cor <- round(cor(select(breast_cancer_dc, -class, -id, -class_char), 
                                  use = "complete.obs",
                                  y=NULL,
                                  method = "pearson"), 2)

breast_cancer_dc.cor


# scatterplot for corrleations
pairs(select(breast_cancer, -class, -class_char),
      main="Scatterplot of Breast Cancer Numeric",
      pch = 20,
      col="goldenrod")





##### PART 1.1 - BUILD NETWORK - LEARN GLOBAL STRUCTURE WITH PC ALGORITHM #####

#list of sufficient statistics for conditional independence tests
suffStat <- list(C = breast_cancer.cor, 
                 n = nrow(breast_cancer))

#list of genes for naming network
gene_list <- c(id=colnames(breast_cancer[,1:50]))

#learn network (constraint based apporach) - using "pcSelect" algorithm by testing Conditional Independence of Gaussians via Fisher's Z 
pc.fit <- pc(suffStat, 
             indepTest = gaussCItest,
             #labels=colnames(select(breast_cancer, -class, -class_char, -id)),
             p=ncol(select(breast_cancer, -class, -class_char, -id)), 
             alpha = 0.01)

#plot network
#plot(pc.fit@graph,
#     cex.lab=0.5)

#dev.off()

##### PART 1.2 - CAUSAL INFERENCE/EFFECT OF NODES ON EBF1 #####

#causal inference Top 5 using IDA 

#empty dataframe
ida_EBF1 <- data.frame()

#loop through 50 genes and calculate the causal effect using IDA on each gene (inc. EBF1)
for (i in 1:50)
{
  
  newIDA <- idaFast(
          i, 
          40, 
          cov(select(breast_cancer, -class, -class_char, -id)), 
          pc.fit@graph)
  
  if(length(newIDA > 1))
  {
    #if there are multiple causal effects the average is is used
    newIDA.df <- data.frame(mean(newIDA))
  }
  else
  {
    newIDA.df <- data.frame(newIDA)
  }
  
  newIDA.df <- cbind(newIDA.df, i)
  
  ida_EBF1 <- rbind(ida_EBF1, newIDA.df)
  
}

#add gene name to the loop created dataframe
ida_EBF1$genes <- colnames(select(breast_cancer, -class, -class_char, -id))

#top 5 genes
EBF1_ida_top5 <- ida_EBF1 %>% 
  arrange(desc(mean.newIDA.)) %>% 
  head(n=5)

EBF1_ida_top5



##### PART 1.3 - VISUALISE NETWORK - EBF1 TOP 5 CAUSAL #####


#pretty graphs
myNet <-  igraph.from.graphNEL(pc.fit@graph, 
                               name = TRUE, 
                               weight = TRUE,
                               unlist.attrs = TRUE)

#add vertex names
V(myNet)$name <- colnames(select(breast_cancer, -class, -class_char, -id))

#convert to "network" object
myNet <- asNetwork(myNet)
#myNet <- asIgraph(myNet)

#create node attribute to colour EBF1 as focal gene
myGene = network.vertex.names(myNet)
myGene = ifelse(myGene %in% "EBF1" , 
                "EBF1",
                ifelse(myGene %in% c("ABCA10","HIF3A","ARHGAP20","TMEM220","CD300LG"),
                       "Top 5 Causal",
                       "Other Gene"))
myNet %v% "EBF1" = myGene


p <- ggnet2(myNet, 
            size = 15, 
            alpha = 0.8,
            #size = "degree",
            #size.legend = "Number of Degrees",
            #size.cut = 3,
            label = TRUE, 
            label.size = 3,
            color = "EBF1", 
            color.legend = "EBF1 Top Causal Genes",
            palette = c("EBF1" = "tomato", "Top 5 Causal" = "goldenrod", "Other Gene" = "steelblue"),
            label.color = "gray15",
            arrow.size = 8, 
            arrow.gap = 0.025,
            edge.color = 'gray50'
) +
  ggtitle('Bayesian Network for BRCA-50 Breast Cancer Genes inc. EBF1 Top 5 Causal') +
  guides(color = FALSE, size = FALSE) 
#  ggsave(paste0('Bayesian Network for BRCA-50 Breast Cancer Genes','.png'))

#interactive
#ggplotly(p)

#static
print(p)




##### PART 2.1 - MARKOV BLANKET FOR EBF1 #####

#genes in Markov Blanket for EBF1
MB.EBF1 <- learn.mb(select(breast_cancer, -class, -class_char, -id), 
                    node='EBF1', 
                    method='iamb', 
                    alpha=0.01)

MB.EBF1


##### PART 2.2 - VISUALISE NETWORK - EBF1 MARKOV BLANKET #####

#create node attribute to colour EBF1 as focal gene
myMB = network.vertex.names(myNet)
myMB = ifelse(myMB %in% "EBF1" , 
              "EBF1",
              ifelse(myMB %in% MB.EBF1,
                     "EBF1 Markov Blanket",
                     "Other Gene"))
myNet %v% "Markov Blanket" = myMB


p <- ggnet2(myNet, 
            size = 15, 
            alpha = 0.8,
            #size = "degree",
            #size.legend = "Number of Degrees",
            #size.cut = 3,
            label = TRUE, 
            label.size = 3,
            color = "Markov Blanket", 
            color.legend = "EBF1 Markov Blanket",
            palette = c("EBF1" = "tomato", "EBF1 Markov Blanket" = "goldenrod", "Other Gene" = "steelblue"),
            label.color = "gray15",
            arrow.size = 8, 
            arrow.gap = 0.025,
            edge.color = 'gray50'
) +
  ggtitle('Bayesian Network for BRCA-50: EBF1 Markov Blanket') 
#guides(color = FALSE, size = FALSE) 
#  ggsave(paste0('Bayesian Network for BRCA-50: EBF1 Markov Blanket','.png'))

#static
print(p)



##### PART 3.1 - BUILD NETWORK USING BINARISED DATASET #####


#list of sufficient statistics for conditional independence tests
suffStat_dc <- list(C = breast_cancer_dc.cor, 
                    n = nrow(breast_cancer_dc))

#list of genes for naming network
gene_list_dc <- c(id=colnames(select(breast_cancer_dc, -class, -class_char, -id)))

#learn network (constraint based apporach) - using "pcSelect" algorithm by testing Conditional Independence of Gaussians via Fisher's Z 
pc.fit_class <- pc(suffStat_dc, 
                   indepTest = gaussCItest,
                   #labels=colnames(select(breast_cancer, -class, -class_char, -id)),
                   p=ncol(select(breast_cancer_dc, -class, -class_char, -id)), 
                   alpha = 0.01)


##### PART 3.2 - IDENTIFY CLASS PARENT/CHILD SET WITH LOCAL STRUCTURE LEARNING#####


#learn causal strucure around gene EBF1 - using pcSimple
pcS_class <- pcSelect(y=select(breast_cancer_dc, class_bin), 
                      dm=select(breast_cancer_dc, -class, -class_char, -id),
                      alpha=0.05)
pcS_class


#find parent/child set using local structure learning algorithm - using SI-HITON
BRCA_class_set <- learn.nbr(select(breast_cancer_dc, -class, -class_char, -id), 
                            'class_bin', 
                            method='si.hiton.pc',
                            #test='x2',      #specify conditional independence test
                            alpha=0.01)
BRCA_class_set


##### PART 3.3 - VISUALISE BINARISED NETWORK - CLASS PARENT/CHILD SET #####

#pretty graphs
myNet_class <-  igraph.from.graphNEL(pc.fit_class@graph, 
                                     name = TRUE, 
                                     weight = TRUE,
                                     unlist.attrs = TRUE)


#add vertex names
V(myNet_class)$name <- colnames(select(breast_cancer_dc, -class, -class_char, -id))

#convert to "network" object
myNet_class <- asNetwork(myNet_class)
#myNet_class <- asIgraph(myNet_class)

#create node attribute to colour EBF1 as focal gene
myClass = network.vertex.names(myNet_class)
myClass = ifelse(myClass %in% "class_bin" , 
                 "class",
                 #ifelse(myClass %in% c("MYOM1","FIGF","SCARA5","CA4","LYVE1","ATP1A2","HIF3A","C2orf40"),
                 #ifelse(myClass %in% c("GPAM", "LEP", "ATP1A2"),
                 ifelse(myClass %in% c("FIGF", "CD300LG", "ARHGAP20", "KLHL29", "SCARA5", "MAMDC2", "CXCL2", "ATP1A2", "TMEM220"),
                                         "parent/child",
                        "other"))
myNet_class %v% "class" = myClass



p2 <- ggnet2(myNet_class, 
             #size = 15, 
             alpha = 0.8,
             size = "degree",
             size.legend = "Number of Degrees",
             size.cut = 3,
             label = TRUE, 
             label.size = 3,
             color = "class", 
             color.legend = "Class Parent/Child Set",
             palette = c("class" = "tomato", "parent/child" = "goldenrod", "other" = "steelblue"),
             label.color = "gray15",
             arrow.size = 8, 
             arrow.gap = 0.025,
             edge.color = 'gray50'
) +
  ggtitle('Bayesian Network with CLASS for Binarised BRCA-50 Breast Cancer Genes') 
#guides(color = FALSE, size = FALSE) 
#  ggsave(paste0('Bayesian Network with CLASS for Binarised BRCA-50 Breast Cancer Genes','.png'))

#interactive
#ggplotly(p)

#static
print(p2)




##### PART 3.4 - NAIVE BAYES CLASSIFICATION USING ALL BINARISED VARIABLES #####

#train naiveBayes model using e1071 package
nb_all <- naiveBayes(as.factor(class_bin)~.,
                     data=select(breast_cancer_dc, -class, -class_char, -id),
                     laplace = 0)

#function to predict
myPredictor <-function(nb_object, newdata) {
  predict(nb_object,
          newdata=newdata)  
}

#prediction object
myPredict <- myPredictor(nb_all,
                         select(breast_cancer_dc, -class, -class_char, -id, -class_bin))


#confusion matrix
table(prediction=myPredict,
      true=breast_cancer_dc$class_bin)

#accuracy
mean(myPredict == breast_cancer_dc$class_bin)



#10 fold cross-validation
breast_cancer_xval <- tune(naiveBayes,
                           as.factor(class_bin)~.,
                           data = select(breast_cancer_dc, -class, -class_char, -id),
                           predict.func = myPredictor
)

#review Cross Validation output
breast_cancer_xval$best.parameters
breast_cancer_xval$best.performance
breast_cancer_xval$performances
breast_cancer_xval$best.model


#use cross validated best model
nb_all_xv <- breast_cancer_xval$best.model


#view model summary
nb_all_xv
summary(nb_all_xv)


#prediction object using Xvalidated model
myPredict <- myPredictor(nb_all_xv,
                         select(breast_cancer_dc, -class, -class_char, -id, -class_bin))


#confusion matrix
table(prediction=myPredict,
      true=breast_cancer_dc$class_bin)


##### PART 3.5 - NAIVE BAYES CLASSIFICATION USING ONLY CLASS PARENT/CHILD BINARISED VARIABLES #####

#train naiveBayes model using e1071 package
nb_set <- naiveBayes(as.factor(class_bin)~.,
#                     data=select(breast_cancer_dc, class_bin, MYOM1,FIGF,SCARA5,CA4,LYVE1,ATP1A2,HIF3A,C2orf40),
#                     data=select(breast_cancer_dc, class_bin, GPAM,LEP,ATP1A2),  
                     data=select(breast_cancer_dc, class_bin, FIGF, CD300LG, ARHGAP20, KLHL29, SCARA5, MAMDC2, CXCL2, ATP1A2, TMEM220),   
                     laplace = 0)



#view model summary
nb_set
summary(nb_set)


#prediction object
myPredict <- myPredictor(nb_set,
                         #select(breast_cancer_dc, class_bin, MYOM1,FIGF,SCARA5,CA4,LYVE1,ATP1A2,HIF3A,C2orf40))
                         select(breast_cancer_dc, class_bin, FIGF, CD300LG, ARHGAP20, KLHL29, SCARA5, MAMDC2, CXCL2, ATP1A2, TMEM220))


#confusion matrix
table(prediction=myPredict,
      true=breast_cancer_dc$class_bin)

#accuracy
mean(myPredict == breast_cancer_dc$class_bin)


#10 fold cross-validation
breast_cancer_set_xval <- tune(naiveBayes,
                               as.factor(class_bin)~.,
#                              data=select(breast_cancer_dc, class_bin, MYOM1,FIGF,SCARA5,CA4,LYVE1,ATP1A2,HIF3A,C2orf40),
                               data=select(breast_cancer_dc, class_bin, FIGF, CD300LG, ARHGAP20, KLHL29, SCARA5, MAMDC2, CXCL2, ATP1A2, TMEM220),          
                               predict.func = myPredictor
)

#review Cross Validation output
breast_cancer_set_xval$best.parameters
breast_cancer_set_xval$best.performance
breast_cancer_set_xval$performances
breast_cancer_set_xval$best.model


#use cross validated best model
nb_set_xv <- breast_cancer_set_xval$best.model

#prediction object using Xvalidated model
myPredict <- myPredictor(nb_set_xv,
                         select(breast_cancer_dc, -class, -class_char, -id, -class_bin))


#confusion matrix
table(prediction=myPredict,
      true=breast_cancer_dc$class_bin)


##### PART 4.1 - LEARN NETWORK OF SELECTED GENES #####

#correlations for dataset of selected variables
breast_cancer_sel.cor <- round(cor(select(breast_cancer_dc, class_bin, BTNL9, CD300LG, IGSF10, ABCA9), 
                                   use = "complete.obs",
                                   y=NULL,
                                   method = "pearson"), 2)

breast_cancer_sel.cor


#list of sufficient statistics for conditional independence tests
suffStat_sel <- list(C = breast_cancer_sel.cor, 
                     n = nrow(breast_cancer_dc))

#list of genes for naming network
gene_list_sel <- c(id=colnames(select(breast_cancer_dc, class_bin, BTNL9, CD300LG, IGSF10, ABCA9)))

#learn network (constraint based apporach) - using "pcSelect" algorithm by testing Conditional Independence of Gaussians via Fisher's Z 
pc.fit_sel <- pc(suffStat_sel, 
                 indepTest = gaussCItest,
                 #labels=colnames(select(breast_cancer, -class, -class_char, -id)),
                 p=ncol(select(breast_cancer_dc, class_bin, BTNL9, CD300LG, IGSF10, ABCA9)), 
                 alpha = 0.01)


##### PART 4.2 - VISUALISE NETWORK - SELECTED BINARISED GENES + CLASS #####


#pretty graphs
myNet_sel <-  igraph.from.graphNEL(pc.fit_sel@graph, 
                                   name = TRUE, 
                                   weight = TRUE,
                                   unlist.attrs = TRUE)

#add vertex names
V(myNet_sel)$name <- colnames(select(breast_cancer_dc, class_bin, BTNL9, CD300LG, IGSF10, ABCA9))

#convert to "network" object
myNet_sel <- asNetwork(myNet_sel)
#myNet <- asIgraph(myNet)

#create node attribute to colour class as focal gene
mySel = network.vertex.names(myNet_sel)
mySel = ifelse(mySel %in% "class_bin", 
               "class",
               "other")
myNet_sel %v% "Selected" = mySel


p <- ggnet2(myNet_sel, 
            #size = 15, 
            alpha = 0.8,
            size = "degree",
            size.legend = "Number of Degrees",
            size.cut = 3,
            label = TRUE, 
            label.size = 3,
            color = "Selected", 
            color.legend = "Select Genes + Class",
            palette = c("class" = "tomato", "other" = "steelblue"),
            label.color = "gray15",
            arrow.size = 8, 
            arrow.gap = 0.025,
            edge.color = 'gray50'
) +
  ggtitle('Bayesian Network for Selected BRCA-50 Breast Cancer Genes - Binarised') +
  guides(color = FALSE, size = FALSE) 
#  ggsave(paste0('Bayesian Network for Selected BRCA-50 Breast Cancer Genes - Binarised','.png'))

#interactive
#ggplotly(p)

#static
print(p)

##### PART 4.3 - BUILD NETWORK MANUALLY + VISUALISE - SELECTED BINARISED GENES + CLASS #####


#construct graph manually
myNet_sel2 <- make_empty_graph(directed = TRUE) %>%
  add_vertices(nv=5) %>%
  set_vertex_attr("name", value = c("class_bin","BTNL9","CD300LG","IGSF10","ABCA9")) %>%
  add_edges(c(2,3,
              2,5,
              3,1,
              1,4,
              4,5   
  )) 

#convert to network for plotting
myNet_sel2 <- asNetwork(myNet_sel2)

#label vertices
mySel2 = network.vertex.names(myNet_sel2)
mySel2 = ifelse(mySel2 %in% "class_bin", 
                "class",
                "other")
myNet_sel2 %v% "Selected2" = mySel2


p <- ggnet2(myNet_sel2, 
            #size = 15, 
            alpha = 0.8,
            #size = "degree",
            #size.legend = "Number of Degrees",
            #size.cut = 3,
            label = TRUE, 
            label.size = 3,
            color = "Selected2", 
            color.legend = "Select Genes + Class",
            palette = c("class" = "tomato", "other" = "steelblue"),
            label.color = "gray15",
            arrow.size = 8, 
            arrow.gap = 0.025,
            edge.color = 'gray50'
) +
  ggtitle('Bayesian Network for Selected BRCA-50 Breast Cancer Genes - Binarised') 

print(p)


##### PART 4.4 - CONDITIONAL PROBABILITY TABLES #####


#dataset for selected gene expressions
breast_cancer_dc_sel <- breast_cancer_dc %>% 
  select(class_bin, BTNL9, CD300LG, IGSF10, ABCA9) 
# %>% lapply(factor)

#build DAG using conditional probability representations
myDAG <- dag(~ class_bin:CD300LG  + 
               IGSF10:class_bin + 
               CD300LG:BTNL9 + 
               ABCA9:BTNL9 +
               ABCA9:IGSF10)

plot(myDAG)

#add data to graoh object
myDAG.data <- grain(myDAG, data = breast_cancer_dc_sel, smooth = 0.1)


#Estimate the probability of having cancer when the expression level of CD300LG is low and the expression level of BTNL9 is high. 
querygrain(myDAG.data, nodes=c("class_bin", "CD300LG", "BTNL9"), type="conditional") 


#Estimate the probability of the four genes in the network having high expression levels
querygrain(myDAG.data, nodes=c("CD300LG", "BTNL9", "IGSF10", "ABCA9"), type="joint") 

