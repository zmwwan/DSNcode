target = as.matrix(read.csv("DDI_target.csv", header = T, row.names = 1))
chem = as.matrix(read.csv("DDI_structure.csv", header = T, row.names = 1))
cell = as.matrix(read.csv("DDI_IC50.csv", header = T, row.names = 1))
colnames(cell)<-rownames(cell)
colnames(chem)<-rownames(chem)
rownames(target)<-rownames(chem)
colnames(target)<-rownames(target)

library(SNFtool)
K = 20; # number of neighbors, usually (10~30)
alpha = 0.5; # hyperparameter, usually (0.3~0.8)
T = 20; # Number of Iterations, usually (10~20)

affini_chem = affinityMatrix(1-chem, K, alpha)
affini_target = affinityMatrix(1-target, K, alpha)
affini_fusion = SNF(list(affini_chem,affini_target), K, T)

## data type contribution
affini_list = NULL
affini_list = c(affini_list, list(affini_fusion))
affini_list = c(affini_list, list(affini_chem))
affini_list = c(affini_list, list(affini_target))
contribution_matrix = affini_cell
contribution_matrix[which(contribution_matrix != 0)] = 0
for (i in 1:nrow(contribution_matrix))
{
  for (j in i:nrow(contribution_matrix))#上三角
  {
    if (j != i)
    {
      temp = matrix(data = 0,nrow = 3,ncol = 1)
      temp[1] = affini_list[[1]][i,j]
      temp[2] = affini_list[[2]][i,j]
      temp[3] = affini_list[[3]][i,j]
      rank_temp = order(-temp)#由大到小降序
      temp = temp[rank_temp]
      diff = (temp[1]-temp[2])/temp[2]
      diff = c(diff,(temp[2]-temp[3])/temp[3])
      if (diff[1]>0.1)# if the first was more than 10% higher than the second, the edge was attributed to the first data type.
      {
        contribution_matrix[i,j] = rank_temp[1]
      }else# if the difference between the first and the second is less than 10%, then look at the difference between the second and the third
      {
        if (diff[2]>0.1)# if the difference between the second and the third is greater than 10%, the edge was attributed to both the first and the second data type
        {
          a = c(rank_temp[1],rank_temp[2])
          if (length(intersect(a,c(1,2)))==2)
          {
            contribution_matrix[i,j] = 4
          }
          if (length(intersect(a,c(1,3)))==2)
          {
            contribution_matrix[i,j] = 5
          }
          if (length(intersect(a,c(2,3)))==2)
          {
            contribution_matrix[i,j] = 6
          }
        }else# if the difference between the second and the third is less than 10%, the edge was attributed to all three data types
        {
          contribution_matrix[i,j] = 7
        }
      }
    }
  }
}
contribution_sta = table(contribution_matrix[which(contribution_matrix!=0)])
pie_col=c("skyblue","green","pink","yellow","orange","red","bisque","burlywood","cadetblue","gray","chocolate","cornsilk","firebrick","khaki","black")
pie(contribution_sta,border = NA,col = pie_col,clockwise = T,cex=0.8)


## NMIs among drug similarity networks
combine = list(affini_cell,affini_fusion,affini_chem, affini_target)
#combine = list(cell,affini_fusion,chem, target)
Comcordance_matrix = concordanceNetworkNMI(combine,32)
rownames(Comcordance_matrix) = c("cell","affini_fusion","chem","target")
colnames(Comcordance_matrix) = c("cell","affini_fusion","chem","target")


drugATC<-read.csv("E:/药物/DSN/Imput_data/drugATC.csv", header = T)
drugATC<-drugATC[-which(drugATC$fusionATC==""),]
drugATC[,"layer"]<-substr(drugATC$fusionATC,1,4)
atc<-unique(drugATC$layer)
feature_ATC<-data.frame(matrix(0,nrow = length(atc),ncol = length(unique(drugATC$drug_name))))
rownames(feature_ATC)<-atc
colnames(feature_ATC)<- unique(drugATC$drug_name)
for (drug in unique(drugATC$drug_name)) {
  feature<-drugATC$layer[which(drugATC$drug_name==drug)]
  feature_ATC[feature,drug]<-1
}
ATC_matrix<-data.frame(matrix(0,nrow = 43,ncol = 43))
rownames(ATC_matrix)<-colnames(feature_ATC);colnames(ATC_matrix)<-colnames(feature_ATC)
for (drug in unique(drugATC$drug_name)) {
  feature<-drugATC$layer[which(drugATC$drug_name==drug)]
  com<-drugATC$drug_name[which(drugATC$layer %in% feature)]
  ATC_matrix[com,drug]<-1
}
affini_ATC = affinityMatrix(1-as.matrix(ATC_matrix), K, alpha)
drugname_ATC<-unique(drugATC$drug_name)
combine = list(affini_ATC,affini_fusion[drugname_ATC,drugname_ATC],affini_chem[drugname_ATC,drugname_ATC], affini_target[drugname_ATC,drugname_ATC])
Comcordance_matrix = concordanceNetworkNMI(combine,32)
rownames(Comcordance_matrix) = c("ATC","affini_fusion","chem","target")
colnames(Comcordance_matrix) = c("ATC","affini_fusion","chem","target")

C = 32 # number of clusters
group = spectralClustering(affini_fusion, C)
drug_name = rownames(chem)
CIDandgroup = cbind(drug_name,group)
N = 276 # number of all drugs

## ATC enrichment analysis
CIDandatc = read.csv("E:/药物/DSN/Imput_data/drugATC.csv", header = T)
CIDandatc<-CIDandatc[-which(CIDandatc$fusionATC==""),]
CIDandatc[,"layer"]<-substr(CIDandatc$fusionATC,1,4)
C<-length(unique(CIDandatc$layer))
#原始数据分四类
group = spectralClustering(affini_fusion, C)
#group = spectralClustering(affini_fusion, 32)
#group = spectralClustering(affini_chem, 32)
#group = spectralClustering(affini_target, 32)
Group_matrix<-data.frame(matrix(nrow = 276,ncol = 2))
colnames(Group_matrix)<-c("drug_name","group")
Group_matrix$drug_name<-rownames(affini_fusion)
Group_matrix$group<-group
Group_matrix<-Group_matrix[which(Group_matrix$drug_name %in% CIDandatc$drug_name),]
Group_matrix[,"ATC"]<-CIDandatc[1:43,"layer"]
#截取片段分四类
drug_atc<-intersect(Group_matrix$drug_name,CIDandatc$drug_name)
group = spectralClustering(affini_fusion[drug_atc,drug_atc], C)
Group_matrix<-data.frame(matrix(nrow = length(group),ncol = 2))
colnames(Group_matrix)<-c("drug_name","group")
Group_matrix$drug_name<-drug_atc
Group_matrix$group<-group
Group_matrix[,"ATC"]<-CIDandatc[1:43,"layer"]
#ARIindex
class<-data.frame(matrix(0,nrow = length(unique(group)),ncol = C))
rownames(class)<-unique(group)
colnames(class)<-unique(CIDandatc$layer)
C1<-0;C2<-0;C3<-0
for (i in 1:C) {
  iGroup<-Group_matrix[which(Group_matrix[,3]==colnames(class)[i]),]
  set<-unique(iGroup$group)
  for (num in set) {
    class[num,i]<-length(iGroup$group[which(iGroup$group==num)])
    C1<-C1+choose(class[num,i],2)
  }
}
class[nrow(class)+1,]<-colSums(class)
rownames(class)[nrow(class)]<-"Colsums"
class[,ncol(class)+1]<-rowSums(class)
colnames(class)[ncol(class)]<-"Rowsums"
for (i in 1:C) {
  C2<-C2+choose(class[nrow(class),i],2)
  C3<-C3+choose(class[i,ncol(class)],2)
}
C4<-choose(class[nrow(class),ncol(class)],2)
ARI_atc<-(C1-(C2*C3)/C4)/((C2+C3)/2-(C2*C3)/C4)

## BPS enrichment analysis
CIDandatc = read.csv("E:/药物/DSN/Imput_data/drugBPS.csv", header = T)
CIDandatc<-CIDandatc[-which(CIDandatc$PHARMACOLOGY_Target_Classification==""),]
C<-length(unique(CIDandatc$PHARMACOLOGY_Target_Classification))
#原始数据分四类
group = spectralClustering(affini_fusion, 32)
#group = spectralClustering(affini_chem, 32)
#group = spectralClustering(affini_target, 32)
Group_matrix<-CIDandatc[which(CIDandatc$drug_name %in% Group_matrix$drug_name),]
for (drug in Group_matrix$drug_name) {
  Group_matrix[which(Group_matrix$drug_name==drug),"group"]<-group[which(rownames(affini_fusion)==drug)]
}
class<-data.frame(matrix(0,nrow = length(unique(group)),ncol = C))
rownames(class)<-unique(group)
colnames(class)<-unique(CIDandatc$PHARMACOLOGY_Target_Classification)
C1<-0;C2<-0;C3<-0
for (i in 1:C) {
  iGroup<-Group_matrix[which(Group_matrix[,2]==colnames(class)[i]),]
  set<-unique(iGroup$group)
  for (num in set) {
    class[num,i]<-length(iGroup$group[which(iGroup$group==num)])
    C1<-C1+choose(class[num,i],2)
  }
}
class[nrow(class)+1,]<-colSums(class)
rownames(class)[nrow(class)]<-"Colsums"
class[,ncol(class)+1]<-rowSums(class)
colnames(class)[ncol(class)]<-"Rowsums"
for (i in 1:C) {
  C2<-C2+choose(class[nrow(class),i],2)
  C3<-C3+choose(class[i,ncol(class)],2)
}
C4<-choose(class[nrow(class),ncol(class)],2)
ARI_atc<-(C1-(C2*C3)/C4)/((C2+C3)/2-(C2*C3)/C4)

## MeSH enrichment analysis
CIDandatc = read.csv("E:/药物/DSN/Imput_data/drugMeSH.csv", header = T)
CIDandatc<-CIDandatc[-which(CIDandatc$ChemicalsandDrugsCategory==""),]
C<-length(unique(CIDandatc$ChemicalsandDrugsCategory))
group = spectralClustering(affini_fusion, 32)
Group_matrix<-CIDandatc[which(CIDandatc$drug_name %in% Group_matrix$drug_name),]
for (drug in Group_matrix$drug_name) {
  Group_matrix[which(Group_matrix$drug_name==drug),"group"]<-group[which(rownames(affini_fusion)==drug)]
}
class<-data.frame(matrix(0,nrow = length(unique(group)),ncol = C))
rownames(class)<-unique(group)
colnames(class)<-unique(CIDandatc$ChemicalsandDrugsCategory)
C1<-0;C2<-0;C3<-0
for (i in 1:C) {
  iGroup<-Group_matrix[which(Group_matrix[,2]==colnames(class)[i]),]
  set<-unique(iGroup$group)
  for (num in set) {
    class[num,i]<-length(iGroup$group[which(iGroup$group==num)])
    C1<-C1+choose(class[num,i],2)
  }
}
class[nrow(class)+1,]<-colSums(class)
rownames(class)[nrow(class)]<-"Colsums"
class[,ncol(class)+1]<-rowSums(class)
colnames(class)[ncol(class)]<-"Rowsums"
for (i in 1:C) {
  C2<-C2+choose(class[nrow(class),i],2)
  C3<-C3+choose(class[i,ncol(class)],2)
}
C4<-choose(class[nrow(class),ncol(class)],2)
ARI_atc<-(C1-(C2*C3)/C4)/((C2+C3)/2-(C2*C3)/C4)

dist_fusion<-1-affini_fusion[intersect(CIDandatc$drug_name,colnames(affini_fusion)),intersect(CIDandatc$drug_name,colnames(affini_fusion))]
CIDandgroup_ATC<-CIDandgroup[which(CIDandgroup[,1] %in% CIDandatc$drug_name),]
dunn<-c()
C=32
for(i in 1:C){
  temp_CID = CIDandgroup_ATC[which(CIDandgroup_ATC[,2] == i),"drug_name"]
  if(length(temp_CID)>0){
    temp_dist<-dist_fusion[temp_CID,temp_CID]
    inclus_max<-max(temp_dist)
    out_dist<-dist_fusion[temp_CID,-which(colnames(dist_fusion) %in% temp_CID)]
    outclus_min<-min(out_dist)
    dunn<-c(dunn,outclus_min/inclus_max)
  }
}
min(dunn) 


                     