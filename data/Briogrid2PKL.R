library(data.table)
I<-fread("C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/BIOGRID-ALL-4.4.222.tab3.txt")
HS.PPI<-I[`Organism Name Interactor A`=="Homo sapiens"][`Organism Name Interactor B`=="Homo sapiens"][`Experimental System Type`=="physical"][,c("Official Symbol Interactor A","Official Symbol Interactor B")]

Genes<-unique(c(unlist(HS.PPI[,"Official Symbol Interactor A"]),unlist(HS.PPI[,"Official Symbol Interactor B"])))
HS.PPI.ADJ_M<-matrix(0,length(Genes),length(Genes))
colnames(HS.PPI.ADJ_M)<-Genes
rownames(HS.PPI.ADJ_M)<-Genes
apply(HS.PPI,1,function(x){HS.PPI.ADJ_M[x[1],x[2]]<<-1;HS.PPI.ADJ_M[x[2],x[1]]<<-1})
library("reticulate")
#py_save_object(HS.PPI.ADJ_M, "C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/BIOGRID-HS-4.4.222.pkl", pickle = "pickle")
#write.table(HS.PPI.ADJ_M, "C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/BIOGRID-HS-4.4.222.csv", 
#            sep=",")


library("rtracklayer")
library("GenomicRanges")
Exons<-rtracklayer::readGFF("C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/Ressources/gencode.v19.annotation.gtf")
Exons<-Exons[Exons$type=="gene" & Exons$gene_type=="protein_coding",]
BIM<-read.table("C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/OUT_Filtered.bim",col.names = c("chr","name","cm","pos","ref","alt"))
Exons.GR<-makeGRangesFromDataFrame(Exons,keep.extra.columns = TRUE)
seqlevels(Exons.GR)<-gsub("chr","",seqlevels(Exons.GR))
names(Exons.GR)<-mcols(Exons.GR)$gene_name

BIM.GR<-makeGRangesFromDataFrame(BIM,seqnames.field = "chr",start.field = "pos",end.field = "pos",keep.extra.columns = TRUE)
write.table(mcols(BIM.GR[unique(to(findOverlaps(Exons.GR,BIM.GR)))])$name,"C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/SNPS_in_Genes.txt",col.names=FALSE)

ALL_GENES<-names(Exons.GR)
#SNPS_IN_GENES<-data.frame(

SNPS_IN_GENES.GR<-BIM.GR[unique(to(findOverlaps(Exons.GR,BIM.GR)))]
Genes.GR_WITH_SNPS<-Exons.GR[unique(from(findOverlaps(Exons.GR,BIM.GR)))]

GenesWithInterectionAndSNPS<-intersect(Genes,names(Genes.GR_WITH_SNPS))
HS.PPI.ADJ_M.SNPS<-HS.PPI.ADJ_M[GenesWithInterectionAndSNPS,GenesWithInterectionAndSNPS]
Genes.GR_WITH_SNPS.WITH_PPI<-Genes.GR_WITH_SNPS[GenesWithInterectionAndSNPS]
SNPS_IN_GENES_WITH_PPI.GR<-BIM.GR[unique(to(findOverlaps(Genes.GR_WITH_SNPS.WITH_PPI,BIM.GR)))]

rm(HS.PPI.ADJ_M)
rm(HS.PPI)
rm(I)
rm(Exons)
rm(Exons.GR)
rm(ALL_GENES)
rm(Genes)
rm(BIM)
rm(BIM.GR)
gc()

#M<-matrix(FALSE,length(SNPS_IN_GENES_WITH_PPI.GR),length(Genes.GR_WITH_SNPS.WITH_PPI))
#
#colnames(M)<-names(Genes.GR_WITH_SNPS.WITH_PPI)
#M<-data.frame(M)

ZE_FILE<-"C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/test.py"

cat("import numpy\nimport pickle\n",file=ZE_FILE)
cat("GENES_x_SNPS = {}\n",file=ZE_FILE,append = TRUE)

for  (x in names(Genes.GR_WITH_SNPS.WITH_PPI)){
   
   x<-"PENK"
   MATCHES<-rep("False",length(SNPS_IN_GENES_WITH_PPI.GR))
   SNPS_INDEX<-to(findOverlaps(
     Genes.GR_WITH_SNPS.WITH_PPI[x],
     SNPS_IN_GENES_WITH_PPI.GR))
   MATCHES[SNPS_INDEX]<-"True"
  if(sum(MATCHES=="True")==0){
    cat(x,"\n")
  }
   
      cat("GENES_x_SNPS['",x,"']=numpy.array([",paste(MATCHES,collapse=","),"])\n",sep="",file=ZE_FILE,append = TRUE)
#   M[,x]<-MATCHES
#   cat(dim(M))
  
  
}
cat("with open('mypicklefile', 'wb') as f1:\n",file=ZE_FILE,append = TRUE)
cat("\tpickle.dump(L, f1)\n\n",file=ZE_FILE,append = TRUE)
   
   
#write.table(M,
#            "C:/Users/Taing/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/taingl/networkGWAS/data/GenesWithSNPS.txt",
#            sep="\t",
#            quote=FALSE,
#            row.names=FALSE)

273462