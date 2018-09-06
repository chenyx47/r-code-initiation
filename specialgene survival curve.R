
#gene survival analysis
#data example(express_matrix;survival_data;genename)
#survival curve
expressdata<-stageIV.rsem                                                                           #输入特定的表达矩阵(first colume is genenames)
genenames<-as.character(stageIV.dwgene$genename)                                                    #选取目标基因名字(input an genename vector (class:character))
analysislabel<-"stageIV.dwgene"                                                                     #input an analysis label
tumorsurvival<-                                                                                     #input a survival data including status(1or2),time,sampleID, 
  
genename1<-sapply(strsplit(genenames,"\\|"),"[",1)     
tumorsurvival<-survivalphenotype[survivalphenotype$tissue=="tumor",]                                #选取肿瘤组织资料（NA的出现？？）
tumorsurvival<-tumorsurvival[!is.na(tumorsurvival$submitter_id.samples),]
tumorsurvival$patient<-substr(tumorsurvival[,1],1,12)                                               #在生存资料中创建病人列表
# tumorsurvival$stage<-tumorsurvival$tumor_stage.diagnoses                                            #选取诊断时的肿瘤分期
survivalpatient<-aggregate(time~patient+status+stage,tumorsurvival,mean)                            #以病人为个体选取生存资料
patient<-unique(substr(names(expressdata),1,12))                                                    #选取表达矩阵的病人资料
survivalnew<-survivalpatient[survivalpatient$patient%in%patient, ]                                  #以表达矩阵中病人资料为基础选出生存资料
result_all<-data.frame()
for(i in 1:length(genenames)){                                                                      #针对每个基因做生存分析的pvalue                                                     
  geneexpression<-expressdata[expressdata[1]==genenames[i],]                                        #有i
  geneexpression<-melt(geneexpression,id=names(geneexpression)[1])                                  #提取某一基因表达量做成矩阵
  geneexpression$level<-NA
  names(geneexpression)[1]<-"gene"
  geneexpression$patient<-substr(geneexpression$variable,1,12)                                      #选取病人资料                     
  patientexpression<-aggregate(value~patient+gene,geneexpression,mean)                              #相同病人组织表达量取均值
  patientexpression$level[patientexpression$value>median(patientexpression$value)]<-"High"          #大于中位数为高表达
  patientexpression$level[patientexpression$value<=median(patientexpression$value)]<-"Low"          #低于中位数者为低表达
  lastsurvival<-merge(survivalnew,patientexpression,id="patient")
  survivaldata<-lastsurvival                                                                        #merge生存资料和表达信息
  gene_surv <- Surv(survivaldata$time, survivaldata$status)~as.character(survivaldata$level)
  gene_survfit <- survfit(Surv(time, status) ~ level, data = survivaldata)                          #生存分析
  gene_survdiff <- tryCatch(survdiff(gene_surv), error = function(e) return(NA))
  p_Chisq <- pchisq(gene_survdiff$chisq, length(gene_survdiff$n) - 1, lower.tail = FALSE)
  gene_coxph <- coxph(gene_surv)
  summary_gene_coxph <- summary(gene_coxph)
  p_sctest = summary_gene_coxph$sctest[[3]]       #Score (logrank) test
  p_waldtest = summary_gene_coxph$waldtest[[3]]   #Wald test
  p_logtest = summary_gene_coxph$logtest[[3]]     #Likelihood ratio test
  HR = summary_gene_coxph$coefficients[1,2]
  result_x = data.frame(genename=genenames[i],p_Chisq = p_Chisq, p_sctest = p_sctest, p_waldtest = p_waldtest, p_logtest = p_logtest,
                        HR = HR, stringsAsFactors = F)                  #有i
  result_all<-rbind(result_all,result_x)
  write.csv(result_all,file=paste(analysislabel,"_pvalue.csv"),quote=F)
  if(p_logtest<0.05){
    pdf(file=paste(genename1[i],"survival curves.pdf"),width = 8,height = 8)
    BestCurve<-ggsurvplot(gene_survfit, data = survivaldata,
                          surv.median.line ="hv",
                          pval = T,
                          title=paste(genename1[i],"Overall survival curve"),
                          legend=c(0.9,0.9),
                          legend.title="",
                          legend.labs=c("High","Low"),     
                          pval.size = 9,
                          font.x = c(28, "plain", "black"), # feature of curve xlab title
                          font.y = c(28, "plain", "black"), # feature of curve ylab title
                          font.tickslab = c(25, "plain", "black"), palette = c("royalblue2","tomato3"),
                          font.title=c(25, "plain", "black")
                          
    )
    BestCurve$plot = BestCurve$plot+ theme(legend.text = element_text(size = 28, color = "black", face = "plain"),
                                           legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=18, face="plain"))
    print(BestCurve)
    dev.off()
  }
}
