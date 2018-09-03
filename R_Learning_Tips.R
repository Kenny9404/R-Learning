1.使用Git指南
ssh -T git@github.com
cd /mnt/d/Scientific\ Research/Research\ Project/Bioinformatics/myGitHub/
git add .
git commit -m "first commit" 
# git remote add origin git@github.com:Kenny9404/R-Learning.git
# git remote rm origin
# git remote add origin git@github.com:Kenny9404/R-Learning.git
git push -u origin master
# git pull origin master
# git push -u origin master

2. 获取芯片Series Matrix数据并做ID转换
# source("https://bioconductor.org/biocLite.R")
# biocLite("GEOquery")
library(GEOquery)
GSE41038 <- getGEO(GEO = "GSE41038", destdir = "")
eset <- exprs(GSE41038[[1]])
eset <- as.data.frame(eset)
eset$probe_id=rownames(eset)
# GPL6883 Illumina HumanRef-8 v3.0 expression beadchip
# A search for GPL6883 on GEO identified the title for this platform as Illumina HumanRef-8 v3.0 expression beadchip. A subsequent search on Bioconductor for illuminahuman identified the appropriate package as illuminaHumanv3.
library('illuminaHumanv3.db')
probe2symbol=toTable(illuminaHumanv3SYMBOL)
dat <- merge(eset,probe2symbol,by='probe_id')
symbol.data <- t(sapply(split(dat,dat$symbol),function(x) colMeans(x[,2:(ncol(x)-1)])))  # 取重复探针平均值

3. 获取芯片原始CEL数据并做ID转换
library(GEOquery)
getGEOSuppFiles("GSE55457", baseDir = "data/") ## 解压数据
library(affy)
library(affyPLM)
library(simpleaffy)
library(RColorBrewer)
cel.files <- list.celfiles("data/GSE55457/", full.name = TRUE)
data <- ReadAffy(filenames = cel.files)
sampleNames(data)
# 质量控制
data.qc <- qc(data) 
plot(data.qc)
第一列是所有样本的名称; 第二列是检测率和平均背景噪声。 第三列蓝色为比例因子，
值（-3,3）），“圆圈”不能超过1.25，否则数据质量不好，“三角形”不能超过3，否则数据质量不好，
bioB表明芯片测试不符合标准
#背景矫正、标准化和汇总
eset.rma <- rma(data)
eset.mas5 <- mas5(data)
eset.gcrma <- gcrma(data)
# 比较算法优劣性
colors<-brewer.pal(12,"Set3")
par(mfrow = c(2,2))
hist(data, main="orignal", col=colors)
hist(eset.mas5, main="MAS 5.0", xlim=c(-150,2^10), col=colors)
hist(eset.rma , main="RMA", col=colors)
hist(eset.gcrma, main="gcRMA", col=colors)
# 提取表达矩阵
eset<-exprs(eset.gcrma)
eset <- as.data.frame(eset)
# 探针ID转换
eset$probe_id=rownames(eset)
# GPL96:Affymetrix Human Genome U133A Array
library('hgu133a.db')
probe2symbol=toTable(hgu133aSYMBOL)
dat <- merge(eset,probe2symbol,by='probe_id')
symbol.data <- t(sapply(split(dat,dat$symbol),function(x) colMeans(x[,2:(ncol(x)-1)]))) ## 取重复探针平均值
