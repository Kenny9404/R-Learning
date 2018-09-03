1.使用Git指南
ssh -T git@github.com
cd /mnt/d/Scientific\ Research/Research\ Project/Bioinformatics/myGitHub/


git add .
git commit -m ‘first commit’ 
git remote add origin git@github.com:Kenny9404/R-Learning.git 
git push -u origin master

2. 获取探针Series Matrix数据并做ID转换
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
