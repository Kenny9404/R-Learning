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
# 读取数据
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
plot(data.qc) ## 第一列是所有样本的名称; 第二列是检测率和平均背景噪声。 第三列蓝色为比例因子，值（-3,3）），“圆圈”不能超过1.25，
# 否则数据质量不好，“三角形”不能超过3，否则数据质量不好，bioB表明芯片测试不符合标准。
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

4. 芯片探针与基因对应关系
# No.    gpl           organism                  bioc_package
# 1     GPL32       Mus musculus                        mgu74a
# 2     GPL33       Mus musculus                        mgu74b
# 3     GPL34       Mus musculus                        mgu74c
# 6     GPL74       Homo sapiens                        hcg110
# 7     GPL75       Mus musculus                     mu11ksuba
# 8     GPL76       Mus musculus                     mu11ksubb
# 9     GPL77       Mus musculus                     mu19ksuba
# 10    GPL78       Mus musculus                     mu19ksubb
# 11    GPL79       Mus musculus                     mu19ksubc
# 12    GPL80       Homo sapiens                        hu6800
# 13    GPL81       Mus musculus                      mgu74av2
# 14    GPL82       Mus musculus                      mgu74bv2
# 15    GPL83       Mus musculus                      mgu74cv2
# 16    GPL85  Rattus norvegicus                        rgu34a
# 17    GPL86  Rattus norvegicus                        rgu34b
# 18    GPL87  Rattus norvegicus                        rgu34c
# 19    GPL88  Rattus norvegicus                         rnu34
# 20    GPL89  Rattus norvegicus                         rtu34
# 22    GPL91       Homo sapiens                      hgu95av2
# 23    GPL92       Homo sapiens                        hgu95b
# 24    GPL93       Homo sapiens                        hgu95c
# 25    GPL94       Homo sapiens                        hgu95d
# 26    GPL95       Homo sapiens                        hgu95e
# 27    GPL96       Homo sapiens                       hgu133a
# 28    GPL97       Homo sapiens                       hgu133b
# 29    GPL98       Homo sapiens                     hu35ksuba
# 30    GPL99       Homo sapiens                     hu35ksubb
# 31   GPL100       Homo sapiens                     hu35ksubc
# 32   GPL101       Homo sapiens                     hu35ksubd
# 36   GPL201       Homo sapiens                       hgfocus
# 37   GPL339       Mus musculus                       moe430a
# 38   GPL340       Mus musculus                     mouse4302
# 39   GPL341  Rattus norvegicus                       rae230a
# 40   GPL342  Rattus norvegicus                       rae230b
# 41   GPL570       Homo sapiens                   hgu133plus2
# 42   GPL571       Homo sapiens                      hgu133a2
# 43   GPL886       Homo sapiens                     hgug4111a
# 44   GPL887       Homo sapiens                     hgug4110b
# 45  GPL1261       Mus musculus                    mouse430a2
# 49  GPL1352       Homo sapiens                       u133x3p
# 50  GPL1355  Rattus norvegicus                       rat2302
# 51  GPL1708       Homo sapiens                     hgug4112a
# 54  GPL2891       Homo sapiens                       h20kcod
# 55  GPL2898  Rattus norvegicus                     adme16cod
# 60  GPL3921       Homo sapiens                     hthgu133a
# 63  GPL4191       Homo sapiens                       h10kcod
# 64  GPL5689       Homo sapiens                     hgug4100a
# 65  GPL6097       Homo sapiens               illuminaHumanv1
# 66  GPL6102       Homo sapiens               illuminaHumanv2
# 67  GPL6244       Homo sapiens   hugene10sttranscriptcluster
# 68  GPL6947       Homo sapiens               illuminaHumanv3
# 69  GPL8300       Homo sapiens                      hgu95av2
# 70  GPL8490       Homo sapiens   IlluminaHumanMethylation27k
# 71 GPL10558       Homo sapiens               illuminaHumanv4
# 72 GPL11532       Homo sapiens   hugene11sttranscriptcluster
# 73 GPL13497       Homo sapiens         HsAgilentDesign026652
# 74 GPL13534       Homo sapiens  IlluminaHumanMethylation450k
# 75 GPL13667       Homo sapiens                        hgu219
# 76 GPL15380       Homo sapiens      GGHumanMethCancerPanelv1
# 77 GPL15396       Homo sapiens                     hthgu133b
# 78 GPL17897       Homo sapiens                     hthgu133a

5. 雷达图的制作
library(devtools)
install_github("ricardo-bion/ggradar")
library("ggplot2")
library("ggradar")
mydata<-matrix(runif(40,0,1),5,8)
rownames(mydata) <- LETTERS[1:5]
colnames(mydata) <- c("Apple","Google","Facebook","Amozon","Tencent","Alibaba","Baidu","Twitter")
mynewdata<-data.frame(mydata)
Name<-c("USA","CHN","UK","RUS","JP")
mynewdata<-data.frame(Name,mynewdata)
ggradar(mynewdata)

6. Rmarkdown中文写作解决方案
header-includes: 
  - \usepackage{xeCJK} 