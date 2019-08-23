### Simulation - add shift on real data
source('~/ktsp/script/function/lowMeans.filter.R')
source('~/ktsp/script/function/cross.platform.norm.R')
source('~/ktsp/script/function/binary.test.R')
source('~/ktsp/script/function/pair.t.filter.R')


bi.simu <- function(data.train, lb.train, lb.test, base.fea.no, toptable, kpair, var) {
  final.res <- data.frame() -> final.roml
  for (v in var) {
    print(paste0('variance=', v))
    res <- c()
    roml <- c()
    
    for (i in 1:5) {
      # Generate new data
      shift <- rnorm(nrow(data.train), mean=0, sd=v)
      new.data <- data.train + shift
      nm.data.test <- cross.platform.norm(data.train=data.train, data.test=new.data, method='QN')
      
      for (j in 1:3) {
        # train models
        baseline.model <- train.bi.model(data.train, lb.train, rownames(toptable[1:base.fea.no, ]))
        roml.model <- train.bi.model(data.train, lb.train, kpair, baseline=FALSE)
        
        # baseline
        print(paste0('==Top', base.fea.no,' DEGs...'))
        tmp <- bi.test(nm.data.test, lb.test, baseline.model, rownames(toptable[1:base.fea.no, ]))
        res <- c(res, tmp[["ACC"]])
        
        # ROML
        tmp2 <- bi.test(nm.data.test, lb.test, roml.model, kpair, baseline=FALSE)
        roml <- c(roml, tmp2[["ACC"]])
      }
    }
    final.res <- rbind(final.res, c(mean(res), sd(res), v))
    final.roml <- rbind(final.roml, c(mean(roml), sd(roml), v))
  }
  
  # test on itself
  res0 <- c() -> roml0
  for (i in 1:5) {
    baseline.model <- train.bi.model(data.train, lb.train, rownames(toptable[1:base.fea.no,]))
    tmp0 <- bi.test(data.train, lb.train, baseline.model, rownames(toptable[1:base.fea.no,]))
    res0 <- c(res0, tmp0[["ACC"]])
    
    roml0.model <- train.bi.model(data.train, lb.train, kpair, baseline=FALSE)
    tmp00 <- bi.test(data.train, lb.train, roml0.model, kpair, baseline=FALSE)
    roml0 <- c(roml0, tmp00[["ACC"]])
  }
  
  final.res <- rbind.data.frame(c(mean(res0), sd(res0), 0), final.res)
  final.roml <- rbind.data.frame(c(mean(roml0), sd(roml0), 0), final.roml)
  colnames(final.roml) <- c('meanACC', 'varACC', 'var') -> colnames(final.res)
  
  final.res$method=rep('baseline', nrow(final.res))
  final.roml$method=rep('ROML', nrow(final.roml))
  
  final.df <- rbind(final.res, final.roml)
  final.df$var <- factor(final.df$var, levels=as.character(unique(final.df$var)))
  final.df
}

draw.res <- function(plot.df) {
  ggplot(plot.df, aes(x=var, y=meanACC, group=method, color=method)) +   
    geom_line(size=1.1) +
    geom_point(size=1.1)+
    geom_errorbar(aes(ymin=meanACC-varACC, ymax=meanACC+varACC), width=.2,
                  position=position_dodge(0.2)) +
    labs(y='ACC', x='Standard Deviation') + ggtitle('RF vs tspRF') +
    scale_y_continuous(breaks = seq(0.5,1, by=0.1), limits = c(0.5, 1)) +
    scale_color_manual(values=c('lightskyblue1','rosybrown1', 'skyblue2', 'lightpink2')) +
    theme_bw()
}


var <- c(0.4, 0.8, 1.2, 1.6, 2)


## Breast cancer 
load("~/ktsp/data/binary/realdata/balanced/TCGA.balance.LumAB.Rdata")
load("~/ktsp/data/binary/realdata/balanced/MB.balance.LumAB.Rdata")


# 1. MB -> MB + shift
data.train <- mb.expr.sub
lb.train <- lb.test <- factor(mb.clin.sub$Pam50Subtype)

load("~/ktsp/data/binary/realdata/balanced/baseline/mb.DEG.Rdata")
tt.train <- tt.mb
data.train <- data.train[match(rownames(tt.mb), rownames(data.train)),] #lowmean filter

Ks <- c(500)
load('~/ktsp/data/binary/realdata/balanced/ROML/mb.TSPscore.Rdata')
sim.kpairs <- mb.tsp[c(1:Ks),]
top.pair=mb.tsp[1:1000,]
tmp.t <- as.data.frame(t(apply(top.pair[,-1], 1, modify_t, data.train=data.train, label=lb.train)))
kpair=cbind(top.pair, tmp.t)
kpair=kpair[order(kpair$min.abs.t, decreasing=T),]
t.kpair <- kpair[1:500,]



plot.df <- bi.simu(data.train, lb.train, lb.test, 500, tt.train, t.kpair, var)

draw.res(plot.df)



# 2. BRCA -> BRCA + shift 
data.train <- brca.expr.sub
lb.train <- lb.test <- factor(brca.clin.sub$final_assign)

load("~/ktsp/data/binary/realdata/balanced/baseline/tcga.DEG.Rdata")
tt.train <- tt.brca
data.train=data.train[match(rownames(tt.brca), rownames(data.train)),]

Ks <- c(500)
load('~/ktsp/data/binary/realdata/balanced/ROML/tcga.TSPscore.Rdata')
sim.kpairs <- tcga.tsp[c(1:Ks),]


plot.df <- bi.simu(data.train, lb.train, lb.test, 500, tt.train, sim.kpairs, var)

draw.res(plot.df)


##### Colorectal cancer

load('~/ktsp/CRC/TCGA.COAD.Rdata')
load('~/ktsp/CRC/kfsyscc.CRC.Rdata')

coad.clin <- coad.clin[coad.clin$cms_label %in% c('CMS2', 'CMS4'), ]
coad.expr <- coad.expr[,match(rownames(coad.clin), colnames(coad.expr))]


rownames(kf.clin) <- unlist(kf.clin$id)
kf.clin <- kf.clin[kf.clin$cms_label %in% c('CMS2', 'CMS4'), ]
kf.expr <- kf.expr[,match(rownames(kf.clin), colnames(kf.expr))]

lb.coad <- factor(coad.clin$cms_label)
lb.kf <- factor(kf.clin$cms_label)

## COAD -> KFSYSCC
load("~/ktsp/CRC/binary/coad.DEG.Rdata")
data.train <- coad.expr[match(rownames(tt.coad), rownames(coad.expr)),]
lb.train=lb.test=lb.coad
tt.train=tt.coad

load("~/ktsp/CRC/multiclass/score/tcga/CMS2_CMS4.score.Rdata")
sim.kpair=tsp[1:500,]

plot.df <- bi.simu(data.train, lb.train, lb.test, 500, tt.train, sim.kpair, var)
draw.res(plot.df)



## KFSYSCC -> COAD
load("~/ktsp/CRC/multiclass/score/kfsyscc/CMS2_CMS4.score.Rdata")
load("~/ktsp/CRC/binary/kfsyscc.DEG.Rdata")

data.train <- kf.expr[match(rownames(tt.kf), rownames(kf.expr)),]
lb.train=lb.test=lb.kf
tt.train=tt.kf
sim.kpair=tsp[1:500,]

plot.df <- bi.simu(data.train, lb.train, lb.test, 500, tt.train, sim.kpair, var)
draw.res(plot.df)
