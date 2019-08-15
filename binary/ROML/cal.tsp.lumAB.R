### Calculate KTSP scores
source('~/ktsp/script/function/cal.TSPscore.R')

# Balanced LumA vs LumB

data.dir <- '~/ktsp/data/binary/realdata/balanced'
load(file.path(data.dir, "TCGA.balance.LumAB.Rdata"))
load(file.path(data.dir, "MB.balance.LumAB.Rdata"))

lb.brca <- factor(brca.clin.sub$final_assign)
lb.mb <- factor(mb.clin.sub$Pam50Subtype)

# TCGA -> MB
load(file.path(data.dir, "baseline/tcga.DEG.Rdata"))
brca.train <- brca.expr.sub[match(rownames(tt.brca), rownames(brca.expr.sub)),]
mb.test <- mb.expr.sub[match(rownames(tt.brca), rownames(mb.expr.sub)),]

tcga.tsp <- cal.TSPscore(brca.train, lb.brca)
save(tcga.tsp, file = file.path(data.dir, 'ROML/tcga.TSPscore.Rdata'))

# MB -> TCGA
load(file.path(data.dir, "baseline/mb.DEG.Rdata"))
mb.train <- mb.expr.sub[match(rownames(tt.mb), rownames(mb.expr.sub)),]
brca.test <- brca.expr.sub[match(rownames(tt.mb), rownames(brca.expr.sub)),]

# Workflow
mb.tsp <- cal.TSPscore(mb.train, lb.mb)
save(mb.tsp, file = file.path(data.dir, 'ROML/mb.TSPscore.Rdata'))



