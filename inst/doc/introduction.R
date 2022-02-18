## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>', fig.retina = 2, message = FALSE)

## -----------------------------------------------------------------------------
library('annotate')
library('data.table')
library('foreach')
library('ggplot2')
library('limma')
library('limorhyde')
library('org.Mm.eg.db')
library('qs')

## -----------------------------------------------------------------------------
period = 24
qvalRhyCutoff = 0.15
qvalDrCutoff = 0.1

## -----------------------------------------------------------------------------
y = qread(system.file('extdata', 'GSE34018_expression_data.qs', package = 'limorhyde'))
y[1:5, 1:5]

metadata = qread(system.file('extdata', 'GSE34018_metadata.qs', package = 'limorhyde'))
metadata

## -----------------------------------------------------------------------------
metadata = cbind(metadata, limorhyde(metadata$time, 'time_'))

## -----------------------------------------------------------------------------
rhyLimma = foreach(condNow = unique(metadata$cond), .combine = rbind) %do% {
  design = model.matrix(~ time_cos + time_sin, data = metadata[cond == condNow])
  fit = lmFit(y[, metadata$cond == condNow], design)
  fit = eBayes(fit, trend = TRUE)
  rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
  setnames(rhyNow, 'rn', 'gene_id')
  rhyNow[, cond := condNow]}

rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(rhyLimmaSummary, 'adj.P.Val')

rhyLimmaSummary[1:5, ]

## -----------------------------------------------------------------------------
design = model.matrix(~ cond * (time_cos + time_sin), data = metadata)

fit = lmFit(y, design)
fit = eBayes(fit, trend = TRUE)
drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
setnames(drLimma, 'rn', 'gene_id')

drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= qvalRhyCutoff]$gene_id]
drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(drLimma, 'adj.P.Val')

drLimma[1:5, ]

## -----------------------------------------------------------------------------
design = model.matrix(~ cond + time_cos + time_sin, data = metadata)

fit = lmFit(y, design)
fit = eBayes(fit, trend = TRUE)
deLimma = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
setnames(deLimma, 'rn', 'gene_id')

deLimma = deLimma[!(gene_id %in% drLimma[adj.P.Val <= qvalDrCutoff]$gene_id)]
deLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(deLimma, 'adj.P.Val')

deLimma[1:5, ]

## ---- fig.width = 4, fig.height = 3-------------------------------------------
ggplot(deLimma) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), size = 0.2, alpha = 0.5) +
  labs(x = expression(log[2]*' fold-change'), y = expression(-log[10]*' '*q[DE]))

## ---- fig.width = 5, fig.height = 4-------------------------------------------
geneIdsNow = c(rhyLimmaSummary$gene_id[1L], drLimma$gene_id[1L], deLimma$gene_id[1L])
geneSymbolsNow = unname(getSYMBOL(geneIdsNow, 'org.Mm.eg.db'))

df = data.table(t(y[geneIdsNow, ]))
setnames(df, geneSymbolsNow)
df[, sample_id := colnames(y[geneIdsNow, ])]

df = merge(df, metadata[, .(sample_id, cond, time)], by = 'sample_id')
df = melt(df, measure.vars = geneSymbolsNow, variable.name = 'symbol',
          value.name = 'expr')
df[, symbol := factor(symbol, geneSymbolsNow)]

ggplot(df) +
  facet_grid(symbol ~ cond, scales = 'free_y') +
  geom_point(aes(x = time, y = expr, shape = cond, color = symbol), size = 2) +
  labs(x = 'Zeitgeber time (h)', y = 'Expression (norm.)') +
  scale_shape(solid = FALSE, guide = 'none') +
  scale_color_brewer(type = 'qual', palette = 'Dark2', guide = 'none') +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 4))

