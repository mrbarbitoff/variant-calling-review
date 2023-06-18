setwd("/media/barbitoff/DATA/Working issues/WES/BRK/NGS_in_rare_disease_Review")
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggsci)

### Figure 2
rma_data = read.table('./rma_counts.tsv', sep='\t', header = T)
head(rma_data)
rma_cnt <- ggplot(rma_data, aes(x=Assembly, y=Count/1000, fill=Concordance)) + 
  geom_bar(col='black', stat='identity') + theme_bw() + 
  scale_fill_simpsons() + ylab('RMA count (x1000)') +
  facet_wrap(~Type, nrow=2, scales='free') +
  theme(legend.position='top', panel.grid = element_blank()) +
  coord_flip()
rma_cnt

mf_data = read.table('MF_assemblies.tsv', sep='\t', header=F)
colnames(mf_data) = c('Sample', 'Assembly', 'MF > 0.4', 'MF > 0.95')
assemblies = c('b37', 'b37 (primary)', 'hg38', 'hg38 (primary)', 'T2T')
names(assemblies) = unique(mf_data$Assembly)
mf_data$Assembly = sapply(mf_data$Assembly, function(x) assemblies[x])
head(mf_data)
mfd = melt(mf_data, measure.vars=c('MF > 0.4', 'MF > 0.95'))
mfplot <- ggplot(mfd, aes(x=Assembly, y=value/1000, fill=Assembly)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(size=2, pch=21, col='black') +
  facet_wrap(~variable, nrow=1) + theme_bw() + 
  ylab('Bases (kbp)') + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid=element_blank()) +
  scale_fill_locuszoom()
mfplot

plot_grid(rma_cnt, mfplot, nrow=1, labels=c('a', 'b'), rel_widths = c(0.7, 1))


### Figure 3

concord = read.table('./giab_allregions.tsv', sep='\t', header=T)
head(concord)

concord$called = rowSums(concord[, c(6:9, 11:14)])
concord = concord[!(concord$HC_2DCNN == 1 & concord$called == 0), ]

caller_num <- ggplot(concord, aes(x=sample, y=called, fill=TRUTH)) + geom_violin(scale='width') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(vars(TRUTH), vars(REGIONS), scales="free_x", space='free_x') + theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_y_continuous(limits=c(0, 8)) +
  scale_fill_brewer(palette = "Accent") +
  ylab('Number of callers')
caller_num

aliases = data.frame(src = colnames(concord)[6:14],
        alias = c('CL', 'DV', 'FB', 'G1', 'G2', 'GH', 'OF', 'OS', 'ST'))

uq_calls <- melt(concord[concord$called == 1, c(1, 4, 6:9, 11:14)], 
                 id.vars=c('sample', 'REGIONS'))
uq_call_df <- aggregate(value~variable+sample+REGIONS, uq_calls, sum)
uq_call_df$variable = sapply(uq_call_df$variable, 
                function(x) aliases[aliases$src == as.character(x), 'alias'])
uq_call_df$total = aggregate(value~variable+sample+REGIONS, 
    melt(concord[, c(1, 4, 6:9, 11:14)], id.vars=c('sample', 'REGIONS')), sum)$value
uq_call_df$pct_unique = uq_call_df$value / uq_call_df$total
head(uq_call_df)

uc <- ggplot(uq_call_df, aes(x=variable, y=pct_unique, fill=REGIONS)) + 
  geom_boxplot(position=position_dodge(1)) +
  theme_bw() + scale_fill_brewer(palette="Accent") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid=element_blank()) +
  ylab('Unique call rate') + guides(fill=F) +
  scale_y_log10() + facet_wrap(~variable, nrow=1, scales='free_x')
uc

uq_noncalls <- melt(concord[concord$called == 7, c(1, 4, 6:9, 11:14)], 
                    id.vars=c('sample', 'REGIONS'))
uq_noncall_df <- aggregate(value~variable+sample+REGIONS, 
                           uq_noncalls, function(x) length(x) - sum(x))
uq_noncall_df$variable = sapply(uq_noncall_df$variable,
                function(x) aliases[aliases$src == as.character(x), 'alias'])
uq_noncall_df$total = aggregate(value~variable+sample+REGIONS, 
          melt(concord[, c(1, 4, 6:9, 11:14)], id.vars=c('sample', 'REGIONS')), sum)$value
uq_noncall_df$pct_unique = uq_noncall_df$value / uq_noncall_df$total
head(uq_call_df)

unc <- ggplot(uq_noncall_df, aes(x=variable, y=pct_unique, fill=REGIONS)) + 
  geom_boxplot(position=position_dodge(1)) +
  theme_bw() + scale_fill_brewer(palette="Accent") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid=element_blank()) +
  ylab('Unique non-call rate') + guides(fill=F) +
  scale_y_log10() + facet_wrap(~variable, nrow=1, scales='free_x')
unc

plot_grid(uc, unc, nrow=2)



### Figure 4
clinvar_stats =  read.table('clinvar_counts.tsv', sep='\t',
                            header=T)
head(clinvar_stats)
clinvar_stats$Version = as.Date(clinvar_stats$Version, "%Y/%m/%d")
clinvar_stats$pct_coding = clinvar_stats$Coding/clinvar_stats$Total
cds_pct <- ggplot(clinvar_stats, aes(x=Version, y=pct_coding, col=Class)) + 
  geom_point(size=2) + geom_line() + theme_bw() + 
  scale_color_brewer(palette='Set1') + 
  scale_y_continuous(limits=c(0, 1)) +
  theme(legend.position = 'top', panel.grid=element_blank()) +
  xlab('ClinVar version (release date)') +
  ylab('% of coding variants')
cds_pct

ggplot(clinvar_stats, aes(x=Version, y=Total, fill=Class)) + 
  geom_area(stat='identity', col='black') + theme_bw() +
  facet_wrap(~Class, scales='free', nrow=2) +
  scale_fill_simpsons()

unmlt <- clinvar_stats[clinvar_stats$Class == 'All', ]
unmlt$pct_patho = clinvar_stats[clinvar_stats$Class == 'Pathogenic', 'Total'] / 
  unmlt$Total
patho_pct <- ggplot(unmlt, aes(x=Version, y=pct_patho)) + geom_point(col='red') +
  geom_line(col='red') + theme_bw() + scale_y_continuous(limits=c(0, 0.2))
patho_pct

plot_grid(cds_pct, patho_pct, nrow=2, rel_heights = c(1, 0.7))


### Figure X - variant counts

counts = read.table('variant_counts.tsv', sep='\t', header=T)
cnts = melt(counts, id.vars='regions')
cnts$expected = ifelse(grepl('_exp', cnts$variable), 'expected', 'observed')
cnts$variable = as.character(cnts$variable)
cnts$variable[grepl('exome', cnts$variable)] = 'exome_variants'
ggplot(cnts, aes(x=regions, y=value, fill=expected)) + 
  geom_bar(stat='identity', position='dodge') + theme_bw() + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~variable, nrow=1, scales='free_y')
