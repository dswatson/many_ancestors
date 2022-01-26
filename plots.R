# Set working directory
setwd('~/Documents/UCL/many_ancestors')

# Load libraries
library(data.table)
library(scales)
library(ggsci)
library(foreach)
library(tidyverse)

################################################################################

### Bivariate ###

# Import bivariate results
df <- fread('./results/lin_biv_benchmark.csv')

# Polish
sims <- expand.grid(n = c(1000, 2000, 4000), g = c('xy', 'ci', 'na')) %>% 
  mutate(sp = 0.5, d_z = 100, rho = 0.25, r2 = 2/3, lin_pr = 1, 
         s_id = row_number()) %>%
  as.data.table(.)
df <- merge(df, sims, by = 's_id') 
orcl <- expand.grid(s_id = 1:9, idx = 1:100) %>%
  mutate(method = 'oracle') %>%
  inner_join(sims, by = 's_id') %>%
  mutate(h = g) %>%
  select(colnames(df))
df <- rbind(df, orcl)
df[, na := sum(h == 'na'), by = .(s_id, method)]
df[, ci := sum(h == 'ci'), by = .(s_id, method)]
df[, xy := sum(h == 'xy'), by = .(s_id, method)]
df[, yx := sum(h == 'yx'), by = .(s_id, method)]
df <- df %>%
  pivot_longer(cols = na:yx, names_to = 'g_hat', values_to = 'prop') %>%
  select(s_id, n, g, method, g_hat, prop) %>%
  unique(.) %>%
  mutate(n = factor(paste('n =', n), 
                    levels = c('n = 1000', 'n = 2000', 'n = 4000')),
         g_hat = factor(g_hat, levels = c('xy', 'yx', 'ci', 'na'))) %>%
  as.data.table(.)
df[g == 'xy', g := '(a)']
df[g == 'ci', g := '(b)']
df[g == 'na', g := '(c)']
df[method == 'cbr', method := 'CBR']
df[method == 'constr', method := 'constraint']
df[, method := factor(method, 
                      levels = c('constraint', 'score', 'CBR', 'oracle'))]
df[, g_hat := recode_factor(g_hat, `xy` = 'X %->% Y', `yx` = 'Y %->% X',
                            `ci` = 'X %~% Y', `na` = 'NA')]

# Plot
p <- ggplot(filter(df, method != 'oracle'), 
            aes(method, prop, fill = g_hat)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_npg(labels = parse_format()) + 
  labs(x = 'Method', y = 'Percentage', title = 'Linear SCM',
       fill = 'Estimated\nStructure') + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(g ~ n) 
ggsave('./plots/biv.pdf', width = 3, height = 6)
  
p + theme(legend.position = 'none')

################################################################################

#### Multivariate ###

# Import data
res1 <- na.omit(readRDS('./results/multivariate.rds'))
res2 <- na.omit(readRDS('./results/rfci_res.rds'))

# Loop through, export summary stats
adj_smry <- function(df, sim_id, m) {
  fn <- function(i, m) {
    tmp <- df %>% filter(s_id == sim_id, idx == i)
    if (m == 'true') {
      out <- tmp$amat[[1]] %>% keep(lower.tri(.))
    } else if (m == 'cbr') {
      out <- tmp$amat_cbr[[1]] %>% keep(lower.tri(.))
    } else if (m == 'ges') {
      out <- tmp$amat_ges[[1]] %>% keep(lower.tri(.))
    } else if (m == 'rfci') {
      out <- tmp$amat_rfci[[1]] %>% keep(lower.tri(.))
    }
    return(out)
  }
  max_idx <- max(df$idx)
  y <- foreach(ii = seq_len(max_idx), .combine = c) %do% fn(ii, 'true')
  y_hat <- foreach(ii = seq_len(max_idx), .combine = c) %do% fn(ii, m)
  tmp <- data.table(y, y_hat)
  pjr <- tmp[, sum(!is.na(y_hat)) / .N]
  acc <- tmp[!is.na(y_hat), sum(y == y_hat) / .N]
  tpr <- tmp[!is.na(y_hat) & y == 1, sum(y_hat) / .N]
  tnr <- tmp[!is.na(y_hat) & y == 0, sum(1 - y_hat) / .N]
  out <- data.table(method = m, s_id = sim_id, pjr, acc, tpr, tnr)
  return(out)
}
df <- foreach(ss = 1:5, .combine = rbind) %:%
  foreach(mm = c('cbr', 'ges'), .combine = rbind) %do% 
  adj_smry(res1, ss, mm)
tmp <- adj_smry(res2, 1, 'rfci')
df <- rbind(df, tmp)

# Polish
sims <- data.table(s_id = 1:5, n = c(500, 1000, 2000, 4000, 8000))
df <- merge(df, sims, by = 's_id')
df[method != 'rfci', B := 15 * 20 * pjr]  # Replicates for SE calculation
df[method == 'rfci', B := 15 * pjr]
df[, se := sqrt(acc * (1 - acc) / B)]
df[method == 'cbr', method := 'CBR']
df[method == 'ges', method := 'GES']
df[method == 'rfci', method := 'RFCI']

# Plot
p <- ggplot(df, aes(n, acc, group = method, color = method)) + 
  geom_point(aes(shape = method), size = 3) + 
  geom_path(size = 0.5) + 
  geom_errorbar(aes(ymin = acc - se, ymax = acc + se), width = 0.05) + 
  scale_color_npg() + 
  scale_x_log10(breaks = c(500, 1000, 2000, 4000, 8000)) + 
  ylim(0, 1) + 
  labs(x = 'Sample Size', y = 'Accuracy') + 
  guides(color = guide_legend(title = 'Method'),
         shape = guide_legend(title = 'Method')) + 
  theme_bw() +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        legend.position = c(0.85, 0.2), 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.box.background = element_rect(colour = 'black'))
ggsave('./plots/multiv.pdf', width = 6, height = 6)




