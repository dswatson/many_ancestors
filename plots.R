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
  mutate(g_hat = g) %>%
  select(colnames(df))
df <- rbind(df, orcl)
df[, na := sum(g_hat == 'na'), by = .(s_id, method)]
df[, ci := sum(g_hat == 'ci'), by = .(s_id, method)]
df[, xy := sum(g_hat == 'xy'), by = .(s_id, method)]
df[, yx := sum(g_hat == 'yx'), by = .(s_id, method)]
df <- df %>%
  select(-g_hat) %>%
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
df[method == 'cbr', method := 'CBL']
df[method == 'constr', method := 'constraint']
df[, method := factor(method, 
                      levels = c('constraint', 'score', 'CBL', 'oracle'))]
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
ggsave('./plots/biv.pdf', width = 4, height = 6)
  
p + theme(legend.position = 'none')

################################################################################

#### Multivariate ###

# Import data
res1 <- na.omit(readRDS('./results/multivariate.rds'))
res1[, d_z := 100]
res2 <- na.omit(readRDS('./results/multivariate2.rds'))
res2[, d_z := 50]
res1 <- rbind(res1, res2)

res1_rfci <- na.omit(readRDS('./results/rfci_res.rds'))
res1_rfci[, d_z := 100]
res2_rfci <- na.omit(readRDS('./results/rfci_res2.rds'))
res2_rfci[, d_z := 50]
res2 <- rbind(res1_rfci, res2_rfci)

# Loop through, export summary stats
adj_smry <- function(df, sim_id, m) {
  fn <- function(id, dz, m) {
    tmp <- df %>% filter(s_id == sim_id, idx == id, d_z == dz)
    tru <- as.numeric(tmp$amat[[1]])
    keep <- !is.na(tru)
    if (m == 'true') {
      out <- tru[keep]
    } else if (m == 'cbr') {
      out <- tmp$amat_cbr[[1]]
      d_x <- ncol(out)
      for (i in seq_len(d_x)) {
        for (j in seq_len(d_x)[-i]) {
          if (!is.na(out[i, j])) {
            if (out[i, j] > 0) {
              out[j, i] <- 0
            }
          } 
        }
      }
      out <- as.numeric(out)[keep]
      # Judgments on partial orders are reflected in the mirror image
      out[out %in% c(-0.5, 0.5)] <- NA_real_
    } else if (m == 'ges') {
      out <- as.numeric(tmp$amat_ges[[1]])[keep]
    } else if (m == 'rfci') {
      out <- as.numeric(tmp$amat_rfci[[1]])[keep]
    }
    return(out)
  }
  max_idx <- max(df$idx)
  y <- foreach(dz = c(50, 100), .combine = c) %:%
    foreach(ii = seq_len(max_idx), .combine = c) %do% fn(ii, dz, 'true')
  y_hat <- foreach(dz = c(50, 100), .combine = c) %:%
    foreach(ii = seq_len(max_idx), .combine = c) %do% fn(ii, dz, m)
  tmp <- data.table(y, y_hat)
  n <- nrow(tmp)
  tmp[, d_z := rep(c(50, 100), each = n/2)]
  pjr <- tmp[, sum(!is.na(y_hat)) / .N, by = d_z]$V1
  acc <- tmp[!is.na(y_hat), sum(y == y_hat) / .N, by = d_z]$V1
  tpr <- tmp[!is.na(y_hat) & y == 1, sum(y_hat) / .N, by = d_z]$V1
  tnr <- tmp[!is.na(y_hat) & y == 0, sum(1 - y_hat) / .N, by = d_z]$V1
  out <- data.table(
    method = m, s_id = sim_id, d_z = c(50, 100), pjr, acc, tpr, tnr
  )
  return(out)
}
df1 <- foreach(ss = 1:5, .combine = rbind) %:%
  foreach(mm = c('cbr', 'ges'), .combine = rbind) %do% 
  adj_smry(res1, ss, mm)

# Need a different function for RFCI
adj_smry <- function(df, sim_id, m) {
  fn <- function(id, dz, m) {
    tmp <- df %>% filter(s_id == sim_id, idx == id, d_z == dz)
    tru <- as.numeric(tmp$amat[[1]])
    keep <- !is.na(tru)
    if (m == 'true') {
      out <- tru[keep]
    } else {
      out <- as.numeric(tmp$amat_rfci[[1]])[keep]
    }
    return(out)
  }
  if (sim_id == 1) {
    y <- foreach(d = c(50, 100), .combine = c) %do% fn(1, d, 'true')
    y_hat <- foreach(d = c(50, 100), .combine = c) %do% fn(1, d, 'rfci')
    tmp <- data.table(y, y_hat)
    n <- nrow(tmp)/2
    tmp[, d_z := c(rep(50, n), rep(100, n))]
    dz <- c(50, 100)
  } else if (sim_id == 2) {
    y <- fn(1, 50, 'true')
    y_hat <- fn(1, 50, 'rfci')
    tmp <- data.table(y, y_hat)
    tmp[, d_z := 50]
    dz <- 50
  }
  pjr <- tmp[, sum(!is.na(y_hat)) / .N, by = d_z]$V1
  acc <- tmp[!is.na(y_hat), sum(y == y_hat) / .N, by = d_z]$V1
  tpr <- tmp[!is.na(y_hat) & y == 1, sum(y_hat) / .N, by = d_z]$V1
  tnr <- tmp[!is.na(y_hat) & y == 0, sum(1 - y_hat) / .N, by = d_z]$V1
  # Export
  out <- data.table(
    method = m, s_id = sim_id, d_z = dz, pjr, acc, tpr, tnr
  )
  return(out)
}
df2 <- foreach(ss = 1:2, .combine = rbind) %do%
  adj_smry(res2, ss, 'rfci')
df <- rbind(df1, df2)

# Polish
sims <- data.table(s_id = 1:5, n = c(500, 1000, 2000, 4000, 8000))
df <- merge(df, sims, by = 's_id')
df[method != 'rfci', B := 30 * 20 * pjr]  # Replicates for SE calculation
df[method == 'rfci', B := 30 * pjr]
df[, se := sqrt(acc * (1 - acc) / B)]
df[method == 'cbr', method := 'CBR']
df[method == 'ges', method := 'GES']
df[method == 'rfci', method := 'RFCI']

# Plot
p <- ggplot(df, aes(n, acc, group = method, color = method)) + 
  geom_point(aes(shape = method), size = 3) + 
  geom_path(size = 0.5) + 
  geom_errorbar(aes(ymin = acc - se, ymax = acc + se), width = 0.05) + 
  scale_shape_manual(values = c(1, 2, 0))+
  scale_color_npg() + 
  scale_x_log10(breaks = c(500, 1000, 2000, 4000, 8000)) + 
  ylim(0.4, 1) + 
  labs(x = 'Sample Size', y = 'Accuracy') + 
  guides(color = guide_legend(title = 'Method'),
         shape = guide_legend(title = 'Method')) + 
  theme_bw() +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        legend.position = c(0.9, 0.2), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.box.background = element_rect(colour = 'black')) + 
  facet_wrap(~d_z, labeller = label_bquote(italic(d[Z])==.(d_z)))
  
ggsave('./plots/multiv.pdf', width = 6, height = 4)




