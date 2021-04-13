library(data.table)
library(ggplot2)
library(ggpubr)

diverge_plot <- function(dd, group, threshold = 4) {
  
  dsum <- copy(dd)
  
  dsum[, prop := N/sum(N), keyby = rx]
  dsum[, cumprop := cumsum(prop), keyby = rx]
  dsum[, cutpoint := qlogis(cumprop)]
  dsum[, cumodds := (1-cumprop)/cumprop]
  dsum[, linept := rx*(2 - 0.3) + (1 - rx)*(1 + 0.3)]
  dsum[, rx := factor(rx, labels = c("control", "treatment"))]
  
  dnames <- dsum[, sum(N), keyby = rx]
  dnames[, axislabel := paste0(rx, "\n(n = ", V1, ")")]
  dnames[, rxlabel := paste0(rx, " (n = ", V1, ")")]
  
  dsum[, axislabel := factor(rx, labels = dnames[, axislabel])]
  dsum[, rxlabel := factor(rx, labels = dnames[, rxlabel])]
  
  dsum[, WHOletter := factor(WHO, levels = 1:11, labels=letters[1:11])]
  dsum_low <- dsum[as.numeric(WHOletter) %in% c(1:threshold)]
  dsum_high <- dsum[as.numeric(WHOletter) %in% c((threshold+1):11)]
  dsum_high[, WHOletter := factor(WHOletter, levels = letters[11:(threshold + 1)])]
  
  cc_low <- scales::seq_gradient_pal("#faa226", "white")(seq(0.2, 0.8, length.out=threshold))
  cc_high <- scales::seq_gradient_pal("white", "#267efa")(seq(0.2, 0.8, length.out=(11-threshold)))
  cc <- c(cc_low, cc_high)
  
  ggplot() +
    geom_bar(
      data = dsum_low,
      aes(x = axislabel, y = -prop, fill = WHOletter),
      width = .5, stat="identity") +
    geom_bar(
      data = dsum_high,
      aes(x = axislabel, y = prop, fill = WHOletter),
      width = .5, stat="identity") +
    scale_fill_manual(values = cc, name = "WHO score", labels = c("virus-free", 1:9, "died")) +
    ylab("proportion") +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9, face = "bold"),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 10)) + 
    scale_y_continuous(limits = c(-.85,.25),
                       breaks = c(-.75, -.50, -.25, 0, 0.25),
                       labels = c("75%","50%", "25%","0%", "25%")) +
    geom_hline(yintercept = 0, color = "grey96") +
    ggtitle(group) +
    coord_flip()  
}
  
### Day 14 - full data set

contain_all_14 <- data.table(rx = rep(c(0,1), each = 11), WHO = rep(c(1:11), 2),
                      N= c(38, 66, 143, 48, 16, 38, 20, 5, 15, 26, 39,
                           39, 67, 156, 49, 11, 25, 31, 7, 17, 15, 34))


### Day 28 - full data set

contain_all_28 <- data.table(rx = rep(c(0,1), each = 11), WHO = rep(c(1:11), 2),
                          N= c(84, 73, 121, 46, 9, 11, 8, 6, 7, 9, 70,
                               99, 64, 126, 39, 10, 12, 7, 4, 8, 12, 57))

### Generate plots

p_all_14 <- diverge_plot(contain_all_14, "WHO score at 14 days", 7)
p_all_28 <- diverge_plot(contain_all_28, "WHO score at 28 days", 7)

### Combine & save

p <- ggarrange(p_all_14, p_all_28, nrow = 2, common.legend = TRUE, legend = "right")
ggsave(file = "diverge_plot_14_28_cutoff_6_7.png", plot = p, width = 6, height = 3, scale = 1.5)

