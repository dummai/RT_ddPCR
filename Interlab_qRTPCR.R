library(reshape2)
library(ggplot2)
library(readr)
library(ggpubr)

theme_set(theme_bw())

dat <- read_csv("Data/InterLab_qRTPCR.csv")
dat$Serotype <- dat$Serotope
dat <- dat[!is.na(dat$Serotype),]

#loLoQ <- log10(700)
loLoQ <- 0.0001

dat <- dat[ dat$si_log >= loLoQ & dat$af_log >= loLoQ, ]
m <- lm( si_log ~ af_log, dat)
sm <- summary(m)
rsq  <- paste0("R^2 == ", round(sm$r.squared,4))
cr <- cor.test(~si_log + af_log, dat)
r <- format(round(cr$estimate,4),nsmall = 4)

#Scatter + Regression plot
g1 <- ggplot(dat, aes( y = si_log, x = af_log, col = Serotype)) +
  geom_abline(intercept = 0,slope = 1, color = "deepskyblue", lty = 2)+
  geom_point() +
  geom_smooth(method='lm', se = F, color = "black")+
  annotate("text", x = 6.5, y = 2, label = paste0("italic(r) == ",r), parse=TRUE)+
  coord_cartesian(ylim=c(0,9), xlim=c(0,9)) +
  scale_x_continuous(breaks = seq(0,9,2)) +
  scale_y_continuous(breaks = seq(0,9,2)) +
  #ggtitle(Sero)+
  xlab("Log10[dengue RNA by AF (copies/ml)]") +
  ylab("Log10[dengue RNA by SI (copies/ml)]") 

#Bland-Altman-Plots (percentage)
tdat <- dat
ftdat <- tdat[ !is.na(tdat$si_log) & !is.na(tdat$af_log),]
ftdat$mean <- (ftdat$si_log + ftdat$af_log)/2
ftdat$diff <- ftdat$si_log - ftdat$af_log
ftdat$perCV <- 100 * ftdat$diff/ftdat$mean

label_shift <-4
label_shift_raw <-0.1

g2 <- ggplot(ftdat, aes(x = mean, y = perCV, col = Serotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(ftdat$perCV), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(ftdat$perCV) - (1.96 * sd(ftdat$perCV)), colour = "red") +
  geom_hline(yintercept = mean(ftdat$perCV) + (1.96 * sd(ftdat$perCV)), colour = "red") +
  #annotate("text", hjust = 0, x = 0, y = mean(ftdat$perCV) - label_shift, label = paste0(round(mean(ftdat$perCV),1),"%"))+
  annotate("text", hjust = 0, x = 0, y = mean(ftdat$perCV) + label_shift, label = paste0("Mean = ",round(mean(ftdat$perCV),1),"%"))+
  #annotate("text", x = 8.5, y = mean(ftdat$perCV) + (1.96 * sd(ftdat$perCV)) - label_shift, label = paste0(round(mean(ftdat$perCV) + (1.96 * sd(ftdat$perCV)),1),"%"))+
  annotate("text", hjust = 0, x = 0, y = mean(ftdat$perCV) + (1.96 * sd(ftdat$perCV)) + label_shift, label = paste0("+1.96SD = ",round(mean(ftdat$perCV) + (1.96 * sd(ftdat$perCV)),1),"%"))+
  #annotate("text", x = 8.5, y = mean(ftdat$perCV) - (1.96 * sd(ftdat$perCV)) - label_shift, label = paste0(round(mean(ftdat$perCV) - (1.96 * sd(ftdat$perCV)),1),"%"))+
  annotate("text", hjust = 0, x = 0 , y = mean(ftdat$perCV) - (1.96 * sd(ftdat$perCV)) + label_shift, label = paste0("-1.96SD = ",round(mean(ftdat$perCV) - (1.96 * sd(ftdat$perCV)),1),"%"))+
  ylab("Percent difference compared to mean\n(SI - AF)")+
  xlab("Log10[Mean count of SI and AF (copies/ml)]") +
  #ggtitle(Sero)+
  coord_cartesian(xlim=c(0,9),ylim=c(-50,50)) +
  scale_x_continuous(breaks = seq(0,9,2)) #+
  
g3 <- ggplot(ftdat, aes(x = mean, y = diff )) +
  geom_smooth(method ="lm", col = "black")+
  geom_point(alpha = 0.5, aes(col = Serotype)) +
  geom_hline(yintercept = mean(ftdat$diff), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(ftdat$diff) - (1.96 * sd(ftdat$diff)), colour = "red") +
  geom_hline(yintercept = mean(ftdat$diff) + (1.96 * sd(ftdat$diff)), colour = "red") +
  annotate("text", hjust = 0, x = 0, y = mean(ftdat$diff) + label_shift_raw, label = paste0("Mean = ",round(mean(ftdat$diff),3),""))+
  annotate("text", hjust = 0, x = 0, y = mean(ftdat$diff) + (1.96 * sd(ftdat$diff)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(ftdat$diff) + (1.96 * sd(ftdat$diff)),3),""))+
  annotate("text", hjust = 0, x = 0 , y = mean(ftdat$diff) - (1.96 * sd(ftdat$diff)) + label_shift_raw, label = paste0("-1.96SD = ",round(mean(ftdat$diff) - (1.96 * sd(ftdat$diff)),3),""))+
  ylab("Log10[Difference count between\nSI and AF (copies/ml)]")+
  xlab("Log10[Average count of SI and AF (copies/ml)]") +
  coord_cartesian(xlim=c(0,9),ylim=c(-2.5,2.5)) +
  scale_x_continuous(breaks = seq(0,9,2))

cor.test(~mean+diff,ftdat)
sd(ftdat$diff, na.rm = T)
mean(ftdat$diff, na.rm = T)

gall <- ggarrange(g1,g2,nrow = 1 , labels = "AUTO")
gall_raw <- ggarrange(g1,g3,nrow = 1 , labels = "AUTO")
ggsave(filename = "Plots/Figure9.jpeg",gall, width = 10, height = 4 ,units = "in", dpi = 600)
ggsave(filename = "Plots/Figure9_raw.jpeg",gall_raw, width = 10, height = 4 ,units = "in", dpi = 600)
