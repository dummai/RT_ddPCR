library(readr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

dat <- read_csv("Data/InterLab_ddPCR_repeated.csv")

#LLOQ cut-off here
log.lloq <- 2.337
log.lloq.spec <-log10(1000*((10^log.lloq)*10)/140)

#ULOQ cut-off here
log.uloq.spec <- log10(100000 * 10 * 1000/140)

#AF vs SI
t_dat <- dat
t_dat$OK <- "Yes"
t_dat$OK[which(t_dat$af_log < log.lloq.spec | 
                 t_dat$af_log > log.uloq.spec |
                 t_dat$si_log < log.lloq.spec |
                 t_dat$si_log > log.uloq.spec)] <- "No"
m <- lm(si_log ~ af_log, t_dat[which(t_dat$OK == "Yes"),])
sm <- summary(m)
rsq <- format(round(sm$r.squared,4), nsmall = 4)
cr <- cor.test(~si_log + af_log, data = t_dat[t_dat$OK == "Yes",] )
r <- format(round(cr$estimate,4),nsmall = 4)

g11 <-ggplot(t_dat[t_dat$OK == "Yes",], aes(x = af_log, y = si_log, col = Serotype))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se = F, col = "black") +
  geom_abline(slope = 1, intercept = 0, lty = 2, col ="deepskyblue") +
  geom_point(data = t_dat[t_dat$OK == "No",],aes(x = af_log, y = si_log, col = Serotype), shape = 1)+
  #annotate("text", x = 6.5, y = 2, label = paste0("R^2 == ",rsq), parse=TRUE)+
  annotate("text", x = 6.5, y = 2, label = paste0("italic(r) == ",r), parse=TRUE)+
  coord_cartesian(ylim=c(0,9), xlim=c(0,9)) +
  scale_x_continuous(breaks = seq(0,9,2)) +
  scale_y_continuous(breaks = seq(0,9,2)) +
  #ggtitle(Sero)+
  xlab("Log10[DENV RNA count by AF (copies/ml)]") +
  ylab("Log10[DENV RNA count by SI (copies/ml)]") 

#Bland-Altman-Plots (percentage)
t_dat <- t_dat[t_dat$OK == "Yes",]
t_dat$mean <- (t_dat$si_log + t_dat$af_log)/2
t_dat$diff <- t_dat$si_log - t_dat$af_log
t_dat$perCV <- 100 * t_dat$diff/t_dat$mean


label_shift <- 2
label_shift_raw <- 0.05
sz <- 4

g21 <- ggplot(t_dat, aes(x = mean, y = perCV, col = Serotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(t_dat$perCV, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat$perCV, na.rm = T) - (1.96 * sd(t_dat$perCV, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat$perCV, na.rm = T) + (1.96 * sd(t_dat$perCV, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 0, y = mean(t_dat$perCV, na.rm = T) + label_shift, label = paste0("Mean = ",round(mean(t_dat$perCV),1),"%"))+
  annotate("text", hjust = 0, size =sz,x = 0, y = mean(t_dat$perCV, na.rm = T) + (1.96 * sd(t_dat$perCV, na.rm = T)) + label_shift, label = paste0("+1.96SD = ",round(mean(t_dat$perCV, na.rm = T) + (1.96 * sd(t_dat$perCV, na.rm = T)),1),"%"))+
  annotate("text",  hjust = 0,size =sz,x = 0, y = mean(t_dat$perCV, na.rm = T) - (1.96 * sd(t_dat$perCV, na.rm = T)) + label_shift, label = paste0("-1.96SD = ",round(mean(t_dat$perCV, na.rm = T) - (1.96 * sd(t_dat$perCV, na.rm = T)),1),"%"))+
  ylab("Percent difference compared to mean\n(SI - AF)")+
  xlab("Log10[Mean count of SI and AF (copies/ml)]") +

  coord_cartesian(xlim=c(0,9),ylim=c(-25,25)) +
  scale_x_continuous(breaks = seq(0,9,2))

g21_raw <- ggplot(t_dat, aes(x = mean, y = diff)) +
  geom_smooth(method ="lm", col = "black")+
  #geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5, aes(col = Serotype)) +
  geom_hline(yintercept = mean(t_dat$diff, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat$diff, na.rm = T) - (1.96 * sd(t_dat$diff, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat$diff, na.rm = T) + (1.96 * sd(t_dat$diff, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 3, y = mean(t_dat$diff, na.rm = T) + label_shift_raw, label = paste0("Mean = ",round(mean(t_dat$diff),3),""))+
  annotate("text", hjust = 0, size =sz,x = 3, y = mean(t_dat$diff, na.rm = T) + (1.96 * sd(t_dat$diff, na.rm = T)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(t_dat$diff, na.rm = T) + (1.96 * sd(t_dat$diff, na.rm = T)),3),""))+
  annotate("text",  hjust = 0,size =sz,x = 3, y = mean(t_dat$diff, na.rm = T) - (1.96 * sd(t_dat$diff, na.rm = T)) + label_shift_raw, label = paste0("-1.96SD = ",round(mean(t_dat$diff, na.rm = T) - (1.96 * sd(t_dat$diff, na.rm = T)),3),""))+
  ylab("Log10[Difference count between\nSI and AF (copies/ml)]")+
  xlab("Log10[Average count of SI and AF (copies/ml)]") +
  
  coord_cartesian(xlim=c(3,7),ylim=c(-1.4,0.6)) +
  #scale_y_continuous(breaks = seq(-1.4,0.6,0.2)) +
  scale_x_continuous(breaks = seq(0,7,1))  

#BIO vs AF
t_dat1 <- dat
t_dat1$OK <- "Yes"
t_dat1$OK[which(t_dat1$af_log < log.lloq.spec | 
                 t_dat1$af_log > log.uloq.spec |
                 t_dat1$bio_log < log.lloq.spec |
                 t_dat1$bio_log > log.uloq.spec)] <- "No"
m <- lm(bio_log ~ af_log, t_dat1[which(t_dat1$OK == "Yes"),])
sm <- summary(m)
rsq <- format(round(sm$r.squared,4), nsmall = 4)
cr <- cor.test(~bio_log + af_log, data = t_dat1[t_dat1$OK == "Yes",] )
r <- format(round(cr$estimate,4),nsmall = 4)


g12 <-ggplot(t_dat1[t_dat1$OK == "Yes",], aes(x = af_log, y = bio_log, col = Serotype))+
  geom_abline(slope = 1, intercept = 0, lty = 2,  col ="deepskyblue") +
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se = F, col = "black") +
  geom_point(data = t_dat1[t_dat1$OK == "No",],aes(x = af_log, y = bio_log, col = Serotype), shape = 1)+
  #annotate("text", x = 6.5, y = 2, label = paste0("R^2 == `",rsq,"`"), parse=TRUE)+
  annotate("text", x = 6.5, y = 2, label = paste0("italic(r) == ",r), parse=TRUE)+
  coord_cartesian(ylim=c(0,9), xlim=c(0,9)) +
  scale_x_continuous(breaks = seq(0,9,2)) +
  scale_y_continuous(breaks = seq(0,9,2)) +
  #ggtitle(Sero)+
  xlab("Log10[DENV RNA count by AF (copies/ml)]") +
  ylab("Log10[DENV RNA count by BIO (copies/ml)]") 

#Bland-Altman-Plots (percentage)
t_dat1 <- t_dat1[t_dat1$OK == "Yes",]
t_dat1$mean <- (t_dat1$bio_log + t_dat1$af_log)/2
t_dat1$diff <- t_dat1$bio_log - t_dat1$af_log
t_dat1$perCV <- 100 * t_dat1$diff/t_dat1$mean


label_shift <- 2
sz <- 4

g22 <- ggplot(t_dat1, aes(x = mean, y = perCV, col = Serotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(t_dat1$perCV, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat1$perCV, na.rm = T) - (1.96 * sd(t_dat1$perCV, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat1$perCV, na.rm = T) + (1.96 * sd(t_dat1$perCV, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 0, y = mean(t_dat1$perCV, na.rm = T) + label_shift, label = paste0("Mean = ",round(mean(t_dat1$perCV),1),"%"))+
  annotate("text", hjust = 0, size =sz,x = 0, y = mean(t_dat1$perCV, na.rm = T) + (1.96 * sd(t_dat1$perCV, na.rm = T)) + label_shift, label = paste0("+1.96SD = ",round(mean(t_dat1$perCV, na.rm = T) + (1.96 * sd(t_dat1$perCV, na.rm = T)),1),"%"))+
  annotate("text",  hjust = 0,size =sz,x = 0, y = mean(t_dat1$perCV, na.rm = T) - (1.96 * sd(t_dat1$perCV, na.rm = T)) + label_shift, label = paste0("-1.96SD = ",round(mean(t_dat1$perCV, na.rm = T) - (1.96 * sd(t_dat1$perCV, na.rm = T)),1),"%"))+
  ylab("Percent difference compared to mean\n(BIO - AF)")+
  xlab("Log10[Mean count of BIO and AF (copies/ml)]") +
  
  coord_cartesian(xlim=c(0,9),ylim=c(-25,25)) +
  scale_x_continuous(breaks = seq(0,9,2))


g22_raw <- ggplot(t_dat1, aes(x = mean, y = diff)) +
  geom_smooth(method ="lm", col = "black")+
  #geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5, aes(col = Serotype)) +
  geom_hline(yintercept = mean(t_dat1$diff, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat1$diff, na.rm = T) - (1.96 * sd(t_dat1$diff, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat1$diff, na.rm = T) + (1.96 * sd(t_dat1$diff, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 3, y = mean(t_dat1$diff, na.rm = T) + label_shift_raw, label = paste0("Mean = ",round(mean(t_dat1$diff),3),""))+
  annotate("text", hjust = 0, size =sz,x = 3, y = mean(t_dat1$diff, na.rm = T) + (1.96 * sd(t_dat1$diff, na.rm = T)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(t_dat1$diff, na.rm = T) + (1.96 * sd(t_dat1$diff, na.rm = T)),3),""))+
  annotate("text",  hjust = 0,size =sz,x = 3, y = mean(t_dat1$diff, na.rm = T) - (1.96 * sd(t_dat1$diff, na.rm = T)) + label_shift_raw, label = paste0("-1.96SD = ",format(round(mean(t_dat1$diff, na.rm = T) - (1.96 * sd(t_dat1$diff, na.rm = T)),3),nsmall = 3),""))+
  ylab("Log10[Difference count between\nBIO and AF (copies/ml)]")+
  xlab("Log10[Average count of BIO and AF (copies/ml)]") +
  
  coord_cartesian(xlim=c(3,7),ylim=c(-1.4,0.6)) +
  #scale_y_continuous(breaks = seq(-1.4,0.6,0.2)) +
  scale_x_continuous(breaks = seq(0,7,1)) 

#SI VS BIO
t_dat2 <- dat
t_dat2$OK <- "Yes"
t_dat2$OK[which(t_dat2$si_log < log.lloq.spec | 
                  t_dat2$si_log > log.uloq.spec |
                  t_dat2$bio_log < log.lloq.spec |
                  t_dat2$bio_log > log.uloq.spec)] <- "No"
m <- lm(si_log ~ bio_log, t_dat2[which(t_dat2$OK == "Yes"),])
sm <- summary(m)
rsq <- format(round(sm$r.squared,4), nsmall = 4)
cr <- cor.test(~bio_log + si_log, data = t_dat2[t_dat2$OK == "Yes",] )
r <- format(round(cr$estimate,4),nsmall = 4)

g13 <-ggplot(t_dat2[t_dat2$OK == "Yes",], aes(x = bio_log, y = si_log, col = Serotype))+
 
  geom_abline(slope = 1, intercept = 0, lty = 2, col ="deepskyblue") +
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se = F, col = "black") +
  geom_point(data = t_dat2[t_dat2$OK == "No",],aes(x = bio_log, y = si_log, col = Serotype), shape = 1)+
  #annotate("text", x = 6.5, y = 2, label = paste0("R^2 == `",rsq,"`"), parse=TRUE)+
  annotate("text", x = 6.5, y = 2, label = paste0("italic(r) == ",r), parse=TRUE)+
  coord_cartesian(ylim=c(0,9), xlim=c(0,9)) +
  scale_x_continuous(breaks = seq(0,9,2)) +
  scale_y_continuous(breaks = seq(0,9,2)) +
  #ggtitle(Sero)+
  xlab("Log10[DENV RNA count by BIO (copies/ml)]") +
  ylab("Log10[DENV RNA count by SI (copies/ml)]") 

#Bland-Altman-Plots (percentage)
t_dat2 <- t_dat2[t_dat2$OK == "Yes",]
t_dat2$mean <- (t_dat2$bio_log + t_dat2$si_log)/2
t_dat2$diff <- t_dat2$si_log - t_dat2$bio_log
t_dat2$perCV <- 100 * t_dat2$diff/t_dat2$mean



label_shift <- 2
sz <- 4

g23 <- ggplot(t_dat2, aes(x = mean, y = perCV, col = Serotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(t_dat2$perCV, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat2$perCV, na.rm = T) - (1.96 * sd(t_dat2$perCV, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat2$perCV, na.rm = T) + (1.96 * sd(t_dat2$perCV, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 0, y = mean(t_dat2$perCV, na.rm = T) + label_shift, label = paste0("Mean = ",round(mean(t_dat2$perCV),1),"%"))+
  annotate("text", hjust = 0, size =sz,x = 0, y = mean(t_dat2$perCV, na.rm = T) + (1.96 * sd(t_dat2$perCV, na.rm = T)) + label_shift, label = paste0("+1.96SD = ",round(mean(t_dat2$perCV, na.rm = T) + (1.96 * sd(t_dat2$perCV, na.rm = T)),1),"%"))+
  annotate("text",  hjust = 0,size =sz,x = 0, y = mean(t_dat2$perCV, na.rm = T) - (1.96 * sd(t_dat2$perCV, na.rm = T)) + label_shift, label = paste0("-1.96SD = ",round(mean(t_dat2$perCV, na.rm = T) - (1.96 * sd(t_dat2$perCV, na.rm = T)),1),"%"))+
  ylab("Percent difference compared to mean\n(SI - BIO)")+
  xlab("Log10[Mean count of SI and BIO (copies/ml)]") +
  coord_cartesian(xlim=c(0,9),ylim=c(-25,25)) +
  scale_x_continuous(breaks = seq(0,9,2)) 

g23_raw <- ggplot(t_dat2, aes(x = mean, y = diff)) +
  geom_smooth(method ="lm", col = "black")+
  #geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5, aes(col = Serotype)) +
  geom_hline(yintercept = mean(t_dat2$diff, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(t_dat2$diff, na.rm = T) - (1.96 * sd(t_dat2$diff, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(t_dat2$diff, na.rm = T) + (1.96 * sd(t_dat2$diff, na.rm = T)), colour = "red") +
  annotate("text",  hjust = 0, size =sz,x = 3, y = mean(t_dat2$diff, na.rm = T) + label_shift_raw, label = paste0("Mean = ",round(mean(t_dat2$diff),3),""))+
  annotate("text", hjust = 0, size =sz,x = 3, y = mean(t_dat2$diff, na.rm = T) + (1.96 * sd(t_dat2$diff, na.rm = T)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(t_dat2$diff, na.rm = T) + (1.96 * sd(t_dat2$diff, na.rm = T)),3),""))+
  annotate("text",  hjust = 0,size =sz,x = 3, y = mean(t_dat2$diff, na.rm = T) - (1.96 * sd(t_dat2$diff, na.rm = T)) + label_shift_raw, label = paste0("-1.96SD = ",round(mean(t_dat2$diff, na.rm = T) - (1.96 * sd(t_dat2$diff, na.rm = T)),3),""))+
  ylab("Log10[Difference count between\nSI and BIO (copies/ml)]")+
  xlab("Log10[Average count of SI and BIO (copies/ml)]") +
  
  coord_cartesian(xlim=c(3,7),ylim=c(-1.4,0.6)) +
  #scale_y_continuous(breaks = seq(-1.4,0.6,0.2)) +
  scale_x_continuous(breaks = seq(0,7,1)) 

fig7 <- ggarrange(g11,g12,g13,nrow = 2, ncol = 2, labels ="AUTO")
fig8 <- ggarrange(g21,g22,g23,nrow = 2, ncol = 2, labels ="AUTO")
fig8_raw <- ggarrange(g21_raw,g22_raw,g23_raw,nrow = 2, ncol = 2, labels ="AUTO")

ggsave("Plots/Figure7.jpeg",fig7, width = 10, height = 8 ,units = "in", dpi = 600)
ggsave("Plots/Figure8.jpeg",fig8, width = 10, height = 8 ,units = "in", dpi = 600)
ggsave("Plots/Figure8_raw.jpeg",fig8_raw, width = 10, height = 8 ,units = "in", dpi = 600)

cor.test(~diff + mean,t_dat)$estimate
sd(t_dat$diff, na.rm = T)
cor.test(~diff + mean,t_dat1)$estimate
sd(t_dat1$diff, na.rm = T)
cor.test(~diff + mean,t_dat2)$estimate
sd(t_dat2$diff, na.rm = T)
mean(t_dat2$diff, na.rm = T)
