library(readr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

dat <- read_csv("Data/SI_ddpcr_v_qpcr.csv")

#LLOQ cut-off here
log.lloq <- 2.337
log.lloq.spec <-log10(1000*((10^log.lloq)*10)/140)

#ULOQ cut-off here
log.uloq.spec <- log10(100000 * 10 * 1000/140)

#Prep data for Bland-Altman
dat$mean <- (dat$qpcr_log + dat$ddpcr_log)/2
dat$diff <- dat$ddpcr_log - dat$qpcr_log
dat$perCV <- 100 * dat$diff/dat$mean

label_shift <- 2
label_shift_raw <- 0.1

gall <- list()
galt <- list()
graw <- list()

for(sero in unique(dat$Serotype)){

fdat_sero <- dat[which(dat$Serotype == sero & 
                         dat$ddpcr_log >= log.lloq.spec &
                         dat$ddpcr_log <= log.uloq.spec),]

udat_sero <- dat[which(dat$Serotype == sero & 
                         (dat$ddpcr_log < log.lloq.spec |
                         dat$ddpcr_log > log.uloq.spec)),]

m <- lm(ddpcr_log ~ qpcr_log, fdat_sero)
sm <- summary(m)
#print(sm)
rsq <- format(round(sm$r.squared,4),nsmall = 4)
cr <- cor.test(~ ddpcr_log + qpcr_log, fdat_sero)
r <- format(round(cr$estimate, 4), nsmall = 4)
print(cr)
print(cor.test(~mean+diff,fdat_sero))
print(sd(fdat_sero$diff, na.rm = T))


g1 <- ggplot(fdat_sero, aes(x = qpcr_log, y = ddpcr_log)) +
  ggtitle(sero)+
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  coord_cartesian(ylim = c(0,9), xlim = c(0,9)) +
  geom_smooth(method='lm', formula= y~x, se = F, col = "orange")+
  geom_point() +
  geom_point(data = udat_sero, aes(x = qpcr_log, y = ddpcr_log), shape = 1) + 
  annotate("text",parse = T, x = 6, y = 2, label = paste0("italic(r) == ",r))+
  scale_x_continuous(name = "Log10[DENV RNA count by qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Log10[DENV RNA count by RT-ddPCR (copies/ml)]", breaks = seq(0,9,2))+
  theme(plot.title = element_text(hjust = 0.5))
gall[[sero]] <- g1

g2 <- ggplot(fdat_sero, aes(x = mean, y = perCV))+
  ggtitle(sero)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(fdat_sero$perCV, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(fdat_sero$perCV, na.rm = T) - (1.96 * sd(fdat_sero$perCV, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(fdat_sero$perCV, na.rm = T) + (1.96 * sd(fdat_sero$perCV, na.rm = T)), colour = "red") +
  coord_cartesian(ylim = c(-40,40), xlim = c(0,9)) +
  annotate("text", hjust = 0, x = 0, y = mean(fdat_sero$perCV, na.rm = T) + label_shift, label = paste0("Mean = ",round(mean(fdat_sero$perCV, na.rm = T),1),"%"))+
  annotate("text", hjust = 0, x = 0, y = mean(fdat_sero$perCV, na.rm = T) + (1.96 * sd(fdat_sero$perCV, na.rm = T)) + label_shift, label = paste0("+1.96SD = ",round(mean(fdat_sero$perCV, na.rm = T) + (1.96 * sd(fdat_sero$perCV, na.rm = T)),1),"%"))+
  annotate("text", hjust = 0, x = 0 , y = mean(fdat_sero$perCV, na.rm = T) - (1.96 * sd(fdat_sero$perCV, na.rm = T)) + label_shift, label = paste0("-1.96SD = ",round(mean(fdat_sero$perCV, na.rm = T) - (1.96 * sd(fdat_sero$perCV, na.rm = T)),1),"%"))+
  scale_x_continuous(name = "Log10[Mean of RT-ddPCR and qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Percent difference compared to mean\n(RT-ddPCR - qRT-PCR)", breaks = seq(-40,40,20))+
  theme(plot.title = element_text(hjust = 0.5))

galt[[sero]] <- g2


g3 <- ggplot(fdat_sero, aes(x = mean, y = diff))+
  ggtitle(sero)+
  geom_smooth(method ="lm", col = "black")+
  #geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(fdat_sero$diff, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(fdat_sero$diff, na.rm = T) - (1.96 * sd(fdat_sero$diff, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(fdat_sero$diff, na.rm = T) + (1.96 * sd(fdat_sero$diff, na.rm = T)), colour = "red") +
  coord_cartesian(ylim = c(-2,1), xlim = c(2,9)) +
  annotate("text", hjust = 0, x = 2, y = mean(fdat_sero$diff, na.rm = T) + label_shift_raw, label = paste0("Mean = ",round(mean(fdat_sero$diff, na.rm = T),3),""))+
  annotate("text", hjust = 0, x = 2, y = mean(fdat_sero$diff, na.rm = T) + (1.96 * sd(fdat_sero$diff, na.rm = T)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(fdat_sero$diff, na.rm = T) + (1.96 * sd(fdat_sero$diff, na.rm = T)),3),""))+
  annotate("text", hjust = 0, x = 2 , y = mean(fdat_sero$diff, na.rm = T) - (1.96 * sd(fdat_sero$diff, na.rm = T)) + label_shift_raw, label = paste0("-1.96SD = ",round(mean(fdat_sero$diff, na.rm = T) - (1.96 * sd(fdat_sero$diff, na.rm = T)),3),""))+
  scale_x_continuous(name = "Log10[Average count of RT-ddPCR and qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Log10[Difference count between\nRT-ddPCR - qRT-PCR (copies/ml)]")



graw[[sero]] <- g3


}

ggall <- ggarrange(plotlist = gall,labels = "AUTO",nrow = 2, ncol = 2)
ggalt <- ggarrange(plotlist = galt,labels = "AUTO",nrow = 2, ncol = 2)
ggraw <- ggarrange(plotlist = graw,labels = "AUTO",nrow = 2, ncol = 2)


ggsave("Plots/Figure4.jpeg",ggall, width = 10, height = 10 ,units = "in", dpi = 600)
ggsave("Plots/Figure5.jpeg",ggalt, width = 10, height = 10 ,units = "in", dpi = 600)
ggsave("Plots/Figure5_raw.jpeg",ggraw, width = 10, height = 10 ,units = "in", dpi = 600)


#All combined
fdat <- dat[which(dat$ddpcr_log >= log.lloq.spec &
                         dat$ddpcr_log <= log.uloq.spec),]

udat <- dat[which(       (dat$ddpcr_log < log.lloq.spec |
                            dat$ddpcr_log > log.uloq.spec)),]

m <- lm(ddpcr_log ~ qpcr_log, fdat)
sm <- summary(m)
#print(sm)
rsq <- format(round(sm$r.squared,4),nsmall = 4)
cr <- cor.test(~ ddpcr_log + qpcr_log, fdat)
r <- format(round(cr$estimate,4),nsmall = 4)

g1 <- ggplot(fdat, aes(x = qpcr_log, y = ddpcr_log)) +
  #ggtitle(sero)+
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  coord_cartesian(ylim = c(0,9), xlim = c(0,9)) +
  geom_smooth(method='lm', formula= y~x, se = F, col = "orange")+
  geom_point() +
  geom_point(data = udat, aes(x = qpcr_log, y = ddpcr_log), shape = 1) + 
  annotate("text",parse = T, x = 6, y = 2, label = paste0("italic(r) == ",r))+
  scale_x_continuous(name = "Log10[DENV RNA count by qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Log10[DENV RNA count by RT-ddPCR (copies/ml)]", breaks = seq(0,9,2))+
  theme(plot.title = element_text(hjust = 0.5))


g2 <- ggplot(fdat, aes(x = mean, y = perCV))+
  #ggtitle(sero)+
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(fdat$perCV, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(fdat$perCV, na.rm = T) - (1.96 * sd(fdat$perCV, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(fdat$perCV, na.rm = T) + (1.96 * sd(fdat$perCV, na.rm = T)), colour = "red") +
  coord_cartesian(ylim = c(-40,40), xlim = c(0,9)) +
  annotate("text", hjust = 0, x = 0, y = mean(fdat$perCV, na.rm = T) + label_shift, label = paste0("Mean = ",round(mean(fdat$perCV, na.rm = T),1),"%"))+
  annotate("text", hjust = 0, x = 0, y = mean(fdat$perCV, na.rm = T) + (1.96 * sd(fdat$perCV, na.rm = T)) + label_shift, label = paste0("+1.96SD = ",round(mean(fdat$perCV, na.rm = T) + (1.96 * sd(fdat$perCV, na.rm = T)),1),"%"))+
  annotate("text", hjust = 0, x = 0 , y = mean(fdat$perCV, na.rm = T) - (1.96 * sd(fdat$perCV, na.rm = T)) + label_shift, label = paste0("-1.96SD = ",round(mean(fdat$perCV, na.rm = T) - (1.96 * sd(fdat$perCV, na.rm = T)),1),"%"))+
  scale_x_continuous(name = "Log10[Mean of RT-ddPCR and qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Percent difference compared to mean\n(RT-ddPCR - qRT-PCR)", breaks = seq(-40,40,20))+
  theme(plot.title = element_text(hjust = 0.5))

g2_raw <- ggplot(fdat, aes(x = mean, y = diff))+
  geom_smooth(method ="lm", col = "black")+
  #ggtitle(sero)+
  #geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(fdat$diff, na.rm = T), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(fdat$diff, na.rm = T) - (1.96 * sd(fdat$diff, na.rm = T)), colour = "red") +
  geom_hline(yintercept = mean(fdat$diff, na.rm = T) + (1.96 * sd(fdat$diff, na.rm = T)), colour = "red") +
  coord_cartesian(ylim = c(-2,1), xlim = c(2,9)) +
  annotate("text", hjust = 0, x = 2, y = mean(fdat$diff, na.rm = T) + label_shift_raw, label = paste0("Mean = ",round(mean(fdat$diff, na.rm = T),3),""))+
  annotate("text", hjust = 0, x = 2, y = mean(fdat$diff, na.rm = T) + (1.96 * sd(fdat$diff, na.rm = T)) + label_shift_raw, label = paste0("+1.96SD = ",round(mean(fdat$diff, na.rm = T) + (1.96 * sd(fdat$diff, na.rm = T)),3),""))+
  annotate("text", hjust = 0, x = 2 , y = mean(fdat$diff, na.rm = T) - (1.96 * sd(fdat$diff, na.rm = T)) + label_shift_raw, label = paste0("-1.96SD = ",round(mean(fdat$diff, na.rm = T) - (1.96 * sd(fdat$diff, na.rm = T)),3),""))+
  scale_x_continuous(name = "Log10[Average count of RT-ddPCR and qRT-PCR (copies/ml)]", breaks = seq(0,9,2))+
  scale_y_continuous(name = "Log10[Difference count between\nRT-ddPCR - qRT-PCR (copies/ml)]")

cor.test(~mean+diff,fdat)
sd(fdat$diff, na.rm = T)
mean(fdat$diff, na.rm = T)

ggcom <- ggarrange(g1,g2,labels = "AUTO",nrow = 1, ncol = 2)
ggsave("Plots/Figure6.jpeg",ggcom, width = 10, height = 5 ,units = "in", dpi = 600)

ggcom_raw <- ggarrange(g1,g2_raw,labels = "AUTO",nrow = 1, ncol = 2)
ggsave("Plots/Figure6_raw.jpeg",ggcom_raw, width = 10, height = 5 ,units = "in", dpi = 600)
