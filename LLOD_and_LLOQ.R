library(readr)
library(reshape2)
library(ggplot2)
library(VFP)
library(VCA)
library(cowplot)
library(ggplotify)
library(ggpubr)



theme_set(theme_bw())

dat <- read_csv("Raw_Data/LoDLoQ_data.csv")

mdat <- melt(dat,id.vars = c("qubit",    "Dilution" ))
mdat$value <- as.numeric(as.character(mdat$value))

mdat$detected <- !is.na(mdat$value)
amdat <- aggregate(detected ~ Dilution, mdat, FUN = mean)
mdat$logdil <- log10(mdat$Dilution)
amdat$logdil <- log10(amdat$Dilution)

log_model <- glm(detected ~ logdil, data = mdat, family = "binomial")

x <- seq(-5,0,0.0001)
y <- predict(log_model, list(logdil = x ),type="response")

pred <- data.frame(x = x , y = y, expx = 10^x)

#find 95% cut-off
cx <- min(pred$x[pred$y > 0.95 & pred$y < 0.96])
dil95 <- 10^(cx)

g1 <- ggplot(pred, aes(x =expx , y = y))+
  geom_line() +
  #geom_hline(yintercept = 0.95, lty = 2) +
  geom_segment(aes(x = 0.000008, xend = 10^cx, y = 0.95, yend =0.95), col  ="blue", lty =5) +
  #geom_vline(xintercept = 10^cx, lty = 2) + 
  geom_segment(aes(x = 10^cx, xend = 10^cx, y = -0.06, yend =0.95), col  ="blue", lty =5) +
  geom_point(data = amdat, aes(x=Dilution, y = detected)) +
  annotate("text", label = "0.95", x = 0.0000055, y = 0.95, col = "blue")+
  annotate("text", label = "1:663", x = dil95*1.4, y = -0.08, col = "blue")+
  ylab("Esitmated probability of detection") +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  xlab("Dilution")+
  scale_x_log10(breaks = amdat$Dilution,labels = c("1:20000", "1:10000", "1:2000", "1:1000", "1:100", "1:10", "Undiluted")) +
  coord_cartesian(xlim = c(0.000015,1), ylim =c(0,1), clip = "off")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
        ,panel.grid.minor = element_blank())



#LLOD95%
ant_conc <- mean(mdat$value[ mdat$Dilution == 1])
raw <- dil95* 10^ant_conc
log.raw <- log10(raw)
print(raw)
print(log.raw)


#Avarage
ave_dat <- aggregate(value ~ Dilution, mdat, mean)
sd_dat <- aggregate(value ~ Dilution, mdat, sd)

ave_dat$perCV <- 100 * sd_dat$value/ave_dat$value
ave_dat$antValue <- log10(ave_dat$Dilution * 10^ave_dat$value[ave_dat$Dilution == 1])

ave_dat$antValue_lo <- 0.8 * ave_dat$antValue 
ave_dat$antValue_up <- 1.2 * ave_dat$antValue 
fave_dat <- ave_dat[ ave_dat$Dilution >= dil95,]

mdat <- merge(mdat, ave_dat[,c("Dilution","antValue")], by = "Dilution", all.x = T)

m <- lm(formula = value ~ log10(Dilution), mdat)

m1 <- lm(value ~  antValue, mdat[which(mdat$Dilution > 10^cx),])


fave_dat$reg <- 1.083*log10(fave_dat$Dilution) + 4.413

x <-  m1$coefficients[1] / (0.8 - m1$coefficients[2])
ant_x <- x*m1$coefficients[2] + m1$coefficients[1]
predict(m1, data.frame(antValue = 2.921), interval = 'confidence')


g2 <- ggplot(ave_dat, aes(x =antValue, y = value)) +
  geom_point(shape = 3, size = 5) +
  geom_abline(slope = 1, intercept = 0)+
  annotate("text", label = round(x,3) , x =2.7, y = -0.5, col = "blue")+
  annotate("text", label = round(ant_x,3) , x =-0.6, y = ant_x , col = "blue")+
  #geom_line(data = ave_dat , aes(x=antValue,y=antValue_lo), lty=3) +
  #geom_line(data = ave_dat , aes(x=antValue,y=antValue_up), lty=3) +
  geom_ribbon(data= ave_dat, 
              aes(ymin=antValue_lo,ymax=antValue_up), fill="grey", alpha=0.5) +
  #geom_smooth(data = ave_dat[which(ave_dat$Dilution > 10^cx),], aes( x = antValue, y = value), method = 'lm', se = F)
  #geom_vline(xintercept = log10((10^cx)*(10^max(ave_dat$antValue))), lty = 2, col = "grey20") + 
  #geom_abline(slope = m1$coefficients[2], intercept = m1$coefficients[1], col = "red", lty = 5) +
  geom_segment(aes(x = 0.9, xend = 5.25, y = 0.9*m1$coefficients[2] + m1$coefficients[1], yend =5.25*m1$coefficients[2] + m1$coefficients[1]), col  ="red", lty =5) +
  geom_segment(aes(x=x,xend =x, y = -0.3, yend=ant_x), col = "blue", lty = 5)+
  geom_segment(aes(x=-0.3,xend =x, y = ant_x, yend=ant_x), col = "blue", lty = 5)+
  geom_point(data = mdat, aes(x=antValue, y=value )) +
  xlab("Log10[Anticipated template count (copies/reaction)]") +
  ylab("Log10[Measured template count (copies/reaction)]") +
  scale_y_continuous(breaks = seq(0,6,1)) +
  scale_x_continuous(breaks = seq(0,6,1)) +
  coord_cartesian(xlim = c(-0.1,5), ylim =c(-0.1,5.5), clip = "off")+
  theme(panel.grid.minor = element_blank()
        ,axis.title.x = element_text(vjust=-1)
        ,axis.title.y = element_text(vjust=10)
        ,plot.margin = unit(c(0.1, 0.1, 0.2, 0.5), "in")
        )


mdat$mea <- as.numeric(as.character(mdat$variable))


mdat$day <- 3
mdat$day[mdat$mea <7] <- 2
mdat$day[mdat$mea <4] <- 1

mdat$rep <- 1 + (mdat$mea %% 3)

mdat$value <- as.numeric(as.character(mdat$value))

mdat$expvalue <- 10^(mdat$value)


VCA.mdat<- anovaVCA(value~day, mdat, by="Dilution")


mat.total <- getMat.VCA(VCA.mdat)

tot.all <- fit.vfp(mat.total, 1:10)


plot(tot.all, type ="cv",mar=c(4, 4, 0.1, 0.1)
     ,ylim=c(0, 100) 
     ,Prediction = list(y=c(20,10,5.9), cex = 1)
     ,Xlabel  = list(text="Log10[Measured template count (copies/reaction)]",cex=1)
     ,Ylabel = list(cex=1)
     ,Grid = list(lty = 1, lwd = 1)
     ,Model=F
     ,Title = NULL
     ,Crit=NULL
)

g3 <- as.ggplot(function(){par(cex = 0.8)
plot(tot.all, type ="cv",mar=c(4, 4, 0.1, 0.1)
     ,ylim=c(0, 100) 
     ,Prediction = list(y=20, cex = 1)
     ,Xlabel  = list(text="Log10[Measured template count (copies/reaction)]",cex=1)
     ,Ylabel = list(cex=1)
     ,Grid = list(lty = 1, lwd = 1)
     ,Model=F
     ,Title = NULL
     ,Crit=NULL
)}

)

g3 <- g3 + theme(plot.margin = unit(c(0.1, 0.1, 0, 0), "in")

)



gall <- ggarrange(g1,g3,g2,nrow = 2,ncol=2, labels = c("A","B","C"))
ggsave("Plots/Figure3_V3.jpeg",gall, width = 11, height = 9, units = "in", dpi = 600)
