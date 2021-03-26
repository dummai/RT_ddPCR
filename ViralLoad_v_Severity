library(tidyverse)

theme_set(theme_bw())

dat  <- read.csv("Data/SenseSpec_formated.csv", stringsAsFactors=FALSE)
dat <- dat[dat$Dx != '',]
dat <- dat[which(dat$ddpcr_log >= 4.190 & dat$ddpcr_log <= 6.854),]

g1 <- ggplot(dat, aes(x = Dx, y = ddpcr_log, col = Dx ))+
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  xlab('') +
  scale_y_continuous(name = "Viral load measured by RT-ddPCR (log10-copies/ml)", breaks = 0:7) +
  scale_color_discrete(name = "Diagnosis") +
  coord_cartesian(ylim = c(4,7)) +
  theme(
    panel.grid = element_blank()
  )

ggsave("Plots/FigureS5.jpeg",g1 , width = 6, height = 5 ,units = "in", dpi = 600)


t.test(ddpcr_log ~ Dx, dat)  
wilcox.test(ddpcr_log ~ Dx, dat)

m1 <- lm(ddpcr_log ~ Dx + Day + serotype, dat)
summary(m1)
