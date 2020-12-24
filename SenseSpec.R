library(readr)
library(exact2x2)

dat <- read_csv("Data/SenseSpec_formatted.csv") #Load RT-ddPCR versus qRT-PCR results

dat$qpcr_detect <- dat$qpcr_log > 0
dat$qpcr_detect[is.na(dat$qpcr_detect)] <- F

log.lloq <- log10(5065) #Detection limit of RT-ddPCR

dat$ddpcr_detect <- dat$ddpcr_log > log.lloq
dat$ddpcr_detect[is.na(dat$ddpcr_detect)] <- F


#Overall sensitivity
tb <- table(dat[dat$simpledx == "DEN",c("ddpcr_detect","qpcr_detect")])
print(tb)
mcnemar.test(tb)
mcnemar.exact(tb)

prop.test(148,156, correct = F)
prop.test(141,156, correct = F)
prop.test(154,156, correct = F)
prop.test(34,40, correct = F)


#Specificity
tb1 <- table(dat[dat$simpledx == "Non-dengue",c("ddpcr_detect","qpcr_detect")])
print(tb1)
#mcnemar.test(tb)
prop.test(40,40)
prop.test(40,40, correct = F)

#Den1
tb <- table(dat[dat$serotype == "DEN1",c("ddpcr_detect","qpcr_detect")])
print(tb)
mcnemar.test(tb)
mcnemar.exact(tb)
prop.test(39,40,correct = F)#ddPCR
prop.test(38,40,correct = F)#qPCR


#Den2
tb <- table(dat[dat$serotype == "DEN2",c("ddpcr_detect","qpcr_detect")])
print(tb)
mcnemar.test(tb)
mcnemar.exact(tb)
prop.test(37,40,correct = F)#ddPCR and qPCR

#Den3
tb <- table(dat[dat$serotype == "DEN3",c("ddpcr_detect","qpcr_detect")])
print(tb)
mcnemar.test(tb)
mcnemar.exact(tb)
prop.test(36,40,correct = F)#ddPCR
prop.test(30,40,correct = F)#qPCR

#Den4
tb <- table(dat[dat$serotype == "DEN4",c("ddpcr_detect","qpcr_detect")])
print(tb)
mcnemar.test(tb)
mcnemar.exact(tb)
prop.test(36,36, correct = F)
