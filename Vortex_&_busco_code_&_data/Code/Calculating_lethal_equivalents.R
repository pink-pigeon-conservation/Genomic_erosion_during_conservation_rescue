
library("ggplot2")

data <- read.csv("../Data/longevity.csv") 



data$logAgedays <- log((data$AgeDays +1))

model1 <- lm(data$logAgedays ~ Inbreeding)
summary(model1)
confint(model1, level = 0.95)

FVlongconfidence <- ggplot(data, aes(Inbreeding, logAgedays)) + geom_point() + geom_smooth(method = lm, size = 1,color = "black") + xlab("Inbreeding Coefficient (F)") + ylab("Longevity") + theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(size = 1.2))
FVlongconfidence
