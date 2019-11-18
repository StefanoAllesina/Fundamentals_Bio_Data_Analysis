library(tidyverse)
genomes <- read_csv("data/genomes.csv")
genomes <- genomes %>% 
  group_by(Year) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(Year)
costs <- read_csv("data/sequencing_cost.csv")
genomes <- genomes %>% 
  inner_join(costs)

# model only intercept
p1 <- glm(n ~ 1, data = genomes, family = "poisson")
plot(genomes$n, p1$fitted.values, log = "xy"); abline(c(0,1))
AIC(p1)

# model Year
p2 <- glm(n ~ Year, data = genomes, family = "poisson")
plot(genomes$n, p2$fitted.values, log = "xy"); abline(c(0,1))
AIC(p2)

# model Cost
p3 <- glm(n ~ Dollars_per_Mb, data = genomes, family = "poisson")
plot(genomes$n, p3$fitted.values); abline(c(0,1))
AIC(p3)


p2 <- glm(n ~ Year, data = genomes, family = "poisson")
plot(genomes$n, p2$fitted.values, log = "xy"); abline(c(0,1))
AIC(p2)

# contrast with negative binomial
p2nb <- MASS::glm.nb(n ~ Year, data = genomes)
plot(genomes$n, p2nb$fitted.values, log = "xy"); abline(c(0,1))
AIC(p2nb)




