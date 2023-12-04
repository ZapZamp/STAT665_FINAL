library(tidyverse)
library(DescTools)
library(boot)
library(lme4)
library(ggplot2)
library(pROC)
set.seed(678)

#preparing data
{
masses <- read.table("masses.data", sep = ",", header = FALSE)
names <- c("bi_rads","age","shape","margin","density","severity")
colnames(masses) <- names

masses[masses == "?"] <- "mis"
masses$bi_rads <- replace(masses$bi_rads, masses$bi_rads == 55, 5)
raw <- masses
masses <- masses %>%
  filter(bi_rads > 0 & bi_rads < 6 & age != "mis") %>%
  mutate(
    severity = as.integer(severity),
    bi_rads = as.factor(bi_rads),
    age = as.numeric(age),
    shape = as.factor(shape),
    margin = as.factor(margin),
    density = as.factor(density)
  )
}

## exploration
fig1 <- table(raw$severity, raw$bi_rads)
fig2 <- table(masses$bi_rads, masses$severity)/nrow(masses)
#Measuring association: Bi-rad score vs Severity

#creating resampling space for bootstrap test for independence
{
counts <- table(masses$bi_rads, masses$severity)
cell.prob <- counts/nrow(masses)

sev_marg <- colSums(cell.prob)
bi_marg <- rowSums(cell.prob)
resamp <- matrix(bi_marg)%*%t(matrix(sev_marg))


resamp <- round(resamp*938,digits=0)
rownames(resamp) <- c("2","3","4","5")
colnames(resamp) <- c("0","1")
obs <- c()

for (i in 1:4) {
  for (j in 1:2) {
  obs <- append(
    obs,
    rep(
      c(rownames(resamp)[i],colnames(resamp)[j]),
      times = resamp[i,j])
    )
  if (i*j == 8) {resamp_long <- matrix(obs,ncol=2, byrow = TRUE)}
  }
}

colnames(resamp_long) <- c("bi_rads","severity")
resamp_long <- as.data.frame(resamp_long)
resamp_long %>%
  mutate(
    bi_rads <- factor(bi_rads, levels = c("2","3","4","5")),
    severity <- factor(severity, levels = c("0","1"), labels = c("Benign","Malignant"))
  )
}

#conducting the bootstrap
{
tauc_boot <- function(d, i){
  d2 <- d[i,]
  freqtab <- table(d2$bi_rads, d2$severity)
  return(StuartTauC(freqtab))
}
B1 <- 1000

results <- boot(resamp_long, statistic = tauc_boot, R = B1)

st_star <- nrow(masses)*sum((results$t-mean(results$t))^2)/(B1-1)
z_boot <- sqrt(nrow(masses))*StuartTauC(counts)/st_star
}

#values for fig 3
asym_res<- StuartTauC(counts,conf.level = .95)
low <- qnorm(.05/2, mean = 0, sd = 1, lower.tail = TRUE)
upp <- qnorm((1-.05/2), mean = 0, sd = 1, lower.tail = TRUE)
rej <- if(z_boot<low | z_boot>upp) {TRUE} else {FALSE}
fig3_values <- rbind(c("bootstrap", asym_res[1],z_boot, rej),c("asymptotic", asym_res[1], paste(asym_res[2],asym_res[3]), (asym_res[3]<0 | asym_res[2]>0)))

#power: testing at variable possible values of proportions within category 4
{
  loop <- 0
  base <- rep(table(masses$bi_rads)/nrow(masses), times = 2) #proportion of sample within each category

power_5 <- c(.995, .98, .95, .05, .005, .02, .05, .95) * base
power_10 <- c(.995, .98, .90, .05, .005, .02, .10, .95) * base
power_50 <- c(.995, .98, .50, .05, .005, .02, .50, .95) * base
power_92 <- c(.995, .98, .10, .05, .005, .02, .90, .95) * base
power_93 <- c(.995, .98, .07, .05, .005, .02, .93, .95) * base


probs <- rbind(power_5,power_10,power_50,power_92,power_93)
ss <- c(50, 100, 150, 300)
power_data <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(power_data) <- c("cat4_prob","ss","tau","asym_reject","boot_sig","z")

for (a in 1:5) {
  for (b in 1:4) {
    sims <- rmultinom(n = 500, size = ss[b], prob = probs[a,]) #each column is a list of frequencies
    for (x in 1:100) {
      #convert sims to frequency table
      sample_freq <- matrix(sims[,x], nrow = 4, byrow = FALSE)
      
      #convert freq tables to resampled space
      {
      cell.prob <- sample_freq/ss[b]
      
      sev_marg <- colSums(cell.prob)
      bi_marg <- rowSums(cell.prob)
      resamp <- matrix(bi_marg)%*%t(matrix(sev_marg))
      resamp <- round(resamp*ss[b],digits=0)
      rownames(resamp) <- c("2","3","4","5")
      colnames(resamp) <- c("0","1")
      obs <- c()
      
      for (i in 1:4) {
        for (j in 1:2) {
          obs <- append(
            obs,
            rep(
              c(rownames(resamp)[i],colnames(resamp)[j]),
              times = resamp[i,j])
          )
          if (i*j == 8) {resamp_long <- matrix(obs,ncol=2, byrow = TRUE)}
        }
      }
      
      colnames(resamp_long) <- c("bi_rads","severity")
      resamp_long <- as.data.frame(resamp_long)
      resamp_long %>%
        mutate(
          bi_rads <- factor(bi_rads, levels = c("2","3","4","5")),
          severity <- factor(severity, levels = c("0","1"), labels = c("Benign","Malignant"))
        )
      }
     
      #conduct bootstrap and save stuart's tau c. Save value and z stat
      {
        B1 <- 500
        results <- boot(resamp_long, statistic = tauc_boot, R = B1)
        
        st_star <- ss[b]*sum((results$t-mean(results$t))^2)/(B1-1)
        z_boot <- sqrt(ss[b])*StuartTauC(sample_freq)/st_star
        if (is.nan(z_boot)) {
          hold <- results
          print(resamp)}
        loop <- loop + 1
        print(loop)
      }
      #save values
      power_data <- add_row(power_data, cat4_prob = a, ss = b, tau = StuartTauC(sample_freq), asym_reject = (StuartTauC(sample_freq, conf.level = .95)[3]<0 | StuartTauC(sample_freq, conf.level = .95)[2]>0), boot_sig = st_star,z = z_boot)
    }
  }
}
}
#
power_data_ammended <- power_data %>% mutate(
  boot_reject = (if(z_boot<low | z_boot>upp) {TRUE} else {FALSE})
)


fig5_ss1 <- power_data_ammended %>% filter(ss == 1) %>%
  group_by(cat4_prob) %>% summarise(asym_reject_rate = mean(asym_reject, na.rm = TRUE),
                                    boot_reject_rate = mean(asym_reject, na.rm = TRUE))

fig5_ss2 <- power_data_ammended %>% filter(ss == 2) %>%
  group_by(cat4_prob) %>% summarise(asym_reject_rate = mean(asym_reject, na.rm = TRUE),
                                    boot_reject_rate = mean(asym_reject, na.rm = TRUE))

fig5_ss3 <- power_data_ammended %>% filter(ss == 3) %>%
  group_by(cat4_prob) %>% summarise(asym_reject_rate = mean(asym_reject, na.rm = TRUE),
                                    boot_reject_rate = mean(asym_reject, na.rm = TRUE))

fig5_ss4 <- power_data_ammended %>% filter(ss == 4) %>%
  group_by(cat4_prob) %>% summarise(asym_reject_rate = mean(asym_reject, na.rm = TRUE),
                                    boot_reject_rate = mean(asym_reject, na.rm = TRUE))


## fitting with "missing"
{
fit_bi <- glm(severity ~ 1 + bi_rads + age,
                 data = masses,
                 family = binomial('logit'))

fit_att <- glm(severity ~ 1 + margin + density + shape + age,
               data = masses,
               family = binomial('logit'))


vis_roc <- roc(masses$severity,
                        predict(fit_bi, masses, 
                                'response'),
                        smoothed = TRUE)

att_roc <- roc(masses$severity,
               predict(fit_att, masses, 
                       'response'),
               smoothed = TRUE)
}
Fig6_bic <- BIC(fit_bi,fit_att)
Fig6_aic <- AIC(fit_bi,fit_att)
Fig7 <- ggroc(list(Bi_rads = vis_roc, Attributes = att_roc))
