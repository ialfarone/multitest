#####################################

rm(list=ls())
library(lavaan)
library(car)

#####################################

# ANOVA

niter = 10

aovsim = function(niter = 10, h1 = FALSE, N = 1000) {
  
  significanceCount = rep(NA,niter)
  ms = as.data.frame(matrix(NA, nrow = niter, ncol = 15))
  
  for(i in 1:niter){
    A = rbinom(N, 1, .5)
    B = rbinom(N, 1, .5)
    C = rbinom(N, 1, .5)
    D = rnorm(N, 0, 1)
    residual = rnorm(N, 0, 1)
    
    if(h1){
      y = 1 + 0.1*C + residual
    } else {
      y = 1 + residual
    }
    
    df = data.frame(A, B, C, D, y)
    
    ps = Anova(aov(y~A*B*C*D, data = df))$"Pr(>F)"
    ps = ps[!is.na(ps)]
    
    significanceCount[i] = sum(ps<0.05)
    ms[i, 1:length(ps)] = ps
  }

  colnames(ms) = rownames(Anova(aov(y~A*B*C*D,data = df)))[-nrow(Anova(aov(y~A*B*C*D,data = df)))]
  
  return(list(
    ms = ms,
    signCount = significanceCount,
    hist = hist(significanceCount),
    percSignCount = mean(significanceCount > 0)
  ))
}


StartAnova0 = Sys.time()
null = aovsim(niter = 5000, h1 = FALSE)
EndAnova0 = Sys.time()
(diffAnova0 = difftime(EndAnova0,StartAnova0,units="mins"))

StartAnova1 = Sys.time()
nonull = aovsim(niter = 5000, h1 = TRUE)
EndAnova1 = Sys.time()
(diffAnova1 = difftime(EndAnova1,StartAnova1,units="mins"))

# 1 - 0.95^15 (15 è il numero dei test fatti)

save.image("multitestResults.RData")

#####################################

# SEM

semsim = function(niter = 10, N = 2000, k = 6, w = 10, loading = 0.7, h1 = F) {
  
  generate_lavaan_model = function(k, w) {
    model_lines = sapply(1:k, function(f) {
      lhs = paste0("F", f, " =~ ")
      rhs = paste0("F", f, "_item", 1:w, collapse = " + ")
      paste0(lhs, rhs)
    })
    paste(model_lines, collapse = "\n")
  }
  
  significanceCount = rep(NA, niter)
  ms.sem = as.data.frame(matrix(NA, nrow = niter, ncol = 15))
  
  for(i in 1:niter){
    residual_sd = sqrt(1 - loading^2)
    latent_vars = matrix(rnorm(N * k), nrow = N, ncol = k)
    colnames(latent_vars) = paste0("F", 1:k)
    
    if (h1) {
      latent_vars[, "F3"] = 0 + 0.25*latent_vars[, "F1"] + rnorm(N, 0, 1)
    }
    
    observed_list = vector("list", k)
    for (j in 1:k) {
      indicators = matrix(
        loading * latent_vars[, j] + residual_sd * matrix(rnorm(N * w), nrow = N),
        nrow = N, ncol = w
      )
      colnames(indicators) = paste0("F", j, "_item", 1:w)
      observed_list[[j]] = indicators
    }
    
    df = as.data.frame(do.call(cbind, observed_list))
    
    modelM = generate_lavaan_model(k = 6,w = 10)
    model = paste(modelM,
                  "\n
              F6 ~ F1 + F2 + F3 + F4 + F5 \n
              F5 ~ F1 + F2 + F3 + F4 \n
              F4 ~ F1 + F2 + F3 \n
              F3 ~ F1 + F2 \n
              F2 ~ F1")
    
    fit = sem(model,df)
    pe = summary(fit)$pe
    ps = pe$pvalue[pe$op=="~"]
    
    significanceCount[i] = sum(ps < 0.05)
    #  print(paste("i:",i))
    #  print(paste("sign.count:",significanceCount[i]))
    ms.sem[i, 1:length(ps)] = ps
  } 
  
  colnames(ms.sem) = paste(pe$lhs[pe$op == "~"], pe$op[pe$op == "~"], pe$rhs[pe$op == "~"])
  
  return(list(
    ms.sem = ms.sem,
    significanceCount = significanceCount,
    hist = hist(significanceCount),
    percSigPaths = mean(significanceCount > 0)
  ))
}

StartSem0 = Sys.time()
null.sem = semsim(niter = 5000, h1 = FALSE)
EndSem0 = Sys.time()
(diffSem0 = difftime(EndSem0,StartSem0,units="mins"))

save.image("multitestResults.RData")

StartSem1 = Sys.time()
nonull.sem = semsim(niter = 5000, h1 = TRUE)
EndSem1 = Sys.time()
(diffSem1 = difftime(EndSem1,StartSem1,units="mins"))

#####################################

save.image("multitestResults.RData")

#####################################

load("multitestResults.RData")

res = nonull$ms
mean(res[,"B:C"]<0.05)
for(i in 1:nrow(res)){
  res[i,] = p.adjust(res[i,],method="fdr")
}
mean(res[,"B:C"]<0.05)

###################################

# New sim (N = 1000, simpler h1)

nonull = aovsim(niter = 1000, h1 = TRUE)

res = nonull$ms
mean(res[,"C"]<0.05)

for(i in 1:nrow(res)){
  res[i,] = p.adjust(res[i,],method="bonferroni")
}


mean(res[,"C"]<0.05)

