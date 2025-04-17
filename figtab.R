library(car)
library(lavaan)
library(semPlot)

N = 100000
A = rbinom(N, 1, .5)
B = rbinom(N, 1, .5)
C = rbinom(N, 1, .5)
D = rnorm(N, 0, 1)
residual = rnorm(N, 0, 1)

y = 1 + 0.20*C + residual

df = data.frame(A, B, C, D, y)
Anova(aov(y~A*B*C*D, data = df), type="II")

## SEM

k = 6 
w = 10
loading = 0.5

generate_lavaan_model = function(k, w) {
  model_lines = sapply(1:k, function(f) {
    lhs = paste0("F", f, " =~ ")
    rhs = paste0("F", f, "_item", 1:w, collapse = " + ")
    paste0(lhs, rhs)
  })
  paste(model_lines, collapse = "\n")
}

residual_sd = sqrt(1 - loading^2)
latent_vars = matrix(rnorm(N * k), nrow = N, ncol = k)
colnames(latent_vars) = paste0("F", 1:k)
  
latent_vars[, "F3"] = 0 + 0.15*latent_vars[, "F1"] + rnorm(N, 0, 1)
  
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
summary(fit)
semPaths(fit, what = "paths", layout = 'circle2')
