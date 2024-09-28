# All distributions together n = 1000, n = 10000

# Each Distribution, Sample size increasing, grid width increasing 

# Simulation Study 

n <- 1000000
B = 100

res_pdf_L1_norm = rep(0,B)
res_pdf_L1_beta = rep(0,B)
res_pdf_L1_gamma = rep(0,B)
res_pdf_L1_logis = rep(0,B)
res_pdf_L1_t = rep(0,B)

res_pdf_L2_norm = rep(0,B)
res_pdf_L2_beta = rep(0,B)
res_pdf_L2_gamma = rep(0,B)
res_pdf_L2_logis = rep(0,B)
res_pdf_L2_t = rep(0,B)

# Parameters 

alpha=2;beta=2;
mu = 10;sigma=1;
scale_val = 0.05;df=3;

# Ranges

range_norm = c(-Inf,Inf)
range_beta = c(0,1)
range_gamma = c(0,Inf)
range_logistic = c(-Inf,Inf)
range_t = c(-Inf,Inf)

# Pdfs 

pdf_beta <- function(x) {
  dbeta(x,alpha,beta)
}

pdf_norm <- function(x) {
  dnorm(x, mean = mu, sd = sigma)
}

pdf_gamma <- function(x) {
  dgamma(x, alpha, beta)
}

pdf_logistic <- function(x) {
  dlogis(x)
}

pdf_t <- function(x) {
  dt(x,df)
}


# Simulation Process 

set.seed(8)

for (i in 1:B) {
  
  # Normal Distribution 
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_pdf_L1_norm[i] = L1_Distance(y,pdf_norm,range_norm)
  res_pdf_L2_norm[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Beta Distribution 
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_pdf_L1_beta[i] = L1_Distance(y,pdf_beta,range_beta)
  res_pdf_L2_beta[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Gamma Distribution 
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_pdf_L1_gamma[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_pdf_L2_gamma[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Logistic Distribution
  
  x = sort(rlogis(n))
  grid  <- seq(-20,20, length.out=101)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_pdf_L1_logis[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_pdf_L2_logis[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # t-Distribution
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 0.5
  y <- floor(x/delta)
  y = y*delta
  
  res_pdf_L1_t[i] = L1_Distance(y,pdf_t,range_t)
  res_pdf_L2_t[i] = L2_Distance(y,pdf_t,range_t)
  
}

# Convergence Limits

pdf_L1_norm_limit = mean(res_pdf_L1_norm)
pdf_L2_norm_limit = mean(res_pdf_L2_norm)
pdf_L1_gamma_limit = mean(res_pdf_L1_gamma)
pdf_L2_gamma_limit = mean(res_pdf_L2_gamma)
pdf_L1_beta_limit = mean(res_pdf_L1_beta)
pdf_L2_beta_limit = mean(res_pdf_L2_beta)
pdf_L1_logis_limit = mean(res_pdf_L1_logis)
pdf_L2_logis_limit = mean(res_pdf_L2_logis)
pdf_L1_t_limit = mean(res_pdf_L1_t)
pdf_L2_t_limit = mean(res_pdf_L2_t)

# Plots

par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,2),cex.lab=1.5,cex.axis=1.5,
    font.axis=1,cex.main=1.5)

colors <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00")

boxplot(res_pdf_L1_norm,res_pdf_L1_beta,res_pdf_L1_gamma,res_pdf_L1_logis,res_pdf_L1_t,main = "L1 Distance vs Distributions\n n = 1000",
        names = c("Normal", "Beta", "Gamma", "Logistic", "Student's t"),ylab="L1 Distance",
        col = colors,ylim = c(0,0.3))

boxplot(res_pdf_L2_norm,res_pdf_L2_beta,res_pdf_L2_gamma,res_pdf_L2_logis,res_pdf_L2_t,main = "L2 Distance vs Distributions\n n = 1000",
        names = c("Normal","Beta","Gamma","Logistic","Student's t"),ylab="L2 Distance",
        col = colors,ylim = c(0,0.3))

# Simulations Within the Distribution 

colors <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00")

# Colors for Normal Distribution 

c11 = darken("#92CAEB",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors1 = c("#92CAEB",c11,c21,c31,c41)

# Colors for Beta Distribution

c11 = darken("#F3CF70",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors2 = c("#F3CF70",c11,c21,c31,c41)

# Colors for Gamma Distribution 

c11 = darken("#66D7A5",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors3 = c("#66D7A5",c11,c21,c31,c41)

# Colors for Logistic Distribution

c11 = darken("#E6A4C6",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors4 = c("#E6A4C6",c11,c21,c31,c41)

# Colors for Student's t Distribution

c11 = darken("#D55E00",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors5 = c("#D55E00",c11,c21,c31,c41)


# Normal Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

res_norm_pdf_n_11 = rep(0,B)
res_norm_pdf_n_21 = rep(0,B)
res_norm_pdf_n_31 = rep(0,B)
res_norm_pdf_n_41 = rep(0,B)
res_norm_pdf_n_51 = rep(0,B)

res_norm_pdf_n_12 = rep(0,B)
res_norm_pdf_n_22 = rep(0,B)
res_norm_pdf_n_32 = rep(0,B)
res_norm_pdf_n_42 = rep(0,B)
res_norm_pdf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Normal Distribution - n1
  
  x <- sort(rnorm(n[1], mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_n_11[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_n_12[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - n2
  
  x <- sort(rnorm(n[2], mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_n_21[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_n_22[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - n3
  
  x <- sort(rnorm(n[3], mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_n_31[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_n_32[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - n4
  
  x <- sort(rnorm(n[4], mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_n_41[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_n_42[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - n5
  
  x <- sort(rnorm(n[5], mean = mu,sd=sigma))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_n_51[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_n_52[i] = L2_Distance(y,pdf_norm,range_norm)
  
}

grDevices::pdf("L1-PDF-1-Normal-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_11,res_norm_pdf_n_21,res_norm_pdf_n_31,res_norm_pdf_n_41,res_norm_pdf_n_51,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L1 Distance",
        col = colors1,ylim=c(0,0.4),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L1_norm_limit, x1 = 10^5, y1 = pdf_L1_norm_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("L2-PDF-1-Normal-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_12,res_norm_pdf_n_22,res_norm_pdf_n_32,res_norm_pdf_n_42,res_norm_pdf_n_52,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L2 Distance",
        col = colors1,ylim=c(0,0.3),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L2_norm_limit, x1 = 10^5, y1 = pdf_L2_norm_limit, col = "red", lty = 2)
dev.off()

# Beta Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

res_beta_pdf_n_11 = rep(0,B)
res_beta_pdf_n_21 = rep(0,B)
res_beta_pdf_n_31 = rep(0,B)
res_beta_pdf_n_41 = rep(0,B)
res_beta_pdf_n_51 = rep(0,B)

res_beta_pdf_n_12 = rep(0,B)
res_beta_pdf_n_22 = rep(0,B)
res_beta_pdf_n_32 = rep(0,B)
res_beta_pdf_n_42 = rep(0,B)
res_beta_pdf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Beta Distribution - n1
  
  x   <- sort(rbeta(n[1],alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_n_11[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_n_12[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - n2
  
  x   <- sort(rbeta(n[2],alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_n_21[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_n_22[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - n3
  
  x   <- sort(rbeta(n[3],alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_n_31[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_n_32[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - n4
  
  x   <- sort(rbeta(n[4],alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_n_41[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_n_42[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - n5
  
  x   <- sort(rbeta(n[5],alpha,beta))
  grid <- seq(0,5, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_n_51[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_n_52[i] = L2_Distance(y,pdf_beta,range_beta)
  
}

grDevices::pdf("L1-PDF-1-Beta-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11,res_beta_pdf_n_21,res_beta_pdf_n_31,res_beta_pdf_n_41,res_beta_pdf_n_51,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L1 Distance",
        col = colors2,ylim=c(0,0.4),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("L2-PDF-1-Beta-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12,res_beta_pdf_n_22,res_beta_pdf_n_32,res_beta_pdf_n_42,res_beta_pdf_n_52,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L2 Distance",
        col = colors2,ylim=c(0,0.3),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Gamma Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

res_gamma_pdf_n_11 = rep(0,B)
res_gamma_pdf_n_21 = rep(0,B)
res_gamma_pdf_n_31 = rep(0,B)
res_gamma_pdf_n_41 = rep(0,B)
res_gamma_pdf_n_51 = rep(0,B)

res_gamma_pdf_n_12 = rep(0,B)
res_gamma_pdf_n_22 = rep(0,B)
res_gamma_pdf_n_32 = rep(0,B)
res_gamma_pdf_n_42 = rep(0,B)
res_gamma_pdf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Gamma Distribution - n1
  
  x   <- sort(rgamma(n[1],alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_n_11[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_n_12[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n2
  
  x   <- sort(rgamma(n[2],alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_n_21[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_n_22[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n3
  
  x   <- sort(rgamma(n[3],alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_n_31[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_n_32[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n4
  
  x   <- sort(rgamma(n[4],alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_n_41[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_n_42[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - n5
  
  x   <- sort(rgamma(n[5],alpha,beta))
  grid <- seq(-10,10, length.out=101)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_n_51[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_n_52[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
}

grDevices::pdf("L1-PDF-1-Gamma-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11,res_gamma_pdf_n_21,res_gamma_pdf_n_31,res_gamma_pdf_n_41,res_gamma_pdf_n_51,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L1 Distance",
        col = colors3,ylim=c(0,0.4),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L1_gamma_limit, x1 = 10^5, y1 = pdf_L1_gamma_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("L2-PDF-1-Gamma-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12,res_gamma_pdf_n_22,res_gamma_pdf_n_32,res_gamma_pdf_n_42,res_gamma_pdf_n_52,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L2 Distance",
        col = colors3,ylim=c(0,0.3),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L2_gamma_limit, x1 = 10^5, y1 = pdf_L2_gamma_limit, col = "red", lty = 2)
dev.off()


# Logistic Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

res_logis_pdf_n_11 = rep(0,B)
res_logis_pdf_n_21 = rep(0,B)
res_logis_pdf_n_31 = rep(0,B)
res_logis_pdf_n_41 = rep(0,B)
res_logis_pdf_n_51 = rep(0,B)

res_logis_pdf_n_12 = rep(0,B)
res_logis_pdf_n_22 = rep(0,B)
res_logis_pdf_n_32 = rep(0,B)
res_logis_pdf_n_42 = rep(0,B)
res_logis_pdf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Logistic Distribution - n1
  
  x = sort(rlogis(n[1]))
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_n_11[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_n_12[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n2
  
  x = sort(rlogis(n[2]))
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_n_21[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_n_22[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n3
  
  x = sort(rlogis(n[3]))
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_n_31[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_n_32[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n4
  
  x = sort(rlogis(n[4]))
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_n_41[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_n_42[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - n5
  
  x = sort(rlogis(n[5]))
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_n_51[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_n_52[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
}

grDevices::pdf("L1-PDF-1-Logistic-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11,res_logis_pdf_n_21,res_logis_pdf_n_31,res_logis_pdf_n_41,res_logis_pdf_n_51,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L1 Distance",
        col = colors4,xlab="Sample Size",ylim=c(0,0.4))
segments(x0 = 0.29, y0 = mean(res_logis_pdf_n_51), x1 = 10^5, y1 = mean(res_logis_pdf_n_51), col = "red", lty = 2)
dev.off()

grDevices::pdf("L2-PDF-1-Logistic-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12,res_logis_pdf_n_22,res_logis_pdf_n_32,res_logis_pdf_n_42,res_logis_pdf_n_52,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L2 Distance",
        col = colors4,ylim=c(0,0.3),xlab="Sample Size")
segments(x0 = 0.29, y0 = mean(res_logis_pdf_n_52), x1 = 10^5, y1 = mean(res_logis_pdf_n_52), col = "red", lty = 2)
dev.off()

# Student's t Distribution - Sample Size Simulation 

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

res_t_pdf_n_11 = rep(0,B)
res_t_pdf_n_21 = rep(0,B)
res_t_pdf_n_31 = rep(0,B)
res_t_pdf_n_41 = rep(0,B)
res_t_pdf_n_51 = rep(0,B)

res_t_pdf_n_12 = rep(0,B)
res_t_pdf_n_22 = rep(0,B)
res_t_pdf_n_32 = rep(0,B)
res_t_pdf_n_42 = rep(0,B)
res_t_pdf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # t Distribution - n1
  
  x = sort(rt(n[1],df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 1
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_n_11[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_n_12[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - n2
  
  x = sort(rt(n[2],df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 1
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_n_21[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_n_22[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - n3
  
  x = sort(rt(n[3],df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 1
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_n_31[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_n_32[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - n4
  
  x = sort(rt(n[4],df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 1
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_n_41[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_n_42[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - n5
  
  x = sort(rt(n[5],df = df))
  grid  <- seq(-20,20, length.out=101)
  delta <- 1
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_n_51[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_n_52[i] = L2_Distance(y,pdf_t,range_t)
  
}

grDevices::pdf("L1-PDF-1-t-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11,res_t_pdf_n_21,res_t_pdf_n_31,res_t_pdf_n_41,res_t_pdf_n_51,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L1 Distance",
        col = colors5,ylim=c(0,0.5),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L1_t_limit, x1 = 10^5, y1 = pdf_L1_t_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("L2-PDF-1-t-n.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12,res_t_pdf_n_22,res_t_pdf_n_32,res_t_pdf_n_42,res_t_pdf_n_52,
        names = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6)),ylab="L2 Distance",
        col = colors5,ylim=c(0,0.3),xlab="Sample Size")
segments(x0 = 0.29, y0 = pdf_L2_t_limit, x1 = 10^5, y1 = pdf_L2_t_limit, col = "red", lty = 2)
dev.off()

# Normal Distribution - Grid-Width Simulation 

n = c(10^3)
B = 100

res_norm_pdf_gw_11 = rep(0,B)
res_norm_pdf_gw_21 = rep(0,B)
res_norm_pdf_gw_31 = rep(0,B)
res_norm_pdf_gw_41 = rep(0,B)
res_norm_pdf_gw_51 = rep(0,B)

res_norm_pdf_gw_12 = rep(0,B)
res_norm_pdf_gw_22 = rep(0,B)
res_norm_pdf_gw_32 = rep(0,B)
res_norm_pdf_gw_42 = rep(0,B)
res_norm_pdf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Normal Distribution - GW1
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, by=1)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_gw_11[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_gw_12[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - GW2
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, by=0.8)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_gw_21[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_gw_22[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - GW3
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, by =0.6)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_gw_31[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_gw_32[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - GW4
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, by =0.4)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_gw_41[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_gw_42[i] = L2_Distance(y,pdf_norm,range_norm)
  
  # Normal Distribution - GW5
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  grid <- seq(-10,10, by =0.2)
  delta <- min(diff(grid))
  
  y <- floor(x/delta)
  y = y*delta
  
  res_norm_pdf_gw_51[i] = L1_Distance(y,pdf_norm,range_norm)
  res_norm_pdf_gw_52[i] = L2_Distance(y,pdf_norm,range_norm)
  
}

grDevices::pdf("L1-PDF-1-Normal-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11,res_norm_pdf_gw_21,res_norm_pdf_gw_31,res_norm_pdf_gw_41,res_norm_pdf_gw_51,
        names = c("1 (9)", "0.8 (11)", "0.6 (14)", "0.4 (21)", "0.2 (41)"),ylab="L1 Distance",
        col = colors1,xlab="Grid Width (k)",ylim=c(0,0.7))
dev.off()

grDevices::pdf("L2-PDF-1-Normal-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12,res_norm_pdf_gw_22,res_norm_pdf_gw_32,res_norm_pdf_gw_42,res_norm_pdf_gw_52,
        names = c("1 (9)", "0.8 (11)", "0.6 (14)", "0.4 (21)", "0.2 (41)"),ylab="L2 Distance",
        col = colors1,xlab="Grid Width (k)",ylim=c(0,0.4))
dev.off()

# Beta Distribution - Grid Width Simulation 

n = 1000
B = 100

res_beta_pdf_gw_11 = rep(0,B)
res_beta_pdf_gw_21 = rep(0,B)
res_beta_pdf_gw_31 = rep(0,B)
res_beta_pdf_gw_41 = rep(0,B)
res_beta_pdf_gw_51 = rep(0,B)

res_beta_pdf_gw_12 = rep(0,B)
res_beta_pdf_gw_22 = rep(0,B)
res_beta_pdf_gw_32 = rep(0,B)
res_beta_pdf_gw_42 = rep(0,B)
res_beta_pdf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Beta Distribution - GW1
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, by = 0.25)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_gw_11[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_gw_12[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - GW2
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, by = 0.2)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_gw_21[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_gw_22[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - GW3
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, by = 0.15)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_gw_31[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_gw_32[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - GW4
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, by = 0.1)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_gw_41[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_gw_42[i] = L2_Distance(y,pdf_beta,range_beta)
  
  # Beta Distribution - GW5
  
  x   <- sort(rbeta(n,alpha,beta))
  grid <- seq(0,5, by = 0.05)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_beta_pdf_gw_51[i] = L1_Distance(y,pdf_beta,range_beta)
  res_beta_pdf_gw_52[i] = L2_Distance(y,pdf_beta,range_beta)
  
}

grDevices::pdf("L1-PDF-1-Beta-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11,res_beta_pdf_gw_21,res_beta_pdf_gw_31,res_beta_pdf_gw_41,res_beta_pdf_gw_51,
        names = c("0.25 (5)", "0.20 (6)", "0.15 (7)", "0.10 (11)", "0.05 (21)"),ylab="L1 Distance",
        col = colors2,xlab="Grid Width (k)",ylim=c(0,0.7))
dev.off()

grDevices::pdf("L2-PDF-1-Beta-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12,res_beta_pdf_gw_22,res_beta_pdf_gw_32,res_beta_pdf_gw_42,res_beta_pdf_gw_52,
        names = c("0.25 (5)", "0.20 (6)", "0.15 (7)", "0.10 (11)", "0.05 (21)"),ylab="L2 Distance",
        col = colors2,xlab="Grid Width (k)",ylim=c(0,0.4))
dev.off()

# Gamma Distribution - Grid Width Simulation 

n = 1000
B = 100

res_gamma_pdf_gw_11 = rep(0,B)
res_gamma_pdf_gw_21 = rep(0,B)
res_gamma_pdf_gw_31 = rep(0,B)
res_gamma_pdf_gw_41 = rep(0,B)
res_gamma_pdf_gw_51 = rep(0,B)

res_gamma_pdf_gw_12 = rep(0,B)
res_gamma_pdf_gw_22 = rep(0,B)
res_gamma_pdf_gw_32 = rep(0,B)
res_gamma_pdf_gw_42 = rep(0,B)
res_gamma_pdf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Gamma Distribution - GW1
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, by=1)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_gw_11[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_gw_12[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW2
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, by = 0.8)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_gw_21[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_gw_22[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW3
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, by = 0.6)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_gw_31[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_gw_32[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW4
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, by = 0.4)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_gw_41[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_gw_42[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
  # Gamma Distribution - GW5
  
  x   <- sort(rgamma(n,alpha,beta))
  grid <- seq(-10,10, by = 0.2)
  delta <- min(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_gamma_pdf_gw_51[i] = L1_Distance(y,pdf_gamma,range_gamma)
  res_gamma_pdf_gw_52[i] = L2_Distance(y,pdf_gamma,range_gamma)
  
}

grDevices::pdf("L1-PDF-1-Gamma-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11,res_gamma_pdf_gw_21,res_gamma_pdf_gw_31,res_gamma_pdf_gw_41,res_gamma_pdf_gw_51,
        names = c("1 (7)", "0.8 (8)", "0.6 (11)", "0.4 (16)", "0.2 (31)"),ylab="L1 Distance",
        col = colors3,xlab="Grid Width (k)",ylim=c(0,0.7))
dev.off()

grDevices::pdf("L2-PDF-1-Gamma-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12,res_gamma_pdf_gw_22,res_gamma_pdf_gw_32,res_gamma_pdf_gw_42,res_gamma_pdf_gw_52,
        names = c("1 (7)", "0.8 (8)", "0.6 (11)", "0.4 (16)", "0.2 (31)"),ylab="L2 Distance",
        col = colors3,xlab="Grid Width (k)",ylim=c(0,0.4))
dev.off()

# Logistic Distribution - Grid Width Simulation 

n = 1000
B = 100

res_logis_pdf_gw_11 = rep(0,B)
res_logis_pdf_gw_21 = rep(0,B)
res_logis_pdf_gw_31 = rep(0,B)
res_logis_pdf_gw_41 = rep(0,B)
res_logis_pdf_gw_51 = rep(0,B)

res_logis_pdf_gw_12 = rep(0,B)
res_logis_pdf_gw_22 = rep(0,B)
res_logis_pdf_gw_32 = rep(0,B)
res_logis_pdf_gw_42 = rep(0,B)
res_logis_pdf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Logistic Distribution - GW1
  
  x = sort(rlogis(n))
  grid  <- seq(-20,20, by = 1.4)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_gw_11[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_gw_12[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW2
  
  x = sort(rlogis(n))
  grid  <- seq(-20,20, by = 1.2)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_gw_21[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_gw_22[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW3
  
  x = sort(rlogis(n))
  grid  <- seq(-20,20, by = 1)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_gw_31[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_gw_32[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW4
  
  x = sort(rlogis(n))
  grid  <- seq(-20,20, by = 0.8)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_gw_41[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_gw_42[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
  # Logistic Distribution - GW5
  
  x = sort(rlogis(n))
  min(x)
  grid  <- seq(-20,20, by = 0.6)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_logis_pdf_gw_51[i] = L1_Distance(y,pdf_logistic,range_logistic)
  res_logis_pdf_gw_52[i] = L2_Distance(y,pdf_logistic,range_logistic)
  
}

grDevices::pdf("L1-PDF-1-Logistic-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11,res_logis_pdf_gw_21,res_logis_pdf_gw_31,res_logis_pdf_gw_41,res_logis_pdf_gw_51,
        names = c("1.4 (20)", "1.2 (23)", "1.0 (28)", "0.8 (34)", "0.6 (46)"),ylab="L1 Distance",
        col = colors4,xlab="Grid Width (k)",ylim=c(0,0.7))
dev.off()

grDevices::pdf("L2-PDF-1-Logistic-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12,res_logis_pdf_gw_22,res_logis_pdf_gw_32,res_logis_pdf_gw_42,res_logis_pdf_gw_52,
        names = c("1.4 (20)", "1.2 (23)", "1.0 (28)", "0.8 (34)", "0.6 (46)"),ylab="L2 Distance",
        col = colors4,xlab="Grid Width (k)",ylim=c(0,0.4))
dev.off()

# Student's t Distribution - Grid Width Simulation 

n = 1000
B = 100

res_t_pdf_gw_11 = rep(0,B)
res_t_pdf_gw_21 = rep(0,B)
res_t_pdf_gw_31 = rep(0,B)
res_t_pdf_gw_41 = rep(0,B)
res_t_pdf_gw_51 = rep(0,B)

res_t_pdf_gw_12 = rep(0,B)
res_t_pdf_gw_22 = rep(0,B)
res_t_pdf_gw_32 = rep(0,B)
res_t_pdf_gw_42 = rep(0,B)
res_t_pdf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # t Distribution - GW1
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, by = 5)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_gw_11[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_gw_12[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - GW2
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, by = 4)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_gw_21[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_gw_22[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - GW3
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, by = 3)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_gw_31[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_gw_32[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - GW4
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, by=2)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_gw_41[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_gw_42[i] = L2_Distance(y,pdf_t,range_t)
  
  # t Distribution - GW5
  
  x = sort(rt(n,df = df))
  grid  <- seq(-20,20, by=1)
  delta <- max(diff(grid))
  y <- floor(x/delta)
  y = y*delta
  
  res_t_pdf_gw_51[i] = L1_Distance(y,pdf_t,range_t)
  res_t_pdf_gw_52[i] = L2_Distance(y,pdf_t,range_t)
  
}

grDevices::pdf("L1-PDF-1-t-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11,res_t_pdf_gw_21,res_t_pdf_gw_31,res_t_pdf_gw_41,res_t_pdf_gw_51,
        names = c("5 (17)", "4 (21)", "3 (27)", "2 (41)", "1 (81)"),ylab="L1 Distance",
        col = colors5,xlab="Grid Width",ylim=c(0,1.5))
dev.off()

grDevices::pdf("L2-PDF-1-t-gw.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12,res_t_pdf_gw_22,res_t_pdf_gw_32,res_t_pdf_gw_42,res_t_pdf_gw_52,
        names = c("5 (17)", "4 (21)", "3 (27)", "2 (41)", "1 (81)"),ylab="L2 Distance",
        col = colors5,xlab="Grid Width",ylim=c(0,0.7))
dev.off()




