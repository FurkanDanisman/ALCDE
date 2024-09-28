library(logcondiscr)
library(logcondens)
library(fdrtool)
library(truncnorm)
library(pracma)

# Plotting the PMF of the Discerete Log-Concave Function 

log_Con_Discr_MLE_plot = function(y,uniform_gw = T,a_i = NA,b_i=NA,k=0){ 
  
  if (uniform_gw == T) {
    
    delta = min(diff(unique(y))) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
    # (Assuming a uniform grid, delta remains the same for any i)
    
    y = floor(y / delta) # Turning the data into integers 
    
    mle <- logConDiscrMLE(y,output = F) # Discrete Log-Concave Function 
    
    psi <- mle$psi # log(p1), log(p2), ...
    
    plot(delta*mle$x, exp(psi), type = "b", col = 2, xlim = c(min(delta*mle$x),max(delta*mle$x)),
         xlab = "x", ylim = c(0, max(exp(psi)) + k), ylab = "PMF", pch = 19)
    
    # delta * mle$x means, we are returning the integers back into the original form. 
    
  }else{
    
    Result_NU = Log_Con_Discr_Algorithm_NU(y,a_i,b_i)
    
    plot(Result_NU$mle_x, exp(Result_NU$mle_psi), type = "b", col = 2, xlim = c(min(Result_NU$mle_x),max(Result_NU$mle_x)),
         xlab = "x", ylim = c(0, max(exp(Result_NU$mle_psi))), ylab = "PMF",
         main = "PMF of Discerete Log-Concave Density", pch = 19)
    
  }
  
  
}

# Log Concave Discrete Algorithm 

Log_Con_Discr_Algorithm = function(y,intercept = F,slope = F,bound_right_to_left=FALSE,bound_left_to_right=FALSE){
  
  delta = min(diff(unique(y))) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
                               # (Assuming a uniform grid, delta remains the same for any i)
  
  mle_y = y/delta # Turning y variables into integers (Way around for logConDiscMLE function to work)
    
  mle  <- logConDiscrMLE(mle_y,output = F) # Discrete Log-Concave Function 
  
  mle_x = mle$x * delta  # Turning y variables back into the original form
  
  mle_psi = mle$psi # log(p1),log(p2),log(p3),...
  
  lcm = gcmlcm(x = mle_x, y = mle$psi,type = "lcm") # Least-Concave Majorant Function 
  
  knots_level = length(lcm$slope.knots) # Taking the length of the knots (dots)
  
  int = rep(0,knots_level) # Intercept
  
  integ = 0
  
  # Integration of the LCM (Recovery Function)
  
  for (i in 1:knots_level) {
    
    int[i] = lcm$y.knots[i] - lcm$x.knots[i]*lcm$slope.knots[i] # Forming the intercept
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) # Function between each two knots
    
    integ = integral(fx,lcm$x.knots[i],lcm$x.knots[i+1]) + integ # Taking the integral of the whole range
    
    coef = 1 / integ # The coefficient that needs to be multiplied by the function so that it integrates into 1.
    
  }
  
  # Taking the Expectation of the LCM (Recovery Function)
  
  LCM_mu = 0
  
  for (i in 1:knots_level) {
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) * coef * x # E(x) = Integral( xf(x) dx )
    
    LCM_mu = integral(fx,lcm$x.knots[i],lcm$x.knots[i+1]) + LCM_mu 
    
  }
  
  EM_Result = EM_Log_Con_Discr_algorithm(y) # EM Algorithm 
  
  if("EM Algorithm did not converge" %in% EM_Result){
    
    mid_point = (mean(y/delta)+0.5)*delta
    delta = 0 # Required Shift Distance 
    
  }else{
    
    mu = EM_Result$`Estimated Mu` # Estimated Mu from EM Algorithm
    
    mid_point = EM_Result$`Mid-point` # Mid-point Assumption 
    
    delta = mu - LCM_mu # Required Shift Distance 
    
  }
  
  x_knots = lcm$x.knots + delta # Shifting the Points 
  
  # If the Function Requires to be bounded to a point. (For example, let minimum point to be 0)
  
  if (is.numeric(bound_right_to_left)) { # If Bound Right to Left = 0 and min(LCM) = 0.2 becomes min(LCM) = 0
    
    x_knots[1] = bound_right_to_left
    
  }
  if (is.numeric(bound_left_to_right)) { # If Bound Left to Right = 0 and min(LCM) = -0.2 becomes min(LCM) = 0
    
    x_knots[x_knots<bound_left_to_right] = bound_left_to_right
    
  }
  
  # New Coefficient Calculation After Shifting the x-points (Integration changes when we move the x's)
  
  integ = 0
  
  for (i in 1:knots_level) {
    
    int[i] = lcm$y.knots[i] - lcm$x.knots[i]*lcm$slope.knots[i] # Intercepts
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) * coef # Recovery Function 
    
    integ = integral(fx,x_knots[i],x_knots[i+1]) + integ 
    
    coef_2 = 1 / integ # The coefficient that needs to be multiplied by the function so that it integrates into 1.
    
  }
  
  Coef_total = coef_2 * coef # New Coefficient
  
  # Output of the Function 
  
  if (intercept == T & slope == T) {
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Intercepts"=int,"Slope"=lcm$slope.knots))
    
  }else if(intercept == T){
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Intercepts"=int))
    
  }else if(slope == T){
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Slope"=lcm$slope.knots))
    
  }else{
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi))
    
  }
  
}


EM_Log_Con_Discr_algorithm <- function(y,round_to=2,tol=1e-10, max_iter = 3000){
  
  # Preparation of the inputs
  
  EM_y = round(y,round_to) 
  EM_y = unique(EM_y) 
  
  delta = min(diff(EM_y))   # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
                            # (Assuming a uniform grid, delta remains the same for any i)
  
  a_i  = EM_y # Lower Bound of the Intervals 
  a_i  = round(seq(min(EM_y),max(EM_y),by = delta),round_to)
  b_i  = a_i[-1] 
  b_i  = c(b_i,max(b_i)+delta) # Upper Bound of the Intervals 
  n_i  = hist(y,breaks = c(min(a_i)-delta,a_i),plot=FALSE)$counts # Frequency of the Intervals
  
  mu_0 = (mean(y/delta)+0.5)*delta # Mid-point as initial mu
  mid_point = mu_0
  
  sd = sd(y)
  
  sum_E = 0
  iter  = 1
  diff = 1
  
  # Algorithm Procedure 
  
  while (diff > tol & iter < max_iter) {
    
    for (i in 1:length(n_i)) {
      
      sum_E = n_i[i] * etruncnorm( a = a_i[i], b = b_i[i], mean = mu_0, sd = sd) + sum_E 
      
    }
    
    mu_hat = sum_E / sum(n_i)
    
    diff = abs(mu_hat - mu_0)
    
    mu_0 = mu_hat
    
    iter = iter + 1
    
    sum_E = 0
    
  }
  
  if (iter == max_iter & diff > tol) {
    return("EM Algorithm did not converge")
    
  }
  return(list("Estimated Mu" = mu_hat,"Iteration" = iter,"Mid-point"=mid_point))
  
}


Log_Con_Discr_Density_Plot = function(y,uniform_gw = T,a_i=NA,b_i=NA, Discrete_LCD_Points = T,LCM_Line = T, Density_LCM = T, Adjusted_Density_LCM = T,
                                      bound_right_to_left = F,bound_left_to_right = F,Bound_Plot = T,legend = T){
  
  Delta_x = min(diff(unique(y))) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
                                 # (Assuming a uniform grid, delta remains the same for any i)
  
  if(uniform_gw == T){
    
    Result_Log_Con_Discr_Alg = Log_Con_Discr_Algorithm(y)
    
  }else{
    
    Result_Log_Con_Discr_Alg = Log_Con_Discr_Algorithm_NU(y,a_i,b_i)
    
    }
  
  mle_psi = Result_Log_Con_Discr_Alg$mle_psi # log(p1), log(p2), ...
  
  # Least Concave Majorant 
  
  lcm = gcmlcm(x = Result_Log_Con_Discr_Alg$mle_x,  y = mle_psi,type = "lcm")
  
  # Initial Plot 
  
  plot(x = 0,y = -1,
       ylim = c(0,max(exp(mle_psi))*Result_Log_Con_Discr_Alg$Total_Coefficient + 0.05),col="white"
       ,xlab="x",ylab = "Probability",bty="l",
       xlim=c(Result_Log_Con_Discr_Alg$x_knots[1]-Delta_x,Result_Log_Con_Discr_Alg$x_knots[length(Result_Log_Con_Discr_Alg$x_knots)]+Delta_x))
  
  # Initial Inputs for Legend
  
  legend_names = NULL
  legend_col = NULL
  legend_pch = NULL
  legend_lty = NULL
  legend_lwd = NULL
  
  # Plotting the Discrete LCD Knots 
  
  if (Discrete_LCD_Points == T) {
    
    plot(x = Result_Log_Con_Discr_Alg$mle_x,y = exp(mle_psi),
         ylim = c(0,max(exp(mle_psi))*Result_Log_Con_Discr_Alg$Total_Coefficient + 0.05),col="black",
         xlab = "x", ylab = "Probability",bty="l",lwd = 2)
    
    legend_names = c(legend_names,"Discrete-LCD")
    legend_col = c(legend_col,"black")
    legend_pch = c(legend_pch,1)
    legend_lty = c(legend_lty,NA)
    legend_lwd = c(legend_lwd,2)
    
  }
  
  # Plotting the Least Concave Majorant Line 
  
  if (LCM_Line == T) {
    
    lines(x = lcm$x.knots,y = exp(lcm$y.knots),col="navy",lwd = 3)
    legend_names = c(legend_names,"LCM")
    legend_col = c(legend_col,"navy")
    legend_pch = c(legend_pch,NA)
    legend_lty = c(legend_lty,1)
    legend_lwd = c(legend_lwd,3)
    
  }
  
  # Plotting the Density of the Least Concave Majorant 
  
  if (Density_LCM == T) {
    
    lines(x = lcm$x.knots,y = exp(lcm$y.knots) * Result_Log_Con_Discr_Alg$First_Coefficient ,col="purple",lwd =3)
    legend_names = c(legend_names,"Density-LCM")
    legend_col = c(legend_col,"purple")
    legend_pch = c(legend_pch,NA)
    legend_lty = c(legend_lty,1)
    legend_lwd = c(legend_lwd,3)
    
  }
  
  # Plotting the Adjusted Density of the Least Concave Majorant (Mean point shifted)
  
  if (Adjusted_Density_LCM == T) {
    
    lines(x = Result_Log_Con_Discr_Alg$x_knots,y = exp(lcm$y.knots) * Result_Log_Con_Discr_Alg$Total_Coefficient,col="brown3", lwd = 3)
    legend_names = c(legend_names,"Adjusted-Density-LCM")
    legend_col = c(legend_col,"brown3")
    legend_pch = c(legend_pch,NA)
    legend_lty = c(legend_lty,1)
    legend_lwd = c(legend_lwd,3)
    
  }
  
  # Plotting the Density of the Adjusted Least Concave Majorant Where Minimum Point is Bounded
  
  if (is.numeric(bound_right_to_left) & Bound_Plot == T) { # If Bound Right to Left = 0 and min(LCM) = 0.2 becomes min(LCM) = 0
    
    Result_Log_Con_Discr_Alg_Bounded = Log_Con_Discr_Algorithm(y,bound_right_to_left = bound_right_to_left)
    lcm = gcmlcm(x = Result_Log_Con_Discr_Alg_Bounded$mle_x,  y = mle_psi,type = "lcm")
    lines(x = Result_Log_Con_Discr_Alg_Bounded$x_knots,y = exp(lcm$y.knots) * Result_Log_Con_Discr_Alg_Bounded$Total_Coefficient,col="navy",lwd = 2)
    legend_names = c(legend_names,paste("Bounded-to",bound_right_to_left))
    legend_col = c(legend_col,"navy")
    legend_pch = c(legend_pch,NA)
    legend_lty = c(legend_lty,1)
    legend_lwd = c(legend_lwd,3)
      
  }else if (is.numeric(bound_left_to_right) & Bound_Plot == T) { # If Bound Left to Right = 0 and min(LCM) = -0.2 becomes min(LCM) = 0
    
    Result_Log_Con_Discr_Alg_Bounded = Log_Con_Discr_Algorithm(y,bound_left_to_right  = bound_left_to_right)
    lcm = gcmlcm(x = Result_Log_Con_Discr_Alg_Bounded$mle_x,  y = mle_psi,type = "lcm")
    lines(x = Result_Log_Con_Discr_Alg_Bounded$x_knots,y = exp(lcm$y.knots) * Result_Log_Con_Discr_Alg_Bounded$Total_Coefficient,col="navy", lwd = 2)
    legend_names = c(legend_names,paste("Bounded-to",bound_left_to_right))
    legend_col = c(legend_col,"navy")
    legend_pch = c(legend_pch,NA)
    legend_lty = c(legend_lty,1)
    legend_lwd = c(legend_lwd,3)
    
  }
  
  # Adding Legend to Plot
  
  if (legend==T) {
    legend("topright", 
           legend = legend_names, 
           col = legend_col, 
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           bty = "n", 
           pt.cex = 2, 
           cex = 2, 
           text.col = "black", 
           horiz = F, 
           inset = c(-0.35, 0),x.intersp = 0.2, y.intersp = 0.5,seg.len = 0.75) # -0.3
  }
  
}



# L1 Distance Function 

L1_Distance = function(y,pdf,range){
  
  # Running the LCD Algorithm 
  
  lcm_results = Log_Con_Discr_Algorithm(y,intercept = T,slope = T) 
  
  Total_coef = lcm_results$Total_Coefficient
  
  x_knots = lcm_results$x_knots
  
  # L1 Algorithm 
  
  L1 = 0 
  
  for (i in 1:(length(x_knots)-1)) {
    
    d_lcm = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
    fx_gx = function(x) abs(  pdf(x) * pdf_coef  - d_lcm(x) )
    L1 = integral(fx_gx,x_knots[i],x_knots[i+1]) + L1
      
  }
  
  L1 = integral(pdf,range[1],x_knots[1]) + L1
  L1 = integral(pdf,x_knots[length(x_knots)],range[2]) + L1
  
  return(L1)
  
}

# L2 Distance Function 

L2_Distance = function(y,pdf,range){
  
  # Running the LCD Algorithm 
  
  lcm_results = Log_Con_Discr_Algorithm(y,intercept = T,slope = T)
  
  Total_coef = lcm_results$Total_Coefficient
  
  x_knots = lcm_results$x_knots
  
  # L2 Algorithm 
  
  L2 = 0 
  
  for (i in 1:(length(x_knots)-1)) {
    
    d_lcm = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
    fx_gx = function(x) ( pdf(x)*pdf_coef - d_lcm(x) )^2
    L2 = integral(fx_gx,x_knots[i],x_knots[i+1]) + L2
    
  }
  
  L2 = integral(pdf,range[1],x_knots[1])^2 + L2
  L2 = integral(pdf,x_knots[length(x_knots)],range[2])^2 + L2

  return(sqrt(L2))
  
}


EM_L1_L2 = function(pdf,y,range){
  
  # EM Algorithm 
  
  EM_format = EM_Algorithm_Format(y)
  
  EM_results = em(bl = EM_format$bl,bu = EM_format$bu
                  ,freq = EM_format$Freq,theta_init = EM_format$Theta)
  
  # EM Algorithm Density Estimation  
  
  pdf_em <- function(x) {
    dnorm(x, mean = EM_results$mu_estimate, sd = EM_results$sigma_estimate)
  }
  
  # L1 Distance 
  
  fx_gx = function(x) abs(  pdf(x) - pdf_em(x) )
  L1 = integral(fx_gx,range[1],range[2])
  
  # L2 Distance 
  
  fx_gx = function(x) (  pdf(x) - pdf_em(x) )^2
  L2 = sqrt(integral(fx_gx,range[1],range[2]))
  
  return(list("L1_Distance" = L1,"L2_Distance" = L2))
  
}


L1_L2_pmf = function(y,cdf){ 
  
  delta = min(diff(unique(y))) # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
                               # (Assuming a uniform grid, delta remains the same for any i)
  
  grid <- seq(min(y), max(y), by = delta)
  grid = c(grid,max(grid)+delta)
  
  p0 <- diff(cdf(grid)) # Extracting the original probability of the interval by taking the differences of the cut points.
  
  grid <- seq(min(y), max(y), by = delta)
  
  y = floor(y / delta) # Turning the data into integers 
  
  mle <- logConDiscrMLE(y,output = F) # Discrete Log-Concave Function 
  
  psi <- mle$psi # log(p1), log(p2), ...
  
  exp_psi = exp(psi) # p1, p2, ... 
  
  x_points = delta*mle$x # delta * mle$x means, we are returning the integers back into the original form. 
  
  L1 = 0; L2 = 0;
  
  for (i in 1:length(x_points)) {
    
    if (length(which(grid == x_points[i])) == 0) {
      
      # L1 Distance 
      
      L1 = abs(  p0[ which(round(grid,8) == round(x_points[i],8)) ] - exp_psi[i] ) + L1
      
      # L2 Distance 
      
      L2 = ( p0[ which(round(grid,8) == round(x_points[i],8)) ] - exp_psi[i] )^2 + L2
      
      next
    }
    
    # L1 Distance 
    
    L1 = abs(  p0[ which(grid == x_points[i]) ]  - exp_psi[i] ) + L1
    
    # L2 Distance 
    
    L2 = ( p0[ which(grid == x_points[i]) ] - exp_psi[i] )^2 + L2
    
  }
  
  L1 = L1 + 1 - sum(p0)
  L2 = L2 + (1-sum(p0))^2
  L2 = sqrt(L2)
  
  return(list("L1_pmf" = L1,"L2_pmf" = L2))
  
}

# Density Function 

d_lcm = function(x,data){
  
  lcm_results = Log_Con_Discr_Algorithm(data,intercept = T,slope = T)
  
  Total_coef = lcm_results$Total_Coefficient
  
  x_knots = lcm_results$x_knots
  
  result = 0 
  
  for (i in 1:(length(x_knots)-1)) {
    
    if(i == 1 & x < x_knots[i]){
      result = 0
      
    }else if(x>=x_knots[i] & x<=x_knots[i+1]){
      
      result = exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
      
      break
      
    }else if(i == (length(x_knots)-1) & x > x_knots[i+1]){
      
      result = 0
      
    }
  }
  
  return(result)
  
}


# Cdf Function

p_lcm = function(q,data){
  
  lcm_results = Log_Con_Discr_Algorithm(data,intercept = T,slope = T)
  
  Total_coef = lcm_results$Total_Coefficient
  
  x_knots = lcm_results$x_knots
  
  cdf = 0 
  
  for (i in 1:(length(x_knots)-1)) {
    
    if(i == 1 & q < x_knots[i]){
      
      cdf = 0
      
      break
      
    }else if(i == (length(x_knots)-1) & q > x_knots[i+1]){
      
      cdf = 1
      
      break
      
    }else if(q>=x_knots[i] & q<=x_knots[i+1]){
      
      fx = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
      
      cdf = integral(fx,x_knots[i],q) + cdf
      
      break
      
    }else if(q>=x_knots[i] & q >= x_knots[i+1]){
      
      fx = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
      
      cdf = integral(fx,x_knots[i],x_knots[i+1]) + cdf
      
    }
  }
  
  return(cdf)
  
}


# Inverse CDF Function 

q_lcm = function(p,data){
  
  lcm_results = Log_Con_Discr_Algorithm(data,intercept = T,slope = T)
  
  Total_coef = lcm_results$Total_Coefficient
  
  x_knots = lcm_results$x_knots
  
  pdf = 0 
  
  for (i in 1:(length(x_knots)-1)) {
    
    fx = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
    
    pdf[i] = integral(fx,x_knots[i],x_knots[i+1])
    
  }
  
  cdf = cumsum(pdf)
  
  stop_f = 0
  
  q = 0
  
  if (p < 0) {
    
    stop(paste(p,"is not a probability value"))
    
  }else if(p == 1){
    
    q = paste("Any real number greater than or equal to",round(x_knots[length(x_knots)],6))
    
    stop_f = 1 
    
  }else if(p == 0){
    
    q = paste("Any real number less than",round(x_knots[1],6))
    
    stop_f = 1
    
  }else if(p > 1){
    
    stop(paste(p,"is not a probability value"))
    
  }
  
  for (i in 1:(length(x_knots)-1)) {
    
    if (stop_f == 1) {
      
      break
      
    }
    
    if(i == 1 & p < cdf[i]){
      
      fx = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
      
      qx = integral(fx,x_knots[i],q)$value - p 
      
      q = uniroot(qx, lower = x_knots[i], upper = x_knots[i+1])$root
      
      break
      
    }else if(p==cdf[i]){
      
      q = cdf[i]
      
      break
      
    }else if(p<cdf[i]){
      
      fx = function(x) exp(lcm_results$Intercepts[i] + lcm_results$Slope[i]*x)*lcm_results$Total_Coefficient
      
      qx = function(k) integral(fx,x_knots[i],k) - p + cdf[i-1]
      
      q = uniroot(qx, lower = x_knots[i], upper = x_knots[i+1])$root
      
      break
      
    }
  }
  
  return(q)
  
}


darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


mest<- function(theta,bl,bu,freq){
  
  Aj<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))
  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
  }
  
  dinom<- NULL
  
  for(i in 1:base::length(bl)){
    dinom[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    if(dinom[i]==0) {
      Aj[i]<- theta[1]-theta[2]*((stats::dnorm(bstar[i])-
                                    stats::dnorm(astar[i]))/0.0001)
    }
    else {
      Aj[i]<- theta[1]-theta[2]*((stats::dnorm(bstar[i])-
                                    stats::dnorm(astar[i]))/dinom[i])
    }
  }
  
  M <- base::sum(Aj*freq)/base::sum(freq)
  return(M)
  
}


ssest<- function(theta,bl,bu,muupdate,freq){
  
  Bj<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))
  
  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
    
  }
  
  astar[1] <- (-1000)
  bstar[base::length(bl)]<- (1000)
  
  dinom<- NULL
  
  for(i in 1:base::length(bl)){
    
    dinom[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    
    if(dinom[i]==0) {
      Bj[i]<- theta[2]^2*(1-(bstar[i]*stats::dnorm(bstar[i])-
                               astar[i]*stats::dnorm(astar[i]))/0.0001)
      +(muupdate-theta[1])^2+
        (2*theta[2]*(muupdate-theta[1])*((stats::dnorm(bstar[i])-
                                            stats::dnorm(astar[i]))/0.0001))
    }
    else{
      Bj[i]<- theta[2]^2*(1-(bstar[i]*stats::dnorm(bstar[i])-
                               astar[i]*stats::dnorm(astar[i]))/dinom[i])
      +(muupdate-theta[1])^2+
        (2*theta[2]*(muupdate-theta[1])*((stats::dnorm(bstar[i])-
                                            stats::dnorm(astar[i]))/dinom[i]))}
    
  }
  
  ss<- base::sum(Bj*freq)/base::sum(freq)
  return(ss)
}

em <- function(bl,bu,freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){
  
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]
  
  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    munew<- mest(theta=cur,bl,bu,freq)
    ssnew<- ssest(theta=cur,bl,bu,
                               muupdate=mest(theta=cur,bl,bu,freq) ,freq)
    
    Mu_new <- munew
    S_new <- base::sqrt(ssnew)
    new_step<- c(Mu_new,S_new)
    
    if(base::abs(cur[1]-new_step[1])<tol1 &
       base::abs(cur[2]-new_step[2])<tol2){
      flag<-1 ;break
    }
    Mu_cur<- Mu_new
    S_cur<- S_new
  }
  
  if(!flag){
    warning("Didn't Converge \n")
  }
  
  updateres<- base::list("mu_estimate" = Mu_cur,
                         "sigma_estimate" = (S_cur))
  
  return(updateres)
}

EM_Algorithm_Format <- function(y,round_to=2){
  
  # Preparation of the inputs
  
  EM_y = round(y,round_to) 
  EM_y = unique(EM_y) 
  
  delta = min(diff(EM_y))   # Taking the interval gap if [a_i, b_i) => delta = b_i - a_i 
  # (Assuming a uniform grid, delta remains the same for any i)
  
  a_i  = EM_y # Lower Bound of the Intervals 
  a_i  = round(seq(min(EM_y),max(EM_y),by = delta),round_to)
  b_i  = a_i[-1] 
  b_i  = c(b_i,max(b_i)+delta) # Upper Bound of the Intervals 
  n_i  = hist(y,breaks = c(min(a_i)-delta,a_i),plot=FALSE)$counts # Frequency of the Intervals
  
  mu_0 = (mean(y/delta)+0.5)*delta # Mid-point as initial mu
  mid_point = mu_0
  
  sd = sd(y)
  
  return(list("Theta" = c(mu_0,sd),"bl" = a_i,"bu"=b_i,"Freq" = n_i))
  
}

EM_L1_L2_sigma_known = function(pdf,y,range,sd){
  
  # EM Algorithm 
  
  EM_format = EM_Algorithm_Format(y)
  EM_results = mest(bl = EM_format$bl,bu = EM_format$bu,freq = EM_format$Freq,theta = c(EM_format$Theta[1],sd))
  
  
  # EM Algorithm Density Estimation  
  
  pdf_em <- function(x) {
    dnorm(x, mean = EM_results, sd = sd)
  }
  
  # L1 Distance 
  
  fx_gx = function(x) abs(  pdf(x) - pdf_em(x) )
  L1 = integral(fx_gx,range[1],range[2])
  
  # L2 Distance 
  
  fx_gx = function(x) (  pdf(x) - pdf_em(x) )^2
  L2 = sqrt(integral(fx_gx,range[1],range[2]))
  
  return(list("L1_Distance" = L1,"L2_Distance" = L2))
  
}

Log_Con_Discr_Algorithm_NU = function(y,a_i, b_i,round_to = 2, intercept = F, slope = F, bound_right_to_left=FALSE, bound_left_to_right=FALSE){
  
  diff = round(b_i - a_i,round_to) # Grid Widths
  delta = c(diff[diff(diff) != 0],diff[length(diff)])
  
  # Calculate the index of the first non-changing number
  index_ai <- c(1,which(diff != dplyr::lag(diff)))
  index_bi <- c(which(diff != dplyr::lag(diff)) - 1,length(diff))
  
  y_i_u = list() # Upper bound for different grid widths
  y_i_l = list() # Lower bound for different grid widths
  y_x = list() # y-values for each unique grid widths
  y_x_mle = list() # To input into mle function
  mle = list() # Storing mles
  xx = list() # Storing y-values for different grid widths
  phat = list() # Storing p-values for different grid widths
  ratios = list() # Ratio of the length of the particular grid width over sample size
  pmle = list() # Storing pmle's for different grid widths
  mu_0 = 0 # Mid-point & Initial mean 
  
  for (i in 1:length(delta)) {
    
    y_i_l[[i]] = a_i[index_ai[i]]
    y_i_u[[i]] = b_i[index_bi[i]]
    y_x[[i]] = y[y_i_u[[i]] > y & y >= y_i_l[[i]]]
    y_x_mle[[i]] = round(y_x[[i]] / delta[i],2)
    mle[[i]] = logConDiscrMLE(y_x_mle[[i]],output = F)
    xx[[i]] = mle[[i]]$xSupp * delta[i]
    phat[[i]] = exp(mle[[i]]$psiSupp)
    ratios[[i]] = length(y_x[[i]]) / length(y)
    pmle[[i]] = phat[[i]] * ratios[[i]] / delta[i]
    mu_0 = (mean(y_x[[i]]/delta[i])+0.5)*delta[i] + mu_0
    
  }
  
  xx_all = sort(unlist(xx)) # Unlisting into a vector
  pmle_all = unlist(pmle) # Unlisting into a vector
  psi_all = log(pmle_all) # Psi values
  
  n_i  = hist(y,breaks = c(min(a_i)-delta[1],b_i),plot=FALSE)$counts # Frequency
  
  mid_point = mu_0 # Mid-point Assumption
  theta_init = c(mu_0,sd(y)) # Initial Points for EM-Algorithm Assuming Normality
  
  EM_result = em(bl = a_i,bu = b_i,freq = n_i,theta_init = theta_init) # EM-Algorithm 
  
  mu = EM_result$mu_estimate # Mu estimation 
  
  lcm = gcmlcm(x = xx_all, y = psi_all, type = "lcm") # Least-Concave Majorant Function 
  
  knots_level = length(lcm$slope.knots) # Taking the length of the knots (dots)
  
  int = rep(0,knots_level) # Intercept
  
  integ = 0
  
  # Integration of the LCM (Recovery Function)
  
  for (i in 1:knots_level) {
    
    int[i] = lcm$y.knots[i] - lcm$x.knots[i]*lcm$slope.knots[i] # Forming the intercept
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) # Function between each two knots
    
    integ = integral(fx,lcm$x.knots[i],lcm$x.knots[i+1]) + integ # Taking the integral of the whole range
    
    coef = 1 / integ # The coefficient that needs to be multiplied by the function so that it integrates into 1.
    
  }
  
  # Taking the Expectation of the LCM (Recovery Function)
  
  LCM_mu = 0
  
  for (i in 1:knots_level) {
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) * coef * x # E(x) = Integral( xf(x) dx )
    
    LCM_mu = integral(fx,lcm$x.knots[i],lcm$x.knots[i+1]) + LCM_mu 
    
  }
  
  delta = mu - LCM_mu # Required Shift Distance 
  
  x_knots = lcm$x.knots - delta # Shifting the Points 
  
  # If the Function Requires to be bounded to a point. (For example, let minimum point to be 0)
  
  if (is.numeric(bound_right_to_left)) { # If Bound Right to Left = 0 and min(LCM) = 0.2 becomes min(LCM) = 0
    
    x_knots[1] = bound_right_to_left
    
  }
  if (is.numeric(bound_left_to_right)) { # If Bound Left to Right = 0 and min(LCM) = -0.2 becomes min(LCM) = 0
    
    x_knots[x_knots<bound_left_to_right] = bound_left_to_right
    
  }
  
  # New Coefficient Calculation After Shifting the x-points (Integration changes when we move the x's)
  
  integ = 0
  
  for (i in 1:knots_level) {
    
    int[i] = lcm$y.knots[i] - lcm$x.knots[i]*lcm$slope.knots[i] # Intercepts
    
    fx = function(x) exp( lcm$slope.knots[i]*x + int[i] ) * coef # Recovery Function 
    
    integ = integral(fx,x_knots[i],x_knots[i+1]) + integ 
    
    coef_2 = 1 / integ # The coefficient that needs to be multiplied by the function so that it integrates into 1.
    
  }
  
  Coef_total = coef_2 * coef # New Coefficient
  
  # Output of the Function 
  
  if (intercept == T & slope == T) {
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Intercepts"=int,"Slope"=lcm$slope.knots))
    
  }else if(intercept == T){
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Intercepts"=int))
    
  }else if(slope == T){
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=mle_x,"mle_psi"=mle_psi,"Slope"=lcm$slope.knots))
    
  }else{
    
    return(list("Total_Coefficient"=Coef_total,"First_Coefficient" = coef,"Expectation" = LCM_mu,"EM_Mean"=mu,"Mid_Point" = mid_point,"Delta_mu" = delta,
                "x_knots" = x_knots,"mle_x"=xx_all,"mle_psi"=psi_all))
    
  }
  
}
