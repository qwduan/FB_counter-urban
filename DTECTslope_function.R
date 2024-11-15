
#### These two functions aim to detect trends and their significance using robust linear regression
#### The 'DTECTslope_ro' function is used for countries where a breakpoint has been detected.
#### The 'DTECTslope_0_ro' function is used for countries where no breakpoint has been detected.

DTECTslope_ro <- function(series, bp) {

  Yt <- series$ts_data
  Vt <- series$deseason
  ti <- series$ts_time
  number <- series$number
  data <- data.frame(Yt = Yt, Vt = Vt, ti = ti, number = number)
  # Initialise variables to store the results
  s1 <- NA
  p1 <- NA
  s10 <-NA
  s2 <- NA
  p2 <- NA
  s20 <- NA
  sall <- NA
  pall <-NA
  sa0 <-NA
  cil_1 <- NA
  cil_2 <-NA
  cil_all <- NA
  # Separate the data into two segments based on the breakpoint
  data_be <- data %>%
    filter(number <= bp) %>%
    arrange(number)
  data_af <- data %>%
    filter(number > bp) %>%
    arrange(number)
  
  # Fit linear models to each segment
  model1 <- rlm(Vt~ ti, data = data_be, maxit = 100)
  model2 <- rlm(Vt~ ti, data = data_af, maxit = 100)
  modelall <- rlm(Vt ~ ti, data = data, maxit = 100)
  # Summarise the models to get slopes and p-values
  sum1 <- summary(model1)
  sum2 <- summary(model2)
  sumall <- summary(modelall)
  
  cil_before <- predict(model1, newdata = data.frame(ti = data_be$ti), interval = "confidence")
  cil_after <- predict(model2, newdata = data.frame(ti = data_af$ti), interval = "confidence")
  cil_all <- predict(modelall, newdata = data.frame(ti = data$ti),interval = "confidence")
  # Extract slopes and p-values from summaries
  s1 <- sum1$coefficients["ti", "Value"] 
  
  s10<- sum1$coefficients["(Intercept)", "Value"]
  s2 <- sum2$coefficients["ti", "Value"] 
  
  s20 <- sum2$coefficients["(Intercept)", "Value"]
  sall <- sumall$coefficients["ti", "Value"] 
  
  sa0 <- sumall$coefficients["(Intercept)", "Value"]
  
  
  # Optional: Calculate p-values assuming normal distribution (not recommended without careful consideration)
  p1 <- f.robftest(model1, var = "ti")[["p.value"]]
  p2 <- f.robftest(model2, var = "ti")[["p.value"]]
  pall <- f.robftest(modelall, var = "ti")[["p.value"]]
  
  cil_1 <-  cil_before
  cil_2 <- cil_after
  cil_all <-  cil_all
  sd1 <- sum1$coefficients["ti", "Std. Error"] 
  sd2 <- sum2$coefficients["ti", "Std. Error"] 
  sdall <- sumall$coefficients["ti", "Std. Error"] 
  # Return a list with the results
  return(list(
    s1 = s1,
    p1 = p1,
    s10 = s10,
    s2 = s2,
    p2 = p2,
    s20 = s20,
    sall = sall,
    sa0 = sa0,
    pall = pall,
    cil_1 = cil_1,
    cil_2 = cil_2,
    cil_all = cil_all,
    sd1 = sd1,
    sd2 = sd2,
    sdall = sdall
  ))
}

DTECTslope_0_ro <- function(series) {
  # Run bfast on the provided time series
  
  Vt <- series$deseason
  ti <- series$ts_time
  data <- data.frame(Vt = Vt, ti = ti, number = series$number)
  # Initialize variables to store the results
  s1 <- NA
  p1 <- NA
  s10 <-NA
  s2 <- NA
  p2 <- NA
  s20 <- NA
  sall <- NA
  pall <-NA
  sa0 <-NA
  cil_1 <- NA
  cil_2 <-NA
  cil_all <- NA
  
  modelall <- rlm(Vt ~ ti, data = data, maxit = 50)
  # Summarize the models to get slopes and p-values
  sumall <- summary(modelall)
  # Extract slopes and p-values from summaries
  sall <- sumall$coefficients["ti", "Value"] 
  
  sa0 <- sumall$coefficients["(Intercept)", "Value"]
  sdall <- sumall$coefficients["ti", "Std. Error"] 
  
  
  pall <- f.robftest(modelall, var = "ti")[["p.value"]]
  # Return a list with the results
  
  cil_all <- predict(modelall, interval = "confidence")
  
  return(list(
    s1 = NA,
    p1 = NA,
    s10 = NA,
    s2 = NA,
    p2 = NA,
    s20 = NA,
    sall = sall,
    sa0 = sa0,
    pall = pall,
    cil_1 = NA,
    cil_2 = NA,
    cil_all = cil_all,
    sd1 = sdall,
    sd2 = sdall,
    sdall = sdall
    
  ))
}
