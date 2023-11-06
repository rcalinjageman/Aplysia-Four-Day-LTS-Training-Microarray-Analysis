#' Estimate a proportion
#' 
#' \code{estimateProportion} estimates a proportion along with a confidence interval for a categorical variable.
#' The function can accept a column of categorical outcomes (a factor).  It will calculate the proportion of 
#' outcomes of a specified level (cases) out of all valid outcomes (n): P = cases/N.  It will then calculate
#' the confidence interval on that proportion.
#' 
#' This function can also accept summary data, where the user pases the number of cases and sample size directly.  
#' 
#' @param cases The number of observations of the specific outcome
#' @param n The total number of observations
#' @param caselabels Optional argument to define labels for outcomes.  Default is "Affected" and "Not Affected".
#' @param conf.level Confidence level for the confidence interval.  Default is 0.95.  Must be >0.5 and <1.
#' 
#' @examples
#' estimateProportion(cases = 10, n = 100)
#' estimateProportion(cases = 10, n = 100, conf.level = 0.90)
#' estimateProportion(cases = 10, n = 100, caselabels = c("Employed", "Unemployed"))


estimateProportion = function(cases, n, caselabels = c("Affected", "Not Affected"), conf.level=0.95) {
  #Function to calculate CI on single proportion from summary data only
  
  # Validate inputs ---------------------------
  # Check case labels
  if (is.null(caselabels)) {
    err_string <- stringr::str_interp(
      "`caselabels` must either be omitted or not null"
    )
    stop(err_string)
  } else {
    if(length(caselabels) < 2) {
      err_string <- stringr::str_interp(
        "`caselabels` must be a vector with two elements.  Currently class is ${class(caselabels)} of length ${length(caselabels)}"
      )
      stop(err_string)
    } else {
      if(caselabels[1] == caselabels[2]) {
        err_string <- stringr::str_interp(
          "caselabels must have 2 distinct elements.  Currently element 1 is ${caselabels[1]} and element 2 is also ${caselabels[2]}."
        )
        caselabels <- c("Affected", "Not Affected")
        warning(err_string)
      }
    }
  }
  
  # Check conf.level is >=.50 and < 1.  Could allow down to near 0, but why bother?
  if (conf.level < 0.50 | conf.level >= 1) {
    err_string <- stringr::str_interp(
      "`conf.level` must be between 0.50 and 1, not ${conf.level}"
    )
    stop(err_string)
  }
  
  # Check cases
  if(is.null(cases)) {
    err_string <- stringr::str_interp(
      "`cases` must be an integer number, not NULL"
    )
    stop(err_string)
  }
  if(cases != as.integer(cases)) {
    err_string <- stringr::str_interp(
      "`cases` must be an integer number, not ${cases}"
    )
    stop(err_string)
  }
  
  #check N
  if(is.null(n)) {
    err_string <- stringr::str_interp(
      "`n` must be an integer number, not NULL"
    )
    stop(err_string)
  }
  if(n != as.integer(n)) {
    err_string <- stringr::str_interp(
      "`n` must be an integer number, not ${n}"
    )
    stop(err_string)
  }
  
  #Check that cases<N
  if(cases>n) {
    err_string <- stringr::str_interp(
      "`cases` must be less than n.  Currently cases is ${cases} and n is ${n}."
    )
    stop(err_string)
  }
  
  
  P= cases/n
  q = 1-P
  z = qnorm(1-(1-conf.level)/2)
  
  A = 2*cases+z^2
  B = z*sqrt(z^2+(4*cases*q))
  C = 2*(n+z^2)
  
  ci.low = (A-B) / C
  ci.high = (A+B) / C
  
  summary_data <- data.frame(
    cases = cases,
    unaffected = n-cases,
    total = n,
    P = P,
    P.ci.low = ci.low,
    P.ci.high = ci.high
  )
  
  names(summary_data)[1] <- caselabels[1]
  names(summary_data)[2] <- caselabels[2]
  names(summary_data)[4] <- paste("P.", caselabels[1], sep="")
  
  res = list(
    summary_data = summary_data,
    P = P, 
    ci.low = ci.low, 
    ci.high = ci.high)
  class(res) <- "estimateProportions"
  return(res)
}


estimateProportionDifference.numeric <- function(cases1, n1, cases2, n2, caselabels = c("Affected", "Not Affected"), grouplabels = c("Group 1", "Group 2"), conf.level=0.95) {
  #Function to compare proportions from two groups using summary data only
  
  # Validate inputs ---------------------------
  # Check labels
  if (is.null(grouplabels)) {
    err_string <- stringr::str_interp(
      "`grouplabels` must either be omitted or not null"
    )
    stop(err_string)
  } else {
    if(length(grouplabels) < 2) {
      err_string <- stringr::str_interp(
        "`grouplabels` must be a vector with two elements.  Currently class is ${class(grouplabels)} of length ${length(grouplabels)}"
      )
      stop(err_string)
    } else {
      if(grouplabels[1] == grouplabels[2]) {
        err_string <- stringr::str_interp(
          "grouplabels must have elemnts 1 and 2 be distinct.  Currently element 1 is ${grouplabels[1]} and element 2 is also ${grouplabels[2]}."
        )
        grouplabels <- c("Group 1", "Group 2")
        warning(err_string)
      }
    }
  }
  
  # Check case labels
  if (is.null(caselabels)) {
    err_string <- stringr::str_interp(
      "`caselabels` must either be omitted or not null"
    )
    stop(err_string)
  } else {
    if(length(caselabels) < 2) {
      err_string <- stringr::str_interp(
        "`caselabels` must be a vector with two elements.  Currently class is ${class(caselabels)} of length ${length(caselabels)}"
      )
      stop(err_string)
    } else {
      if(caselabels[1] == caselabels[2]) {
        err_string <- stringr::str_interp(
          "caselabels must have elemtns 1 and 2 be distinct.  Currently element 1 is ${caselabels[1]} and element 2 is also ${caselabels[2]}."
        )
        caselabels <- c("Affected", "Not Affected")
        warning(err_string)
      }
    }
  }
  
  # Check conf.level is >=.50 and < 1.  Could allow down to near 0, but why bother?
  if (conf.level < 0.50 | conf.level >= 1) {
    err_string <- stringr::str_interp(
      "`conf.level` must be between 0.50 and 1, not ${conf.level}"
    )
    stop(err_string)
  }

  
  # Calculations ---------------------------
  P1 = cases1/n1
  P2 = cases2/n2
  q1 = 1-P1
  q2 = 1-P2
  z = 1.96
  
  P1_res = estimateProportion(cases1, n1, conf.level = conf.level)
  P2_res = estimateProportion(cases2, n2, conf.level = conf.level)
  
  P1.ci.low = P1_res$ci.low
  P1.ci.high = P1_res$ci.high
  
  P2.ci.low = P2_res$ci.low
  P2.ci.high = P2_res$ci.high
  
  
  Pdiff = P1 - P2
  
  ci.low = Pdiff - sqrt( (P1-P1.ci.low)^2 + (P2.ci.high - P2)^2  )
  ci.high = Pdiff + sqrt( (P2 - P2.ci.low)^2 + (P1.ci.high - P1)^2 )
  
  
  # Prep Ouput ---------------------------
  summary_data <- rbind(P1_res$summary_data, P2_res$summary_data)
  summary_data[nrow(summary_data) + 1,] = list(NA, NA, NA, Pdiff, ci.low, ci.high)
  summary_data <- cbind(data.frame(Condition = c(grouplabels, "Difference")), summary_data)
  
  names(summary_data)[2] <- caselabels[1]
  names(summary_data)[3] <- caselabels[2]
  names(summary_data)[5] <- paste("P.", caselabels[1], sep="")
  
  
  # Gather output in a list
  out <- list(
    conf.level = conf.level,
    summary_data = summary_data,
    level1 = grouplabels[1],
    level2 = grouplabels[2],
    cases1 = cases1,
    cases2 = cases2,
    n1 = n1,
    n2 = n2,
    P1 = P1,
    P1.ci.low = P1.ci.low,
    P1.ci.high = P1.ci.high,
    P2 = P2,
    P2.ci.low = P2.ci.low,
    P2.ci.high = P2.ci.high,
    Pdiff = Pdiff,
    ci.low = ci.low,
    ci.high = ci.high,
    type = "fromSummary"
  )
  
  # Set class as estimate for printing purposes
  class(out) <- "estimateProportions"
  return(out)
}