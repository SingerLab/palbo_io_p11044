#' http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
#' Gives count, mean, standard deviation, standard error of the mean, and confidence
#' interval (default 95%).
#'   data: a data frame.
#'   measurevar: the name of a column that contains the variable to be summariezed
#'   groupvars: a vector containing names of columns that contain grouping variables
#'   na.rm: a boolean that indicates whether to ignore NA's
#'   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          median = median(xx[[col]], na.rm=na.rm), ## added by RGM
          var = var     (xx[[col]], na.rm=na.rm), ## added by RGM
          sd   = sd     (xx[[col]], na.rm=na.rm),
          min  = min(xx[[col]], na.rm=na.rm),  ## added by RGM
          max  = max(xx[[col]], na.rm=na.rm),   ## added by RGM
          q1  =  as.numeric(quantile(xx[[col]], .5)),
          q3  =  as.numeric(quantile(xx[[col]], .95))
        )
      },
      measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

#' Norms the data within specified groups in a data frame; it normalizes each
#' subject (identified by idvar) so that they have the same mean, within each group
#' specified by betweenvars.
#'   data: a data frame.
#'   idvar: the name of a column that identifies each subject (or matched subjects)
#'   measurevar: the name of a column that contains the variable to be summariezed
#'   betweenvars: a vector containing names of columns that are between-subjects variables
#'   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}

#' Summarizes data, handling within-subjects variables by removing inter-subject variability.
#' It will still work if there are no within-S variables.
#' Gives count, un-normed mean, normed mean (with same between-group mean),
#'   standard deviation, standard error of the mean, and confidence interval.
#' If there are within-subject variables, calculate adjusted values using method from Morey (2008).
#'   data: a data frame.
#'   measurevar: the name of a column that contains the variable to be summariezed
#'   betweenvars: a vector containing names of columns that are between-subjects variables
#'   withinvars: a vector containing names of columns that are within-subjects variables
#'   idvar: the name of a column that identifies each subject (or matched subjects)
#'   na.rm: a boolean that indicates whether to ignore NA's
#'   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}


#' function for number of observations
#' n placed at y = -1
#'
#' @examples
#' geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
#' 
#' @references https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp

#' @export
give.n <- function(x) {
  return(c(y = -1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#' function for number of observations
#' n placed at y = 0
#' 
#' @examples
#' geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.n0, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.n0, geom = "text", fun.y = mean, colour = "red")
#' 
#' @export
give.n0 <- function(x) {
  return(c(y = -0.02, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#' function for mean(x) labels
#' 
#' mean placed at y = -1
#' 
#' @examples
#' ggplot(mtcars, aes(factor(cyl), mpg, label=rownames(mtcars))) +
#'   geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.mean, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.mean, geom = "text", fun.y = mean, colour = "red")
#' 
#' @references https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
#' @export
give.mean <- function(x) {
  return(c(y = -1, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}


#' function for mean(x) labels
#' 
#' mean placed at y ~ 1
#' 
#' @examples
#' ggplot(mtcars, aes(factor(cyl), mpg, label=rownames(mtcars))) +
#'   geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.mean, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.mean, geom = "text", fun.y = mean, colour = "red")
#' 
#' @references https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
#' @export
give.n1 <- function(x) {
    return(c(y = 0.34, label = length(x)))
}


#' function for mean(x) labels
#' 
#' mean placed at y ~ 2
#' 
#' @examples
#' ggplot(mtcars, aes(factor(cyl), mpg, label=rownames(mtcars))) +
#'   geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.mean, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.mean, geom = "text", fun.y = mean, colour = "red")
#' 
#' @references https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
#' @export
give.n2 <- function(x) {
    return(c(y = .528, label = length(x)))
}

#' mode
#' @param x vector
#' @param na.rm if remove na.rm
#' 
#' @references https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
estimate_mode <- function(x) {
    d <- density(x)
    out <- d$x[which.max(d$y)]

    return(out)
}
