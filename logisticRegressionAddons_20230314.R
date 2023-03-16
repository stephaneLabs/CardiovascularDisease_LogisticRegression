#' Tools for Logit regression
#'
#' Copyright Stephane Lassalvy 2023 03 14
#' 
#' Licence GPL-3
#' Disclaimer of Warranty.
#'
#' THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
#' EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
#' PROVIDE THE PROGRAM AS IS WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
#' INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#' FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
#' SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
#' REPAIR OR CORRECTION.
#' 
#' 
#' 
library(ggplot2)
library(ggpubr)
library(DescTools)
library(ResourceSelection)
library(moments)
library(arules)
library(car)
library(GGally)
library(InformationValue)
library(pROC)

#' Function : quantitativeVariableDescription
#' Histogram of a quantitative variable, with skewness and kurtosis information
#' Inputs :
#' @param df data frame including the variable of interest
#' @param variableName name of the variable of interest
#' 
#' Outputs :
#' Descriptive graphics for the variable of interest
#' 
quantitativeVariableDescription <- function(df, variableName){
  
  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }
  
  # Extract the variable of interest
  extractVariable <- function(df, variableName) {variable <- eval(parse(text = paste("df$", variableName, sep="")));
                                                 cat("Lenght of the variable :\n")
                                                 cat(length(variable))
                                                 cat("\n")
                                                 cat("Variable summary :\n")
                                                 print(summary(variable))
                                                 return(variable)
                                                }

  variable <- extractVariable(df, variableName)

  # Data frame from the variable of interest only, NA removed
  df_Graphics <- data.frame(variable = variable)
  df_Graphics <- subset(df_Graphics, !is.na(variable))
  df_Graphics <- mutate(df_Graphics, outlier = is_outlier(variable))

  # Graphics
  figureHistogram <- ggplot(df_Graphics, aes(x = variable)) + 
                            geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "lightblue") +
                            geom_density(color="darkgreen") +
                            xlab(paste("Values, Skewness = ", round(skewness(df_Graphics$variable),2), "Kurstosis = ", round(kurtosis(df_Graphics$variable), 2))) +
                            ylab("Prob.") +
                            ggtitle(paste("Histogram of ", toupper(variableName), sep = "")) +
                            geom_vline(xintercept = mean(df_Graphics$variable), color = "blue", linewidth = 1) +
                            geom_vline(xintercept = median(df_Graphics$variable), color = "darkgreen", linewidth = 1)  +
                            stat_function(fun = dnorm, args = list(mean = mean(df_Graphics$variable), sd = sd(df_Graphics$variable)), col = "blue", linewidth = 1)

  figureBoxplot <- ggplot(df_Graphics, aes(x=eval(variableName), y=variable)) + 
                          geom_boxplot(outlier.size=2, outlier.colour="black") + 
                          xlab(toupper(variableName)) + ylab(paste("Values of ", toupper(variableName), sep="")) +
                          geom_text(data = subset(df_Graphics, outlier == TRUE), aes(label = variable), na.rm = TRUE, hjust = -0.3) + 
                          ggtitle(paste("Boxplot and outliers of ", toupper(variableName), sep=""))
  
  figure <- ggarrange(figureHistogram, figureBoxplot, nrow = 1, ncol = 2)
  figure
}


#' Function : removeOutliers
#' Removes the rows of a dataframe having outliers for a given variable
#' Inputs :
#' @param df data frame including the variable of interest
#' @param variableName name of the variable of interest
#' 
#' Outputs :
#' 
#' a list with the initial dataframe and the cleaned dataframe
#' 
removeOutliers <- function(df, variableName ="name"){
  #find Q1, Q3, and interquartile range for values in column A
  cat("Number of rows of initial dataframe : ")
  cat(nrow(df))
  cat("\n")
  
  variable <- eval(parse(text = paste("df$",variableName,sep=""))) 
  Q1 <- quantile(variable, .25)
  Q3 <- quantile(variable, .75)
  IQR <- IQR(variable)
  
  #only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
  df_without_outliers <- subset(df, variable > (Q1 - 1.5*IQR) & variable < (Q3 + 1.5*IQR))
  
  cat("Number of rows of toiletted dataframe : ")
  cat(nrow(df_without_outliers))
  cat("\n")
  
  cat("Percent of rows deleted : ")
  cat(paste(round(nrow(df_without_outliers) / nrow(df), 2), "%"))
  cat("\n")
  return(list(initial_df = df, without_df = df_without_outliers))
}


#' Function : logOddsVsQuantitativePredictor
#' Relationship between the log odds and the quantitative predictor
#' Inputs :
#' @param df data frame including the variable of interest
#' @param binaryResponse name of the binary response"
#' @param quantitativePredictor name of the quantitative predictor
#' @param method method to discretize the quantitative predictor ("round" rounds the quantitative variable to "digits" digits,
#'        "interval" (equal interval width), "frequency" (equal frequency), "cluster" (k-means clustering) 
#'         and "fixed" (categories specifies interval boundaries).
#'         The function then use the rounded variable ("round" option) or compute the mean for each level of the discretization.
#' @param digits number of digits for the quantitative predictor
#' @param breaks option for the discretize function from package arules
#' @param variableName name of the variable of interest
#' 
#' Outputs :
#' Descriptive graphics for the variable of interest
#' 
logOddsVsQuantitativePredictor <- function(df, binaryResponse = "dm", quantitativePredictor, method = "round", digits = 2, breaks = 3){
  # Get the variables of interest
  binaryResp     <- as.factor(eval(parse(text=paste("df$", binaryResponse, sep = ""))))
  xQuantitative  <- as.numeric(eval(parse(text=paste("df$", quantitativePredictor, sep = ""))))
  
  # Discretize the quantitative predictor to observe the relationship between means of predictor values and log(odds) of the response
  if(method=="round"){
    quantPredictor  <- as.factor(round(xQuantitative, digits = digits))
  } else {
    xDiscretization <- as.factor(as.character(discretize(xQuantitative, method = method, breaks = breaks)))
    ddf <- data.frame(xQuantitative = xQuantitative, xDiscretization  = xDiscretization)
    quantPredictor <- aggregate(xQuantitative ~ xDiscretization, data = ddf, FUN = "mean", na.action = "na.omit")
    quantPredictor <- rename(quantPredictor, Xmeans = xQuantitative)
    ddf            <- full_join(ddf, quantPredictor, by = "xDiscretization")
    quantPredictor <- as.numeric(ddf$Xmeans)
  }
  if(any(binaryResp == "")){binaryResp[binaryResp == ""] <- NA}
  if(any(quantPredictor == "")){quantPredictor[quantPredictor == ""] <- NA}

  # Compute the odds for each value of predictor
  frequencyTable <- prop.table(table(quantPredictor , binaryResp), margin = 1)
  odds <- frequencyTable[, "1"] / frequencyTable[, "0"]

  # Make a table with log(odds) and predictor
  oddsTable <- data.frame(predictor = as.numeric(row.names(frequencyTable)),
                          "_1" = frequencyTable[, "1"],
                          "_0" = frequencyTable[, "0"],
                          odds = odds,
                          logodds = log(odds))

  
  oddsTable <- oddsTable[oddsTable$logodds != -Inf & oddsTable$logodds != Inf, ]
  
  # regression and correlations between log(odds) and predictor
  coefOddsModel <- coef(lm(logodds ~ predictor, data = oddsTable))
  correlationPredictorVsLogOdds <- with(oddsTable, round(cor(predictor, odds, use = "complete.obs"), 3))

  # print relationship between log(odds) and predictor
  if(method=="round"){
    with(oddsTable, {figure <- ggplot(data = oddsTable, aes(x = predictor, y = logodds)) +
      geom_point() +
      geom_abline(intercept = coefOddsModel[1], slope = coefOddsModel[2], linetype = 2, color = "red") +
      xlab(paste("Predictor : ", toupper(quantitativePredictor), ", digits = ", digits)) +
      ylab("Log(Odds)") +
      ggtitle(paste("Predictor :", toupper(quantitativePredictor)," vs Log(Odds),\n corr = ", correlationPredictorVsLogOdds))
    figure
    })
  } else {
    with(oddsTable, {figure <- ggplot(data = oddsTable, aes(x = predictor, y = logodds)) +
                               geom_point() +
                               geom_abline(intercept = coefOddsModel[1], slope = coefOddsModel[2], linetype = 2, color = "red") +
                               xlab(paste("Predictor : ", toupper(quantitativePredictor), ", means of ", nrow(oddsTable), " groups")) +
                               ylab("Log(Odds)") +
                               ggtitle(paste("Predictor :", toupper(quantitativePredictor)," vs Log(Odds),\n corr = ", correlationPredictorVsLogOdds))
                     figure
    })
  }
}


#' Function : binaryResponseVSCategoricalPredictor
#' Percentage table and conditionnal percentage tables crossing the binary response with a categorical variable
#' Inputs :
#' @param binaryResponse : name of the binary response variable
#' @param categoricalPredictor : name of the categorical predictor
#' @param digits : number of digits in the text output
#' 
#' Output :
#' Only text outputs
#' 
binaryResponseVSCategoricalPredictor <- function(df, binaryResponse = "dm", categoricalPredictor, digits = 2){
  # Cross table of the two categorical variables
  binaryResp <- as.character(eval(parse(text=paste("df$", binaryResponse, sep = ""))))
  categoPred <- as.character(eval(parse(text=paste("df$", categoricalPredictor, sep = ""))))
  
  binaryResp[binaryResp == ""] <- NA
  categoPred[categoPred == ""] <- NA
  
  binaryResp <- as.factor(binaryResp)
  categoPred <- as.factor(categoPred)
  
  # Discard lines with missing data
  df <- na.omit(data.frame(binaryResp, categoPred))

  responseVsPredictorFrequencies <- with(df, table(binaryResp, categoPred))
  
  # Proportion of observation in each cell of the table
  responseVsPredictorProportions <- round(addmargins(100*prop.table(responseVsPredictorFrequencies)), digits = digits)
  
  # Chi-square test of independance based on the observed frequencies
  responseVsPredictorChi2 <- chisq.test(responseVsPredictorFrequencies)
  
  # Conditional proportions 
  responseVsPredictorConditionnalProportion1 <- round(addmargins(100*prop.table(responseVsPredictorFrequencies, margin = 1), margin = 2), digits = digits)
  responseVsPredictorConditionnalProportion2 <- round(addmargins(100*prop.table(responseVsPredictorFrequencies, margin = 2), margin = 1), digits = digits)
  
  return(invisible(list(frequencies = responseVsPredictorFrequencies,
                        proportions = responseVsPredictorProportions,
                        independanceChi2 = responseVsPredictorChi2,
                        ConditionnalProportions1 = responseVsPredictorConditionnalProportion1,
                        ConditionnalProportions2 = responseVsPredictorConditionnalProportion2
                        )
                   )
         )
}


#' Function plotProportionsResponseConditionnedByPredictor
#' Plot stacked barplot of proportions of binary outcome conditioned by the qualitative predictor levels
#' Inputs :
#' @param tableConditionnalProportions2 conditioned proportions of type 2 given by the function binaryResponseVSCategoricalPredictor
#' 
#' Outputs :
#' Stacked bar plot of the conditioned proportions
#' 
plotProportionsResponseConditionnedByPredictor <- function(tableConditionnalProportions2, VariableName = "Variable"){
  df_graph <- as.data.frame(tableConditionnalProportions2[1:2,])
  figure <- ggplot(df_graph, aes(x = categoPred, y = Freq, fill = binaryResp)) +
            geom_bar(position = "stack", stat = "identity") +
            ylab("Proportions in %") +
            xlab(VariableName) +
            scale_fill_discrete(name="Response") +
            ggtitle(paste("Proportion of each group of the response within each levels of ", VariableName))
  print(figure)
}


#' plotVIF
#' Plot Variable Inflation Factor for variables in a linear of logit model to assess multicolinearity among the predictors
#' Inputs :
#' @param model the considered model
#' 
#' Outputs :
#' Diagram of the VIF value withe 3 thresholds as a rule of thumb : 
#' 2.5 warning, 5 some concerns, 10 serious problem
#' 
plotVIF <- function(model){
  vifDf <- as.data.frame(vif(model))
  vifDf <- rownames_to_column(vifDf, var = "rowName")
  barplotFigure <- ggplot(data = vifDf, aes(x = rowName, y = GVIF)) +
                   geom_col(color ="blue", fill = "lightblue") +
                   xlab("Variable/Factor Predictor") +
                   ylab("GVIF value") +
                   ggtitle("VIF Values for the considered model") +
                   geom_hline(yintercept = c(2.5, 5, 10), color = c("darkgreen", "blue", "red")) +
                   theme(axis.text.x=element_text(angle  = 90, vjust = 0.5, hjust = 1))
  print(barplotFigure)
  return(vifDf)
}


#' Function mcFaddenR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' Mc Fadden R2
#'  
mcFaddenR2 <- function(model, nullModel){
  R2 <- as.numeric(1 - logLik(model) / logLik(nullModel))
  return(R2)
}


#' Function adjMcFaddenR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' R2 of Mc Fadden value 
#' 
adjMcFaddenR2 <- function(model, nullModel){
  R2 <-  as.numeric(1 - (logLik(model) - length(coef(model))) / (logLik(nullModel) - 1))
  return(R2)
}


#' Function sasR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' SAS R2
#' 
sasR2 <- function(model, nullModel){n  <- nrow(na.omit(model$data))
                                    R2 <- as.numeric(1 - exp(-2 * (logLik(model) - logLik(nullModel)) / n))
                                    return(R2)
}

#' Function adjSasR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' Ajusted SAS R2
#' 
adjSasR2 <- function(model, nullModel){
                                       n  <- nrow(na.omit(model$data))
                                       R2 <- as.numeric(1 - exp(-2 * (logLik(model) - logLik(nullModel)) / n))
                                       adjR2 <- as.numeric(R2 / (1 - exp(2 * logLik(nullModel) / n)))
                                       return(adjR2)
}

#' Function devR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' dev R2
#' 
devR2 <- function(model, nullModel){
  n <- nrow(na.omit(model$data))
  nullDeviance  <- model$null.deviance
  modelDeviance <- model$deviance
  deltaDeviance <- nullDeviance - modelDeviance
  # it's strange that 2*(logL(M) - logL(M0)) does not give the same delta-deviance
  R2 <- deltaDeviance / nullDeviance
  return(R2)
}


#' Function R2s
#' Sum up the different R2 measure in a data frame
#' 
#' Inputs :
#' @param model : the model to be assessed
#' @param nullModel : the null model, having the Intercept as the only predictor
#' 
#' Output:
#' result : data frame with the different outputs
#' 
computeR2s <- function(model, nullModel){
  results <- data.frame(mcFadden_R2    = mcFaddenR2(model, nullModel),
                        adjMcFadden_R2 = adjMcFaddenR2(model, nullModel),
                        sas_R2         = sasR2(model, nullModel),
                        adjSas_R2      = adjSasR2(model, nullModel),
                        dev_R2         = devR2(model, nullModel)
  )
}


#' Function HoslemTest
#' Performs the Hoslem Lemeshow test and the related charts
#' Inputs :
#' @param model   : the current model M
#' @param ngroups : number of groups to perform de test
#' 
#' Outputs :
#' Print the results of the test and displays the related graphics 
#' 
HoslemTest <- function(model, ngroups = 10){
  # Compute and display Hoslem Ledeshow test results
  hoslemTest <- hoslem.test(x = model$y, y = fitted(model), g = ngroups)
  print(hoslemTest)

  df <- as_tibble(cbind(hoslemTest$observed, hoslemTest$expected))
  
  # Display the related graphics
  # Observed cases vs Predicted cases for each of the 10 groups
  observed1VsPredicted1 <- ggplot(data = df, aes(x = y1, y = yhat1))  +
                           geom_point() +
                           geom_abline(intercept = 0, slope = 1, color="red") +
                           geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                           xlab("Observed cases") + 
                           ylab("Predicted cases") +
                           ggtitle("Obs. cases Vs Pred. cases")

  # print(observed1VsPredicted1)
  observed0VsPredicted0 <- ggplot(data = df, aes(x = y0, y = yhat0))  +
                           geom_point() +
                           geom_abline(intercept = 0, slope = 1, color="red") +
                           geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                           xlab("Observed cases") + 
                           ylab("Predicted cases") +
                           ggtitle("Obs. non cases Vs Pred. non cases")
  
  df <- df %>% mutate(total = y1 + y0, 
                      totalexpected = yhat1 + yhat0,
                      y1pct = y1 / total * 100,
                      yhat1pct= yhat1 / totalexpected * 100)
  
  # Observed cases vs Predicted cases for each of the 10 groups in %
  observed1VsPredicted1pct <- ggplot(data = df, aes(x = y1pct, y = yhat1pct))  +
                              geom_point() +
                              geom_abline(intercept = 0, slope = 1, color="red") +
                              geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                              xlab("Observed cases") + 
                              ylab("Predicted cases") +
                              ggtitle("Obs. cases Vs Pred. cases in %")
  
  outputGraph <- ggarrange(observed1VsPredicted1, observed0VsPredicted0, observed1VsPredicted1, nrow = 1, ncol = 3)
  print(outputGraph)
  }

#' Residuals plots
#' 
residualPlots <- function(model, binaryresponse = "DiabeteDiagostic", quantitativePredictor = "none", label.size = 3){
 df <- model$data
 df$response <- as.factor(eval(parse(text = paste("df$", binaryresponse))))
 
 df$cooksDistances    <- cooks.distance(model)
 df$cooksIndex        <- seq(along.with = as.numeric(df$cooksDistances))
 
 df$leverages        <- influence(model)$hat
 df$indexLeverages   <- seq(along.with = as.numeric(df$leverages))
 
 df$pearsonResiduals  <- residuals(model, type = "pearson")
 df$pearsonIndex      <- seq(along.with = as.numeric(df$pearsonResiduals))
 
 df$devianceResiduals <- residuals(model, type = "deviance")
 df$devianceIndex     <- seq(along.with = as.numeric(df$devianceResiduals))

 df$devianceResiduals <- residuals(model, type = "deviance")
 df$devianceIndex     <- seq(along.with = as.numeric(df$devianceResiduals))

 df <- df %>% mutate(stdPearsonResiduals  = pearsonResiduals  / sqrt(1 - leverages),
                     stdDevianceResiduals = devianceResiduals / sqrt(1 - leverages))
 
 if(quantitativePredictor == "none"){
   # Cook's distances figure
   cooksTreshold <- 4 / nrow(df)
   CooksFigure  <- ggplot(data = df, aes(ymax = cooksDistances, ymin = 0, x=cooksIndex, color = response)) +
                          geom_linerange() +
                          geom_text(data = subset(df, cooksDistances > cooksTreshold),
                                    aes(y = cooksDistances, x=cooksIndex, label = cooksIndex), 
                                    hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                          geom_hline(yintercept = cooksTreshold, color = "red") +
                          xlab("Observation index") +
                          ylab("Cook's distances") +
                          ggtitle("Cook's distances") +
                          scale_colour_discrete(name = "Response")
   
   # Pearson residuals figure
   pearsonResFigure  <- ggplot(data = df, aes(y = pearsonResiduals, x=pearsonIndex, color = response)) +
                        geom_point() +
                        geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                  aes(y = pearsonResiduals, x=pearsonIndex, label = pearsonIndex),
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        geom_hline(yintercept = 0, color = "black") +
                        xlab("Observation index") +
                        ylab("Pearson residuals") +
                        ggtitle("Pearson residuals") +
                        scale_colour_discrete(name = "Response") 
   
   # Pearson standardized residuals figure
   stdPearsonResFigure  <- ggplot(data = df, aes(y = stdPearsonResiduals, x=pearsonIndex, color = response)) +
                           geom_point() +
                           geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                     aes(y = stdPearsonResiduals, x=pearsonIndex, label = pearsonIndex), 
                                     hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                           geom_hline(yintercept = 0, color = "black") +
                           xlab("Observation index") +
                           ylab("Standardized Pearson residuals") +
                           ggtitle("Standardized Pearson residuals") +
                           scale_colour_discrete(name = "Response")
   
   # Leverages figures
   nParameters <- sum(df$leverages)
   leveragesTreshold <- 3 * nParameters / nrow(df)
   leveragesFigure   <- ggplot(data = df, aes(ymax = leverages, ymin = 0, x=indexLeverages, color = response)) +
                        geom_linerange() +
                        geom_hline(yintercept = leveragesTreshold, color = "red") +
                        geom_text(data = subset(df, leverages > leveragesTreshold), 
                                  aes(y = leverages, x=indexLeverages, label = indexLeverages), 
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        xlab("Observation index") +
                        ylab("Leverages") +
                        ggtitle("Leverages")  +
                        scale_colour_discrete(name = "Response")
   
     
   # Deviance residuals figure
   devianceResFigure <- ggplot(data = df, aes(y = devianceResiduals, x=devianceIndex, color = response)) +
                        geom_point() +
                        geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                  aes(y = devianceResiduals, x=devianceIndex, label = devianceIndex), 
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        geom_hline(yintercept = 0, color = "black") +
                        xlab("Observation index") +
                        ylab("Deviance residuals") +
                        ggtitle("Deviance residuals") +
                        scale_colour_discrete(name = "Response")       
   
   # Deviance standardized residuals figure
   stdDevianceResFigure <- ggplot(data = df, aes(y = stdDevianceResiduals, x=devianceIndex, color = response)) +
                           geom_point() +
                           geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                     aes(y = stdDevianceResiduals, x=devianceIndex, label = devianceIndex), 
                                     hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                           geom_hline(yintercept = 0, color = "black") +
                           xlab("Observation index") +
                           ylab("Standardized Deviance residuals") +
                           ggtitle("Standardized Deviance residuals")  +
                           scale_colour_discrete(name = "Response")
   
   figure <- ggarrange(CooksFigure, pearsonResFigure, stdPearsonResFigure, leveragesFigure, devianceResFigure, stdDevianceResFigure, nrow = 2, ncol = 3)
   figure
 }
}

#' Function plotOddsRatios
#'
#' Inputs :
#' @param oddsRatiosDf odds ratios dataframe with confidence intervals boundaries
#' 
#' Outputs : 
#' Display Odds Ratios with confidence intervals figure
#' 
plotOddsRatios <- function(oddsRatiosDf){
  figure <- ggplot(oddsRatiosDf, aes(x = factor(rowname), y = oddsRatios, ymin = Lower2.5pct, ymax = Upper97.5pct)) +
    geom_pointrange() + geom_hline(yintercept = 1, color = "red", ) +
    ggtitle("Odds ratios and confidence intervals") + xlab("Variables/Factor Predictor") + ylab("Odds Ratios") +
    theme(axis.text.x=element_text(angle  = 90, vjust = 0.5, hjust = 1))
  print(figure)
}


#' Function rocCurve
#' Inputs
#' @param trainedModel   : model trained on the trained data set
#' @param testDf         : test data set to assess the predictive ability of the model
#' @param optimizeMethod : Calculus of the optimal treshold for a decision rule : Youden" calculus based on the maximum Youden metric, 
#'                         "topleftcorner" calculus based on the minimum distance between the ROC curve and the point of coordinates (0; 1),
#'                         "MCER" is based on the misclassification error rate which have to be minimal. These indicators are geometric and heuristic,
#'                         their computation may be enhanced by more complicated calculus, but it's a good starting point.
#' @param optimizeDigits : number of decimals for threshold accuracy
#' Outputs :
#' Graphics of ROC curves and misclassification error rate
#' 
rocCurve <- function(trainedModel, df_train = NULL, df_test = NULL, outcomeVariable = "Outcome", optimizeMethod = "Youden", optimizeDigits = 2){
  if(is.null(df_train)){stop("You must specify a training data set.")}
  
  cat("Training data set\n")
  cat("-----------------\n")
  
  cat("Dimensions of the testing data set : ")
  cat(dim(df_train))
  cat("\n")
  
  # Observed outcomes
  trainObservedOutcome <- eval(parse(text = paste("df_train$", outcomeVariable, sep = "")))
  cat("length of the outcome variable : ")
  cat(length(trainObservedOutcome))
  cat("\n")
  
  # Predictions with the model for the training data set
  if("glm" %in% class(trainedModel) && trainedModel$family$family == "binomial" && trainedModel$family$link=="logit" ){
    trainPredictedProbs  <- try(predict(trainedModel, type = "response"))
  }
  
  if("train" %in% class(trainedModel) && trainedModel$method == "glm"){
    trainPredictedProbs  <- try(predict(trainedModel, type = "prob")[,"X1"])
  }
  
  if(model$method == "LogitBoost"){
    trainPredictedProbs  <- try(predict(trainedModel, type = "prob")[,"X1"])
  }
  
  cat("length of the predicted variable : ")
  cat(length(trainPredictedProbs))
  cat("\n")
  
  #' Function sensibilityAndSpecificity
  #' 
  #' Inputs 
  #' @param observedOutcome observed binary outcome 1 = "yes"
  #' @param predictedProbs  predicted probabilities for binary outcome being 1 = "yes"
  #' 
  #' Outputs
  #' 
  #' @return result a data frame containing False predictive rate, True Positive Rate and Misclassification error rate
  sensibilityAndSpecificity <- function(observedOutcome, predictedProbs, threshold){
    # Compute predicted binary outcome
    predictedOutcome  <- factor(x = ifelse(predictedProbs > threshold, "X1", "X0"), levels = c("X0", "X1"))

    # Confusion matrix between observed and predicted binary outcome
    conf_matrix <- table(predictedOutcome, observedOutcome)
    
    # False predictive rate, True predictive rate and Misclassification error rate
    fac <- factor(c("X0","X1"), levels = c("X0", "X1"))
    truePositiveRate  <- caret::sensitivity(conf_matrix, reference = fac, positive = fac[1])
    falsePositiveRate <- 1 - caret::specificity(conf_matrix, reference = fac, negative = fac[2])
    misClassErrRate     <- (conf_matrix[1,2] + conf_matrix[2,1]) / sum(conf_matrix)
    result              <- data.frame(FPR = falsePositiveRate, TPR = truePositiveRate, MCER = misClassErrRate, threshold = threshold) 
    return(result)
  }

  # Thresholds for ROC curves (for both training data set and testing data set)
  thresholds <- seq(0.00, 1.00, by = 1 / (10^optimizeDigits))

  # Results for training data set
  results_train.list <- lapply(X = thresholds, FUN = sensibilityAndSpecificity, observedOutcome = trainObservedOutcome, predictedProbs = trainPredictedProbs)
  results_train      <- do.call(what = rbind, args = results_train.list)
  results_train      <- results_train %>% mutate(YoudenMetric = TPR - FPR, TLCMetric = sqrt(FPR^2 + (TPR-1)^2))

  # Compute the threshold for an optimal decision rule
  cat("Train data set\n")
  if(optimizeMethod == "topleftcorner"){optim_train_rank <- which.min(results_train$TLCMetric);    cat("Optimization method : topleftcorner\n")}
  if(optimizeMethod == "Youden"       ){optim_train_rank <- which.max(results_train$YoudenMetric); cat("Optimization method : Youden metric\n")}
  if(optimizeMethod == "MCER"         ){optim_train_rank <- which.min(results_train$MCER);         cat("Optimization method : MCER metric\n")}
  if(optimizeMethod != "topleftcorner" && optimizeMethod != "Youden" && optimizeMethod != "MCER"){stop("Use topleftcorner, Youden or MCER as optimization method\n")}

  optim_train_rank  <- min(optim_train_rank)
  optim_train <- results_train[optim_train_rank,]

  # ROC Curve
  rocCurve <- ggplot(data = results_train, aes(x = FPR, y = TPR)) +
              geom_line(color = "#009E73") +
              geom_abline(slope = 1, intercept = 0, color ="blue") +
              geom_text(aes(x = 0.4, y = 0.5, label = paste("AUC (train) = ", round(pROC::auc(trainObservedOutcome, trainPredictedProbs), 4))), color = "#009E73") +
              geom_vline(xintercept = optim_train$FPR, color = "#009E73") +
              geom_point(data = optim_train, aes(x = FPR, y = TPR), color = "red") +
              geom_text(data  = optim_train, aes(x = FPR, y = TPR, label = paste("thr. = ", round(optim_train$threshold, 2))), nudge_y = 0.025, color = "red") +
              xlab("False Positive Rate = 1 - Specificity") + ylab("True Positive Rate = Sensitivity") + ggtitle(paste("ROC Curve, optimization of cutoff with option", optimizeMethod))
  
  # Misclassification Error curve
  MisClassifErrCurve <- ggplot(data = results_train, aes(x = threshold, y = MCER)) +
                        geom_line(color = "#009E73") +
                        geom_vline(xintercept = optim_train$threshold, color = "#009E73") +
                        xlab("Thresholds") + ylab("Misclassification error rate") + ggtitle("Misclassification error rate")
  
  # Results for testing data set if there is a testing data set
  if(!is.null(df_test) == TRUE){
    cat("Training data set\n")
    cat("-----------------\n")
    
    cat("Dimensions of the testing data set : ")
    cat(dim(df_test))
    cat("\n")
    
    if(any(names(df_test) != names(df_train))){stop("df_train and df_test must have the same columns")}
    testObservedOutcome <- eval(parse(text = paste("df_test$", outcomeVariable, sep = "")))
    
    cat("length of the outcome variable : ")
    cat(length(testObservedOutcome))
    cat("\n")
    
    # Predictions with the same model fit on the test data set
    if("glm" %in% class(trainedModel) && trainedModel$family$family == "binomial" && trainedModel$family$link=="logit" ){
      testPredictedProbs  <- try(predict(trainedModel, newdata  = df_test, type = "response"))
    }
    
    if("train" %in% class(trainedModel) && trainedModel$method == "glm"){
      testPredictedProbs  <- try(predict(trainedModel, newdata  = df_test, type = "prob")[,"X1"])
    }
    
    if(model$method == "LogitBoost"){
      testPredictedProbs  <- try(predict(trainedModel, newdata  = df_test, type = "prob")[,"X1"])
    }
    
    cat("length of the predicted variable : ")
    cat(length(testPredictedProbs))
    cat("\n")

    results_test.list   <- lapply(X = thresholds, FUN = sensibilityAndSpecificity, observedOutcome = testObservedOutcome, predictedProbs = testPredictedProbs)
    results_test        <- do.call(what = rbind, args = results_test.list)
    results_test        <- results_test %>% mutate(YoudenMetric = TPR - FPR, TLCMetric = sqrt(FPR^2 + (TPR-1)^2))
    
    # Compute the threshold for an optimal decision rule
    cat("Test data set\n")
    if(optimizeMethod == "topleftcorner"){optim_test_rank <- which.min(results_test$TLCMetric);    cat("Optimization method : topleftcorner\n")}
    if(optimizeMethod == "Youden"       ){optim_test_rank <- which.max(results_test$YoudenMetric); cat("Optimization method : Youden metric\n")}
    if(optimizeMethod == "MCER"         ){optim_test_rank <- which.min(results_test$MCER);         cat("Optimization method : MCER  metric\n")}
    if(optimizeMethod != "topleftcorner" && optimizeMethod != "Youden" && optimizeMethod != "MCER"){stop("Use topleftcorner, Youden or MCER as optimization method")}
    
    optim_test_rank  <- min(optim_test_rank)
    optim_test <- results_test[optim_test_rank,]
    
    # ROC Curve
    rocCurve <- rocCurve +
      geom_line(data = results_test, aes(x = FPR, y = TPR), color = "#E69F00") +
      geom_abline(slope = 1, intercept = 0, color ="blue") +
      geom_text(aes(x = 0.4, y = 0.4, label = paste("AUC (test) = ", round(pROC::auc(testObservedOutcome, testPredictedProbs), 4))), color = "#E69F00") +
      geom_vline(xintercept = optim_test$FPR, color = "#E69F00") +
      geom_point(data = optim_test, aes(x = FPR, y = TPR), color = "red") +
      geom_text(data  = optim_test, aes(x = FPR, y = TPR, label = paste("thr. = ", round(optim_test$threshold, 2))), nudge_y = 0.025, color = "red")
    
    # Misclassification Error curve
    MisClassifErrCurve <- MisClassifErrCurve +
      geom_line(data = results_test, aes(x = threshold, y = MCER), color = "#E69F00") +
      geom_vline(xintercept = optim_test$threshold, color = "#E69F00")
  }
  
  # Final plot
  rocDisplay <- ggarrange(rocCurve, MisClassifErrCurve, nrow = 1, ncol = 2)
  print(rocDisplay)
  return(invisible(list(results_train = results_train, optim_train = optim_train, results_test = results_test, optim_test = optim_test)))
}
