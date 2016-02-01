# runApp("C:/Users/L/Desktop/App/diProject")

library(doBy)
library(shiny)
library(ggvis)
library(dplyr)
library(survival)
library(googleVis)

load("./data/data648.RData")

model <- coxph(Surv(start, stop, sepshock) ~ hr + sbp + wbc + bun + creatinine + platelets + temperature + resp_rate, data=data648)

bestmodels <- NULL

########################
##### Shiny Server #####
########################

shinyServer(function(input, output) {

  # Drop-down selection box for patient list
  output$choose_patient <- renderUI({
      patientList <- renewList()
      selectInput("patientid", "Choose A Patient", patientList)
  })

  renewList <- eventReactive(input$renew_list,{
    patientList <- sample(data648$icustay_id, 10, replace=F)
    return(patientList)
  }, ignoreNULL = FALSE)


  # Time Frame for chosen patient
  output$time_slider <- renderUI({
      if(is.null(input$patientid))
        return()

      patientid <- input$patientid
      patienttime <- which(data648$icustay_id == patientid)
      timemax <- data648$stop[patienttime[length(patienttime)]]
      ticks <- data648$stop[which(data648$icustay_id == patientid)]
      sliderInput(inputId = "time",
                label = "Time",
                min = 1,
                max = timemax,
                value = 1,
                step = 10,
                animate = animationOptions(interval=200, loop=FALSE)
      )
  })


############################################
### Exhaustive Search for the Besr Model ###
############################################

  findBest <- eventReactive(input$varSelect,{
    # Create vectors for outcome and predictors
    outcome    <- c("survival.vector")
    predictors <- names(data648)[5:12]
    survival.vector    <- Surv(data648$start, data648$stop, data648$sepshock)
    dataset    <- data648

    # Create list of models
    list.of.models <- lapply(seq_along((predictors)), function(n) {
      left.hand.side  <- outcome
      right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

      paste(left.hand.side, right.hand.side, sep = "  ~  ")
    })

    # Convert to a vector
    vector.of.models <- unlist(list.of.models)

    # Fit coxph to all models
    list.of.fits <- list()

    withProgress(message = 'Finding the Best Model', value = 0, {
      # Number of times we'll go through the loop
      n <- length(vector.of.models)
      n <- 100

      for (i in 1:length(vector.of.models)){
        formula    <- as.formula(vector.of.models[i])
        fit        <- coxph(formula, data = dataset)
        result.AIC <- extractAIC(fit)

        list.of.fits[[i]] <- data.frame(num.predictors = result.AIC[1], AIC = result.AIC[2], model = vector.of.models[i])
      
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Fitting Model ", i, ", we're ", round(i / n * 100), "% done.", sep=""))

        # Pause for 0.0001 seconds to simulate a long computation.
        Sys.sleep(0.0001)
      }
    })
 
    # Collapse to a data frame
    result <- do.call(rbind, list.of.fits)

    # Sort and print
    bestmodels <- orderBy(~ AIC, result)

    return(bestmodels)
  })


######################################
### Fit Individual Survival Models ###
######################################

  getSurv <- eventReactive(input$plotSurv,{
    if(is.null(input$patientid))
      return()

    patientid <- input$patientid
    datatemp <- data648[which(data648$icustay_id == patientid),]
  
    starts <- data648$start[which(data648$icustay_id == patientid)]
    nstart <-length(starts)
    endtime <- starts[length(starts)]

    survList <- list()

    withProgress(message = 'Fitting Survival Models', value = 0, {
      # Number of times we'll go through the loop is nstart

      for (i in 1:nstart){
        nd <- datatemp[i, 5:12]
        temp <- survfit(model, newdata = nd)

        survprobs <- summary(temp)$surv
        uppers <- summary(temp)$upper
        lowers <- summary(temp)$lower

        survList[[i]] <- data.frame(survprobs=survprobs, uppers=uppers, lowers=lowers)  

        # Increment the progress bar, and update the detail text.
        incProgress(1/nstart, detail = paste("There are ", nstart, ", measurements for this parient. ",
                                             "It's gonna take about ", floor(nstart * 2 / 60), ", minutes. ",
                                             "Fitting Model ", i,", we're ", round(i / nstart * 100), "% done.", sep=""))

        # Pause for 0.0001 seconds to simulate a long computation.
        Sys.sleep(0.0001)     
      }
    })

    times <- summary(survfit(model))$time
 
    list(survList=survList, starts=starts, endtime=endtime, times=times)
  })


#########################################
### Calculate Individual Hazard Rates ###
#########################################
  # now get baseline curve
  baseline <- basehaz(model, centered=FALSE)

  getHaz <- reactive({
    if(is.null(input$patientid))
      return()

    patientid <- input$patientid
    datatemp <- data648[which(data648$icustay_id == patientid),]

    starts <- data648$start[which(data648$icustay_id == patientid)]
    nstart <-length(starts)
    endtime <- starts[length(starts)]

    hazPred <- matrix(0, nrow(baseline), nstart)

    for (j in 1:nstart){
      hazPred[,j] <- exp(model$coef %*% t(datatemp[j, 5:12])) * baseline$hazard
    }

    haztimes <- baseline$time

    list(starts=starts, endtime=endtime, haztimes=haztimes, hazPred=hazPred)
  })


################################################
### Calculate Fitted Cumulative Hazard Rates ###
################################################
  # compute fitted hazards for FIRST m subjects
  m <- 10
  subjectTime <- data648 %>%
                   group_by(icustay_id) %>% 
                   summarise_each(funs(max(.,na.rm=TRUE)), hosp_time) 
  subjectTime <- subjectTime$hosp_time[1:m]
  hazardHat <- matrix(max(baseline$hazard), length(baseline$time), m)
    
  plotHazard <- eventReactive(input$plotHazard,{
    withProgress(message = 'Calculating Fitted Cumulative Hazards', value = 0, {
      # Number of times we'll go through the loop is m 
      for (j in 1:m){
        counter <- 1
        tmpIndx <- which(baseline$time==subjectTime[j])
        for (i in 1:tmpIndx){
          hazardHat[i,j] <- exp(model$coef %*% t(data648[counter,c(5:12)])) * baseline$hazard[i]
          counter <- counter + 1
        }
        # Increment the progress bar, and update the detail text.
        incProgress(1/m, detail = paste("There are ", m, ", paitients. ", "Calculating fitted hazards for patient ", j,
                                             ", we're ", round(j / m * 100), "% done.", sep=""))
      }
    })
    
    subIndx <- floor(seq(1, 24860,length.out=50))
    hazardHat <- hazardHat[subIndx, ]

    nsubs <- unique(data648$icustay_id)[1:m]
    hazards <- c()
    ids <- c()
    times <- c()

    for (i in 1:m){
      hazards <- c(hazards, hazardHat[,i])
      ids <- c(ids, paste("subject ", rep(nsubs[i], nrow(hazardHat)), sep=""))
      times <- c(times, baseline$time[subIndx])
    }

    plotData <- data.frame(ids, times, cumHaz=hazards)
    plotData$times <- plotData$times%/%60+1

    for (k in 2:nrow(plotData)){
      if (plotData[k,3] == plotData[k-1,3]) plotData[k,3] <- plotData[k,3] + 0.0001
    }
 
    return(plotData)
  })



##############
### Output ###
##############
  
  output$haz_predicted <- renderPlot({ 
    hazResult <- getHaz()

    if(is.null(hazResult))
      return() 
    
    starts   <- hazResult$starts
    endtime  <- hazResult$endtime
    haztimes <- hazResult$haztimes
    hazPred  <- hazResult$hazPred

    currenttime <- input$time
    starttemp <- c(currenttime, starts)
    indxtemp <- which(order(starttemp)==1) - 1 
    if (indxtemp == 0) {indxtemp <- 1}

    xmax <- endtime + haztimes[length(haztimes)]
    ymax <- max(hazPred)

    if(length(haztimes)==0 | is.null(currenttime))
      return()

    plot(haztimes+currenttime, hazPred[,indxtemp], type="l", xlim=c(1,xmax), ylim=c(0, ymax), xlab="Time", ylab="Predicted Hazards")
    abline(v = haztimes[1]+currenttime, col = "darkgray", lty = 3)
  })

  output$plot_predicted <- renderPlot({ 
    survResult <- getSurv()

    if(is.null(survResult))
      return() 
   
    survList <- survResult$survList
    starts  <- survResult$starts
    endtime <- survResult$endtime
    times <- survResult$times

    currenttime <- input$time
    starttemp <- c(currenttime, starts)
    indxtemp <- which(order(starttemp)==1) - 1 
    if (indxtemp == 0) {indxtemp <- 1}

    survprobs <- survList[[indxtemp]]$survprobs
    uppers <- survList[[indxtemp]]$uppers
    lowers <- survList[[indxtemp]]$lowers

    indxymin <- which(times >= endtime)[1]
    ymin <- round(min(sapply(survList, function(x){min(x[1:indxymin,2])})), 4)
    xmax <- endtime + times[length(times)]
    
    if(length(times)==0)
      return() 

    plot(times+currenttime, survprobs, type="s", xlim=c(1,xmax), ylim=c(ymin, 1), xlab="Time", ylab="Proportion of No Septic Shock")
    lines(times+currenttime, uppers, type="s", lty = 2)
    lines(times+currenttime, lowers, type="s", lty = 2)
    abline(v = times[1]+currenttime, col = "darkgray", lty = 3)
  })
  
  output$summary <- renderPrint({
    summary(model)
  })
  
  output$raw_data <- renderTable({
    head(data648, 100)
  }, include.rownames = FALSE)
   

  output$best_model <- renderTable({
    bestmodels <- findBest()
 
    if(is.null(bestmodels))
      return()

    head(bestmodels, 10)
  }, include.rownames = FALSE)

  
  output$test <- renderText({
    hazResult <- getHaz()

    if(is.null(hazResult))
      return() 
    
    starts   <- hazResult$starts
    endtime  <- hazResult$endtime
    haztimes <- hazResult$haztimes
    hazPred  <- hazResult$hazPred

    currenttime <- input$time
    starttemp <- c(currenttime, starts)
    indxtemp <- which(order(starttemp)==1) - 1 
    if (indxtemp == 0) {indxtemp <- 1}

    xmax <- endtime + haztimes[length(haztimes)]
    ymax <- max(hazPred)

    if(length(haztimes)==0)
      return()
    
    c(length(haztimes+currenttime), length(hazPred[,indxtemp]))
  })


  output$haz_fitted <- renderGvis({
    plotData <- plotHazard()
  
    if(is.null(plotData))
      return()
    
    return (gvisMotionChart(plotData, idvar="ids", timevar="times", sizevar="cumHaz",colorvar="cumHaz"))
  })
 
})
