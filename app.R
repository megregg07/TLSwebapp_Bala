
##
## 12/01/2025
## TLS WEBAPP - BALA'S VERSION
##
## This app is similar to the one that was created summer 2024 
## (that will be polished later into the final practitioner version). 
## The differences are:
##  1) it only performs the 'Part II' analysis (does not include 'Part I')
##  2) the data input will be two files, with each file containing Cartesian coordinates
##     for all positions collected at that timepoint. 
##     e.g., It no longer expects one file per position
##
##  updated 03/02/2026
##  3) User's choice to either test just the angular residuals or angles + ranging 
##     (this will change the dimension of the covariance matrices)
##  update 03/03/2026
##  4) added residual plots to 'data check' tab
##


library(shiny)
library(shinythemes)
library(ICSNP)     ## for Hotellings T2 test
library(heplots)   ## for visualizing data ellipses
library(dplyr)     ## %>% transformations
library(geometry)  ## for 'cart2sph' function
library(MASS)      ## for mvrnorm() 
library(ggplot2)
library(tidyr)
library(tools)     ## for determining file extension
library(stringr)   ## for naming things


## define where to get functions from
##**USE THIS WHEN POSTING/DISTRIBUTING*
source('R/utils.R')  

##**USE THIS WHEN ON MY COMPUTER
## Import the functions from Bala's verison of the webapp
#source("C:/Users/meg3/OneDrive - NIST/Work/OSAC/CSIR subcommitee/TLS/code/WebApp/TLSwebapp_BalaVersion/R/utils.R")


##############################
##
## DEFINE THE UI
##
##############################
ui <- fluidPage(theme = shinytheme("spacelab"),
                
                # Application title
                titlePanel("TLS Part II Analysis - Bala's version"),
                
                # Data input box, with dropdown for units
                sidebarLayout(
                  sidebarPanel(
                    radioButtons('upload_or_default', "Upload Files or Run Example?", 
                                 choices = c('Upload Files', 'Run Example Files'), 
                                 selected = 'Upload Files'),
                    conditionalPanel(
                      condition = "input.upload_or_default == 'Upload Files'", 
                      h2('TLS Data'),
                      fileInput('file1',"Baseline data", accept = c(".csv", ".txt")),
                      div(style = "margin-top: -20px"),
                      textInput(inputId = "dataname_base",
                                label = "Baseline Data Nickname:",
                                value = "Baseline"),
                      fileInput('file2', "Test data", accept = c(".csv", ".txt")),
                      div(style = "margin-top: -20px"),
                      textInput(inputId = "dataname_test",
                                label = "Test Data Nickname:",
                                value = "Test"), 
                      selectInput('units','TLS units',choices=c('mm', 'meters', 'inches')),
                      ),
                    #fileInput('filename_baseline_p1', "Position 1", accept = c(".csv", ".txt")), 
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_baseline_p2', "Position 2", accept = c(".csv", ".txt")),
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_baseline_p3', "Position 3", accept = c(".csv", ".txt")),
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_baseline_p4', "Position 4", accept = c(".csv", ".txt")),
                    #h2('TLS - Test Data'),
                    #fileInput('filename_testdata_p1', "Position 1", accept = c(".csv", ".txt")), 
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_testdata_p2', "Position 2", accept = c(".csv", ".txt")), 
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_testdata_p3', "Position 3", accept = c(".csv", ".txt")), 
                    #div(style = "margin-top: -25px"),
                    #fileInput('filename_testdata_p4', "Position 4", accept = c(".csv", ".txt")),
                    #div(style = "margin-top: -25px"),
                    #selectInput('units','TLS units',choices=c('mm', 'meters', 'inches')),
                    #h2("Tape Measurements"),
                    #fileInput('filename_tapedata', "Dataset of Tape Measurements", accept = c(".csv", ".txt")), 
                    #div(style = "margin-top: -25px"),
                    #selectInput('units_tape','Tape measure units',choices=c('mm', 'cm', 'inches')),
                    #selectInput('includeThreshold', "Set error threshold?", choices = c('No', 'Yes')),
                    ## conditionl numeric input that appears if 'includeThreshold' == Yes
                    #conditionalPanel(
                    #  condition = "input.includeThreshold == 'Yes'",
                    #  numericInput("threshold","Error threshold",value = 1, step = 0.1)
                    #),
                    #fluidRow(
                    #  column(3,
                    #         numericInput("tapeA", label = "Length A", value = 1000, min = 0.01, step = 0.01)),
                    #  column(3, 
                    #         selectInput("t1_A", label = "Target 1", choices=1:20)), 
                    #  column(3, 
                    #         selectInput("t2_A", label = "Target 2", choices=1:20, selected = 2)),
                    #),
                    
                    actionButton('runAnalysis','Run Analysis', class = "btn-primary"), 
                  ),
                  
                  
                  # Prepare space for all the output
                  mainPanel(
                    tabsetPanel(selected = "Part II - Analysis Results", 
                                br(),
                                tabPanel("Information", 
                                         uiOutput("info")),
                                ## prepare the output space for "Part II - Analysis" tab 
                                tabPanel("Part II - Analysis Results",
                                         h1("Conclusion - angles + ranging"),
                                         textOutput("p2_conclusion"),
                                         verbatimTextOutput("p2_analysis"), 
                                         br(),
                                         #verbatimTextOutput("p2_interpretation"),
                                         h1("Conclusion - angles only"),
                                         textOutput("p2_conclusion_bi"),
                                         verbatimTextOutput("p2_analysis_bi"),
                                         br(),
                                         #textOutput("p2_interpretation"),
                                         #br(),
                                         #br(),
                                         h3("Estimated covariance matrices"),
                                         uiOutput("covMat_info"),
                                         h4("Sigma hat - Baseline data"),
                                         verbatimTextOutput("SigmaHat1"),
                                         h4("Sigma hat - Test data"),
                                         verbatimTextOutput("SigmaHat2"), 
                                         #br(),
                                         #h4("Data Summary"),
                                         #textOutput("basedata_summary"),
                                         #textOutput("testdata_summary"),
                                         #textOutput("units")
                                         #br(),
                                         #h4("Statistical Details for Robust Omnibus Test"),
                                         #verbatimTextOutput("p2_analysis"), 
                                         #verbatimTextOutput("p2_analysis_bi")
                                ), 
                                tabPanel("Part II - Data Ellipses", 
                                         plotOutput("ellipsePlot", width = "600px", height = "600px")
                                ), 
                                tabPanel("Data Summary", 
                                         h3("Data Summary"),
                                         textOutput("basedata_summary"),
                                         textOutput("testdata_summary"),
                                         textOutput("units"), 
                                         br(),
                                         h3("Check for problem targets"), 
                                         uiOutput("residualStatement"),
                                         plotOutput("resPlot_base", width = "600px", height = "600px"), 
                                         plotOutput("resPlot_test", width = "600px", height = "600px")
                                         )
                    )
                  )
                )
)


##############################
##
## DEFINE THE SERVER LOGIC
##
##############################
server <- function(input, output) {
  
  ##IMPORT BASELINE DATA 
  ## (Shendong's 2022 data set, if using example file)
  data_p1base = reactive({
    if(input$upload_or_default== 'Run Example Files') {
      #data_p1base = read.table('ExampleData/Surphaser_Shendong_2022.txt', 
      #                         header = FALSE)
      data_p1base = read.csv('ExampleData2/Surphaser_Shendong_2022_FF.csv')
      return(data_p1base)
    } else if(is.null(input$file1)) {
      return(NULL)
    } else {
      ext = tools::file_ext(input$file1$name)
      data_p1base = switch(ext, csv = read.csv(input$file1$datapath), 
                           txt = read.table(input$file1$datapath, header = FALSE), 
                           validate("Invalid file; Please upload a .csv or .txt file"))
      return(data_p1base)
    }
  })
  

  ## IMPORT TEST DATA
  ## (Bala's 2025 data set, if using example file)
  data_p1test = reactive({
    if(input$upload_or_default== 'Run Example Files') {
      #data_p1test = read.table('ExampleData/Surphaser_Bala_2025.txt', 
      #                         header = FALSE)
      data_p1test = read.csv('ExampleData2/Surphaser_Bala_2025_FF.csv')
      return(data_p1test)
    } else if(is.null(input$file2)) {
      return(NULL)
    } else {
      ext = tools::file_ext(input$file2$name)
      data_p1test = switch(ext, csv = read.csv(input$file2$datapath), 
                           txt = read.table(input$file2$datapath, header = FALSE), 
                           validate("Invalid file; Please upload a .csv or .txt file"))
      return(data_p1test)
    }
  })
  
  ## 
  ## NOW DO STUFF, WHEN THE BUTTON IS CLICKED
  ##
  all_results <- eventReactive(input$runAnalysis, {
    
    ## ESTABLISH PROGRESS BAR
    withProgress(message = 'Performing', value = 0, {
      
      results_out = list()
      
      ## USE THIS TO DEBUG
      #browser()
      
      # GET UNITS FOR TLS AND TAPE
      #results_out$units = input$units
      #results_out$units_tape = input$units_tape
      
      ## GET NAMES FOR THE DATASETS
      dataname_base <- input$dataname_base
      dataname_test <- input$dataname_test
      
      ## UNIT STATEMENTS
      results_out$TLS_UnitStatement = if(input$units=='mm') {
        paste("Coordinates were reported in mm.")
      } else {paste("Coordinates have been converted from ", input$units, "to mm.")
      }
      
      ## READ IN THE TWO DATA FILES AND CONVERT TO mm, IF NECESSARY
      data1 = Zc_to_mm(data_p1base(), units = input$units)$Zc
      data2 = Zc_to_mm(data_p1test(), units = input$units)$Zc
      
      ## DATA SUMMARIES
      nTP_base <- TLS_getTandP(data1)
      nTP_test <- TLS_getTandP(data2)
      results_out$basedata_summary <- paste("The", dataname_base, "data are coordinates from", nTP_base$nTargets, "targets 
                                            measured from", nTP_base$nPositions, "TLS positions.")
      results_out$testdata_summary <- paste("The", dataname_test, "data are coordinates from", nTP_test$nTargets, "targets 
                                            measured from", nTP_test$nPositions, "TLS positions.")
      
      
      ## READ IN THE FOUR FILES AND CREATE THE Zc DATA SETS
      #data1 = create_Zc(data_p1base(), data_p2base(), data_p3base(), data_p4base(), units = input$units)$Zc
      #data2_list = create_Zc(data_p1test(), data_p2test(), data_p3test(), data_p4test(), units = input$units)
      #data2 = data2_list$Zc
      #data2_targetList = data2_list$target_list
      
      ## READ IN THE TAPE MEASURE DATA
      #Dtape = Dtape()
      ## If tape_units are not already mm, convert to mm (unit options are mm, cm or inches)
      #Dtape = Dtape_ToMM(dat = Dtape, units = input$units_tape)
      
      ## save the data in 'results_out'
      #results_out$data1 = data1
      #results_out$data2 = data2
      

      ##*********************
      ## PART II CALCULATIONS
      ##*********************
      ## Step 1:
      ## transform the Cartesian coordinates into spherical residuals
      incProgress(1/4, detail = paste("rigid body transformation on baseline data"))
      results_out$R1 = Zc_to_R(data1)
      incProgress(2/4, detail = "rigid body transformation on test data")
      results_out$R2 = Zc_to_R(data2)
      
      ## Step 2:
      ## Transform the residuals from wide to long
      results_out$X1 = RtoX(results_out$R1)
      results_out$X2 = RtoX(results_out$R2)
      
      ## BIVARIATE DATA (NO RANGING RESIDUALS)
      results_out$X1_bi = results_out$X1[,1:2]
      results_out$X2_bi = results_out$X2[,1:2]
      

      ## Step 3:
      ## Perform the hypothesis test
      ## Step 3: test equality of the covariance matrices
      incProgress(3/4, detail = "Part II analysis")
      p2_results = TLS_cov_check(results_out$X1, results_out$X2, conf.level = .99)
      results_out$test_results = p2_results$results
      results_out$p2_conclusion = p2_results$conclusion
      results_out$p2_interpretation = p2_results$interpretation
      
      ##
      ## BIVARIATE RESULTS
      p2_results_bi = TLS_cov_check(results_out$X1_bi, results_out$X2_bi, conf.level = .99)
      results_out$test_results_bi = p2_results_bi$results
      results_out$p2_conclusion_bi = p2_results_bi$conclusion
      results_out$p2_interpretation_bi = p2_results_bi$interpretation
      
      ## Step 4:
      ## Combine the data in preparation for plotting the data ellipses
      Xcomb <- rbind(data.frame(results_out$X1, dataset = dataname_base), 
                     data.frame(results_out$X2, dataset = dataname_test))
      ## Need to remove any NA values
      if(sum(is.na(Xcomb[,1]))>0) {
        Xcomb <- Xcomb[-c(which(is.na(Xcomb[,1]))),]
      }
      results_out$Xcomb = Xcomb
      
      ## BIVARIATE
      Xcomb_bi <- rbind(data.frame(results_out$X1_bi, dataset = dataname_base), 
                     data.frame(results_out$X2_bi, dataset = dataname_test))
      ## Need to remove any NA values
      if(sum(is.na(Xcomb_bi[,1]))>0) {
        Xcomb_bi <- Xcomb_bi[-c(which(is.na(Xcomb_bi[,1]))),]
      }
      results_out$Xcomb_bi = Xcomb_bi
      
      
      ## Now estimate the covariance matrices of the spherical residuals
      results_out$SigmaHat1 = cov(results_out$X1)
      results_out$SigmaHat2 = cov(results_out$X2)
      
      ##
      ## CHECK FOR BAD TARGETS -- punt over the 'R_pretty' datasets, and render the plots in the output section
      results_out$Rpretty_base <- create_pretty_R(results_out$R1, data_name = dataname_base)
      results_out$Rpretty_test <- create_pretty_R(results_out$R2, data_name = dataname_test)

      return(results_out)
    })
    #browser()
  })
  

  #############################
  ##** PART II ANALYSIS
  ##*
  ## The string of text stating the conclusion
  output$p2_conclusion <- renderText({
    if(is.null(all_results()$p2_conclusion)) {
      return(NULL)
    }
    return(all_results()$p2_conclusion)
  })
  
  output$p2_conclusion_bi <- renderText({
    if(is.null(all_results()$p2_conclusion_bi)) {
      return(NULL)
    }
    return(all_results()$p2_conclusion_bi)
  })
  
  output$p2_interpretation <- renderText({
    if(is.null(all_results()$p2_interpretation)) {
      return(NULL)
    }
    return(all_results()$p2_interpretation)
  })
  
  ## The data summary statements
  output$units <- renderText({
    if(is.null(all_results()$TLS_UnitStatement)) {
      return(NULL)
    }
    return(paste(all_results()$TLS_UnitStatement))
  })
  
  output$basedata_summary <- renderText({
    if(is.null(all_results()$basedata_summary)) {
      return(NULL)
    }
    return(paste(all_results()$basedata_summary))
  })
  
  output$testdata_summary <- renderText({
    if(is.null(all_results()$testdata_summary)) {
      return(NULL)
    }
    return(paste(all_results()$testdata_summary))
  })
  
  ## The information from the Hotellings T2 test
  output$p2_analysis <- renderPrint({
    if(is.null(all_results()$test_results)) {
      return(NULL)
    }
    return(all_results()$test_results)
  })
  
  output$p2_analysis_bi <- renderPrint({
    if(is.null(all_results()$test_results_bi)) {
      return(NULL)
    }
    return(all_results()$test_results_bi)
  })
  
  ## The estimated variance/covariance matrices
  output$covMat_info <- renderUI({
    if(is.null(all_results()$SigmaHat1)) {
      return(NULL)
    }
    ## state some info on the estimated covariance matrices 
    HTML(paste("theta = azimuth angle (arcsec)<br>", 
               "phi = polar angle (arcsec) <br>",
               "r = ranging direction (mm)"))
  })
  
  output$SigmaHat1 <- renderPrint({
    if(is.null(all_results()$SigmaHat1)) {
      return(NULL)
    }
    return(round(all_results()$SigmaHat1,3))
  })
  
  output$SigmaHat2 <- renderPrint({
    if(is.null(all_results()$SigmaHat2)) {
      return(NULL)
    }
    return(round(all_results()$SigmaHat2,3))
  })
  
  ###############################
  ##**PART II - DATA ELLIPSES
  output$ellipsePlot <- renderPlot({
    if(is.null(all_results()$Xcomb)) {
      return(NULL)
    }
    #ggplot(all_results()$plot_data, aes(x=x1, y=x2))+ 
    #  geom_point()
    covEllipses(all_results()$Xcomb[,1:3], as.factor(all_results()$Xcomb[,4]), variables=1:3, pooled = FALSE, 
                label.pos = c("center", "top"), col = c("darkred", "dodgerblue"))
  }, res = 96)
  
  ################################
  ## CHECK FOR BAD TARGETS
  ## 1) generate the lsit of plots
  ## 2) plot the plots
  output$resPlot_base <- renderPlot({
    if(is.null(all_results()$Rpretty_base)) {
      return(NULL)
    }
    plot_my_residuals(all_results()$Rpretty_base)
  }, res = 96)
  
  output$resPlot_test <- renderPlot({
    if(is.null(all_results()$Rpretty_test)) {
      return(NULL)
    }
    plot_my_residuals(all_results()$Rpretty_test)
  })
  
  output$residualStatement <- renderUI({
    HTML(paste("Angular residuals are reported in arcseconds. Ranging residuals are reported in mm.<br>", 
               "Target names are derived from row position in the respective dataset."
    ))
  })
  
  ###############################
  ##**Information
  ##*
  output$info <- renderUI({
    HTML(paste("This web application tool implements the required calculations 
    to perform the OSAC CSIR `Terrestrial Laser Scanner Interim Performance Assessment 
    Test Procedure - Part II. This procedure tests if a TLS instrument's precision has 
    significantly changed since a baseline state.<br><br>",
    "This tool requires specific data input to function correctly:<br><br>", 
    "-- Two sets of data are required: ‘Baseline’ and ‘Test’.<br><br>",
    "-- Each set of data must be Cartesian coordinates [X,Y,Z] corresponding 
    to target centers collected from a target array from multiple TLS positions.<br><br>",
    "-- Rows in the data files must correspond to target.<br><br>",
    "-- Each triplicate set of columns must correspond to the Cartesian coordinates 
    [X,Y,Z] for the given target.<br><br>",
    "-- The dimension of the data file must correspond to the number of targets and number of 
    TLS positions, e.g., if 20 targets are measured from four TLS positions, 
    the corresponding data file must have 20 rows and 12 columns.<br><br>",
    "--	Data must be input as .txt. or .csv files."
    ))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
