shinyaltrap<-function(){
##if (!require("DT")) install.packages('DT')
#xfun::session_info('DT')
#library(DT)
library(shiny)
library(shinybusy)
library(plyr)
library(poibin)
library(rstanarm)
library(ResourceSelection)
library(ggplot2)
library(DescTools)
library(Rcpp)
  #load the direct transfer data located in library file inst/directfile
pgkPath<-path.package("altrap",quiet = FALSE)
load(paste0(pgkPath,"/directfile"))
#PGfile<-read.csv("DirectTransfer.csv")
#PGfile<-read.csv("G:/Dropbox/ElidaBayesNetCorrectData/DirectTransfer.csv")
#OffsetData<-read.csv("Pareto.csv")
#OffsetData<-read.csv("G:/Dropbox/ElidaBayesNetCorrectData/Pareto.csv")
# Define UI for slider demo app ----
#OffsetData<-as.data.frame(OffsetData)
ui <- fluidPage(
  #Busy indicator
  add_busy_spinner(spin = "cube-grid",color = "#1D6FB3",position=c("full-page"),height = "120px", width = "120px"),

  # App title ----
  titlePanel("Activity Level Recovery, Transfer and Persistence (ALTRaP) Program"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar to demonstrate various slider options ----
    sidebarPanel(

      # Input: Simple integer interval ----
      sliderInput("integer", "Direct Transfer Time:",
                  min = 0, max = 30,
                  value = 10),


      # Input: Specification of range within an interval ----
      sliderInput("range", "Secondary Transfer Time",
                  min = 0, max = 30,
                  value = c(10,15)),

      # Input: Specification of range within an interval ----
      sliderInput("bckgrnd", "Probability of background",
                  min = 0.01, max = 0.99,
                  value = c(0.36)),
      helpText("Default value=0.36"),

      radioButtons("radio", h3("Select Innocent Direct"),
                   choices = list("Include Direct Innocent Transfer" = 1, "no Innocent Direct Transfer"=2),selected = 2),


      conditionalPanel(  #Innocent transfer time


        condition = "input.radio == '1'",
        sliderInput("innocent", "Innocent Transfer", min=0, max=30, value=0),
      ),

      fluidRow(

        column(5,offset = 1,
               checkboxInput("checkbox2", label = ("Select random seed 101 to reproduce calculations exactly"), value = TRUE)),###set seed to 101 if selected
        hr(),
        column(6, align="center",
               actionButton("goButton", "Go!")),

        hr(),
        helpText("If sliders are altered, the Go!' button must be pressed to update changes"),
        hr(),
        hr(),

      ),
      ##############
      #################
      fluidRow(
        column(6,
               selectInput("select_n", label = ("Select no. of contacts per hour (n) for sec Transfer"),
                           choices = list("n=1" = 1, "n=2" = 2, "n=3" = 3, "Choice 3" = 3, "n=4"=4),
                           selected = 1),
               helpText("Can be used to update n in all tables and plots without pressing 'Go!'"),
               hr(),
        ),


        column(6,
               selectInput("select_x", label = ("Select logistic regression decision threshold (x)"),
                           choices = list("x=1" = 1, "x=2" = 2, "x=3" = 3, "x=4" = 4, "x=6" = 6, "x=8"=8),
                           selected = 3),
               helpText("Can be used to update x in all tables and plots of 'Sensitivity tab' without pressing 'Go!'.
                        If sliders are altered, then 'Go!' button must be pressed first"),
               hr(),
        ),

        #verbatimTextOutput("value"))),
        # helpText("Select number of contacts per hour (n)"),
        #  ),
        #################
        ######################

        #verbatimTextOutput("valuec"))
      )
    ),
    ###########
    ##########
    ###### Main panel for displaying outputs ----
    mainPanel(
      ##########################################################
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("POI only",
                           ###################################################
                           # Output: Table summarizing the values entered ----

                           fluidRow(
                             tableOutput("values"),
                           ),
                           fluidRow(
                             column(3,offset=1,
                                    h4("Median LRs for values of x (logistic regression decision threshold) and n (no of contacts per hour", align="center"),

                                    tableOutput("Results")
                             ),

                             ###################
                             # column(4,offset=2,
                             #       h4("Median LRs for values of x (logistic regression decision threshold) and n (no of contacts per hour", align="center"),

                             #tableOutput("ResultsT")
                             #),



                             #######################
                             column(3,offset=1,
                                    h4("Tabulated Quantiles of log10LR illustrated in the figure below", align="center"),
                                    tableOutput("Results2")
                             ),

                             ###########
                             column(3,offset=1,
                                    checkboxInput("checkbox", label = ("Select log-scale for quantile table"), value = TRUE),
                                    hr(),
                                    verbatimTextOutput("input$select_n")
                             )
                             ############
                           ),
                           fluidRow(
                             plotOutput("plot2"),
                             p("A total of 4000 log_{10}LRs simulated from logistic regression coefficients
        using T, S and n. For each value
        of logistic regression decision threshold x, a density (violin) plot is shown in
        white. Superimposed is a box-whiskers plot in green, and behind, the blue
        rectangle delimits 0.05 and 0.95 percentiles, whereas the red rectangle
        delimits 0.025 and 0.975 percentiles.")
                             # DT::dataTableOutput("Results")
                           ),
                  ),
                  ###Insert tab here
                  tabPanel("POI and Unknown" ,


                           fluidRow(
                             tableOutput("valuesT"),####CHANGED from pA
                           ),

                           fluidRow(
                             column(3,offset=1,
                                    h4("Median LRs for values of x (logistic regression decision threshold) and n (no of contacts per hour", align="center"),

                                    tableOutput("ResultsT")
                             ),

                             column(3, offset=1,
                                    h4("Tabulated Quantiles of log10LR illustrated in the figure below", align="center"),
                                    tableOutput("Results2T")
                             ),

                             column(3, offset=1,
                                    checkboxInput("checkbox3", label = ("Select log-scale for quantile table"), value = TRUE),
                                    hr(),
                                    #verbatimTextOutput("input$select_n")
                             )
                             ),

                             fluidRow(
                               plotOutput("plot2T"),
                               p("A total of 4000 log_{10}LRs simulated from logistic regression coefficients
        using T, S and n. For each value
        of logistic regression decision threshold x, a density (violin) plot is shown in
        white. Superimposed is a box-whiskers plot in green, and behind, the blue
        rectangle delimits 0.05 and 0.95 percentiles, whereas the red rectangle
        delimits 0.025 and 0.975 percentiles.")
                               # DT::dataTableOutput("Results")
                             ),



                           #),
                           ################################Leave brackets below
                  ),
                  tabPanel("Sensitivity" , ####Tables
                           fluidRow(
                             tableOutput("valuesX"),
                           ),

                           fluidRow(
                             column(5,
                                    h4("Tabulated Quantiles of LRs (POI only)", align="center"),

                                    tableOutput("Results2sen"),
                                    p("Columns are quantiles corresponding to LRs shown in figure below. Rows are quantiles corresponding to probability of dropout based on the distribution of
                                      bootstrapped samples - represented as Pr(S) in the figure below, corresponding to 0.5, 0.75, 0.9, 0.95, 0.975, 0.99 quantiles ")
                             ),

                             column(5,offset=2,
                                    h4(" Tabulated Quantiles of LRs (POI and Unknown(s))", align="center"),

                                    tableOutput("Results2Tsen"),
                                    p("Columns are quantiles corresponding to LRs shown in figure below. Rows are quantiles corresponding to probability of dropout based on the distribution of
                                      bootstrapped samples - represented as Pr(S) in the figure below, corresponding to 0.5, 0.75, 0.9, 0.95, 0.975, 0.99 quantiles ")
                             )#column bracket

                             ),#fluid row bracket
##checkbox############
fluidRow(
  column(6,offset=6,
       checkboxInput("checkbox4", label = ("Select log-scale for quantile table"), value = TRUE),
       hr(),
       #verbatimTextOutput("input$select_n")
)
                  ),


####Plots sensitivity tab
                           fluidRow(
                             column(5,offset=1,
                                    h4("POI only", align="center"),

                                    plotOutput("plot2sen")
                              ),

                              column(5,offset=1,
                              h4("POI and Unknown", align="center"),

                              plotOutput("plot2Tsen")

                              )
                             ),


###############################################################################
                  )#tabpanel sen bracket
      )#tabpanel set bracket
    )
  )
)
# Define server logic for slider examples ----
server <- function(input, output,session) {

  output$valuec <- renderPrint({ input$checkbox })
  output$value <- renderPrint({ input$select_n })
  #output$valueSen <- renderPrint({ input$select_x })

  # Reactive expression to create data frame of all input values ----

  sliderValues <- eventReactive(input$goButton,{

    data.frame(
      Event = c("Direct Transfer Time",
                "Secondary Transfer Time","Innocent Transfer Time","Pr. Background"),
      Time  = as.character(c(input$integer,
                             paste(input$range, collapse = " "),input$innocent,input$bckgrnd
      )),
      stringsAsFactors = FALSE)

  })

  ##Can be used to update
  #observeEvent(input$innocent,{
  #updateNumericInput(session, "innocent", value = 5)

  # })




  # Show the values in an HTML table ----
  output$values <- renderTable({
    sliderValues()
  })

  output$valuesT <- renderTable({
    sliderValues()
  })

  output$valuesX <- renderTable({
    sliderValues()
  })

  #run both StanRes and StanSens main functions
  observeEvent(input$goButton,{
    x<-as.integer(input$select_x)
    if (input$checkbox2==TRUE){set.seed(101)} #set seed 101 if checkbox is true
    if (input$radio==1){
      b<<-StanRes(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd,input$innocent)
      sen<<-StanSens(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd,x,input$innocent)
    } else {
      b<<-StanRes(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd)#,input$innocent>0)
      sen<<-StanSens(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd,x)#,input$innocent>0)
    }
  })

 #senitivity page only for Stansen function
  observeEvent(input$select_x,{
    if (input$goButton==0)
      return()
    x<-as.integer(input$select_x)
    if (input$checkbox2==TRUE){set.seed(101)} #set seed 101 if checkbox is true
    if (input$radio==1){
      #b<<-StanRes(PGfile,OffsetData,input$range[1],input$range[2],input$integer,.36,input$innocent)
       sen<<-StanSens(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd,x,input$innocent)
    } else {
      #b<<-StanRes(PGfile,OffsetData,input$range[1],input$range[2],input$integer,.36)#,input$innocent>0)
      sen<<-StanSens(PGfile,input$range[1],input$range[2],input$integer,input$bckgrnd,x)#,input$innocent>0)
    }
  })
######################
#Make table for the LRs POI tab
  Results<-eventReactive(input$goButton,{
    Results<-b$LRmod
    Results<-as.data.frame(Results)
    LR<-c("x=1","x=2","x=3","x=4","x=6","x=8")
    Results<-cbind(LR,Results)
    rename(Results, c("V1"="n=1", "V2"="n=2", "V3"="n=3", "V4"="n=4"))
  })

  #Make table for LRs POI and Unknown tab
  ResultsT<-eventReactive(input$goButton,{
    ResultsT<-b$LRmodT
    ResultsT<-as.data.frame(ResultsT)
    LR<-c("x=1","x=2","x=3","x=4","x=6","x=8")
    ResultsT<-cbind(LR,ResultsT)
    rename(ResultsT, c("V1"="n=1", "V2"="n=2", "V3"="n=3", "V4"="n=4"))
  })


  #####Quantile tables for POI (Results2) and POI and U (Results2T)
  ######print log quantiles for the ggplot(plot2)
  Results2<-eventReactive(input$goButton,{
    Results2<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    if (input$checkbox==TRUE){
      for(z in 1:6){
        Results2[z,]=quantile(log10(b$LRtot[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2[z,]=quantile(b$LRtot[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2<-as.data.frame(Results2)
    LR<-c("x=1","x=2","x=3","x=4","x=6","n=8")
    Results2<-cbind(LR,Results2)
    rename(Results2, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })

  Results2<-reactive({
    if (input$goButton==0)
      return()
    Results2<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    if (input$checkbox==TRUE){
      for(z in 1:6){
        Results2[z,]=quantile(log10(b$LRtot[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2[z,]=quantile(b$LRtot[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2<-as.data.frame(Results2)
    LR<-c("x=1","x=2","x=3","x=4","x=6","n=8")
    Results2<-cbind(LR,Results2)
    rename(Results2, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })

  Results2T<-eventReactive(input$goButton,{ #second LR table
    Results2T<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    if (input$checkbox3==TRUE){
      for(z in 1:6){
        Results2T[z,]=quantile(log10(b$LRtog[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2T[z,]=quantile(b$LRtog[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2T<-as.data.frame(Results2T)
    LR<-c("x=1","x=2","x=3","x=4","x=6","n=8")
    Results2T<-cbind(LR,Results2)
    rename(Results2T, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })

  Results2T<-reactive({
    if (input$goButton==0)
      return()
    Results2T<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    if (input$checkbox3==TRUE){
      for(z in 1:6){
        Results2T[z,]=quantile(log10(b$LRtog[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2T[z,]=quantile(b$LRtog[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2T<-as.data.frame(Results2T)
    LR<-c("x=1","x=2","x=3","x=4","x=6","n=8")
    Results2T<-cbind(LR,Results2T)
    rename(Results2T, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })

  #SENSITIVITY TAB########################################

  #####Sensitivity table for POI tab
  ######print log quantiles for the ggplot(plot2)
  #Reactive with Go button
  Results2sen<-eventReactive(input$goButton,{
    Results2sen<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    # x<- as.integer(input$select_x)
    if (input$checkbox4==TRUE){#checkbox to decide if the data are logscale or not
      for(z in 1:6){
        Results2sen[z,]=quantile(log10(sen$LRtot[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2sen[z,]=quantile(sen$LRtot[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2sen<-as.data.frame(Results2sen)
    Quantile<-c("0.5","0.75","0.9","0.95","0.975","0.99")
    Results2sen<-cbind(Quantile,Results2sen)
    rename(Results2sen, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })
  ##Sensitivity table 'POI only' reactive to 'n' choice selection only
  Results2sen<-reactive({
    if (input$goButton==0)
      return()
    Results2sen<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    if (input$checkbox4==TRUE){
      for(z in 1:6){
        Results2sen[z,]=quantile(log10(sen$LRtot[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2sen[z,]=quantile(sen$LRtot[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2sen<-as.data.frame(Results2sen)
    Quantile<-c("0.5","0.75","0.9","0.95","0.975","0.99")
    Results2sen<-cbind(Quantile,Results2sen)
    rename(Results2sen, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })
  #Sensitivity table 'POI and Unknown' reactive to Go button
  Results2Tsen<-eventReactive(input$goButton,{ #second LR table
    Results2Tsen<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    if (input$checkbox4==TRUE){#checkbox for log scale
      for(z in 1:6){
        Results2Tsen[z,]=quantile(log10(sen$LRtog[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2Tsen[z,]=quantile(sen$LRtog[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2Tsen<-as.data.frame(Results2Tsen)
    Quantile<-c("0.5","0.75","0.9","0.95","0.975","0.99")
    Results2Tsen<-cbind(Quantile,Results2Tsen)
    rename(Results2Tsen, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })
  ##Sensitivity table 'POI and unknown' reactive to 'n' choice selection only
  Results2Tsen<-reactive({
    if (input$goButton==0)
      return()
    Results2Tsen<-matrix(0,6,8)
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    if (input$checkbox4==TRUE){
      for(z in 1:6){
        Results2Tsen[z,]=quantile(log10(sen$LRtog[z,n,]),probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }else{

      for(z in 1:6){
        Results2Tsen[z,]=quantile(sen$LRtog[z,n,],probs=c(.01,.025,.25,.5,.75,.95,.975,.99))
      }
    }
    Results2Tsen<-as.data.frame(Results2Tsen)
    Quantile<-c("0.5","0.75","0.9","0.95","0.975","0.99")
    Results2Tsen<-cbind(Quantile,Results2Tsen)
    rename(Results2Tsen, c("V1"="0.01","V2"="0.025", "V3"="0.25", "V4"="0.5", "V5"="0.75","V6"="0.95","V7"="0.975","V8"="0.99"))

  })


  #############################################################################################################################
  #OUTPUTS###################
  ##SENSITIVITY TAB
  output$Results2sen <- renderTable({
    Results2sen()
  })
  ##POI + unknown

  output$Results2Tsen <- renderTable({
    Results2Tsen()
  })

  #OUTPUTS###################



  output$Results <- renderTable({
    Results()
  })


  output$Results2 <- renderTable({
    Results2()
  })
  ##POI + unknown
  output$ResultsT <- renderTable({
    ResultsT()
  })

  output$Results2T <- renderTable({
    Results2T()
  })

  ##ggplot2 Tab 1
  plot2<-eventReactive(input$goButton,{
    n<- as.integer(input$select_n)
    violin(b$LRtot,n)
  })

  plot2<-reactive({
    if (input$goButton==0)
      return()
    n<- as.integer(input$select_n)
    violin(b$LRtot,n)
  })


  output$plot2 <- renderPlot({
    plot2()
  })


  ##ggplot2 commands Tab 2
  plot2T<-eventReactive(input$goButton,{
    n<- as.integer(input$select_n)
    violin(b$LRtog,n)
  })

  plot2T<-reactive({
    if (input$goButton==0)
      return()
    n<- as.integer(input$select_n)
    violin(b$LRtog,n)
  })

  output$plot2T <- renderPlot({
    plot2T()
  })


  #SENSITIVITY TAB #ggplot2 commands Tab3 SENS TAB
  #PoI only
  plot2sen<-eventReactive(input$goButton,{
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    violinSens(sen$LRtot,n,x,sen$emp.data)
  })

  plot2sen<-reactive({
    if (input$goButton==0)
      return()
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    violinSens(sen$LRtot,n,x,sen$emp.data)
  })

  output$plot2sen <- renderPlot({
    plot2sen()
  })

 #SENS TAB #POI and U ##################

  plot2Tsen<-eventReactive(input$goButton,{
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    violinSens(sen$LRtog,n)
  })

  plot2Tsen<-reactive({
    if (input$goButton==0)
      return()
    n<- as.integer(input$select_n)
    x<- as.integer(input$select_x)
    violinSens(sen$LRtog,n,x,sen$emp.data)
  })

  output$plot2Tsen <- renderPlot({
    plot2Tsen()
  })

}

# Create Shiny app ----
shinyApp(ui, server)

}


