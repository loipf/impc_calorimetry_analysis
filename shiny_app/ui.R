# shows table and plot of KO mouse

library(shiny)
library(DT)
library(plotly)

ui = fluidPage(
  
  wellPanel(
    fluidRow(
      column(4,selectInput("selectInst", label = "select institute", 
                           choices = list('HMGU' = 'hmgu_', 'UC Davis'='ucd_', 'BCM'='bcm_','ICS'='ics_','KMPC'='kmpc_','MARC'='marc_','MRC Harwell'='mrc_', 'RBRC'='rbrc_','TCP'='tcp_'), 
                           selected = 1) ),
    column(4,selectInput("selectMeas", label = "select measurement", 
                           choices = list('RER smoothed'='RER_smooth_', 'RER all' = 'RER_derived_', "O2" = 'O2_', "CO2" = 'CO2_','heat'='heat_'),  #,'activity'='total_activity_','food intake'='cumulative_food_intake_'), 
                           selected = 1) ),
    column(4,selectInput("selectFeature", label = "select max/min/mean/amp", 
                choices = list('max' = 'max', "min" = 'min', "mean" = 'mean','amp'='amp' ), 
                selected = 1) )
    ),
      DT::dataTableOutput('tableKOmice')
  ),
  
  wellPanel(
    fluidRow(
      column(4, selectInput("selectIC", label = "select measurement", 
                            choices = list('RER all' = 'RER_derived_', "O2" = 'O2_', "CO2" = 'CO2_','heat'='heat_','RER smoothed'='RER_smooth_'),  #,'activity'='total_activity_','food intake'='cumulative_food_intake_'), 
                            selected = 1)),
      column(6,
        p('black: control of experiment,  red: KO mouse,  green: regression line of both groups'),
        p('p-vales of linear regression model') )
    ),
    fluidRow(
      column(3, plotOutput('plotKOmice_min')),
      column(3, plotOutput('plotKOmice_max')),
      column(3, plotOutput('plotKOmice_mean')),
      column(3, plotOutput('plotKOmice_amp'))
    )
  )
  ,
  
  wellPanel(
    h4('raw data only night overall overview:'),
    fluidRow(
      column(3, selectInput('selectAllSex', label = 'select sex', choices = list('male'='male', 'female'='female','both'='both'), selected=1)   ),
      column(3, selectInput("selectAllFeature", label = "select max/min/mean/amp", choices = list('max' = 'max', "min" = 'min', "mean" = 'mean','amp'='amp' ), selected = 1) ),
      column(3, selectInput('selectAllControl', label = 'select control', choices = list('only KO'='ko', 'only control'='wt','both'='both'), selected=1)   )
    ),
    fluidRow(
      column(6, plotlyOutput('plotlyAll')),
      column(6, plotOutput('plot_singleMouse') )
      # column(6, tags$iframe(style="height:400px; width:100%; scrolling=yes",
      #                       src="hmgu_RER_dataPlots.pdf"))
      )
  ),
  
  wellPanel(
    actionButton("buttonShowGene", label = "show gene and influence"),
    verbatimTextOutput("textShowGene")
  )

  
  
  ) # page end
  




















