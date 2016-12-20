library(shiny)

# Define UI
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("Dynamics Simulation of Honey Bee Population and Food Storage"),
  sidebarPanel(
    
    wellPanel(tags$h3("Axis"),
              sliderInput("tmaxday", "Duration (Days):", 
                          value=1000, min=365, max=4000, step=5),
              sliderInput("ymax", "Population (#):", 
                          value=50000, min=2000, max=100000, step=2000)),
    
    wellPanel(tags$h3("Toxicity"),
              sliderInput("d", "Chronic dose of imidacloprid(ug/L):",
                          min=100, max=5000, value=2300, step=100)),
    
    wellPanel(tags$h3("Initial State"),
              sliderInput("f", "Food storage:", 
                          min=1000, max=3000, value=3000, step=100),
              sliderInput("B", "Brood:", 
                          min=1000, max=10000, value=2000, step=500),
              sliderInput("H", "Hive Bees:", 
                          min=1000, max=10000, value=2000, step=500),
              sliderInput("F", "Foragers:",
                          min=100, max=5000, value=1000, step=500)),
    
    wellPanel(tags$h3("Seasonal Rate Parameters"),
              sliderInput("c0", "Food collection rate:", 
                          min=0.05, max=0.5, value=0.1, step=0.02),
              sliderInput("l0", "Laying rate:", 
                          min=1000, max=3000, value=2000, step=100),
              sliderInput("m0", "Mortality rate:", 
                          min=0.033, max=0.3, value=0.1, step=.02)),
              
    wellPanel(tags$h3("Rate Parameters"),
              sliderInput("ga", "Food consumption rate for adult bee:", 
                          min=0.005, max=0.015, value=0.007, step=0.001),
              sliderInput("gb", "Food consumption rate for brood:", 
                          min=0.015, max=0.025, value=0.018, step=0.001),
              sliderInput("amin", "Minimun transition rate:", 
                          min=0.1, max=0.4, value=0.25, step=.05),
              sliderInput("amax", "Maximum transition rate:", 
                          min=0.1, max=0.4, value=0.25, step=.05),
              sliderInput("fi", "Pupation rate:", 
                          min=0.05, max=0.5, value=0.1, step=0.05),
              sliderInput("s", "Social inhibition rate:", 
                          min=0.6, max=0.9, value=0.75, step=.05)), 
    
    wellPanel(tags$h3("Control constants"),
              sliderInput("b", "Food stores:", 
                          min=100, max=2000, value=500, step=100),
              sliderInput("v", "Hive bees:", 
                          min=1000, max=9000, value=5000, step=1000))),
  
  mainPanel(plotOutput("guessPlot", height = 600, width = 800)
  )))