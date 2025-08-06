# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

###############################
###### 1 Variable Module ######
###############################

setwd("/path/to/your/simde")
source("plot_functions.R")
source("model_formulas.R")
library(bslib)
library(bsicons)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)

# install.packages(c("stringr", "bslib", "bsicons", "shinyWidgets", "shinyjs", "shinycssloaders"))

one_var_UI <- function(id) {
  ns <- NS(id) # we are making a namespace to make ID names unique for every time this module is called so we are making a ns function that defines each input for this id
  # so e.g., if I call one_var_UI(id = "one_var"), then in this I'll make a ns("one_var") function that is now used for each of the input below.
  
  # tags$head(
  #   tags$style(HTML('
  #     /* CSS to change the color of the slider input */
  #     .irs-bar {
  #       background-color: green; /* Change the color of the slider track */
  #     }
  #     .irs-handle {
  #       background-color: blue; /* Change the color of the slider handle */
  #       border-color: black; /* Optional: Change the border color of the handle */
  #     }
  #   '))
  # )
  
  # taglist() to put all the input items here?
  tabPanel("Univariate",
           value = ns("one_var_system"),
           br(),
           fluidRow(
             column(1),
             column(3,
                    fluidRow(column(12, class = "panel-custom",
                          selectInput(inputId = ns("modeltype"), label = HTML("<font size='4'> Type of Model</font>"),
                                      choices = list("1st Order Models" = list("Linear ODE" = "l_ode",
                                                                               "Polynomial ODE" = "p_ode",
                                                                               "Linear SDE (OU Model)" = "ou",
                                                                               "Polynomial SDE" = "p_sde"),
                                                     "2nd Order Models" = list("Harmonic Oscillator (Undamped) ODE" = "harmonic_ode",
                                                                               "Damped Oscillator ODE" = "damp_ode",
                                                                               "Harmonic Oscillator (Undamped) SDE" = "harmonic_sde",
                                                                               "Damped Oscillator SDE" = "damp_sde")), 
                                      selected = "l_ode", width = "400px"), 
                          div(id = ns("param_input_1var"),
                              br(),
                              h3("Define Parameters", align = "center"),
                                numericInput(inputId = ns("mean"), label = HTML("<font size='4'>&#956;</font>")  |>
                                                                              tooltip(HTML("Name: Mean/Intercept <br> Range: (-&#8734;, +&#8734;)")), # If you want to include an icon, then you just do this: span(HTML("<font size='4'>&#956;</font>"), bs_icon("info-circle")) |> tooltip("Some message") 
                                                   value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                              conditionalPanel("input.modeltype == 'l_ode' | input.modeltype == 'ou'", ns = ns,
                                               
                                               numericInput(inputId = ns("beta"), label = HTML("<font size='4'>&#946;</font>") |>
                                                                                            tooltip(HTML("Name: Auto-Effect <br> Range: (-&#8734;, 0]")),
                                                            value = -0.8, min = -Inf, max = 0, step = 0.1, width = '300px')),
                              
                              conditionalPanel("input.modeltype == 'p_ode' | input.modeltype == 'p_sde'", ns = ns,
                                               
                                               numericInput(inputId = ns("beta1"), label = HTML("<font size='4'>&#946;<sub>1</sub></font>") |> 
                                                                                            tooltip(HTML("Name: Linear Auto-Effect <br> Range: (-&#8734;, 0]")),
                                                            value = -0.8, min = -Inf, max = 0, step = 0.1, width = '300px'),
                                               
                                               numericInput(inputId = ns("beta2"), label = HTML("<font size='4'>&#946;<sub>2</sub></font>") |>
                                                                                            tooltip(HTML("Name: Quadratic Auto-Effect <br> Range: (-&#8734;, 0]")),
                                                            value = 0, min = -Inf, max = 0, step = 0.1, width = '300px'),
                                               
                                               numericInput(inputId = ns("beta3"), label = HTML("<font size='4'>&#946;<sub>3</sub></font>") |> 
                                                                                            tooltip(HTML("Name: Cubic Auto-Effect <br> Range: (-&#8734;, 0]")),
                                                            value = 0, min = -Inf, max = 0, step = 0.1, width = '300px'),
                                               
                                               numericInput(inputId = ns("beta4"), label = HTML("<font size='4'>&#946;<sub>4</sub></font>") |> 
                                                                                            tooltip(HTML("Name: Bi-quadratic Auto-Effect <br> Range: (-&#8734;, 0]")),
                                                            value = 0, min = -Inf, max = 0, step = 0.1, width = '300px')),
                              
                              conditionalPanel("input.modeltype == 'harmonic_ode' | input.modeltype == 'damp_ode' | input.modeltype == 'harmonic_sde' | input.modeltype == 'damp_sde'", ns = ns,
                                               
                                               numericInput(inputId = ns("omega"), label = HTML("<font size='4'>&#969;</font>") |>
                                                                                            tooltip(HTML("Name: Frequency <br> Range: (-&#8734;, 0]")),
                                                            value = -2, min = -Inf, max = 0, step = 1, width = '300px')),
                              
                              conditionalPanel("input.modeltype == 'damp_ode' | input.modeltype == 'damp_sde'", ns = ns,
                                               
                                               numericInput(inputId = ns("zeta"), label = HTML("<font size='4'>&#950;</font>") |>
                                                                                            tooltip(HTML("Name: Damping <br> Range: (-&#8734;, +&#8734;)")),
                                                            value = 0.1, min = -Inf, max = +Inf, step = 0.1, width = '300px')),
                              
                              conditionalPanel("input.modeltype == 'ou' | input.modeltype == 'p_sde' | input.modeltype == 'harmonic_sde' | input.modeltype == 'damp_sde'", ns = ns,
                                               
                                               numericInput(inputId = ns("sigma"), label = HTML("<font size='4'>&#963;</font>") |>
                                                                                            tooltip(HTML("Name: White Noise Scaling SD <br> Range: [0, +&#8734;)")),
                                                            value = 0, min = 0, max = +Inf, step = 1, width = '300px')),
                              
                              numericInput(inputId = ns("start_value"), label = HTML("<font size='4'>y<sub>0</sub></font>")
                                                                                        |> tooltip(HTML("Name: Initial Value of y at time t<sub>0</sub> <br> Range: (-&#8734;, +&#8734;)")),
                                           value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                              
                              actionButton(ns("reset_params"), "Reset Parameters", width = '300px')
                              ),
                          conditionalPanel("input.modeltype == 'ou' | input.modeltype == 'p_sde' | input.modeltype == 'harmonic_sde' | input.modeltype == 'damp_sde'", ns = ns,
                                           br(),
                                           numericInput(inputId = ns("seed"), label = "Set Your Seed" |>
                                                                                            tooltip(HTML("Name: Seed for Random Number Generator <br> Range: {1, 2, ...}")),
                                                        value = 1, min = 1, max = .Machine$integer.max, step = 1, width = '300px'),
                                           
                                           actionButton(ns("random_seed"), "Choose a Random Seed", width = '300px')),
                    )),
                    fluidRow(column(12, class = "panel-custom",
                            h3("Define Time", align = "center"),
                            sliderInput(ns("period"), label = HTML("<font size='4'> Time Period </font>") |>
                                                                  tooltip(HTML("Use the slider to choose the time period of simulation: (t<sub>0</sub>, T)")), min = 0, max = 300,
                                        value = c(0, 20)),
                            numericInput(inputId = ns("step_size"), label = HTML("<font size='4'>&#916;t</font>") |>
                                                                              tooltip(HTML("Name: Time Step <br> Choose an appropriate time step to iterate the process in the specified Time Period.")),
                                         value = 0.1, min = 0, max = 10, step = 0.1, width = '300px'),
                            h5(uiOutput(ns("n_timepoints")), style = "font-weight:bold")
                           )),
             ),
             column(1),
             column(6, 
                    h4("Equation of the System"),
                    uiOutput(ns('formula')),
                    br(),
                    h4("Plot of Y Over Time", align = "center", style = "font-weight:bold"),
                    br(),
                    plotOutput(outputId = ns("ode_plot")) |> withSpinner(type = 5, color = "#0044A4", size = 1, hide.ui = FALSE),
                    actionButton(ns("plot_button"), "Generate Plot", class = "plot-btn")),
             column(1)
           )
  )
}

one_var_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      observeEvent(input$order1, {
        updateSelectInput(session, 'modeltype', choices = list("Linear ODE" = "l_ode",
                                                               "Polynomial ODE" = "p_ode",
                                                               "Linear SDE (OU Model)" = "ou",
                                                               "Polynomial SDE" = "p_sde"))
        
      })
      
      observeEvent(input$order2, {
        updateSelectInput(session, 'modeltype', choices = list("Harmonic Oscillator (Undamped) ODE" = "harmonic_ode",
                                                               "Damped Oscillator ODE" = "damp_ode",
                                                               "Harmonic Oscillator (Undamped) SDE" = "harmonic_sde",
                                                               "Damped Oscillator SDE" = "damp_sde"))
        
        
      })
      
      observeEvent(list(input$beta, input$omega, input$sigma), {
        
        if(input$beta > 0 & is.numeric(input$beta)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#946; is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
        
        if(input$omega > 0 & is.numeric(input$omega)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#969; is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
        
        if(input$sigma < 0 & is.numeric(input$sigma)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#963; is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
      })
    
      
      output$formula <- renderUI({
        withMathJax(model_formula_univariate(input$modeltype))
      })
      
      reactive_period <- reactive({input$period})
      reactive_step_size <- reactive({input$step_size})
      
      output$n_timepoints <- renderUI({
        n_obs <- floor((reactive_period()[2] - reactive_period()[1]) / reactive_step_size()) + 1
      
        if(is.finite(n_obs)){
          sprintf("Number of Timepoints = %i", n_obs)
        } else{
          HTML("Choose an appropriate value for &#916;t")
        }
      })
      
      observeEvent(input$step_size,{
        n_obs <- floor((reactive_period()[2] - reactive_period()[1]) / reactive_step_size()) + 1
        
        if(n_obs <  35){
          showNotification("Your number of timepoints might be quite low. Adjust time parameters to increase this!",
                           duration = 10)
        }
      })
      
      observeEvent(input$period,{
        if(input$period[1] != 0){
          updateSliderInput(
            session,
            "period", value = c(0, input$period[2])
          )
        }
      })
      
      
      observeEvent(input$random_seed, {
        new_seed <- sample.int(1e6, 1)
        updateNumericInput(session, "seed", value = new_seed)
      })
      
      observeEvent(input$reset_params,{
        
        shinyjs::reset("param_input_1var")
        
      })
      
        
      plot_data <- eventReactive(input$plot_button, {
        model_name <- input$modeltype
        
        validate(need(input$step_size > 0, "Please check your inputs!"))
        
        if(model_name == "l_ode"){
          
          plot_ode(input)
          
        }
        else if(model_name == "p_ode"){
          
          plot_pde(input)
          
        } else if(model_name == "ou"){
          
          plot_ou(input)
          
        } else if(model_name == "p_sde"){
          
          plot_p_sde(input)
          
        } else if(model_name == "harmonic_ode"){
          
          plot_harmonic_ode(input)
          
        } else if(model_name == "damp_ode"){
          
          plot_damp_ode(input)
          
        } else if(model_name == "harmonic_sde"){
          
          plot_harmonic_sde(input)
          
        } else if(model_name == "damp_sde"){
          
          plot_damp_sde(input)
        }
      })
      
      output$ode_plot <- renderPlot({
        plot_data()
      })
    }
  )
}


# This is not working!
# observeEvent(input$modeltype, {
#   model_name <- input$modeltype
#   if(model_name == "l_ode"){
#     text = HTML("Gist:1st order changes in Y are bring predicted by mean centered Y with an auto-effect coefficient of &#946;")
#     
#     } else if(model_name == "p_ode"){
#     text = "some"
#     } else if(model_name == "ou"){
#     text = "Something"
#     } else if(model_name == "p_sde"){
#       text = "Something"
#     }
#   update_tooltip("formula_info", text)
# })