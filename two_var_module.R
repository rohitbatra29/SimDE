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

################################
###### 2 Variables Module ######
################################


setwd("/path/to/your/simde")

source("plot_functions.R")
source("model_formulas.R")

two_var_UI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Bivariate",
           value = "2_var",
           br(),
           fluidRow(
             column(6, style = "padding-left:2em",
                    fluidRow(
                      column(5, class = "panel-custom",
                             h4("Variable X", style = "font-weight:bold"),
                             selectInput(inputId = ns("modeltype_x"), label = "Type of model",
                                         choices = list("1st Order Models" = list("Linear ODE" = "l_ode",
                                                                                  "Linear SDE (OU Model)" = "ou"),
                                                        "2nd Order Models" = list("Harmonic Oscillator (Undamped) ODE" = "harmonic_ode",
                                                                                  "Damped Oscillator ODE" = "damp_ode",
                                                                                  "Harmonic Oscillator (Undamped) SDE" = "harmonic_sde",
                                                                                  "Damped Oscillator SDE" = "damp_sde")), 
                                         selected = "l_ode"),
                             div(id = ns("x_params_inputs"),
                                 br(),
                                 h4("Define Parameters of X"),
                                 numericInput(inputId = ns("mean_x"), label = HTML("&#956;<sub>x</sub>") |>
                                                                                tooltip(HTML("Name: Mean of X <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                                 
                                 conditionalPanel("input.modeltype_x == 'l_ode' | input.modeltype_x == 'ou'", ns = ns,
                                                  numericInput(inputId = ns("beta_x"), label = HTML("&#946;<sub>x</sub>") |>
                                                                                                    tooltip(HTML("Name: Auto-Effect of X <br> Range: (-&#8734;, 0]")),
                                                               value = -0.8, min = -Inf, max = 0, step = 0.1, width = '300px')),
                                 
                                 conditionalPanel("input.modeltype_x == 'harmonic_ode' | input.modeltype_x == 'damp_ode' | input.modeltype_x == 'harmonic_sde' | input.modeltype_x == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("omega_x"), label = HTML("&#969;<sub>x</sub>") |>
                                                                                                  tooltip(HTML("Name: Frequency of X <br> Range: (-&#8734;, 0]")),
                                                               value = -2, min = -Inf, max = 0, step = 1, width = '300px')),
                                 
                                 numericInput(inputId = ns("gamma_x"), label = HTML("&#947;<sub>x</sub>") |>
                                                                                    tooltip(HTML("Name: Coupling Effect of Y on X <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 0.1, width = '300px'),
                                 
                                 conditionalPanel("input.modeltype_x == 'damp_ode' | input.modeltype_x == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("zeta_x"), label = HTML("&#950;<sub>x</sub>") |>
                                                                                                   tooltip(HTML("Name: damping of X <br> Range: (-&#8734;, +&#8734;)")),
                                                               value = 0.1, min = -Inf, max = +Inf, step = 0.1, width = '300px')),
                                 
                                 conditionalPanel("input.modeltype_x == 'ou' | input.modeltype_x == 'harmonic_sde' | input.modeltype_x == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("sigma_x"), label = HTML("&#963;<sub>x</sub>") |>
                                                                                                    tooltip(HTML("Name: White Noise Scaling SD of X <br> Range: [0, +&#8734;)")),
                                                               value = 0, min = 0, max = +Inf, step = 1, width = '300px')),
                                 
                                 numericInput(inputId = ns("start_x"), label = HTML("x<sub>0</sub>") |>
                                                                                tooltip(HTML("Name: Initial Value of X at time t<sub>0</sub> <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                                 
                                 actionButton(ns("reset_x_params"), "Reset X's Parameters")
                                 )
                             ),
                      column(1),
                      column(5, class = "panel-custom",
                             h4("Variable Y", style = "font-weight:bold"), 
                             selectInput(inputId = ns("modeltype_y"), label = "Type of model",
                                         choices = list("1st Order Models" = list("Linear ODE" = "l_ode",
                                                                                  "Linear SDE (OU Model)" = "ou"),
                                                        "2nd Order Models" = list("Harmonic Oscillator (Undamped) ODE" = "harmonic_ode",
                                                                                  "Damped Oscillator ODE" = "damp_ode",
                                                                                  "Harmonic Oscillator (Undamped) SDE" = "harmonic_sde",
                                                                                  "Damped Oscillator SDE" = "damp_sde")), 
                                         selected = "l_ode"),
                             div(id = ns("y_params_inputs"),
                                 br(),
                                 h4("Define Parameters of Y"),
                                 numericInput(inputId = ns("mean_y"), label = HTML("&#956;<sub>y</sub>") |>
                                                                                 tooltip(HTML("Name: Mean of Y <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                                 
                                 conditionalPanel("input.modeltype_y == 'l_ode' | input.modeltype_y == 'ou'", ns = ns,
                                                  numericInput(inputId = ns("beta_y"), label = HTML("&#946;<sub>y</sub>") |>
                                                                                                tooltip(HTML("Name: Auto-Effect of Y <br> Range: (-&#8734;, 0]")),
                                                               value = -0.8, min = -Inf, max = 0, step = 0.1, width = '300px')),
                                 
                                 conditionalPanel("input.modeltype_y == 'harmonic_ode' | input.modeltype_y == 'damp_ode' | input.modeltype_y == 'harmonic_sde' | input.modeltype_y == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("omega_y"), label = HTML("&#969;<sub>y</sub>") |>
                                                                                                    tooltip(HTML("Name: Frequency of Y <br> Range: (-&#8734;, 0]")),
                                                               value = -2, min = -Inf, max = 0, step = 1, width = '300px')),
                                 
                                 numericInput(inputId = ns("gamma_y"), label = HTML("&#947;<sub>y</sub>") |>
                                                                                  tooltip(HTML("Name: Coupling Effect of X on Y <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 0.1, width = '300px'),
                                 
                                 conditionalPanel("input.modeltype_y == 'damp_ode' | input.modeltype_y == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("zeta_y"), label = HTML("&#950;<sub>y</sub>") |>
                                                                                                  tooltip(HTML("Name: damping of Y <br> Range: (-&#8734;, +&#8734;)")),
                                                               value = 0.1, min = -Inf, max = +Inf, step = 0.1, width = '300px')),
                                 
                                 conditionalPanel("input.modeltype_y == 'ou' | input.modeltype_y == 'harmonic_sde' | input.modeltype_y == 'damp_sde'", ns = ns,
                                                  numericInput(inputId = ns("sigma_y"), label = HTML("&#963;<sub>y</sub>") |>
                                                                                                   tooltip(HTML("Name: White Noise Scaling SD of Y <br> Range: [0, +&#8734;)")),
                                                               value = 0, min = 0, max = +Inf, step = 1, width = '300px')),
                                 
                                 numericInput(inputId = ns("start_y"), label = HTML("y<sub>0</sub>") |>
                                                                                        tooltip(HTML("Name: Initial Value of Y at time t<sub>0</sub> <br> Range: (-&#8734;, +&#8734;)")),
                                              value = 0, min = -Inf, max = +Inf, step = 1, width = '300px'),
                                 
                                 actionButton(ns("reset_y_params"), "Reset Y's Parameters")
                                 )
                             ),
                      column(1)),
                    fluidRow(
                      column(12, 
                             conditionalPanel("input.modeltype_x == 'ou' | input.modeltype_x == 'harmonic_sde' | input.modeltype_x == 'damp_sde' |
                                              input.modeltype_y == 'ou' | input.modeltype_y == 'harmonic_sde' | input.modeltype_y == 'damp_sde'", ns = ns,
                                              numericInput(inputId = ns("seed"), label = "Set Your Seed" |>
                                                                                            tooltip(HTML("Name: Seed for Random Number Generator <br> Range: {1, 2, ...}")),
                                                           value = 1, min = 1, max = .Machine$integer.max, step = 1, width = '300px'),
                                              actionButton(ns("random_seed"), "Choose a Random Seed", width = '300px')),
                             br(),
                             h4("Define Time"), align = "center",
                             sliderInput(ns("period"), label = "Time Period" |>
                                                                  tooltip(HTML("Use the slider to choose the time period of simulation: (t<sub>0</sub>, T)")), 
                                         min = 0, max = 300, value = c(0, 20)),
                             numericInput(inputId = ns("step_size"), label = HTML("&#916;t") |>
                                                                                      tooltip(HTML("Name: Time Step <br> Choose an appropriate time step to iterate the process in the specified Time Period.")),
                                          value = 0.1, min = 0, max = 10, step = 0.1, width = '300px'),
                             h5(uiOutput(ns("n_timepoints")), style = "font-weight:bold"))) 
                    
             ),
             column(6, style = "padding-right:2em",
                    br(),
                    selectInput(ns("xy_coupling"), "Define the Coupling Between X and Y",
                                list("Variable-Centered Coupling" = "coup_diff",
                                     "Mean-Centered Coupling" = "coup_mean"), selected = "coup_diff"),
                    h4("Equation of the System"),
                    uiOutput(ns('formula_xy')),
                    br(),
                    h4("Plot of X and Y Over Time", style = "font-weight:bold", align = "center"),
                    plotOutput(outputId = ns("bivariate_plot")) |> withSpinner(type = 5, color = "#0044A4", size = 1, hide.ui = FALSE),
                    actionButton(ns("plot_button"), "Generate Plot", class = "plot-btn")
                    )
           )
  )
}

two_var_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      listen_coupling <- reactive({
        list(input$modeltype_x, input$modeltype_y)
      })
      
      observeEvent(list(input$beta_x, input$omega_x, input$sigma_x, input$beta_y, input$omega_y, input$sigma_y), {
        if(input$beta_x > 0 & is.numeric(input$beta_x)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#946;<sub>x</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning", 
            html = TRUE
          )
        }
        
        if(input$beta_y > 0 & is.numeric(input$beta_y)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#946;<sub>y</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning", 
            html = TRUE
          )
        }
        
        if(input$omega_x > 0 & is.numeric(input$omega_x)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#969;<sub>x</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
        
        if(input$omega_y > 0 & is.numeric(input$omega_y)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#969;<sub>y</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
        
        if(input$sigma_x < 0 & is.numeric(input$sigma_x)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#963;<sub>x</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
        
        if(input$sigma_y < 0 & is.numeric(input$sigma_y)){
          sendSweetAlert(
            session = session,
            title = "Uh-oh!",
            text = HTML("It looks like &#963;<sub>y</sub> is out of range. I wouldn't trust that plot!"),
            type = "warning",
            html = TRUE
          )
        }
      })
      
      
      reactive_period <- reactive({input$period})
      reactive_step_size <- reactive({input$step_size})
      
      # output$n_timepoints <- renderText({
      #   sprintf("Number of Timepoints = %i", (reactive_period()[2] - reactive_period()[1]) / reactive_step_size())
      # })
      
      output$n_timepoints <- renderUI({
        n_obs <- floor((reactive_period()[2] - reactive_period()[1]) / reactive_step_size()) + 1
        
        if(is.finite(n_obs)){
          sprintf("Number of Timepoints = %i", n_obs)
        } else{
          HTML("Choose an appropriate value for &#916;t")
        }
      })
      
      output$formula_xy <- renderUI({
        withMathJax(model_formula_bivariate(input$modeltype_x, input$modeltype_y, input$xy_coupling))
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
      
      observeEvent(input$reset_x_params,{
        
        shinyjs::reset("x_params_inputs")
        
      })
      
      observeEvent(input$reset_y_params,{
        
        shinyjs::reset("y_params_inputs")
        
      })


      
      plot_data <- eventReactive(input$plot_button, {
        model_name_x <- input$modeltype_x
        model_name_y <- input$modeltype_y
        if(model_name_x == "l_ode" & model_name_y == "l_ode"){
          
          # 1. Both are Linear ODEs
          plot_xy_ode2(input)
          
        } else if( (model_name_x == "l_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "l_ode") ){
          
          # 2. Linear SDE and Linear ODE
          plot_xy_l2_ode_sde(input)
          
        } else if ( model_name_x == "ou" & model_name_y == "ou" ){
          
          # 3. Both are Linear SDEs
          plot_xy_l2_ode_sde(input)
          
        } else if ( (model_name_x == "l_ode" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "l_ode") ){
          
          # 4. Linear ODE and Harmonic Oscillator ODE
          plot_xy_l_harmonic_odes(input)
          
        } else if ( (model_name_x == "harmonic_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "harmonic_ode") ){
          
          # 5. Linear SDE and Harmonic Oscillator ODE
          plot_xy_l_sde_harmonic(input)
          
        } else if( model_name_x == "harmonic_ode" & model_name_y == "harmonic_ode" ){
          
          # 6. Both are Haromonic Oscillators ODE
          plot_xy_harmonics_ode(input)
          
        } else if( (model_name_x == "harmonic_sde" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "harmonic_sde") ){
          
          # 7. Harmonic Oscillator SDE and Linear ODE
          plot_xy_harmonic_sde_l_ode(input)
          
        } else if( (model_name_x == "harmonic_sde" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "harmonic_sde") ){
          
          # 8. Harmonic Oscillator SDE and Linear SDE
          plot_xy_harmonic_l_sdes(input)
          
        } else if( (model_name_x == "damp_ode" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "damp_ode") ){
          
          # 9. Damped Oscillator ODE and Linear ODE
          plot_xy_damped_linear_odes(input)
          
        } else if( (model_name_x == "damp_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "damp_ode") ){
          
          # 10. Damped Oscillator ODE and Linear SDE
          plot_xy_l_sde_damped(input)
          
        } else if( (model_name_x == "damp_ode" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "damp_ode") ){
          
          # 11. Damped Oscillator ODE and Harmonic ODE
          plot_xy_damped_harmonic_odes(input)
          
        } else if( model_name_x == "damp_ode" & model_name_y == "damp_ode" ){
          
          # 12. Both are Damped Oscillators ODE
          plot_xy_damped_odes(input)
          
        } else if( (model_name_x == "damp_sde" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "damp_sde") ){
          
          # 13. Damped Oscillator SDE and Linear ODE
          plot_xy_damp_sde_l_ode(input)
          
        } else if( (model_name_x == "damp_sde" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "damp_sde") ){
          
          # 14. Damped Oscilator SDE and Linear SDE
          plot_xy_damp_sde_l_sde(input)
          
        } else if( (model_name_x == "harmonic_sde" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "harmonic_sde") ){
          
          # 15. Harmonic Oscillator SDE and ODE
          plot_xy_harmonic_sde_ode(input)
          
        } else if( model_name_x == "harmonic_sde" & model_name_y == "harmonic_sde"){
          
          # 16. Both are Harmonic Oscillator SDEs
          plot_xy_harmonic_sdes(input)
          
        } else if( (model_name_x == "damp_ode" & model_name_y == "harmonic_sde") | (model_name_x == "harmonic_sde" & model_name_y == "damp_ode") ){
          
          # 17. Damped Oscillator ODE and Harmonic Oscillator SDE
          plot_xy_damp_ode_harmonic_sde(input)
          
        } else if( (model_name_x == "damp_sde" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "damp_sde") ){
          
          # 18. Damped Oscillator SDE and Harmonic Oscillator ODE
          plot_xy_damp_sde_harmonic_ode(input)
          
        } else if( (model_name_x == "damp_sde" & model_name_y == "harmonic_sde") | (model_name_x == "harmonic_sde" & model_name_y == "damp_sde") ){
          
          # 19. Damped Oscillator SDE and Harmonic Oscillator SDE
          plot_xy_damp_harmonic_sdes(input)
          
        } else if( (model_name_x == "damp_sde" & model_name_y == "damp_ode") | (model_name_x == "damp_ode" & model_name_y == "damp_sde") ){
          
          # 20. Damped Oscillator SDE and ODE
          plot_xy_damp_sde_ode(input)
          
        } else if( model_name_x == "damp_sde" & model_name_y == "damp_sde"){
          
          # 21. Both are Damped Oscillator SDEs
          plot_xy_damp_sdes(input)
          
        }
      })
      
      output$bivariate_plot <- renderPlot({
        plot_data()
      })
      
    }
  )
}

