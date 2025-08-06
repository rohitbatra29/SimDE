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

#########################
######## ui file ########
#########################


setwd("/path/to/your/simde")
library(shiny)
library(bslib)
library(shinyjs)
library(bsicons)

source("one_var_module.R")
source("two_var_module.R")

# work on aligning the tab names and their icons!
SimDE_ui <- tagList(
  tags$head(
    tags$style(HTML("
          .navbar-nav { display: flex !important;
                        justify-content: flex-end !important;
                      }
          .navbar-header {
                          display: flex;
                          align-items: center;
                    }
          .nav-pills .nav-link.active {
                                        background-color: #0044A4;} 
        
      .nav-pills .nav-link {color: #0044A4;}
      
      .panel-custom {
                      border: 1px solid #ccc;
                      border-radius: 10px;
                      margin-bottom: 20px;
                      padding-top:2em;
                      padding-bottom:2em;
                      padding-left:2em;
      }
      .plot-btn{background-color:white} 
      .plot-btn:hover{background-color:#0044A4}
                   "
                  
                    )
               )
    ),
  navbarPage(div(tags$img(src = "PNG.png", height='110', width='110', style = "margin-right: 20px;"),
                 tags$text(HTML("Simulate and Visualize Verbal Theories as Differential Equations")),
                 style = "display: flex; align-items: center; font-size: 30px; font-family: 'Palatino', 'Palatino Linotype', 'Palatino LT STD'"),
             id = "topnavbar",
             selected = "simulate",
             useShinyjs(),
             tabPanel(div(
                          style = "align-items: baseline",
                          img(src='graph-time-series.svg', height='30', width='30'), 
                          tags$span("Simulate", style = "font-size: 25px; font-family: 'Palatino', 'Palatino Linotype', 'Palatino LT STD';")
                          ),
                      value = "simulate",
                      h2("Define Your Target System"),
                      withMathJax(),
                      br(),
                      h4("Number of Variables"),
                      tabsetPanel(id = "variables",
                                  type = "pill",
                                  one_var_UI(id = "onevar"),
                                  two_var_UI(id = "twovars")
                                  )
                      
                      ),
             tabPanel(div(
                          img(src='https://www.svgrepo.com/show/501685/question.svg', height='30', width='30'),
                          tags$span("User Guide", style = "font-size: 25px; vertical-align: middle; font-family: 'Palatino', 'Palatino Linotype', 'Palatino LT STD'"),
                          style = "align-items: baseline"
                          ),
                      value = "user_guide",
                      uiOutput("userguide_window")), 
             tabPanel(div(img(src = "https://www.svgrepo.com/show/317909/contact-mail.svg", height = "30", width = "30"), 
                          tags$span("Contact", style = "font-size: 25px; vertical-align: middle; font-family: 'Palatino', 'Palatino Linotype', 'Palatino LT STD'")),
                      value = "info",
                      uiOutput("info_window")),
             collapsible = T,
             theme =  bs_theme(version = 5, bg = "white", fg = "#0d0c0c")
  )
)

# Things to do!
# add a small description along with equation so that people who don't know the particular model can get a small jist
# tool tips - pop ups for parameter names (https://educationshinyappteam.github.io/Style_Guide/commonElements.html#popovers-and-tooltips)
# how can we add a random button beside set your seed, so that whenever someone clicks on it, it fills a random value in the seed input. Would be neat if this is possible!
# require/validate
# Maybe remove the word Undamped from Harmonic because it might not be undamped in stochastic systems.
# the plots of some linear ones have the mu as a dashed line while some plots don't. Decide on the consistent look!
# Will have to say that there might be a lag with the oscillatory SDE functions in the bivariate and it is better to write your number, instead of increasing/decreasing in the input, because then it tries to simulate each of those possibilities. You can also show this as a popup in the cases where Julia is used.

# Other Stuff (not a priority)
# maybe when you make the github repository public later, change the svg icon links to be from files in the github repository instead of the website.
# check if you can progress indicator about how much time is left for the plot to be done so that the user knows to wait. Julia has a progress bar thing https://docs.sciml.ai/DiffEqDocs/stable/features/progress_bar/ but Idk how to integrate it into Shiny's progress bar https://shiny.posit.co/r/articles/build/progress/
# add a run button so they don't run it agian and again?


# Helpful tips
# conditonal goes for every variable that changes the type of output panel you have
# observe in server goes for every variable that changes the type of options you have
# switch - function for calling functions?

# always input "unevaluated" reactive values into a module server!!
#? you can pass an unevaluated reactive value to the module and then evaluate it inside the module??
# For example (refer to screenshot in old folder) but basically, put the reactive value (just the name) in the module server function and then within the function, put it as an evaluation with () paranethesis.
# similarly you can do this with multiple reactive values if you have them in a reactiveValues object.
# similarly you can return reactive values from a module serve using a return function, so that they can be used in the main app.
