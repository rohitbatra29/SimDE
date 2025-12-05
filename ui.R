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
