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
###### server file ######
#########################

library(shiny)
library(bslib)
setwd("/path/to/your/simde")

source("plot_functions.R")
source("model_formulas.R")
source("one_var_module.R")

SimDE_server <- function(input, output, session) {
  
  one_var_server(id = "onevar")
  
  two_var_server(id = "twovars")
  
  shiny::addResourcePath("library", "/opt/shiny-server/simde")
                  
  output$userguide_window <- renderUI({
    tags$iframe(
      src = "library/user_guide.html",
      style="border:0; width: 100%; height:2500px;"
    )
  })
  
  output$info_window <- renderUI({
    tags$iframe(
      src = "library/contact.html",
      style="border:0; width: 100%; height:2500px;"
    )
  })
  
  session$onSessionEnded(stopApp)
  
  #bs_themer()
  
}




