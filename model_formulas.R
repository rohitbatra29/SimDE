# This program is free softWare: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free SoftWare Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it Will be useful,
# but WITHOUT ANY WARRANTY; Without even the implied Warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along With this program. If not, see <https://WWW.gnu.org/licenses/>.

# This file contains the formulas for the different models (and their combinations) in the Server

library(stringr)


# Formulas for Univariate Models ------------------------------------------


model_formula_univariate <- function(modeltype) {
  if(modeltype == "l_ode") {
    eq <- "$$\\frac{dy(t)}{dt} \\: = \\: \\beta \\cdot (y(t) \\: - \\: \\mu)$$"
    
  } else if (modeltype == "p_ode") {
    eq <- "$$\\frac{dy(t)}{dt} \\: = \\: \\beta_{1} \\cdot y(t) \\: + \\: \\beta_{2} \\cdot y(t)^{2} \\: 
    + \\: \\beta_{3} \\cdot y(t)^{3} \\: + \\: \\beta_{4} \\cdot y(t)^{4} \\: - \\: \\mu$$"
    
  } else if (modeltype == "ou") {
    eq <- "$$dy(t) \\: = \\: \\beta \\cdot (y(t) - \\mu) \\cdot dt \\: + \\: \\sigma \\cdot dW(t)$$"
    
  } else if (modeltype == "p_sde") {
    eq <- "$$dy(t) \\: = \\: (\\beta_{1} \\cdot y(t) \\: + \\: \\beta_{2} \\cdot y(t)^{2} \\: 
    + \\: \\beta_{3} \\cdot y(t)^{3} \\: + \\: \\beta_{4} \\cdot y(t)^{4} \\: - \\: \\mu) \\cdot dt \\: + \\: \\sigma \\cdot dW(t)$$"
    
  } else if (modeltype == "harmonic_ode") {
    eq <- "$$\\frac{d^{2}y(t)}{dt^{2}} \\: = \\: \\omega \\cdot (y(t) \\: - \\: \\mu)$$"
    
  } else if(modeltype == "damp_ode") {
    eq <- "$$\\frac{d^{2}y(t)}{dt^{2}} \\: = \\: \\omega \\cdot (y(t) \\: - \\: \\mu) \\: + \\: \\zeta \\cdot \\frac{dy(t)}{dt}$$"
    
  } else if(modeltype == "harmonic_sde") {
    eq <- "$$ d^{2}y(t) \\: = \\: \\left( \\omega \\cdot (y(t) \\: - \\: \\mu) \\right) dt^{2} \\: + \\: \\sigma \\cdot dW(t) $$"
    
  } else if(modeltype == "damp_sde") {
    eq <- "$$d^{2}y(t) \\: = \\: \\left( \\omega \\cdot (y(t) \\: - \\: \\mu) \\: + \\: \\zeta \\cdot \\frac{dy(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma \\cdot dW(t)$$"
  }
  return(eq)
}


# Formulas for Bivariate Models -------------------------------------------


model_formula_bivariate <- function(model_name_x, model_name_y, xy_coupling) {
  if(model_name_x == "l_ode" & model_name_y == "l_ode") {
    

## 1. Both are Linear ODEs -----------------------------------------------
    
            if(xy_coupling == "coup_diff"){
              eq <- "$$
                    \\begin{align}
                        \\frac{dx(t)}{dt} \\: &= \\: \\beta_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\\\
                        \\frac{dy(t)}{dt} \\: &= \\: \\beta_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t))
                    \\end{align}
                    $$"
            } else if(xy_coupling == "coup_mean") {
              eq <- "$$
                    \\begin{align}
                        \\frac{dx(t)}{dt} \\: &= \\: \\beta_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\\\
                        \\frac{dy(t)}{dt} \\: &= \\: \\beta_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x})
                    \\end{align}
                    $$"
            }
    
  } else if ( (model_name_x == "l_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "l_ode") ) {
    

## 2. Linear SDE and Linear ODE ------------------------------------------

    
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_ou <- "dv1(t) \\: &= \\: \\left(\\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\right) dt \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_ou <- "dv1(t) \\: &= \\: \\left(\\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\right) dt \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt"
              
             }
    
          if(model_name_x == "ou"){
            
            # Assign ou equation to X and l_ode equation to Y
            
            eq_ou <- str_replace_all(eq_ou, c("v1" = "x", "v2" = "y"))
            eq_l_ode <- str_replace_all(eq_l_ode, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_ou, "\\\\", eq_l_ode, "\\end{align}$$")
            
          } else if(model_name_y == "ou"){
            
            # Assign l_ode equation to X and ou equation to Y
            
            eq_l_ode <- str_replace_all(eq_l_ode, c("v2" = "x", "v1" = "y"))
            eq_ou <- str_replace_all(eq_ou, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_l_ode, "\\\\", eq_ou, "\\end{align}$$")
            }
    
  } else if (model_name_x == "ou" & model_name_y == "ou"){
    

## 3. Both are Linear SDEs -----------------------------------------------

          if(xy_coupling == "coup_diff"){
            eq <- "$$
                    \\begin{align}
                        dx(t) \\: &= \\: \\left(\\beta_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\right) dt \\: + \\: \\sigma_{x} \\cdot dW_{v1}(t)\\\\
                        dy(t) \\: &= \\: \\left(\\beta_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t)) \\right) dt \\: + \\: \\sigma_{y} \\cdot dW_{v2}(t)
                    \\end{align}
                    $$"
          } else if(xy_coupling == "coup_mean") {
            eq <- "$$
                    \\begin{align}
                        dx(t) \\: &= \\: \\left(\\beta_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\right) dt \\: + \\: \\sigma_{x} \\cdot dW_{v1}(t)\\\\
                        dy(t) \\: &= \\: \\left(\\beta_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x}) \\right) dt \\: + \\: \\sigma_{y} \\cdot dW_{v2}(t)
                    \\end{align}
                    $$"
          }
    
  } else if ( (model_name_x == "harmonic_ode" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "harmonic_ode") ){
    

## 4. Linear ODE and Harmonic Oscillator ODE -----------------------------

    
          # Assign the equation based on the coupling
          if(xy_coupling == "coup_diff"){
            
            eq_l_ode <- "\\frac{dv1(t)}{dt} \\: &= \\: \\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t))"
            
            eq_harmonic_ode <- "\\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t))"
            
          } else if(xy_coupling == "coup_mean") {
            
            eq_l_ode <- "\\frac{dv1(t)}{dt} \\: &= \\: \\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2})"
            
            eq_harmonic_ode <- "\\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1})"
            
          } else if(xy_coupling == "coup_diff_osc"){
          
            eq_l_ode <- "\\frac{dv1(t)}{dt} \\: &= \\: \\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t))"
            
            eq_harmonic_ode <- "\\frac{dv2(t)}{dt} \\: &= \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\\\ \\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2})"
            
          } else if(xy_coupling == "coup_mean_osc"){
          
            eq_l_ode <- "\\frac{dv1(t)}{dt} \\: &= \\: \\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2})"
            
            eq_harmonic_ode <- "\\frac{dv2(t)}{dt} \\: &= \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\\\ \\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2})"
            
          }
    
          if(model_name_x == "l_ode"){
            
            # Assign l_ode equation to X and harmonic_ode equation to Y
            
            eq_l_ode <- str_replace_all(eq_l_ode, c("v1" = "x", "v2" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_l_ode, " \\\\ ", eq_harmonic_ode, "\\end{align}$$")
            
          } else if(model_name_y == "l_ode"){
            
            # Assign harmonic_ode equation to X and l_ode equation to Y
            eq_l_ode <- str_replace_all(eq_l_ode, c("v2" = "x", "v1" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_ode, " \\\\ ", eq_l_ode, "\\end{align}$$")
          }
    
  } else if ( (model_name_x == "harmonic_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "harmonic_ode") ){
  

## 5. Linear SDE and Harmonic Oscillator ODE -----------------------------

    
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_ou <- "dv1(t) \\: &= \\: \\left(\\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\right) dt \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt^{2}"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_ou <- "dv1(t) \\: &= \\: \\left(\\beta_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\right) dt \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt^{2}"
              
            }
    
          if(model_name_x == "ou"){
            
            # Assign ou equation to X and harmonic_ode equation to Y
            
            eq_ou <- str_replace_all(eq_ou, c("v1" = "x", "v2" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_ou, " \\\\ ", eq_harmonic_ode, "\\end{align}$$")
            
          } else if(model_name_y == "ou"){
            
            # Assign harmonic_ode equation to X and ou equation to Y
            eq_ou <- str_replace_all(eq_ou, c("v2" = "x", "v1" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_ode, " \\\\ ", eq_ou, "\\end{align}$$")
          }
    
  } else if( model_name_x == "harmonic_ode" & model_name_y == "harmonic_ode" ){
    

## 6. Both are Harmonic Oscillators ODE ----------------------------------

            if(xy_coupling == "coup_diff"){
              eq <- "$$
                      \\begin{align}
                          \\frac{d^{2}x(t)}{dt^{2}} \\: &= \\: \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\\\
                          \\frac{d^{2}y(t)}{dt^{2}} \\: &= \\: \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t))
                      \\end{align}
                      $$"
            } else if(xy_coupling == "coup_mean") {
              eq <- "$$
                      \\begin{align}
                          \\frac{d^{2}x(t)}{dt^{2}} \\: &= \\: \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\\\
                          \\frac{d^{2}y(t)}{dt^{2}} \\: &= \\: \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x})
                      \\end{align}
                      $$"
            }
    
  } else if( (model_name_x == "harmonic_sde" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "harmonic_sde") ){
    

## 7. Harmonic Oscillator SDE and Linear ODE -----------------------------

            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right)"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right)"
              
            }
    
      if(model_name_x == "harmonic_sde"){
        
        # Assign harmonic_sde equation to X and l_ode equation to Y
        eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v1" = "x", "v2" = "y"))
        eq_l_ode <- str_replace_all(eq_l_ode, c("v1" = "x", "v2" = "y"))
        
        eq <- paste("$$\\begin{align}", eq_harmonic_sde, "\\\\", eq_l_ode, "\\end{align}$$")
        
      } else if(model_name_y == "harmonic_sde"){
        
        # Assign l_ode equation to X and harmonic_sde equation to Y
        
        eq_harmonic_sde  <- str_replace_all(eq_harmonic_sde , c("v2" = "x", "v1" = "y"))
        eq_l_ode <- str_replace_all(eq_l_ode, c("v2" = "x", "v1" = "y"))
        
        eq <- paste("$$\\begin{align}", eq_l_ode, "\\\\", eq_harmonic_sde, "\\end{align}$$")
      }
    
  } else if( (model_name_x == "harmonic_sde" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "harmonic_sde") ){
    

## 8. Harmonic Oscillator SDE and Linear SDE -----------------------------
    
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            }
            
            if(model_name_x == "harmonic_sde"){
              
              # Assign harmonic_sde equation to X and ou equation to Y
              eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v1" = "x", "v2" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_harmonic_sde, "\\\\", eq_ou, "\\end{align}$$")
              
            } else if(model_name_y == "harmonic_sde"){
              
              # Assign ou equation to X and harmonic_sde equation to Y
              
              eq_harmonic_sde  <- str_replace_all(eq_harmonic_sde , c("v2" = "x", "v1" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_ou, "\\\\", eq_harmonic_sde, "\\end{align}$$")
            }
    
  } else if ( (model_name_x == "damp_ode" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "damp_ode") ){
    

## 9. Damped Oscillator ODE and Linear ODE -------------------------------
        
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_ode <- "\\frac{d^{2}v1(t)}{dt^{2}} \\: &= \\: \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt}"
              
              eq_l_ode <- "\\frac{dv2(t)}{dt} \\: &= \\: \\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t))"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_ode <- "\\frac{d^{2}v1(t)}{dt^{2}} \\: &= \\: \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt}"
              
              eq_l_ode <- "\\frac{dv2(t)}{dt} \\: &= \\: \\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1})"
              
            }
            
            if(model_name_x == "damp_ode"){
              
              # Assign damp_ode equation to X and l_ode equation to Y
              
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v1" = "x", "v2" = "y"))
              eq_l_ode <- str_replace_all(eq_l_ode, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_ode, " \\\\ ", eq_l_ode, "\\end{align}$$")
              
            } else if(model_name_y == "damp_ode"){
              
              # Assign l_ode equation to X and damp_ode equation to Y
              
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v2" = "x", "v1" = "y"))
              eq_l_ode <- str_replace_all(eq_l_ode, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_l_ode, " \\\\ ",  eq_damp_ode, "\\end{align}$$")
            }
    
  } else if( (model_name_x == "damp_ode" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "damp_ode") ){
    

## 10. Damped Oscillator ODE and Linear SDE ------------------------------

    
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_ode <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2}"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left( \\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_ode <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2}"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left( \\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            }
            
            if(model_name_x == "damp_ode"){
              
              # Assign damp_ode equation to X and ou equation to Y
              
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v1" = "x", "v2" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_ode, " \\\\ ", eq_ou, "\\end{align}$$")
              
            } else if(model_name_y == "damp_ode"){
              
              # Assign ou equation to X and damp_ode equation to Y
              
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v2" = "x", "v1" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_ou, " \\\\ ",  eq_damp_ode, "\\end{align}$$")
            }
    
  } else if ( (model_name_x == "damp_ode" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "damp_ode") ){
    

## 11. Damped Oscillator ODE and Harmonic Oscillator ODE -----------------

    
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_ode <- "\\frac{d^{2}v1(t)}{dt^{2}} \\: &= \\: \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt}"
              
              eq_harmonic_ode <- "\\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t))"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_ode <- "\\frac{d^{2}v1(t)}{dt^{2}} \\: &= \\: \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt}"
              
              eq_harmonic_ode <- "\\frac{d^{2}v2(t)}{dt^{2}} \\: &= \\: \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1})"
              
            }
            
            if(model_name_x == "damp_ode"){
              
              # Assign damp_ode equation to X and harmonic_ode equation to Y
              
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v1" = "x", "v2" = "y"))
              eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_ode, " \\\\ ", eq_harmonic_ode, "\\end{align}$$")
              
            } else if(model_name_y == "damp_ode"){
              
              # Assign harmonic_ode equation to X and l_ode equation to Y
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v2" = "x", "v1" = "y"))
              eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_harmonic_ode, " \\\\ ", eq_damp_ode, "\\end{align}$$")
            }
    
  } else if( model_name_x == "damp_ode" & model_name_y == "damp_ode" ){
    

## 12. Both are Damped Oscillator ODEs -----------------------------------

    if(xy_coupling == "coup_diff"){
      eq <- "$$ \\begin{align}
                    \\frac{d^{2}x(t)}{dt^{2}} \\: &= \\: \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\: + \\: \\zeta_{x} \\cdot \\frac{dx(t)}{dt} \\\\
                    \\frac{d^{2}y(t)}{dt^{2}} \\: &= \\: \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t)) \\: + \\: \\zeta_{y} \\cdot \\frac{dy(t)}{dt}
                \\end{align} 
                $$"
    } else if(xy_coupling == "coup_mean") {
      eq <- "$$ \\begin{align}
                    \\frac{d^{2}x(t)}{dt^{2}} \\: &= \\: \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\zeta_{x} \\cdot \\frac{dx(t)}{dt} \\\\
                    \\frac{d^{2}y(t)}{dt^{2}} \\: &= \\: \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\zeta_{y} \\cdot \\frac{dy(t)}{dt}
                \\end{align} 
                $$"
    }
    
  } else if( (model_name_x == "damp_sde" & model_name_y == "l_ode") | (model_name_x == "l_ode" & model_name_y == "damp_sde") ){
    

## 13. Damped Oscillator SDE and Linear ODE ------------------------------

            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right)"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_l_ode <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right)"
              
            }
            
            if(model_name_x == "damp_sde"){
              
              # Assign damp_sde equation to X and l_ode equation to Y
              eq_damp_sde <- str_replace_all(eq_damp_sde, c("v1" = "x", "v2" = "y"))
              eq_l_ode <- str_replace_all(eq_l_ode, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_sde, "\\\\", eq_l_ode, "\\end{align}$$")
              
            } else if(model_name_y == "damp_sde"){
              
              # Assign l_ode equation to X and damp_sde equation to Y
              
              eq_damp_sde  <- str_replace_all(eq_damp_sde , c("v2" = "x", "v1" = "y"))
              eq_l_ode <- str_replace_all(eq_l_ode, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_l_ode, "\\\\", eq_damp_sde, "\\end{align}$$")
            }
    
  } else if( (model_name_x == "damp_sde" & model_name_y == "ou") | (model_name_x == "ou" & model_name_y == "damp_sde") ){
    

## 14. Damped Oscillator SDE and Linear SDE ------------------------------

            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_ou <- "dv2(t) \\: &= \\: \\left(\\beta_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
              
            }
            
            if(model_name_x == "damp_sde"){
              
              # Assign damp_sde equation to X and ou equation to Y
              eq_damp_sde <- str_replace_all(eq_damp_sde, c("v1" = "x", "v2" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_sde, "\\\\", eq_ou, "\\end{align}$$")
              
            } else if(model_name_y == "damp_sde"){
              
              # Assign ou equation to X and damp_sde equation to Y
              
              eq_damp_sde  <- str_replace_all(eq_damp_sde , c("v2" = "x", "v1" = "y"))
              eq_ou <- str_replace_all(eq_ou, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_ou, "\\\\", eq_damp_sde, "\\end{align}$$")
            }
    
  } else if( (model_name_x == "harmonic_sde" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "harmonic_sde") ){
    

## 15. Harmonic Oscillator SDE and ODE -----------------------------------

          # Assign the equation based on the coupling
          if(xy_coupling == "coup_diff"){
            
            eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt^{2}"
            
          } else if(xy_coupling == "coup_mean") {
            
            eq_harmonic_sde <- "d^{2}v1(t) \\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt^{2}"
            
          }
          
          if(model_name_x == "harmonic_sde"){
            
            # Assign harmonic_sde equation to X and harmonic_ode equation to Y
            eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v1" = "x", "v2" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_sde, "\\\\", eq_harmonic_ode, "\\end{align}$$")
            
          } else if(model_name_y == "harmonic_sde"){
            
            # Assign harmonic_ode equation to X and harmonic_sde equation to Y
            
            eq_harmonic_sde  <- str_replace_all(eq_harmonic_sde , c("v2" = "x", "v1" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_ode, "\\\\", eq_harmonic_sde, "\\end{align}$$")
          }
    
  } else if( model_name_x == "harmonic_sde" & model_name_y == "harmonic_sde"){
    

## 16. Both are Harmonic Oscillator SDE ------------------------------------

          if(xy_coupling == "coup_diff"){
            eq <- "$$
                            \\begin{align}
                                d^{2}x(t) \\: &= \\: \\left( \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\right) dt^{2} \\: + \\: \\sigma_{x} \\cdot dW_{x}(t) \\\\
                                d^{2}y(t) \\: &= \\: \\left( \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t)) \\right) dt^{2} \\: + \\: \\sigma_{y} \\cdot dW_{y}(t)
                            \\end{align}
                            $$"
          } else if(xy_coupling == "coup_mean") {
            eq <- "$$
                            \\begin{align}
                                d^{2}x(t) \\: &= \\: \\left( \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\right) dt^{2} \\: + \\: \\sigma_{x} \\cdot dW_{x}(t) \\\\
                                d^{2}y(t) \\: &= \\: \\left( \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x}) \\right) dt^{2} \\: + \\: \\sigma_{y} \\cdot dW_{y}(t)
                            \\end{align}
                            $$"
          }
  } else if( (model_name_x == "damp_ode" & model_name_y == "harmonic_sde") | (model_name_x == "harmonic_sde" & model_name_y == "damp_ode") ){
    

## 17. Damped Oscillator ODE and Harmonic Oscillator SDE -----------------

        # Assign the equation based on the coupling
        if(xy_coupling == "coup_diff"){
          
          eq_damp_ode <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2}"
          
          eq_harmonic_sde <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt^{2} \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
          
        } else if(xy_coupling == "coup_mean") {
          
          eq_damp_ode <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2}"
          
          eq_harmonic_sde <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt^{2} \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
          
        }
        
        if(model_name_x == "damp_ode"){
          
          # Assign damp_ode equation to X and harmonic_sde equation to Y
          eq_damp_ode <- str_replace_all(eq_damp_ode, c("v1" = "x", "v2" = "y"))
          eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v1" = "x", "v2" = "y"))
          
          eq <- paste("$$\\begin{align}", eq_damp_ode, "\\\\", eq_harmonic_sde, "\\end{align}$$")
          
        } else if(model_name_y == "damp_ode"){
          
          # Assign harmonic_sde equation to X and damp_ode equation to Y
          
          eq_damp_ode  <- str_replace_all(eq_damp_ode, c("v2" = "x", "v1" = "y"))
          eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v2" = "x", "v1" = "y"))
          
          eq <- paste("$$\\begin{align}", eq_harmonic_sde, "\\\\", eq_damp_ode, "\\end{align}$$")
        }
    
  } else if( (model_name_x == "damp_sde" & model_name_y == "harmonic_ode") | (model_name_x == "harmonic_ode" & model_name_y == "damp_sde") ){
    

## 18. Damped Oscillator SDE and Harmonic Oscillator ODE -----------------
          # Assign the equation based on the coupling
          if(xy_coupling == "coup_diff"){
            
            eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt^{2}"
            
          } else if(xy_coupling == "coup_mean") {
            
            eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt^{2}"
            
          }
          
          if(model_name_x == "damp_sde"){
            
            # Assign damp_sde equation to X and harmonic_ode equation to Y
            eq_damp_sde <- str_replace_all(eq_damp_sde, c("v1" = "x", "v2" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_damp_sde, "\\\\", eq_harmonic_ode, "\\end{align}$$")
            
          } else if(model_name_y == "damp_sde"){
            
            # Assign harmonic_ode equation to X and damp_sde equation to Y
            
            eq_damp_sde  <- str_replace_all(eq_damp_sde, c("v2" = "x", "v1" = "y"))
            eq_harmonic_ode <- str_replace_all(eq_harmonic_ode, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_ode, "\\\\", eq_damp_sde, "\\end{align}$$")
          }
    
  } else if( (model_name_x == "damp_sde" & model_name_y == "harmonic_sde") | (model_name_x == "harmonic_sde" & model_name_y == "damp_sde") ){
    

## 19. Damped Oscillator SDE and Harmonic Oscillator SDE -----------------
          # Assign the equation based on the coupling
          if(xy_coupling == "coup_diff"){
            
            eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_sde <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\right) dt^{2} \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
            
          } else if(xy_coupling == "coup_mean") {
            
            eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
            
            eq_harmonic_sde <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\right) dt^{2}  \\: + \\: \\sigma_{v2} \\cdot dW_{v2}(t)"
            
          }
          
          if(model_name_x == "damp_sde"){
            
            # Assign damp_sde equation to X and harmonic_sde equation to Y
            eq_damp_sde <- str_replace_all(eq_damp_sde, c("v1" = "x", "v2" = "y"))
            eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v1" = "x", "v2" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_damp_sde, "\\\\", eq_harmonic_sde, "\\end{align}$$")
            
          } else if(model_name_y == "damp_sde"){
            
            # Assign harmonic_sde equation to X and damp_sde equation to Y
            
            eq_damp_sde  <- str_replace_all(eq_damp_sde, c("v2" = "x", "v1" = "y"))
            eq_harmonic_sde <- str_replace_all(eq_harmonic_sde, c("v2" = "x", "v1" = "y"))
            
            eq <- paste("$$\\begin{align}", eq_harmonic_sde, "\\\\", eq_damp_sde, "\\end{align}$$")
          }
    
    
  } else if( (model_name_x == "damp_sde" & model_name_y == "damp_ode") | (model_name_x == "damp_ode" & model_name_y == "damp_sde") ){
    

## 20. Damped Oscillator SDE and ODE -------------------------------------
            # Assign the equation based on the coupling
            if(xy_coupling == "coup_diff"){
              
              eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v1(t) - v2(t)) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_damp_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v2(t) - v1(t)) \\: + \\: \\zeta_{v2} \\cdot \\frac{dv2(t)}{dt} \\right) dt^{2}"
              
            } else if(xy_coupling == "coup_mean") {
              
              eq_damp_sde <- "d^{2}v1(t)\\: &= \\: \\left( \\omega_{v1} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\gamma_{v1} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\zeta_{v1} \\cdot \\frac{dv1(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{v1} \\cdot dW_{v1}(t)"
              
              eq_damp_ode <- "d^{2}v2(t) \\: &= \\: \\left( \\omega_{v2} \\cdot (v2(t) \\: - \\: \\mu_{v2}) \\: + \\: \\gamma_{v2} \\cdot (v1(t) \\: - \\: \\mu_{v1}) \\: + \\: \\zeta_{v2} \\cdot \\frac{dv2(t)}{dt} \\right) dt^{2}"
              
            }
            
            if(model_name_x == "damp_sde"){
              
              # Assign damp_sde equation to X and damp_ode equation to Y
              eq_damp_sde <- str_replace_all(eq_damp_sde, c("v1" = "x", "v2" = "y"))
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v1" = "x", "v2" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_sde, "\\\\", eq_damp_ode, "\\end{align}$$")
              
            } else if(model_name_y == "damp_sde"){
              
              # Assign damp_ode equation to X and damp_sde equation to Y
              
              eq_damp_sde  <- str_replace_all(eq_damp_sde, c("v2" = "x", "v1" = "y"))
              eq_damp_ode <- str_replace_all(eq_damp_ode, c("v2" = "x", "v1" = "y"))
              
              eq <- paste("$$\\begin{align}", eq_damp_ode, "\\\\", eq_damp_sde, "\\end{align}$$")
            }
    
  } else if( model_name_x == "damp_sde" & model_name_y == "damp_sde"){
    

## 21. Both are Damped Oscillator SDEs -----------------------------------
          # Assign the equation based on the coupling
          if(xy_coupling == "coup_diff"){
            
            eq <- "$$
                    \\begin{align}
                       d^{2}x(t)\\: &= \\: \\left( \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (x(t) - y(t)) \\: + \\: \\zeta_{x} \\cdot \\frac{dx(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{x} \\cdot dW_{x}(t) \\\\
                        d^{2}y(t)\\: &= \\: \\left( \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (y(t) - x(t)) \\: + \\: \\zeta_{y} \\cdot \\frac{dy(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{y} \\cdot dW_{y}(t)
                    \\end{align}
                    $$"
          } else if(xy_coupling == "coup_mean") {
            
            eq <- "$$
                    \\begin{align}
                       d^{2}x(t)\\: &= \\: \\left( \\omega_{x} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\gamma_{x} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\zeta_{x} \\cdot \\frac{dx(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{x} \\cdot dW_{x}(t) \\\\
                        d^{2}y(t)\\: &= \\: \\left( \\omega_{y} \\cdot (y(t) \\: - \\: \\mu_{y}) \\: + \\: \\gamma_{y} \\cdot (x(t) \\: - \\: \\mu_{x}) \\: + \\: \\zeta_{y} \\cdot \\frac{dy(t)}{dt} \\right) dt^{2} \\: + \\: \\sigma_{y} \\cdot dW_{y}(t)
                    \\end{align}
                    $$"
          }
    
  }
  
  return(eq)
}









