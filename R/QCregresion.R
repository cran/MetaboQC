#' Equation to be used internally to predict values from a regression curve of grade 3
#' @param b coefficient from order 0 part of the equation
#' @param c coefficient from order 1 part of the equation
#' @param d coefficient from order 2 part of the equation
#' @param e coefficient from order 3 part of the equation
#' @param x the x-axis value from which the y-axis value wanted to be predicted for the equation given by the coefficients
#' @examples
#' \dontrun{
#' prediction<-QCregression(b,c,d,e,x)
#' }
#' #' @export
#' @import plyr
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @return A y-value calculated for the x-value especified, taking into account the curve described by the coefficients given 

QCregression<-function(b,c,d,e,x){if (is.na(b)){b=0};
                                  if (is.na(c)){c=0};
                                  if (is.na(d)){d=0};
                                  if (is.na(e)){e=0};
                                  b+c*x+d*x^2+e*x^3}
