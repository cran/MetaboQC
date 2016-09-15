#' Equation to be used internally to predict values from a regression curve of grade 6
#' @param b coefficient from order 0 part of the equation
#' @param c coefficient from order 1 part of the equation
#' @param d coefficient from order 2 part of the equation
#' @param e coefficient from order 3 part of the equation
#' @param f coefficient from order 4 part of the equation
#' @param g coefficient from order 5 part of the equation
#' @param h coefficient from order 6 part of the equation
#' @param x the x-axis value from which the y-axis value wanted to be predicted for the equation given by the coefficients
#' @examples
#' \dontrun{
#' prediction<-QCregression4(b,c,d,e,f,g,h,x)
#' }
#' @export
#' @import plyr
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @return A y-value calculated for the x-value especified, taking into account the curve described by the coefficients given 

QCregression6<-function(b,c,d,e,f,g,h,x){if (is.na(b)){b=0};
                                  if (is.na(c)){c=0};
                                  if (is.na(d)){d=0};
                                  if (is.na(e)){e=0};
                                  if (is.na(f)){f=0};
                                  if (is.na(g)){g=0};
                                  if (is.na(h)){h=0};
                                  b+c*x+d*x^2+e*x^3+f*x^4+g*x^5+h*x^6}
