tiesort <- function (xx=NULL, yy=NULL)
# sorts y values in increasing x-order
# within blocks of tied x-values, sorts y values in increasing y-order
{
  bignumber=100
  x=round(xx*bignumber)/bignumber
  y=round(yy*bignumber)/bignumber
  xtemp = x[order(x,y)]
  ytemp = y[order(x,y)]
  return(list(x=xtemp, y=ytemp))
}
