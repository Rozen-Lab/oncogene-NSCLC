

# Formula from e.g. 
# https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module8-CategoricalData/PH717-Module8-CategoricalData5.html
my.risk.ratio <- function(mm, zz = 1.96) {
  n <- rowSums(mm)
  pp <- (mm[1,1] / n[1]) / (mm[2,1] / n[2])
  xx <- zz * 
    sqrt(
    (mm[1,2]/mm[1,1]) / n[1] +
      mm[2,2]/mm[2,1] / n [2]
  )
  yy <- c(log(pp) - xx, log(pp) + xx)
  
  return(list(RR = pp, conf.int = exp(yy), plus.minus = xx))
  
}

# Test example from the web page above, without the rounding error

test <- matrix(c(992, 165, 2260, 1017), nrow = 2)

m <- matrix(c(792, 1106 - 792, 274, 679 - 274), byrow = TRUE, nrow = 2)

tt <- matrix(c(9, 20, 41, 29), nrow = 2)

t2 <- matrix(c(23, 50 - 23, 11, 50 - 11), byrow = T, nrow = 2)

risk.ratio(m)
