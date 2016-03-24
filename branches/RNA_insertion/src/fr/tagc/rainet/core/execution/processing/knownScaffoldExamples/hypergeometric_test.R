
# Script to produce p-value of hypergeometric test

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 4){
  stop("Rscript: Bad argument number")
}

x = as.numeric(args[1])
m = as.numeric(args[2])
n = as.numeric(args[3])
k = as.numeric(args[4])

# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# 
# x, q vector of quantiles representing the number of white balls drawn
# without replacement from an urn which contains both black and white
# balls.
# 
# m the number of white balls in the urn.
# 
# n the number of black balls in the urn.
# 
# k the number of balls drawn from the urn.

#phyper(10,11,578,578, lower.tail = FALSE) 


print (phyper(x,m,n,k, lower.tail = FALSE) )





