
args <- commandArgs(TRUE)

x = as.numeric(args[1])
m = as.numeric(args[2])
n = as.numeric(args[3])
k = as.numeric(args[4])

phyper(x,m,n,k, lower.tail = FALSE) 
