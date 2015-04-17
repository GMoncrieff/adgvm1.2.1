# 
# 
b <- c(seq(0,20,0.05),seq(20,500,1))
b <- c(seq(0,1.5,0.001))
par(mar=c(4,4,1,1))


b1 <- c( 1.295720, 1.2, 1. )
b2 <- c( 0.392157, 0.45, 0.3 )


# plot(  b, b1[1]*b^b2[1],type="l")
# lines( b, b1[2]*b^b2[2],col="darkred")
# lines( b, b1[3]*b^b2[3],col="darkgreen")


g1 <- c( 0.6604492, 0.1604492, 0.9604492 )
g2 <- c( 2.5499324, 2.3499324, 4.9499324 )


plot(  b, exp( (log(b)+g1[1])/g2[1] ), col="black", type="l")
lines( b, exp( (log(b)+g1[2])/g2[2] ), col="darkred")
lines( b, exp( (log(b)+g1[3])/g2[3] ), col="darkgreen")


# 0.6604492 2.5499324
# 
# p <- c( 3.5413, 3.5337 )
# yy <- b1[1]*b^b2[1]
# 
# f <- function(p)
# {
# 	y <- exp( (log(b)+p[1])/p[2] )
# 	return( sum((y-yy)^2) )
# }
# 
# optim(p,f)


# grasse
# 
# b1 <- c(3.5)
# b2 <- c(0.5)
# 
# 
# p <- c( 3.5, 0.5 )
# yy <- b1[1]*b^b2[1]
# 
# f <- function(p)
# {
# 	y <- exp( (log(b)+p[1])/p[2] )
# 	return( sum((y-yy)^2) )
# }
# 
# optim(p,f)
# 
# g1 <- c( 2.505441 )
# g2 <- c( 1.999963 )
# 
# plot(  b, b1[1]*b^b2[1],type="l")
# lines( b, exp( (log(b)+2.505441)/1.999963 ), col="darkgreen")
# 










