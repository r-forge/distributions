#############################################################################
#   Copyright (c) 2009 Christophe Dutang                                                                                                  #
#                                                                                                                                                                        #
#   This program is free software; you can redistribute it and/or modify                                               #
#   it under the terms of the GNU General Public License as published by                                         #
#   the Free Software Foundation; either version 2 of the License, or                                                   #
#   (at your option) any later version.                                                                                                            #
#                                                                                                                                                                         #
#   This program is distributed in the hope that it will be useful,                                                             #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 #
#   GNU General Public License for more details.                                                                                    #
#                                                                                                                                                                         #
#   You should have received a copy of the GNU General Public License                                           #
#   along with this program; if not, write to the                                                                                           #
#   Free Software Foundation, Inc.,                                                                                                              #
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             #
#                                                                                                                                                                         #
#############################################################################

example.continuous <- function(distr="unif")
{
#density definitions
    dtransgamma <- function(x, alpha, lambda, tau)
    {
        tau*x^(alpha*tau)*exp(-(x/lambda)^tau)/( lambda^(alpha*tau)*x*gamma(alpha) )
    }
    mydinvexp <- function(x, rate) rate/x^2*exp(-rate/x) 
    generlang <- function(x, lam1, lam2, lam3)
    {
        return(
               lam2/(lam2-lam1)*lam3/(lam3-lam1)*lam1*exp(-lam1*x)+	
               lam1/(lam1-lam2)*lam3/(lam3-lam2)*lam2*exp(-lam2*x)+	
               lam1/(lam1-lam3)*lam2/(lam2-lam3)*lam3*exp(-lam3*x)
               )
    }
    dchi <- function(x, df, noncentral=0)
    {
        x^(df-1)*exp(-x^2/2)/(2^(df/2-1)*gamma(df/2))
    }
    dinvchisq<-function(x,k)
    {
        2^(-k/2)/gamma(k/2)*x^(-k/2-1)*exp(-1/(2*x))
    }
    
        
    
    if(distr=="unif")
    {   
        x <- seq(0,3, .01)            
        y <- dunif(x)
        y2 <- dunif(x, 0, 2)
        y3 <- dunif(x, 0, 3)
        leg.txt <- c("U(0,1)","U(0,2)","U(0,3)")
    }
    if(distr=="triangle")
    {
        x <- seq(0,2, .01)            
        y <- dtriangle(x, 0, 2, 1)      
        y2 <- dtriangle(x, 0, 2, 1/2)
        y3 <- dtriangle(x, 0, 2, 4/3)
        leg.txt <- c("T(0,2,1)","T(0,2,1/2)","T(0,2,4/3)")
    }
    if(distr=="beta")
    {
        x <- seq(0,1, .01)         
        y <- dbeta(x, 2, 2)
        y2 <- dbeta(x, 3, 1)
        y3 <- dbeta(x, 1, 5)
        y4 <- dbeta(x, 1/2, 1/2)
        leg.txt <- c("B(2,2)","B(3,1)","B(1,5)","B(1/2,1/2)")
    }
    if(distr=="kumar")
    {
        x <- seq(0,1, .01)         
        y <- dkumar(x, 5, 2)
        y2 <- dkumar(x, 2, 2.5)  
        y3 <- dkumar(head(x, -1), 1/2, 1/2)
        y4 <- dkumar(x, 1, 3)
        leg.txt <- c("K(5,2)","K(2,2.5)","K(1/2,1/2)","K(1,3)")
    }
    if(distr=="genbeta")
    {
        x <- seq(0,2, .01)         
        y <- dgenbeta(x, 2, 2, 2, 1/2)
        y2 <- dgenbeta(x, 3, 1, 2, 1/2)
        y3 <- dgenbeta(x, 3, 1, 1/2, 1/2)
        y4 <- dgenbeta(x, 1/2, 2, 1/3, 1/2)
        leg.txt <- c("GB(2,2,2,1/2)","GB(3,1,2,1/2)","GB(3,1,1/2,1/2)","GB(1/2,2,1/3,1/2)")
    }
    if(distr=="exp")
    {
        x <- seq(0,5,.01)
        y <- dexp(x)  
        y2 <- dexp(x, 2)
        y3 <- dexp(x, 1/2)
        leg.txt <- c("E(1)","E(2)","E(1/2)")
    }    
    if(distr=="shiftexp")
    {
        x <- seq(0,5,.01)
        y <- dexp(x)  
        y2 <- dexp(x-1, 1/2)
        y3 <- dexp(x-2, 1/2)
        leg.txt <- c("E(0,1)","E(1,1/2)","E(2,1/2)")
    }        
    if(distr=="invexp")
    {
        x <- seq(0,5,.01)
        y <- mydinvexp(x, 1)
        y2 <- mydinvexp(x, 2)
        y3 <- mydinvexp(x, 3)
        leg.txt <- c("IE(1)","IE(2)","IE(3)")
    }
    
        
    if(distr=="shiftlnorm")
    {
        x <- seq(0,5,.01)
        y <- dlnorm(x)  
        y2 <- dlnorm(x-1, 0, 1)
        y3 <- dlnorm(x-1/2, 0, 1)
        y4 <- dnorm(x, -1, 1)
        leg.txt <- c("LN(0,0,1)","LN(1,0,1)","LN(1/2,0,1)","LN(0,-1,1)")
    }    
    if(distr=="gamma")
    {
        x <- seq(0,5,.01)
        y <- dgamma(x,1)  
        y2 <- dgamma(x, 2, 1)
        y3 <- dgamma(x, 2, 2)
        y4 <- dgamma(x, 1/2, 1)
        leg.txt <- c("G(1,1)","G(2,1)","G(2,2)","G(1/2,1)")
    }    
    if(distr=="gig")
    {
        x <- seq(0,5,.01)
        y <- dgig(x , lambda = -1/2, chi = 5, psi = 1)
        y2 <- dgig(x , lambda = -1, chi = 5, psi = 1)
        y3 <- dgig(x , lambda = -1, chi = 2, psi = 3)
        y4 <- dgig(x , lambda = 1, chi = 5, psi = 1)
        leg.txt <- c("GIG(-1/2,5,1)","GIG(-1,5,1)","GIG(-1,2,3)","GIG(1,5,1)")
    }    
    if(distr=="generlang")
    {
        x <- seq(0,5,.01)
        y<- generlang(x, 1,2,3)
        y2<- generlang(x, 1,2,4)
        y3<- generlang(x, 1,3,5)
        y4 <- generlang(x, 2,3,4)
        leg.txt <- c("GenEr(1,2,3)","GenEr(1,2,4)","GenEr(1,3,5)","GenEr(2,3,4)")
    }    
    if(distr=="invgamma")
    {
        x <- seq(0,5,.01)
        y <- dinvgamma(x, 3/2, 1)
        y2 <- dinvgamma(x, 3/2, 3/2)
        y3 <- dinvgamma(x, 1, 3)
        leg.txt <- c("IG(3/2,1)","IG(3/2,3/2)","IG(1,3)")
    }        
    if(distr=="invtrgamma")
    {
        x <- seq(0,5,.01)
        y <- dinvtrgamma(x, 3, 2, 1)
        y2 <- dinvtrgamma(x, 3, 2, 1/2)
        y3 <- dinvtrgamma(x, 3, 2, 4/3)
        leg.txt <- c("ITG(3,2,1)","ITG(3,2,1/2)","ITG(3,2,4/3)")
    }        
    if(distr=="transgamma")
    {
        x <- seq(0,5,.01)
        y <- dtransgamma(x, 3, 2, 1)
        y2 <- dtransgamma(x, 3, 1/2, 1/3)
        y3 <- dtransgamma(x, 3, 1/2, 4/3)
        leg.txt <- c("TG(3,2,1)","TG(3,1/2,1/3)","TG(3,1/2,4/3)")
    }        
    if(distr=="lgamma")
    {
        x <- seq(0,5,.01)
        y <- dlgamma(x, 0, 2, 1)
        y2 <- dlgamma(x, 1, 2, 1)
        y3 <- dlgamma(x, 0, 2, 2)
        leg.txt <- c("LG(0,2,1)","LG(1,2,1)","LG(0,2,2)")
    }        
    if(distr=="laplace")
    {
        x <- seq(0,5,.01)
        y <- dlaplace(x, 0,1)
        y2 <- dlaplace(x, 0, 2)
        y3 <- dlaplace(x, 0, 3)
        leg.txt <- c("L(0,1)","L(0,2)","L(0,3)")
    }        
    if(distr=="chisq")
    {
        x <- seq(0,5,.01)
        y <- dchisq(x, 2)
        y2 <- dchisq(x, 3)  
        y3 <- dchisq(x, 4)
        y4 <- dchisq(x, 5)
        leg.txt <- c("Chisq(2)","Chisq(3)","Chisq(4)","Chisq(5)")
    }
    if(distr=="ncchisq")
    {
        x <- seq(0,5,.01)
        y <- dchisq(x, 2)
        y2 <- dchisq(x, 2, 1)  
        y3 <- dchisq(x, 4)
        y4 <- dchisq(x, 4, 1)
        leg.txt <- c("Chisq(2)","Chisq(2,1)","Chisq(4)","Chisq(4,1)")
    }
    if(distr=="chi")
    {
        x <- seq(0,5,.01)
        y <- dchi(x, 2)
        y2 <- dchi(x, 3)
        y3 <- dchi(x, 4)
        y4 <- dchi(x, 5)
        leg.txt <- c("Chi(2)","Chi(3)","Chi(4)","Chi(5)")
    }        
            
    if(distr=="weibull")
    {
        x <- seq(0,5,.01)
        y <- dweibull(x, 3, 1)
        y2 <- dweibull(x, 3, 2)
        y3 <- dweibull(x, 4, 2)
        y4 <- dweibull(x, 4, 3)
        leg.txt <- c("W(3,1)","W(3,2)","W(4,2)","W(4,3)")
    }        
    if(distr=="invweibull")
    {
        x <- seq(0,5,.01)
        y <- dinvweibull(x,3,1)
        y2 <- dinvweibull(x, 3, 2)
        y3 <- dinvweibull(x, 4, 2)
        y4 <- dinvweibull(x, 4, 3)
        leg.txt <- c("IW(3,1)","IW(3,2)","IW(4,2)","IW(4,3)")
    }        
    if(distr=="invchisq")
    {
        x <- seq(0,5,.01)
        y <- dinvchisq(x, 2)
        y2 <- dinvchisq(x, 3)
        y3 <- dinvchisq(x, 4)
        y4 <- dinvchisq(x, 2.5)
        leg.txt <- c("IChisq(2)","IChisq(3)","IChisq(4)","IChisq(2.5)")
    }        
    if(distr=="student")
    {
        x <- seq(0,5,.01)
        y <- dt(x, 1, 0)
        y2 <- dt(x, 2, 0)
        y3 <- dt(x, 3, 0)
        y4 <- dt(x, 4, 0)
        leg.txt <- c("T(1)","T(2)","T(3)","T(4)")
    }        
    if(distr=="cauchy")
    {
        x <- seq(0,5,.01)
        y <- dcauchy(x)
        y2 <- dcauchy(x, 1, 1)
        y3 <- dcauchy(x, 1, 1/2)
        y4 <- dcauchy(x, 1, 2)
        leg.txt <- c("Cau(1,1)","Cau(1,1/2)","Cau(1,2)")
    }        
    if(distr=="pareto1")
    {
        x <- seq(0,5,.01)
        y <- dparetoI(x)
        y2 <- dparetoI(x,2,1)
        y3 <- dparetoI(x,2,2)
        y4 <- dparetoI(x,2,3)
        leg.txt  <- c("Pareto1(1,1)","Pareto1(2,1)","Pareto1(2,2)","Pareto1(2,3)")
    }        
    if(distr=="pareto2")
    {
        x <- seq(0,5,.01)
        y <- dpareto2(x,2, 1)
        y2 <- dpareto2(x,2, 2)
        y3 <- dpareto2(x,2, 3)
        y4 <- dpareto2(x,3, 2)
        leg.txt  <- c("Pareto2(2,1)","Pareto2(2,2)","Pareto2(2,3)","Pareto2(3,2)")
    }        
    if(distr=="pareto3")
    {
        x <- seq(0,5,.01)
        y <- dparetoIII(x, 0,1,1)
        y2 <- dparetoIII(x, 1,1,1)
        y3 <- dparetoIII(x, 0,2,1)
        y4 <- dparetoIII(x, 0,1,1.5)        
        leg.txt  <- c("Pareto3(0,1,1)","Pareto3(1,1,1)","Pareto3(0,2,1)","Pareto3(0,1,3/2)")
    }        
    if(distr=="pareto4")
    {
        x <- seq(0,5,.01)
        y <- dparetoIV(x, 0, 1, 1, 1)
        y2 <- dparetoIV(x, 0, 2, 1, 1)
        y3 <- dparetoIV(x, 0, 1, 3/2, 1)
        y4 <- dparetoIV(x, 0, 1, 1, 2)
        leg.txt  <- c("Pareto4(0,1,1,1)","Pareto4(0,2,1,1)","Pareto4(0,1,3/2,1)","Pareto4(0,1,1,2)")
    }        
    if(distr=="invpareto")
    {
        x <- seq(0,5,.01)
        y <- dinvpareto(x, 1, 1)
        y2 <- dinvpareto(x, 2, 1)
        y3 <- dinvpareto(x, 2, 2)
        y4 <- dinvpareto(x, 1, 2)
        leg.txt  <- c("InvPareto(1,1)","InvPareto(2,1)","InvPareto(2,2)","InvPareto(1,2)")
    }        
    if(distr=="burr")
    {
        x <- seq(0,5,.01)
        y <- dburr(x,1,1,1)
        y2 <- dburr(x,2,1,1)
        y3 <- dburr(x,2,2,1)
        y4 <- dburr(x,2,2,2)
        leg.txt <- c("Burr(1,1,1)","Burr(2,1,1)","Burr(2,2,1)","Burr(2,2,2)")
    }        
    if(distr=="gpd")
    {
        x <- seq(0,5,.01)
        y <- dgpd(x,0,1, 0)
        y2 <- dgpd(x,0,1, 1/2)
    y22 <- dgpd(x,0,1, 1)
    y23 <- dgpd(x,0,1, 2)
    y24 <- dgpd(x,0,1, 3)

        y3 <- dgpd(x,0,1, -1/3)
    y32 <- dgpd(x,0,1, -2/3)    
    y33 <- dgpd(x,0,1, -1)    
    y34 <- dgpd(x,0,1, -5/4)    
        leg.txt  <- c("GPD(0)","GPD(1/2)","GPD(1)","GPD(2)","GPD(3)","GPD(-1/3)","GPD(-2/3)","GPD(-1)","GPD(-5/4)")
    }        
    if(distr=="invburr")
    {
        x <- seq(0,5,.01)
        y <- dinvburr(x,1,1,1)
        y2 <- dinvburr(x,1,2,1)
        y3 <- dinvburr(x,2,2,1)
        y4 <- dinvburr(x,1,2,2)    
        leg.txt <- c("InvBurr(1,1,1)","InvBurr(2,1,1)","InvBurr(2,2,1)","InvBurr(2,2,2)")
    }        
    if(distr=="gumbel")
    {
        x <- seq(0,5,.01)
        y <- dgumbel(x, 0, 1)
        y2 <- dgumbel(x, 1/2, 1)
        y3 <- dgumbel(x, 0, 1/2)
        y4 <- dgumbel(x, -1, 2)
        leg.txt <- c("Gum(0,1)","Gum(1/2,1)","Gum(0,1/2)","Gum(-1,2)")
    }        
    


plot(x, y, xlab="x", ylab="f(x)", main="density function", ylim=c(0,.71), col="black", type="l")    

if(length(leg.txt) >= 2)
    lines(x,y2, col="blue")
    if(distr=="GPD")
    {
    lines(x,y22, col="blue", lty=2)
    lines(x,y23, col="blue", lty=3)
    lines(x,y24, col="blue", lty=4)    
    }
    
if(length(leg.txt) >= 3)
    lines(x,y3, col="red")

    if(distr=="GPD")
    {
    lines(x,y32, col="red", lty=2)
    lines(x,y33, col="red", lty=3)
    lines(x,y34, col="red", lty=4)
    }

if(length(leg.txt) >= 4 && distr !="GPD")    
    lines(x,y4, col="green")


if(distr != "GPD")
    legend("topright",leg=leg.txt, col=c("black","blue","red","green"),lty=1)
else
    legend("topright",leg=leg.txt,col=c("black","blue","blue","blue","blue","red","red","red","red"),lty=c(1,1:4,1:4))

    if(distr == "GPD")
        return(cbind(x,y,y2,y22,y23,y24,y3,y32,y33,y34))

    if(length(leg.txt) == 2)
        return(cbind(x,y,y2))
    if(length(leg.txt) == 3)
        return(cbind(x,y,y2,y3))
    if(length(leg.txt) == 4)
        return(cbind(x,y,y2,y3,y4))

}

example.discrete <- function(distr="binom")
{
    x <- 0:10
    if(distr == "binom")
    {
        y <- dbinom(x, 4, 2/3)
        y2 <- dzibinom(x, 0, 1, 1/2)
        y3 <- y
        y3[1] <- .25
        K <- (1-y3[1])/(1-(1/3)^4)
        y3 <- dzibinom(x, 1/10, 1, 1/3)
    }
    if(distr == "pois")
    {
        y <- dpois(x, .5)
        K <- 1/(1-(1/3)^4)
        y2 <- y
        y2[1] <- 0
        y2[2:6] <- 1/(exp(1/2)-1)*.5^(x[-1])/factorial(x[-1])
        
    }
    if(distr == "geom")
    {
        y <- dgeom(x, 1/3)
        y2 <- dgeom(x, 1/3)
        y2 <- c(0, 1/3*(2/3)^(x[-1]-1))
        y3 <- dgeom(x, 1/4)
    }
    if(distr == "nbinom")
    {
        y <- dnbinom(x, 4, 1/2)
        y2 <- dnbinom(x, 4, 1/3)
        y3 <- dnbinom(x, 3, 1/2)
    }
    
    
    if(!missing(y))
        plot(x, y, xlab="k", ylab="P(X=k)", main="mass probability function", ylim=c(0,.3), col="black")
    if(!missing(y2))
        points(x,y2, col="blue")
    if(!missing(y3))
        points(x,y3, col="red")
    
    
    y3[2:5] <- K*y[2:5]
    y3 <- dpois(x, 1)
    y3[1] <- .25
    y3[2:6] <- .9/(exp(1/2)-1)*.5^(x[-1])/factorial(x[-1])
    
    
    if(!missing(y))
    legend("topleft",leg=c("NB(4,1/2)","NB(4,1/3)","NB(3,1/2)"), col=c("black","blue","red"),lty=1)

    if(missing(y3))
        return(cbind(x,y,y2))
    if(!missing(y3))
        return(cbind(x,y,y2,y3))

    
}