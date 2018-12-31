library(Ryacas)
library(tm)

a.leg.poly = function(l,m){
  
  x <- Sym("x")
  
  easy.part = ((-1)^m)/(2^l * factorial(l)) * (1 - x^2)^(m/2)
  deriv.part = (x^2 - 1)^l
  
  for(i in 1:(l+m)){
    deriv.part = deriv(deriv.part, x)
  } 
  
  out = easy.part*deriv.part
  return(out)
}

plm.positive <- function(l,m, x.entered ){
  
  if(l==0 & m==0){
    return(1)
  } 
  
  easy.part <- ((-1)^m)/(2^l * factorial(l)) * (1 - x.entered^2)^(m/2)
  x <- Sym("x")
  
  deriv.part <- (x^2 - 1)^l
  
  for(i in 1:(l+m)){
    deriv.part <- deriv(deriv.part, x)
  }  
  
  eval.deriv.part <- capture.output(deriv.part)
  eval.deriv.part <- gsub("expression", "", eval.deriv.part)
  eval.deriv.part <- gsub("x", x.entered, eval.deriv.part)
  eval.deriv.part <- eval(parse(text=eval.deriv.part))
  
  out <- easy.part*eval.deriv.part
  
  return(out)
  
}



plm <- function(l,m, x.entered ){
  
  if(l==0 & m==0){
    return(1)
  }
  
  if( m>=0 ){
    out <- plm.positive(l=l,m=m, x.entered=x.entered)
    return(out)
  } 
  
  if( m<0 ){
    
    m.pos <- abs(m)
    
    plm.part <- plm.positive(l=l, m=m.pos, x.entered=x.entered)
    
    other.part <- (-1)^( m.pos ) * ( factorial(l - m.pos) / factorial(l + m.pos) )
    
    out <- other.part * plm.part
    
    return(out)
  }
  
}



spherical.harmonic <- function(l, m, theta, phi){
  
  ylm <- 0
  
  if( m < 0){
    ylm <- sqrt(2) *
      sqrt( ((2*l + 1)/(4*pi)) * ((factorial(l - abs(m))) / (factorial(l + abs(m)))) ) *
      plm(l, abs(m), cos(theta)) * sin( abs(m) * phi )
    
  } else if (m == 0){
    ylm <-sqrt((2*l + 1) / (4*pi)) * plm(l, m, cos(theta))
    
  } else if( m > 0){
    ylm <- sqrt(2) *
      sqrt( ((2*l + 1)/(4*pi)) * ((factorial(l - m)) / (factorial(l + m))) ) *
      plm(l, m, cos(theta)) * cos( m*phi )
  }
  
  if(m %% 2 == 0){
    ylm <- ylm
  } else {
    ylm <- -1*ylm
  }
  
  return(ylm)
}
