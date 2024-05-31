#Auxiliar de implementación de los esquemas Roboustos Bootstrap
# funciont <- (muestrasBotss, nivConfianza,tipoIntervalo)
#intervalos => (1:Percentil, 2:Bootstrap-t, 3:BootstrapIterativo, 4:liu1, 5:liu2, 6:normal)
# return => c(1:2)
ImplementarIntervalosConfianza <- function(data,originalR2,B=100,muestrasR2Boot,nivConfianza=0.95,tipoIntervalo){
    switch(tipoIntervalo,
       'Perc' <- CalcularIntervaloConfianzaPercentil(originalR2,muestrasR2Boot,nivConfianza=0.95),
       'BootT' <- CalcularIntervaloConfianzaBootstrapT(originalR2,muestrasR2Boot,nivConfianza=0.95),
       'Iter' <- CalcularIntervaloConfianzaBootstrapIt(originalR2,muestrasR2Boot,nivConfianza=0.95),
       
       'BCa' <- CalcularIntervaloConfianzaBootstrapBCa(data,originalR2,B,muestrasR2Boot,nivConfianza=0.95),
       'ABC' <- CalcularIntervaloConfianzaBootstrapABC(muestrasR2Boot,nivConfianza=0.95),
       
       'Pond' <- CalcularIntervaloConfianzaBootstrapPond(muestrasR2Bootstrap,nivConfianza),
       
       'Simt'<-CalcularIntervaloConfianzaBootstrapS(originalR2,muestrasR2Boot,nivConfianza=0.95),
       
       stop("Intervalo no válido")
    )
}


#Funcion para obtener el intervalo de confianza-Percentil
CalcularIntervaloConfianzaPercentil<- function(originalR2,muestrasR2Boot,nivConfianza=0.95){
  alpha <- 1-nivConfianza
  vectorR2Bootstrap <- muestrasR2Boot
  n <- length(vectorR2Bootstrap)
  media <- mean(vectorR2Bootstrap)
  puntosCriticos <- quantile(vectorR2Bootstrap, c(alpha/2, 1 - alpha/2)) # Aproximación bootstrap de los puntos críticos
  ICInfBootP <- media - puntosCriticos[2] / sqrt(n)
  ICSupBootP <- media - puntosCriticos[1] / sqrt(n)
  intervaloConfianzaPercentil <-as.vector(c(ICInfBootP, ICSupBootP))
  return(intervaloConfianzaPercentil)
}


#Funcion para calcular el intervalo de confianza-Bootstrap-T
CalcularIntervaloConfianzaBootstrapT <- function(originalR2,muestrasR2Boot,nivConfianza=0.95){
  alpha <- 1-nivConfianza
  vectorR2Bootstrap <- muestrasR2Boot
  n <- length(vectorR2Bootstrap)
  media <-mean(vectorR2Bootstrap)
  puntosCriticos <- quantile(vectorR2Bootstrap, c(alpha/2, 1 - alpha/2)) # Aproximación bootstrap de los puntos críticos
  ICInfBootT <- media - puntosCriticos[2] * sd(vectorR2Bootstrap)/sqrt(n)
  ICSupBootT <- media - puntosCriticos[1] * sd(vectorR2Bootstrap)/sqrt(n)
  intervaloConfianzaBootT <- as.vector(c(ICInfBootT, ICSupBootT))
  return(intervaloConfianzaBootT)
}

#reconsiderar
#Funcion para calcular el intervalo de confianza con bootstrap iterado
CalcularIntervaloConfianzaBootstrapIt <- function(originalR2,B,muestrasR2Boot,remuestrasBoot,nivConfianza=0.95){
  #Paso 1: encuentra theta^ usando x = (X1,...Xn) = orginalR2
  theta_hat <- originalR2
  
  #Paso 2: extraer remuestras x^*1...x^*B de x . Para cada remuestra calcular
  #theta^*1...theta^*B. aqui tendriamos los R^*1 y los x^*B = (tt*residualesRobustosPonderados)/(sqrt(1-hii))
  x_star <- remuestrasBoot #sean los residuales
  theta_star <-muestrasR2Boot
  
  
  #Paso 3: Para cada uno de las B remuestras x^*1...x^*B de x, remuestrear B1 o C veces
  #para obtener  remuestras x^** 11 ... x^** 1C ... x^**B1 ... x^** BC. Para cada conjunto
  #de remuestras internas C  x^** b1...x^** bC, b= 1...B calcular sus versiones 
  #theta^** b1... theta^** bC, donde el estimado bootstrap de P(theta^** <= theta^ | x^* b)
  #es la proporción de los C valores theta^** b1... theta^** bC que son menores o iguales a 
  #theta^. 
  #Elija varios niveles l1, l2, l3... cerca del nivel alfa deseado y cuente el número 
  #de veces que se cumple la condición i para las remuestras B. 
  #La proporción de remuestras B para las cuales se cumple la condición i es una aproximación a
  # pi^(li), i=1,2..
  C <- 100
  for(b in 1:B){
    
  }
  
  #Paso 4: obten un aproximado valor de delta^alpha que cumple pi^(delta^alpha) = alpha para la
  #interpilación entre (li,pi^(delta^alpha))
  
  #Paso 5: el intervalo de confianza del bootstrap interativo es el metodo percentil del 
  #nivel nominal  delta^alpha construido usando los theta^*1...theta^*B del paso 2
  
  #I1(a) = [ theta^* B,[(1 - delta_alpha)B/2] + 1, theta^* B,[(1 + delta_alpha)B/2] + 1] 
}

#Que bondad tendria 
#Funcion para calcular el intervalo de confianza con pecercentil-simetrizado
CalcularIntervaloConfianzaBootstrapS <- function(originalR2,muestrasR2Boot,nivConfianza=0.95){
  alpha <- 1-nivConfianza
  vectorR2Bootstrap <- muestrasR2Boot
  cuasi_dt <- sd(vectorR2Bootstrap)
  n <- length(vectorR2Bootstrap)
  media <-mean(vectorR2Bootstrap)
  puntosCriticos <- quantile(vectorR2Bootstrap, 1 - alpha)# Aproximación bootstrap de los puntos críticos
  # Construcción del IC
  ICInfBootS <- media - puntosCriticos * cuasi_dt/sqrt(n)
  ICSupBootS <- media + puntosCriticos * cuasi_dt/sqrt(n)
  intervaloConfianzaBootS <- as.vector(c(ICInfBootS, ICSupBootS))
  return(intervaloConfianzaBootS)
}




#####################################

#Intervalo de confianza-BCa
CalcularIntervaloConfianzaBootstrapBCa <- function(data,originalR2,B,muestrasR2Boot,nivConfianza=0.95){
  alfa<-1-nivConfianza
  z <- data[[1]]
  y <- data[[2]]
  vectorR2Bootstrap <- muestrasR2Boot
  n <- length(y)
  p0 <- length(which(vectorR2Bootstrap >= originalR2))/B #proporcion de e Rˆ2i’s ≥ Rˆ2 duda
  vectorR2i <- numeric(n)

  #Obteniendo las Rˆ2−i
  for (i in 1:n){
    vectorR2i[i] <- summary(lm(y[-i]~z[-i]))$r.squared
  }
  
  #Calculo de Rˆ2−iprom
  R2iPromedio <- mean(vectorR2i)
  
  #Calculando las diferencias Rˆ2−iprom - Rˆ2−i
  diferenciasR2iPi <- R2iPromedio - vectorR2i
  
  #suma de las diferencias elevadas a cudrado y al cubo
  sumaDiferenciasCuad <- sum(diferenciasR2iPi**2)
  sumaDiferenciasCubo <- sum(diferenciasR2iPi**3)
  
  #Constante de aceleracion
  a <- sumaDiferenciasCubo/(6*(sumaDiferenciasCuad**1.5))
  
  ZL <- qnorm(1-p0) + (qnorm(1-p0) - qnorm(1-alfa/2) )/(1-a *(qnorm(1-p0) - qnorm(1-alfa/2)))
  ZU <- qnorm(1-p0) + (qnorm(1-p0) + qnorm(1-alfa/2) )/(1-a *(qnorm(1-p0) + qnorm(1-alfa/2))) 

  ICInfBootBCa <- quantile(vectorR2Bootstrap,probs= pnorm(ZL))
  ICSupBootBCa <- quantile(vectorR2Bootstrap,probs= pnorm(ZU))
  intervaloConfianzaBootBCa <- as.vector(c(ICInfBootBCa, ICSupBootBCa))
  return(intervaloConfianzaBootBCa)
}

#Funcion para calcular el intervalo de confianza con ABC
CalcularIntervaloConfianzaBootstrapABC <- function(muestrasR2Boot,nivConfianza=0.95){
  vectorR2Bootstrap <- muestrasR2Boot
  fabc <-function(x, w) w%*%x #ABC usando sus pesos
  intervalo<- abc.ci(vectorR2Bootstrap, fabc, conf = nivConfianza) #ABC method C.I.
  intervaloConfianzaBootABC <- c(intervalo[2],intervalo[3])
  return(intervaloConfianzaBootABC)
}



#Funcion para calcualr el intervalo de confianza con Ponderado
CalcularIntervaloConfianzaBootstrapPond <- function(muestrasR2Boot,nivConfianza=0.95){
#
}
