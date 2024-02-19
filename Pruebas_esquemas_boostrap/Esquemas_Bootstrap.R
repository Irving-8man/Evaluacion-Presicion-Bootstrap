
# #Funcion de boostrap simple
# CalcularMuestrasBootstrapSimple <- function(z, yAjustados, e ,n, B=100){
#   vectorModelosBoostrap <- c()
#   vectorR2Bootstrap <-c()
#   
#   #Fabricar las muestras
#   for (i in 1:B) {
#     residualesBootstrap <- sample(e, n, replace =TRUE) #remplazo
#     ySimuladas <- yAjustados + residualesBootstrap
#     nuevoModeloBootstrap <- lm(ySimuladas ~ z)
#     nuevoRCuadradoBootstrap <- summary(nuevoModeloBootstrap)$r.squared
#     
#     vectorModelosBoostrap[i] <- i
#     vectorR2Bootstrap[i] <- nuevoRCuadradoBootstrap
#   }
#   muestrasR2Bootstrap <- data.frame(vectorModelosBoostrap,vectorR2Bootstrap)
#   return(muestrasR2Bootstrap)
# }

#Auxiliar de implementación de los esquemas Roboustos Bootstrap
#Numero de esquemas, parte 1
#return=> c(1:B)
#esquemas => (1:wu1, 2:wu2, 3:wu3, 4:liu1, 5:liu2, 6:normal)
ImplementarRemuestreosBootsR <- function(z,y,B,tipoEsquema){
  constantePeso <- 3
  
  #Proceso regresion lineal simple
  modeloLineal <- lm(y~z)
  residuales <- residuals(modeloLineal)
  coeficientes <- coef(modeloLineal)
  yAjustados <- fitted(modeloLineal)
  nResiduales <- length(residuales)
  originalR2 <- summary(modeloLineal)$r.squared
  #CME <- summary(modeloLineal)$sigma**2
  hii <- hatvalues(modeloLineal)
  
  #Proceso regresión robusta
  modeloLinealRobusto <- lmrob(y ~ z, method = "MM")
  residualesRobustos <- modeloLinealRobusto$residuals
  nResidualesRobustos <- length(residualesRobustos)
  yAjustadosRobustos <- modeloLinealRobusto$fitted.values
  CMERobusto <- modeloLinealRobusto$scale**2
  
  #Ponderacion de los residulales
  #Propuesto
  x=abs(residualesRobustos)/sqrt(CMERobusto)
  w <- rep(1,nResidualesRobustos)
  xx <- which(x > constantePeso) #indices paso 2
  w[xx] <- (constantePeso / w[xx])# Actualizar los valores de W para los casos donde x > 3,paso3
  residualesRobustosPonderados <- w*residualesRobustos
  
  Remuestra <- function(residualesRobustosPonderados,tipoEsquema){
      switch(tipoEsquema,
        'Wu1' = CalcularMuestrasBootstrapWu1(residualesRobustosPonderados),
        'Wu2' = CalcularMuestrasBootstrapWu2(residualesRobustosPonderados),
        'Wu3' = CalcularMuestrasBootstrapWu3(residualesRobustosPonderados),
        'Liu1' = CalcularMuestrasBootstrapLiu1(residualesRobustosPonderados),
        'Liu2' = CalcularMuestrasBootstrapLiu2(residualesRobustosPonderados),
        'Wild' = CalcularMuestrasBootstrapWild(residualesRobustos),
        
        stop()
        )
  }
  
  return(Remuestra(residualesRobustosPonderados,tipoEsquema))
}








#Funcion con el esquema Wu 1
CalcularMuestrasBootstrapWu1 <- function(residualesRobustosPonderados){
  muestrasBootstrapWu1 <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    tt <- rnorm(nResidualesRobustos)
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <- summary(modeloBoots)$r.squared
    muestrasBootstrapWu1[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapWu1)
}


#Funcion con el esquema Wu 2
CalcularMuestrasBootstrapWu2 <- function(residualesRobustosPonderado){
  muestrasBootstrapWu2 <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    ai <- (residuales-mean(residuales))/sd(residuales)
    tt <- sample(ai,replace=T)
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWu2[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapWu2)
}

#Funcion con el esquema Wu 3
CalcularMuestrasBootstrapWu3 <- function(residualesRobustosPonderados){
  muestrasBootstrapWu3 <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    mediana <- median(residualesRobustosPonderados)
    NMAD <- (1/0.6745)*median( abs(residualesRobustosPonderados-mediana) )
    Rai <- (residualesRobustosPonderados-mediana)/NMAD
    tt <- sample(Rai,replace=T)
    
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWu3[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapWu3)
}
#Funcion con el esquema Liu 1
CalcularMuestrasBootstrapLiu1 <- function(residualesRobustosPonderados){
  muestrasBootstrapLui1 <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    tt <- rgamma(nResidualesRobustos,2,4)
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapLui1[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapLui1)
}

#Funcion con el esquema Liu 2
CalcularMuestrasBootstrapLiu2 <- function(residualesRobustosPonderados){
  muestrasBootstrapLui2 <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    media1 <- 0.5*sqrt(17/6)+sqrt(1/6)
    media2 <- 0.5*sqrt(17/6)-sqrt(1/6)
    H <- rnorm(nResidualesRobustos,media1,sqrt(0.5))
    D <- rnorm(nResidualesRobustos,media2,sqrt(0.5))
    tt <- H*D-media1*media2
   
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapLui2[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapLui2)
}

#Funcion con el esquema Wild 
CalcularMuestrasBootstrapWild <- function(residualesRobustos){
  muestrasBootstrapWild <- matrix(nrow=B,ncol=2)
  for (i in 1:B) {
    residualesBootstrap <- sample(residualesRobustos, replace =TRUE) #remplazo
    yBoots <- yAjustadosRobustos + residualesBootstrap
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWild[i,] <- c(i,r2Boots)
  }
  return(muestrasBootstrapWild)
}



