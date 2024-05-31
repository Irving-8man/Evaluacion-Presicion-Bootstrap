#Funcion auxiliar de implementación de los esquemas Roboustos Bootstrap
#Numero de esquemas, parte 1
#esquemas => (1:wu1, 2:wu2, 3:wu3, 4:liu1, 5:liu2, 6:normal)
#return=> c(1:B)
ImplementarRemuestreosBootsR <- function(z,residuales,residualesRobustos,residualesRP,yAjRob,hii,B,tipo){
  switch(tipo,
  'Wu1' <- CalcularMuestrasBootstrapWu1(z,residualesRP,yAjRob,hii,B),
  'Wu2' <- CalcularMuestrasBootstrapWu2(z,residuales,residualesRP,yAjRob,hii,B),
  'Wu3' <- CalcularMuestrasBootstrapWu3(z,residualesRP,yAjRob,hii,B),
  'Liu1' <- CalcularMuestrasBootstrapLiu1(z,residualesRP,yAjRob,hii,B),
  'Liu2' <- CalcularMuestrasBootstrapLiu2(z,residualesRP,yAjRob,hii,B),
  'Wild' <- CalcularMuestrasBootstrapWild(z,residualesRobustos,yAjRob,B),
  stop("Esquema de remuestreo no válido")
  )
}


#Funcion con el esquema Wu 1
CalcularMuestrasBootstrapWu1 <- function(z,residualesRP,yAjRob,hii,B=100,caso=1){
  muestrasBootstrapWu1 <- numeric(B)
  nRRP <-length(residualesRP)
  for (i in 1:B) {
    tt <- rnorm(nRRP)
    yBoots <- yAjRob + (tt*residualesRP)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <- summary(modeloBoots)$r.squared
    muestrasBootstrapWu1[i] <- c(r2Boots)
  }
  return(muestrasBootstrapWu1)
}


#Funcion con el esquema Wu 2
CalcularMuestrasBootstrapWu2 <- function(z,residuales,residualesRP,yAjRob,hii,B=100){
  muestrasBootstrapWu2 <- numeric(B)
  for (i in 1:B) {
    ai <- (residuales-mean(residuales))/sd(residuales)
    tt <- sample(ai,replace=T)
    yBoots <- yAjRob + (tt*residualesRP)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWu2[i] <- c(r2Boots)
  }
  return(muestrasBootstrapWu2)
}

#Funcion con el esquema Wu 3
CalcularMuestrasBootstrapWu3 <- function(z,residualesRP,yAjRob,hii,B=100){
  muestrasBootstrapWu3 <- numeric(B)
  for (i in 1:B) {
    mediana <- median(residualesRP)
    NMAD <- (1/0.6745)*median( abs(residualesRP-mediana) )
    Rai <- (residualesRP-mediana)/NMAD
    tt <- sample(Rai,replace=T)
    
    yBoots <- yAjRob + (tt*residualesRP)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWu3[i] <- c(r2Boots)
  }
  return(muestrasBootstrapWu3)
}

#Funcion con el esquema Liu 1
CalcularMuestrasBootstrapLiu1 <- function(z,residualesRP,yAjRob,hii,B=100){
  muestrasBootstrapLui1 <- numeric(B)
  nRRP <-length(residualesRP)
  for (i in 1:B) {
    tt <- rgamma(nRRP,2,4)
    yBoots <- yAjRob + (tt*residualesRP)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapLui1[i] <- c(r2Boots)
  }
  return(muestrasBootstrapLui1)
}

#Funcion con el esquema Liu 2
CalcularMuestrasBootstrapLiu2 <- function(z,residualesRP,yAjRob,hii,B=100){
  muestrasBootstrapLui2 <- numeric(B)
  nRRP <-length(residualesRP)
  for (i in 1:B) {
    media1 <- 0.5*sqrt(17/6)+sqrt(1/6)
    media2 <- 0.5*sqrt(17/6)-sqrt(1/6)
    H <- rnorm(nRRP,media1,sqrt(0.5))
    D <- rnorm(nRRP,media2,sqrt(0.5))
    tt <- H*D-media1*media2
   
    yBoots <- yAjRob + (tt*residualesRP)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapLui2[i] <- c(r2Boots)
  }
  return(muestrasBootstrapLui2)
}

#Funcion con el esquema Wild 
CalcularMuestrasBootstrapWild <- function(z,residualesRobustos,yAjRob,B=100){
  muestrasBootstrapWild <- numeric(B)
  for (i in 1:B) {
    residualesBootstrap <- sample(residualesRobustos, replace =TRUE) #remplazo
    yBoots <- yAjRob + residualesBootstrap
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    muestrasBootstrapWild[i] <- c(r2Boots)
  }
  return(muestrasBootstrapWild)
}


#Funcion de apoyo para ponderar residuales
PoderarResidualesRobustos <-function(residualesRobustos,CMERobusto){
  constantePeso <- 3
  n <- length(residualesRobustos)
  x <- abs(residualesRobustos)/sqrt(CMERobusto)
  w <- rep(1,n)
  xx <- which(x > constantePeso) 
  w[xx] <- (constantePeso / w[xx])
  residualesRobustosPonderados <- w*residualesRobustos
  return(residualesRobustosPonderados)
}

