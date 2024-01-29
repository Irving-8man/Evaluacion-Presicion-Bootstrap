
#Funcion de boostrap simple
CalcularMuestrasBootstrapSimple <- function(z, yAjustados, e ,n, B=100){
  vectorModelosBoostrap <- c()
  vectorR2Bootstrap <-c()
  
  #Fabricar las muestras
  for (i in 1:B) {
    residualesBootstrap <- sample(e, n, replace =TRUE) #remplazo
    ySimuladas <- yAjustados + residualesBootstrap
    nuevoModeloBootstrap <- lm(ySimuladas ~ z)
    nuevoRCuadradoBootstrap <- summary(nuevoModeloBootstrap)$adj.r.squared
    
    vectorModelosBoostrap[i] <- i
    vectorR2Bootstrap[i] <- nuevoRCuadradoBootstrap
  }
  muestrasR2Bootstrap <- data.frame(vectorModelosBoostrap,vectorR2Bootstrap)
  return(muestrasR2Bootstrap)
}



#Funcion con el esquema Wu 1
CalcularMuestrasBootstrapWu1 <- function(z,B=100,nResidualesRobustos,yAjustadosRobustos,residualesRobustosPonderados,modeloLineal){
  muestrasBootstrapWu1 <- matrix(nrow=B,ncol=3)
  hii <- hatvalues(modeloLineal)
  for (i in 1:B) {
    tt <- rnorm(nResidualesRobustos)
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <- summary(modeloBoots)$r.squared
    CMEBoots <- summary(modeloBoots)$sigma**2
    muestrasBootstrapWu1[i,] <- c(i,r2Boots,CMEBoots)
  }
  return(muestrasBootstrapWu1)
}


#Funcion con el esquema Wu 2
CalcularMuestrasBootstrapWu2 <- function(z,B=100,yAjustadosRobustos,residualesRobustosPonderados,modeloLineal,residuales){
  muestrasBootstrapWu2 <- matrix(nrow=B,ncol=3)
  hii <- hatvalues(modeloLineal)
  for (i in 1:B) {
    ai <- (residuales-mean(residuales))/sd(residuales)
    tt <- sample(ai,replace=T)
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    CMEBoots <- summary(modeloBoots)$sigma**2
    muestrasBootstrapWu2[i,] <- c(i,r2Boots,CMEBoots)
  }
  return(muestrasBootstrapWu2)
}

#Funcion con el esquema Wu 3
CalcularMuestrasBootstrapWu3 <- function(z,B=100,yAjustadosRobustos,residualesRobustosPonderados,modeloLineal){
  muestrasBootstrapWu3 <- matrix(nrow=B,ncol=3)
  hii <- hatvalues(modeloLineal)
  for (i in 1:B) {
    mediana <- median(residualesRobustosPonderados)
    NMAD <- (1/0.6745)*median( abs(residualesRobustosPonderados-mediana) )
    Rai=(residualesRobustosPonderados-mediana)/NMAD
    tt=sample(Rai,replace=T)
    
    yBoots <- yAjustadosRobustos + (tt*residualesRobustosPonderados)/(sqrt(1-hii))
    modeloBoots <- lm(yBoots~z)
    r2Boots <-summary(modeloBoots)$r.squared
    CMEBoots <-summary(modeloBoots)$sigma**2
    muestrasBootstrapWu3[i,] <- c(i,r2Boots,CMEBoots)
  }
  return(muestrasBootstrapWu3)
}
#Funcion con el esquema Liu 1
CalcularMuestrasBootstrapLiu1 <- function(){
  
}

#Funcion con el esquema Liu 2
CalcularMuestrasBootstrapLiu1 <- function(){
  
}

#Funcion con el esquema Wild 
CalcularMuestrasBootstrapWild <- function(){
  
}




