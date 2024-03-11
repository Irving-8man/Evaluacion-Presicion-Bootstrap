library(nortest)
library(lmtest)

#Funcion para la comprobacion de normalidad en los residuales
ComprobarSupuestoNormalidad <- function (residuales, nivConfianza = 0.95){
  normalidad <- FALSE
  alpha <- 1-nivConfianza
  #Creacion de las pruebas
  pValShap <- shapiro.test(residuales)$p.value  #Shapiro-Wilk
  ValCShap = shapiro.test(residuales)$statistic
  
  pValLILLIE <- lillie.test(residuales)$p.value   #Lilliefort
  ValCLILLIE = lillie.test(residuales)$statistic
  
  
  #Resultados
  pValores <- c(pValShap, pValLILLIE)
  ValorCal <- c(ValCShap, ValCLILLIE)
  
  #Construccion de la tabla
  estadistica <- c("Shapiro-Wilk","Liliefort")
  tablaNormal <- data.frame(estadistica, pValores, ValorCal)
  cat("Resultados para las pruebas de normalidad\n")
  print(tablaNormal)
  
  pValminimo <- min(pValores)
  if(alpha < pValminimo) normalidad <- TRUE
  return(normalidad)
}


#Funcion para la comprobacion de varianza constante
ComprobarVarianzaConstante <- function(z,modeloLineal, nivConfianza = 0.95){
  varianzaConstante <- FALSE
  alpha <- 1-0.95
  #Pruebas
  datosBrush <- bptest(modeloLineal)$p.value#Brush pagan
  ValBrush <- bptest(modeloLineal)$statistic
  
  datosWhite <- bptest(modeloLineal, varformula = ~ I(z^2))$p.value #White
  ValWhite <- bptest(modeloLineal, varformula = ~ I(z^2))$statistic
  
  #Resultados
  pValores <- c(datosBrush,datosWhite)
  ValorCal <- c(ValBrush,ValWhite)
  
  #Construccion de la tabla
  estadistica = c("Brush pagan","White")
  tablaVar <- data.frame(estadistica, pValores, ValorCal)
  cat("\nResultados para las pruebas de homocedasticidad\n")
  print(tablaVar)
  
  pValminimo <- min(pValores)
  if(alpha < pValminimo) varianzaConstante <- TRUE
  
  return (varianzaConstante)
}


#Funcion para comprobar la media cero
ComprobarMediaCero <- function(residuales, nivConfianza = 0.95){
  mediaCero <- FALSE
  alpha <- 1-nivConfianza
  #Prueba
  pValor <-t.test(residuales)$p.value #test T-student
  ValorCal <- t.test(residuales)$statistic

  #Construccion de la tabla
  estadistica = c("T-Student")
  tablaVar <- data.frame(estadistica, pValor, ValorCal)
  cat("\nResultados para la prueba de T-Student para media cero\n")
  print(tablaVar)
  
  if(alpha < pValor) mediaCero <- TRUE
  
  return (mediaCero)
}


#Funcion para comprobar independencia
ComprobarIndependencia <- function(modeloLineal,nivConfianza= 0.95){
  independencia <- FALSE
  alpha <- 1-nivConfianza
  #Pruebas
  pValor <- dwtest(modeloLineal)$p.value #Durbin-Watson test
  ValorCal <- dwtest(modeloLineal)$statistic
  
  #Construccion de la tabla
  estadistica = c("Durbin-Watson test")
  tablaVar <- data.frame(estadistica, pValor, ValorCal)
  cat("\nResultados para la prueba de Durbin-Watson test para independencia\n")
  print(tablaVar)
  
  if(alpha < pValor) independencia <- TRUE
  
  return (independencia)
}

