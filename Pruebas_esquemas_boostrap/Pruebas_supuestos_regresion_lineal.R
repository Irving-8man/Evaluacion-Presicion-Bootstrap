library(nortest)
library(lmtest)
library(car)



#Función para la comprobacion de normalidad en los residuales
#Sea el vector de los residuales
#Se retorna si existe o no normalidad
ComprobarSupuestoNormalidad <- function (residuales,NivSignicancia = 0.95){
  normalidad <- FALSE
  alfa <- 1-NivSignicancia
  
  #Creacion de las pruebas
  pValShap = shapiro.test(residuales)$p.value  #Shapiro-Wilk
  pValCVM = cvm.test(residuales)$p.value        #Cramer-von-misses
  pValAD = ad.test(residuales)$p.value         #Anderson-Darling
  pValSF = sf.test(residuales)$p.value        #Shapiro-francia
  pValLILLIE = lillie.test(residuales)$p.value   #Lilliefort
  
  #Resultados
  pValores = c(pValShap, pValCVM, pValAD, pValSF, pValLILLIE)
  pValminimo = min(pValores)
  
  if(alfa<pValminimo) normalidad <- TRUE
    
  return(normalidad)
}


#Función para la comprobacion de varianza constante
#Sea modelo_lineal el modelo creado
#Se retorna si existe o no varianza constante
ComprobarVarianzaConstante <- function(modeloLineal,NivSignicancia = 0.95){
  varianzaConstante = FALSE
  alfa <- 1-NivSignicancia
  
  #Pruebas
  datosBrush <- bptest(modeloLineal)#Brush pagan
  datosScore <- ncvTest(modeloLineal)#con score test
  
  #Resultados
  pValores = c(datosBrush[[4]], datosScore[[5]])
  pValminimo = min(pValores)
  
  if(alfa<pValminimo) varianzaConstante <- TRUE
  
  return (varianzaConstante)
}