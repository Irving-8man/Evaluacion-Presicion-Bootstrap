library(nortest)
library(lmtest)

#Ajustar y modificar

#Función para la comprobacion de normalidad en los residuales
#Sea el vector de los residuales
#Se retorna si existe o no normalidad
ComprobarSupuestoNormalidad <- function (residuales,NivSignicancia = 0.95){
  normalidad <- FALSE
  alfa <- 1-NivSignicancia
  #Creacion de las pruebas
  pValShap = shapiro.test(residuales)$p.value  #Shapiro-Wilk
  pValLILLIE = lillie.test(residuales)$p.value   #Lilliefort
  #Resultados
  pValores = c(pValShap, pValLILLIE)
  pValminimo = min(pValores)
  
  if(alfa<pValminimo) normalidad <- TRUE
  return(normalidad)
}


#Función para la comprobacion de varianza constante
#Sea modelo_lineal el modelo creado
#Se retorna si existe o no varianza constante
ComprobarVarianzaConstante <- function(z,modeloLineal,NivSignicancia = 0.95){
  varianzaConstante = FALSE
  alfa <- 1-NivSignicancia
  
  #Pruebas
  datosBrush <- bptest(modeloLineal)$p.value#Brush pagan
  datosWhite <- bptest(modeloLineal, varformula = ~ I(z^2))$p.value #White
    
  #Resultados
  pValores = c(datosBrush,datosWhite)
  pValminimo = min(pValores)
  
  if(alfa<pValminimo) varianzaConstante <- TRUE
  
  return (varianzaConstante)
}