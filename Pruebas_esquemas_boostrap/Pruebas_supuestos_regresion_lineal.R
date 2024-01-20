library(nortest)

#Función para la comprobacion de normalidad en los residuales
#Sea e el vector de los residuales
#Se retorna si existe o no normalidad
ComprobarSupuestoNormalidad <- function (e,NivSignicancia = 0.95){
  normalidad <- FALSE
  alfa <- 1-NivSignicancia
  
  #Creacion de las pruebas
  #Shapiro-Wilk
  pValShap = shapiro.test(e)$p.value
  
  #Cramer-von-misses
  pValCVM = cvm.test(e)$p.value
  
  #Anderson-Darling
  pValAD = ad.test(e)$p.value
  
  #Shapiro-francia
  pValSF = sf.test(e)$p.value
  
  #Lilliefort
  pValLILLIE = lillie.test(e)$p.value
  
  #Resultados
  pValores = c(pValShap, pValCVM, pValAD, pValSF, pValLILLIE)
  pValminimo = min(pValores)
  
  if(alfa<pValminimo) normalidad <- TRUE
    
  return(normalidad)
}



ComprobarVarianzaConstante <- function(){
  varianzaConstante = FALSE
  
  return (varianzaConstante)
}