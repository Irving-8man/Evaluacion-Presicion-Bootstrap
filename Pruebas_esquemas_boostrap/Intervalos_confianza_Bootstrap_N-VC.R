

#Funcion para obtener el intervalo de confianza de R2 de un modelo
#con N-VC

CalcularIntervalosConfianzaBootstrapNVC<- function(muestrasR2Bootstrap,B, nivSignicancia=0.95){
  alfa <- 1-nivSignicancia
  vectorR2Bootstrap <- muestrasR2Bootstrap[[2]]
  #Metodo percentil
  intervaloConfianzaPercentil <- mean(vectorR2Bootstrap) + qt(nivSignicancia, B - 1) * sd(vectorR2Bootstrap) * c( - 1, +1)/sqrt(B) 
  
  #Metodo Bootstrap-T
  return(intervaloConfianzaPercentil)
}



#Funcion para calcular el intervalo de confianza con Bootstrap-T
CalcularIntervalosConfianzaBootstrapT <- function(vectorR2Bootstrap,B,originalR2,nivSignicancia){
}