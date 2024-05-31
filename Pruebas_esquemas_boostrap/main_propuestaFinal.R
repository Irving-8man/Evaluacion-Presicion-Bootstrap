data0 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_0.csv")
data1 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_1.csv")
data2 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_2.csv")
data3 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_3.csv")
data4 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/N-VC.csv")
data <- data0

#Cargando librerias
library("robustbase")
library("boot")
library("bootstrap")

#Carga de las funciones de pruebas de supuestos
source("Pruebas_esquemas_boostrap/Pruebas_supuestos_regresion_lineal.R")
source("Pruebas_esquemas_boostrap/Esquemas_Bootstrap.R")
source("Pruebas_esquemas_boostrap/Intervalos_confianza_Bootstrap.R")

#Funcion principal
CalcularPrecision <- function(data, alpha=0.05, nivConfianza=0.95){
  z <- data[[1]]
  y <- data[[2]]
  
  #Remuestras bootstrap necesarias
  B <- 100 
  
  modeloLineal <- lm(y~z)
  residuales <- residuals(modeloLineal)
  originalR2 <- summary(modeloLineal)$r.squared
  hii <- hatvalues(modeloLineal)
 
  #Comprobar normalidad, homocedasticidad, media cero e idependendia de los residuales
  hayNormalidad <- ComprobarSupuestoNormalidad(residuales, nivConfianza)
  hayHomocedasticidad <- ComprobarVarianzaConstante(z, modeloLineal, nivConfianza)
  esMediaCero <- ComprobarMediaCero(residuales, nivConfianza)
  esIndependiente <- ComprobarIndependencia(modeloLineal, nivConfianza)
  print(paste(hayNormalidad," ",hayHomocedasticidad," ",esMediaCero," ",esIndependiente))
  caso <- 0
  
  #
  #Primer caso si cumple los 4 supuestos 
  if(hayNormalidad && hayHomocedasticidad && esMediaCero && esIndependiente){
    # se utiliza el estimador de minimos cuadrados
    #implmentar los 6 esquemas de remuestreo
    #calcular y retornar los intervalos a este caso;percentil, boots-t e iterativo
    caso <- 1
    yAju <- fitted(modeloLineal)
    residualesRP <- residuales
    residualesRobustos <- residuales
    yAjRob <- yAju
    
  }else{
    
    modeloLinealRobusto <- lmrob(y ~ z, method = "MM")
    residualesRobustos <- modeloLinealRobusto$residuals
    yAjRob <- modeloLinealRobusto$fitted.values
    CMERobusto <- modeloLinealRobusto$scale**2
    
    #Segundo caso, cuando no cumple normalidad
    if( !hayNormalidad && esIndependiente && hayHomocedasticidad){
      # se utiliza el estimador robusto MM  
      # los residuales robustos sin ponderacion
      #implmentar los 6 esquemas de remuestreo
      #calcular y retornar los intervalos a este caso; BCa, ABC y ponderado
      caso <- 2
      residualesRP <- residualesRobustos
    }
    
    #Tercer caso, cuando no hay homocedasticidad
    if(!hayHomocedasticidad){
      # se utiliza el estimador robusto MM  
      # los residuales robustos con ponderacion
      #implmentar los 6 esquemas de remuestreo
      #calcular y retornar los intervalos a este caso; BCa, ABC y ponderado
      caso <- 3
      residualesRP <- PoderarResidualesRobustos(residualesRobustos,CMERobusto)
      
    }
    
  }
  
  remuestrasBootR <- matrix(0,nrow = B, ncol = 6) #todas las remuestras
  
  for(i in 1:6){
    rems <- ImplementarRemuestreosBootsR(z,residuales,residualesRobustos,residualesRP,yAjRob,hii,B=100,tipo=i)
    remuestrasBootR[,i] <- rems
  }
  print(remuestrasBootR)
  
  resultados <- vector("list", 6) #todas las remuestras
  
  #Calculo de los intervalos
  if(caso == 1){
    print("hola en 1")
    
    for(i in 1:6){
      muestrasR2Boot <- remuestrasBootR[,i]
      #Percentil
      perc <- ImplementarIntervalosConfianza(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=1)
      #Bootstrap-t
      bootT <- ImplementarIntervalosConfianza(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=2)
      #Percentil-simetrico
      simt <- ImplementarIntervalosConfianza(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=7)
      
      resultados[[i]] <- list( perc, bootT,simt)
    }
  
    print(resultados)
  }
  
  if(caso==2 || caso==3){
    print(paste("hola en ", caso))
    for(i in 1:6){
      
      muestrasR2Boot <- remuestrasBootR[,i]
      #BCa
      bca <- ImplementarIntervalosConfianza(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=4)
      #ABC
      abc <- ImplementarIntervalosConfianza(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=5)
      #Ponderado ? Iterativo
      
      resultados[[i]] <- list( bca,abc)
      
    }
    print(resultados)
  }
}







# CalcularIntervaloConfianzaPercentil(originalR2,muestrasR2Bootstrap=muestrasBootstrapLiu2[,2])
# CalcularIntervaloConfianzaBootstrapT(originalR2,muestrasR2Bootstrap=muestrasBootstrapLiu2[,2])
# CalcularIntervaloConfianzaBootstrapBCa(y,z,originalR2,B,muestrasR2Bootstrap=muestrasBootstrapLiu2[,2])
# CalcularIntervaloConfianzaBootstrapABC(muestrasR2Bootstrap=muestrasBootstrapLiu2[,2])


