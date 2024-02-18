#cargando librerias
library("robustbase")

#Carga de las funciones de pruebas de supuestos
source("Pruebas_esquemas_boostrap/Pruebas_supuestos_regresion_lineal.R")
source("Pruebas_esquemas_boostrap/Esquemas_Bootstrap.R")
source("Pruebas_esquemas_boostrap/Intervalos_confianza_Bootstrap_N-VC.R")

data0 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_0.csv")
data1 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_1.csv")
data2 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_2.csv")
data3 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_3.csv")
data4 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/N-VC.csv")

#Casos de analisis
caso0 <- 0
caso1 <- 1
caso2 <- 2
caso3 <- 3
casoRealizar <- 0


#Carga de los datos
y <- data0[[1]]
z <- data0[[2]]

B <- 100 #numero de muestras Booststrap
nivSignicancia <- 0.95
alpha <- 1-nivSignicancia
constantePeso <- 3

#Proceso regresion lineal simple
modeloLineal <- lm(y~z)
residuales <- residuals(modeloLineal)
coeficientes <- coef(modeloLineal)
yAjustados <- fitted(modeloLineal)
nResiduales <- length(residuales)
originalR2 <- summary(modeloLineal)$r.squared
CME <- summary(modeloLineal)$sigma**2
hii <- hatvalues(modeloLineal)

#Proceso regresión robusta, paso 1
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

#Muestras boostrap
muestrasBootstrapWu1 <- CalcularMuestrasBootstrapWu1(z,B=100,nResidualesRobustos,yAjustadosRobustos,residualesRobustosPonderados,hii)
muestrasBootstrapWu2 <-CalcularMuestrasBootstrapWu2(z,B=100,yAjustadosRobustos,residualesRobustosPonderados,hii,residuales)
muestrasBootstrapWu3 <-CalcularMuestrasBootstrapWu3(z,B=100,yAjustadosRobustos,residualesRobustosPonderados,hii)

muestrasBootstrapLiu1 <-CalcularMuestrasBootstrapLiu1(z,B=100,yAjustadosRobustos,nResidualesRobustos,residualesRobustosPonderados,hii)
muestrasBootstrapLiu2 <-CalcularMuestrasBootstrapLiu2(z,B=100,yAjustadosRobustos,nResidualesRobustos,residualesRobustosPonderados,hii)

muestrasBootstrapWild <-CalcularMuestrasBootstrapWild(z,B=100,yAjustadosRobustos,nResidualesRobustos,residualesRobustos)




#comprobar la normalidad y homoestecidad de los residuales
hayNormalidad <- ComprobarSupuestoNormalidad(residuales)
hayHomoestacidad <- ComprobarVarianzaConstante(z,modeloLineal)




