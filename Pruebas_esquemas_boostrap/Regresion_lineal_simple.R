#Carga de las funciones de pruebas de supuestos
source("Pruebas_esquemas_boostrap/Pruebas_supuestos_regresion_lineal.R")


data <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_prubas/N-VC.csv")


z = data$ValPreKgZ
y = data$ValObsKgY
modelo_lineal = lm(y~z)
e = residuals(modelo_lineal)

#comprobar la normalidad de los residuales
normalidad <- ComprobarSupuestoNormalidad(e)