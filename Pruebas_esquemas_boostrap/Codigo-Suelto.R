#Datos modelo de regresion lineal simple de los datos 
modeloLineal <- lm(y~z)
residuales <- residuals(modeloLineal)
coeficientes <- coef(modeloLineal)
yAjustados <- fitted(modeloLineal)
nResiduales <- length(residuales)
originalR2 <- summary(modeloLineal)$r.squared



for(i in 1:B){
  if(TipoBoot==1) {tt=rnorm(n)}
  if(TipoBoot==2){t1=(Res-mean(Res))/sd(Res)
  tt=sample(t1,replace=T)
  }
  if(TipoBoot==3){Mediana=median(ResPond)
  NMAD=(1/0.6745)*median(abs(ResPond-Mediana))
  t1=(ResPond-Mediana)/NMAD
  tt=sample(t1,replace=T)
  }
  
  if(TipoBoot==4) {tt=rgamma(n,4,2)}
  if(TipoBoot==5){Media1=0.5*sqrt(17/6)+sqrt(1/6)
  Media2=0.5*sqrt(17/6)-sqrt(1/6)
  H=rnorm(n,Media1,sqrt(0.5))
  D=rnorm(n,Media2,sqrt(0.5))
  tt=H*D-Media1*Media2
  }
  
  
  
  
  
  IntBootsBCaR2=function(b,x,y,B,NivConf)
  {
    alfa=1-NivConf
    if (b==1) M=residualesbalanceado(x,y,B)
    else M=pareadobalanceado(x,y,B)
    
  R2=summary(lm(y~x)) $r.squared
  p0=length(M[,3] [M[,3]>R2])/B
  suma0=0
  suma1=0
  suma2=0
  for (i in 1:length(y))
  {
    R2MI=summary(lm(y[-i]~x[-i])) $r.squared
    suma0=R2MI+suma0
  }
  
  R2PMI=suma0/length(y)
  
  for (i in 1:length(y)){
    DifR2=R2PMI-summary(lm(y[-i]~x[-i])) $r.squared
    suma1=DifR2^2+suma1
    suma2=DifR2^3+suma2
  }
  
  aR2=suma2/(6*(suma1^1.5))
  ZR2L=(qnorm(1-p0)-qnorm(1-alfa/2))/
    (1-aR2*(qnorm(1-p0)-qnorm(1-alfa/2)))+qnorm(1-p0)
  ZR2U=(qnorm(1-p0)+qnorm(1-alfa/2))/
    (1-aR2*(qnorm(1-p0)+qnorm(1-alfa/2)))+qnorm(1-p0)
  q00=pnorm(ZR2L)
  q01=pnorm(ZR2U)
  LIR2=quantile(M[,3],probs=q00)
  LSR2=quantile(M[,3],probs=q01)
  Tabla=matrix(nrow=1,ncol=2)
  dimnames(Tabla)=list(c("R2"),c("LI","LS"))
  Tabla[1,]=c(LIR2,LSR2)
  return(Tabla)
  }
  
  
  
  
  x <- c(41.28, 45.16, 34.75, 40.76, 43.61, 39.05, 41.20, 
         41.02, 41.33, 40.61, 40.49, 41.77, 42.07, 
         + 44.83, 29.12, 45.59, 41.95, 45.78, 42.89, 40.42, 49.31, 
         44.01, 34.87, 38.60, 39.63, 38.52, 38.52, 
         + 43.95, 49.08, 50.52, 43.85, 40.64, 45.86, 41.25, 50.35, 
         45.18, 39.67, 43.89, 43.89, 42.16) 

n <-length(x) #sample size 

mean(x) + qt(0.975, n - 1) * sd(x) * c( - 1, +1)/sqrt(n)
  

#Otra version
alfa<- 0.05
x_barra <- mean(x)
pto_crit <- quantile(x, c(alfa/2, 1 - alfa/2))
ic_inf_boot <- x_barra - pto_crit[2]/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1]/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot







data0 <- read.csv("C:/Users/irving/Downloads/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_0.csv")
z <- data0[[1]]
y <- data0[[2]]

B <- 100
constantePeso <- 3
#Proceso regresion lineal simple
modeloLineal <- lm(y~z)
residuales <- residuals(modeloLineal)
coeficientes <- coef(modeloLineal)
yAjustados <- fitted(modeloLineal)
nResiduales <- length(residuales)
originalR2 <- summary(modeloLineal)$r.squared
#CME <- summary(modeloLineal)$sigma**2
hii <- hatvalues(modeloLineal)

#Proceso regresión robusta
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









ImplementarRemuestreosBootsR <- function(z,y,B,residuales,tipoEsquema){
  
  Remuestra <- function(residualesRobustosPonderados,tipoEsquema){
    switch(tipoEsquema,
           'Wu1' <- CalcularMuestrasBootstrapWu1(residualesRobustosPonderados),
           'Wu2' <- CalcularMuestrasBootstrapWu2(residualesRobustosPonderados),
           'Wu3' <- CalcularMuestrasBootstrapWu3(residualesRobustosPonderados),
           'Liu1' <- CalcularMuestrasBootstrapLiu1(residualesRobustosPonderados),
           'Liu2' <- CalcularMuestrasBootstrapLiu2(residualesRobustosPonderados),
           'Wild' <- CalcularMuestrasBootstrapWild(residualesRobustos),
           
           stop()
    )
  }
  
  return(Remuestra(residualesRobustosPonderados,tipoEsquema))
  
}

CME <- summary(modeloLineal)$sigma**2


CalcularIntervaloConfianzaPercentil<- function(originalR2,muestrasR2Bootstrap,nivSignicancia=0.95){
  alfa <- 1-nivSignicancia
  vectorR2Bootstrap <- muestrasR2Bootstrap
  n <- length(vectorR2Bootstrap)
  puntosCriticos <- quantile(vectorR2Bootstrap, c(alfa/2, 1 - alfa/2)) # Aproximación bootstrap de los puntos críticos
  ICInfBootP <- originalR2 - puntosCriticos[2] / sqrt(n)
  ICSupBootP <- originalR2 - puntosCriticos[1] / sqrt(n)
  intervaloConfianzaPercentil <-as.vector(c(ICInfBootP, ICSupBootP))
  return(intervaloConfianzaPercentil)
}








resultados <- vector("list", 6) #todas las remuestras

#Calculo de los intervalos
if(NVC == caso){
  #percentil, boots-t e iterativo
  print("hola en 1")
  
  for(i in 1:6){
    muestrasR2Boot <- remuestrasBootR[,i]
    perc <- ContruirIntervBoot(data,R2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=1)
    bootT <- ContruirIntervBoot(data,R2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=2)
    simt <- ContruirIntervBoot(data,R2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=7)
    resultados[[i]] <- list( perc, bootT,simt)
  }
  
  print(resultados)
}else{
  print(paste("hola en ", caso))
  for(i in 1:6){
    
    muestrasR2Boot <- remuestrasBootR[,i]
    bca <- ContruirIntervBoot(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=4)
    abc <- ContruirIntervBoot(data,originalR2,B,muestrasR2Boot,nivConfianza,tipoIntervalo=5)
    resultados[[i]] <- list( bca,abc)
  }
  print(resultados)
}





#Funcion para aplicar esquema de bootstrap
AplicarEsqBoot <- function(z,residuales,residualesRob,residualesRP,yAjRob,hii,B,tipo){
  n <- length(residualesRP)
  muestrasBoot <- numeric(B)
  sqrt_hii <- sqrt(1 - hii)
  
  for (i in 1:B) {
    
    residualBT <- switch(
      tipo,
      {
        tt <- rnorm(n)
        (tt * residualesRP) / sqrt_hii
      },
      {
        ai <- (residuales-mean(residuales))/sd(residuales)
        tt <- sample(ai, replace = TRUE)
        (tt * residualesRP) / sqrt_hii
      },
      {
        mediana <- median(residualesRP)
        NMAD <- (1 / 0.6745) * median(abs(residualesRP - mediana))
        Rai <- (residualesRP - mediana) / NMAD
        tt <- sample(Rai, replace = TRUE)
        (tt * residualesRP) / sqrt_hii
      },
      {
        tt <- rgamma(n, 2, 4)
        (tt * residualesRP) / sqrt_hii
      },
      {
        media1 <- 0.5 * sqrt(17 / 6) + sqrt(1 / 6)
        media2 <- 0.5 * sqrt(17 / 6) - sqrt(1 / 6)
        H <- rnorm(n, media1, sqrt(0.5))
        D <- rnorm(n, media2, sqrt(0.5))
        tt <- H * D - media1 * media2
        (tt * residualesRP) / sqrt_hii
      },
      {
        sample(residualesRob, replace = TRUE)
      }
    )
    
    yBoots <- yAjRob + residualBT
    modeloBoots <- lm(yBoots ~ z)
    muestrasBoot[i] <- summary(modeloBoots)$r.squared
  }
  
  muestrasBoot
  
}





{#percentil-t fuente(book_remuestreo_RicardoCao-RubenFernandez_05sep2022.pdf/)
  # desvEst <- sd(vectorR2Bootstrap)
  # puntosCriticos <- quantile(vectorR2Bootstrap, c(alpha/2, 1 - alpha/2)) # Aproximación bootstrap de los puntos críticos
  # ICInfBootT <- media - puntosCriticos[2] * desvEst/raiz_n
  # ICSupBootT <- media - puntosCriticos[1] * desvEst/raiz_n
  # as.vector(c(ICInfBootT, ICSupBootT))
},{#percentil-t simetrizado fuente(book_remuestreo_RicardoCao-RubenFernandez_05sep2022.pdf/pag74)
  #Checar
  # desvEst <- sd(vectorR2Bootstrap)
  # puntosCriticos <- quantile(vectorR2Bootstrap, 1 - alpha)# Aproximación bootstrap de los puntos críticos
  # # Construcción del IC
  # ICInfBootS <- media - puntosCriticos * desvEst/raiz_n
  # ICSupBootS <- media + puntosCriticos * desvEst/raiz_n
  # as.vector(c(ICInfBootS, ICSupBootS))
}



for (i in 1:B){
  pos1=(i-1)*n+1
  pos2=i*n
  z=N[pos1:pos2]
  for (j in 1:length(x)){
    Bootsresiduos[i,j]=residuos[z[j]]
  }
    YBoots=Bootsresiduos[i,]+ajustados
    BetasBoots=coefficients(lm(YBoots~x))
    modelo=lm(YBoots~x)
    summary(modelo)
    R2=summary(modelo)$r.squared
    MuestraBoots[i,]=c(BetasBoots[1],BetasBoots[2],R2) 
}
return(MuestraBoots)}



carpetas <- sort(list.dirs(directorio, full.names = FALSE, recursive = FALSE))
archivos <- sort(list.files(directorio))







# Crear un data frame para almacenar los resultados
resultados_df <- data.frame(
  Esquema = integer(),
  Intervalo = integer(),
  R2_original = numeric(),
  Inferior = numeric(),
  Superior = numeric(),
  Longitud = numeric(),
  R2_en_intervalo = integer()
)

# Llenar el data frame con los intervalos y calcular sus longitudes
for (es in 1:length(resultadosInter)) {
  esquema <- resultadosInter[[es]]
  
  for (esj in 1:length(esquema)) {
    intervalo <- esquema[[esj]]
    longitud <- intervalo[2] - intervalo[1]
    R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
    resultados_df <- rbind(resultados_df, data.frame(
      Esquema = es,
      Intervalo = esj,
      R2_original = R2_modelo,
      Inferior = intervalo[1],
      Superior = intervalo[2],
      Longitud = longitud,
      R2_en_intervalo = R2_intervalo
    ))
  }
}








######################################3
#Función para implementar el calculo de las precision de las muestras
# Dado archivos_encontrados <-  list(muestra,R2) con rutas de archivos
# para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
# para replicas sea el número de replicas
# con nivConfinza entre 0 - 1
# y sea para N el tamaño de las muestras
CalPrecMuestras <- function(archivos_encontrados, caso, replicas, nivConfianza, N) {
  archivo_muestra <- archivos_encontrados$muestra
  archivo_R2 <- archivos_encontrados$R2
  library(readxl)
  
  # Leer los archivos .xlsx considerando que tienen cabecera
  data_muestra <- read_excel(archivo_muestra, col_names = TRUE)
  data_R2 <- read_excel(archivo_R2, col_names = TRUE)
  
  block_size <- 1000
  cols_per_model <- 2
  
  # Convertir a matriz para procesamiento
  data_matrix <- as.matrix(data_muestra)
  R2_matrix <- as.matrix(data_R2)
  limit_model <- 3
  # Ajustando por replicas
  for (replica in 1:replicas) {
    print(paste("Replica #", replica))
    # Reiniciar el contador de modelos procesados para cada réplica
    m <- 0
    
    # Extraer el bloque de datos para la réplica actual
    filaInicio <- (replica - 1) * N + 1
    filaFin <- replica * N
    replica_data <- data_matrix[filaInicio:filaFin, ]
    R2_replica <- R2_matrix[replica, ]
    
    # Procesar cada bloque de 1000 columnas
    for (i in seq(1, ncol(replica_data), by = block_size)) {
      block_end <- min(i + block_size - 1, ncol(replica_data))
      blockCaso <- replica_data[, i:block_end]  # Extraer el bloque actual de columnas
      R2_block <- R2_replica[i:block_end]  # Extraer el bloque actual de R2
      
      for (j in seq(1, ncol(blockCaso), by = cols_per_model)) {
        model_end <- min(j + cols_per_model - 1, ncol(blockCaso))
        modeloActual <- blockCaso[, j:model_end]  # Extraer las columnas para el modelo actual
        R2_modelo <- R2_block[ceiling(j / cols_per_model)]  # Extraer la R2 del modelo actual, redondeo
        
        # Aquí procesas los datos del modelo actual
        print(paste("Procesando modelo con columnas:", j, "-", model_end))
        resultadosInter <- CalPrecicion(modeloActual, alpha, nivConfianza, caso)
        print("procesado")
        print(R2_modelo)
        
        # Procesando resultados
        # Crear un data frame para almacenar los resultados
        resultados_df <- data.frame(
          Esquema = integer(),
          Intervalo = integer(),
          R2_original = numeric(),
          Inferior = numeric(),
          Superior = numeric(),
          Longitud = numeric(),
          R2_en_intervalo = integer()
        )
        
        # Llenar el data frame con los intervalos y calcular sus longitudes
        for (es in 1:length(resultadosInter)) {
          esquema <- resultadosInter[[es]]
          
          for (esj in 1:length(esquema)) {
            intervalo <- esquema[[esj]]
            longitud <- intervalo[2] - intervalo[1]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            resultados_df <- rbind(resultados_df, data.frame(
              Esquema = es,
              Intervalo = esj,
              R2_original = R2_modelo,
              Inferior = intervalo[1],
              Superior = intervalo[2],
              Longitud = longitud,
              R2_en_intervalo = R2_intervalo
            ))
          }
        }
        
        # Mostrar el data frame
        rownames(resultados_df) <- NULL  # Eliminar nombres automáticos de filas
        print(resultados_df)
        
        # Procesando ganadores
        resultados_ganadores <- data.frame(
          Esquema = integer(),
          IntervaloGanador = integer(),
          GanadorPor = integer(),
          TamanoIntervalo = numeric()
        )
        
        for (es in 1:length(resultadosInter)) {
          esquema <- resultadosInter[[es]]
          intervalos_ganadores <- data.frame(
            Intervalo = integer(),
            Inferior = numeric(),
            Superior = numeric(),
            Longitud = numeric(),
            R2_en_intervalo = integer()
          )
          
          for (esj in 1:length(esquema)) {
            intervalo <- esquema[[esj]]
            longitud <- intervalo[2] - intervalo[1]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            
            if (R2_intervalo == 1) {
              intervalos_ganadores <- rbind(intervalos_ganadores, data.frame(
                Intervalo = esj,
                Inferior = intervalo[1],
                Superior = intervalo[2],
                Longitud = longitud,
                R2_en_intervalo = R2_intervalo
              ))
            }
          }
          
          if (nrow(intervalos_ganadores) == 0) {
            # No hay ganadores, seleccionamos el mejor de los primeros tres intervalos
            if (length(esquema) >= 3) {
              mejores_intervalos <- esquema[1:3]
              longitudes <- sapply(mejores_intervalos, function(x) x[2] - x[1])
              mejor_intervalo_idx <- which.min(longitudes)
              intervalo_ganador <- mejor_intervalo_idx
              ganador_por <- 0
              tamano_intervalo <- longitudes[mejor_intervalo_idx]
            } else {
              mejores_intervalos <- esquema
              longitudes <- sapply(mejores_intervalos, function(x) x[2] - x[1])
              mejor_intervalo_idx <- which.min(longitudes)
              intervalo_ganador <- mejor_intervalo_idx
              ganador_por <- 0
              tamano_intervalo <- longitudes[mejor_intervalo_idx]
            }
          } else if (nrow(intervalos_ganadores) == 1) {
            intervalo_ganador <- intervalos_ganadores$Intervalo
            ganador_por <- 1
            tamano_intervalo <- intervalos_ganadores$Longitud
          } else if (nrow(intervalos_ganadores) == 2) {
            intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
            ganador_por <- 2
            tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
          } else if (nrow(intervalos_ganadores) == 3) {
            intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
            ganador_por <- 3
            tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
          } else {
            intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
            ganador_por <- 3
            tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
          }
          
          resultados_ganadores <- rbind(resultados_ganadores, data.frame(
            Esquema = es,
            IntervaloGanador = intervalo_ganador,
            GanadorPor = ganador_por,
            TamanoIntervalo = tamano_intervalo
          ))
        }
        
        # Mostrar el data frame
        rownames(resultados_ganadores) <- NULL  # Eliminar nombres automáticos de filas
        print(resultados_ganadores)
        
        m <- m + 1
        if (m == limit_model) {
          break
        }
      } # Fin procesando modelo
      if (m == limit_model) {
        break
      }
    } # Fin proceso bloques de 1000
  } # Fin replicas 
}



# Determinar los ganadores únicos, empates y triple empate
if (length(intervalos_ganadores) == 1) {
  if (intervalos_ganadores[[1]]$Intervalo == 1) FrecEficIB1Unico <- FrecEficIB1Unico + 1
  if (intervalos_ganadores[[1]]$Intervalo == 2) FrecEficIB2Unico <- FrecEficIB2Unico + 1
  if (intervalos_ganadores[[1]]$Intervalo == 3) FrecEficIB3Unico <- FrecEficIB3Unico + 1
} else if (length(intervalos_ganadores) == 2) {
  FrecEficIB1Emp2 <- FrecEficIB1Emp2 + (1 %in% sapply(intervalos_ganadores, function(x) x$Intervalo))
  FrecEficIB2Emp2 <- FrecEficIB2Emp2 + (2 %in% sapply(intervalos_ganadores, function(x) x$Intervalo))
  FrecEficIB3Emp2 <- FrecEficIB3Emp2 + (3 %in% sapply(intervalos_ganadores, function(x) x$Intervalo))
} else if (length(intervalos_ganadores) == 3) {
  FrecEficIB1Emp3 <- FrecEficIB1Emp3 + 1
  FrecEficIB2Emp3 <- FrecEficIB2Emp3 + 1
  FrecEficIB3Emp3 <- FrecEficIB3Emp3 + 1
}else{
  no_entro_ninguno <-no_entro_ninguno+1
}

#Ganador absoluto
ganador_absoluto <- resultados_ganadores[which.min(resultados_ganadores$TamanoIntervalo),]





####################################################3
#Bloque entero de busqueda del mejor por data.frame
# Procesando ganadores
resultados_ganadores <- data.frame(
  Esquema = integer(),
  IntervaloGanador = integer(),
  GanadorPor = integer(),
  TamanoIntervalo = numeric()
)

for (es in 1:length(resultadosInter)) {
  esquema <- resultadosInter[[es]]
  intervalos_ganadores <- data.frame(
    Intervalo = integer(),
    Inferior = numeric(),
    Superior = numeric(),
    Longitud = numeric(),
    R2_en_intervalo = integer()
  )
  
  for (esj in 1:length(esquema)) {
    intervalo <- esquema[[esj]]
    longitud <- intervalo[2] - intervalo[1]
    R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
    
    if (R2_intervalo == 1) {
      intervalos_ganadores <- rbind(intervalos_ganadores, data.frame(
        Intervalo = esj,
        Inferior = intervalo[1],
        Superior = intervalo[2],
        Longitud = longitud,
        R2_en_intervalo = R2_intervalo
      ))
    }
  }
  
  if (nrow(intervalos_ganadores) == 0) {
    intervalo_ganador <- 0
    ganador_por <- 0
    tamano_intervalo <- 0
  } else if (nrow(intervalos_ganadores) == 1) {
    intervalo_ganador <- intervalos_ganadores$Intervalo
    ganador_por <- 1
    tamano_intervalo <- intervalos_ganadores$Longitud
  } else if (nrow(intervalos_ganadores) == 2) {
    intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
    ganador_por <- 2
    tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
  } else if (nrow(intervalos_ganadores) == 3) {
    intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
    ganador_por <- 3
    tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
  } else {
    intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
    ganador_por <- 3
    tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
  }
  
  resultados_ganadores <- rbind(resultados_ganadores, data.frame(
    Esquema = es,
    IntervaloGanador = intervalo_ganador,
    GanadorPor = ganador_por,
    TamanoIntervalo = tamano_intervalo
  ))
  
  
}

# Mostrar el data frame
rownames(resultados_ganadores) <- NULL  # Eliminar nombres automáticos de filas
print(resultados_ganadores)


#Filtrar ganadores
resultados_ganadores_filtrados <- subset(resultados_ganadores, IntervaloGanador!=0)

#Seleccionar al mejor
if(nrow(resultados_ganadores_filtrados)>0){
  ganador_absoluto <-resultados_ganadores_filtrados[which.min(resultados_ganadores_filtrados$TamanoIntervalo),]
}else{
  anador_absoluto <- resultados_ganadores[1,0] #fallbac en caso de que sean 0
}

#datos
datos_ganador <- ganador_absoluto[,c("Esquema","IntervaloGanador","GanadorPor")]
print(datos_ganador)



# CalPrecMuestras <- function(archivos_encontrados, caso, replicas, nivConfianza, N) {
#   archivo_muestra <- archivos_encontrados$muestra
#   archivo_R2 <- archivos_encontrados$R2
#   library(readxl)
#   
#   # Leer los archivos .xlsx considerando que tienen cabecera
#   data_muestra <- read_excel(archivo_muestra, col_names = TRUE)
#   data_R2 <- read_excel(archivo_R2, col_names = TRUE)
#   
#   block_size <- 1000
#   cols_per_model <- 2
#   
#   # Convertir a matriz para procesamiento
#   data_matrix <- as.matrix(data_muestra)
#   R2_matrix <- as.matrix(data_R2)
#   limit_model <- 3
#   
#   # Ajustando por replicas
#   for (replica in 1:replicas) {
#     print(paste("Replica #", replica))
#     # Reiniciar el contador de modelos procesados para cada réplica
#     m <- 0
#     
#     #Datos de excel
#     
#     
#     # Extraer el bloque de datos para la réplica actual
#     filaInicio <- (replica - 1) * N + 1
#     filaFin <- replica * N
#     replica_data <- data_matrix[filaInicio:filaFin, ]
#     R2_replica <- R2_matrix[replica, ]
#     
#     # Procesar cada bloque de 1000 columnas
#     for (i in seq(1, ncol(replica_data), by = block_size)) {
#       block_end <- min(i + block_size - 1, ncol(replica_data))
#       blockCaso <- replica_data[, i:block_end]  # Extraer el bloque actual de columnas
#       R2_block <- R2_replica[i:block_end]  # Extraer el bloque actual de R2
#       
#       for (j in seq(1, ncol(blockCaso), by = cols_per_model)) {
#         model_end <- min(j + cols_per_model - 1, ncol(blockCaso))
#         modeloActual <- blockCaso[, j:model_end]  # Extraer las columnas para el modelo actual
#         R2_modelo <- R2_block[ceiling(j / cols_per_model)]  # Extraer la R2 del modelo actual, redondeo
#         
#          # Aquí procesas los datos del modelo actual
#          print(paste("Procesando modelo con columnas:", j, "-", model_end))
#          resultadosInter <- CalPrecicion(modeloActual, alpha, nivConfianza, caso)
#          print("procesado")
#          print(R2_modelo)
#          
#          #Procesando resultados
#          
#          # Crear un data frame para almacenar los resultados
#          resultados_df <- data.frame(
#            Esquema = integer(),
#            Intervalo = integer(),
#            R2_original = numeric(),
#            Inferior = numeric(),
#            Superior = numeric(),
#            Longitud = numeric(),
#            R2_en_intervalo = integer()
#          )
#          
#          # Llenar el data frame con los intervalos y calcular sus longitudes
#          for (es in 1:length(resultadosInter)) {
#            esquema <- resultadosInter[[es]]
#            
#            for (esj in 1:length(esquema)) {
#              intervalo <- esquema[[esj]]
#              longitud <- intervalo[2] - intervalo[1]
#              R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
#              resultados_df <- rbind(resultados_df, data.frame(
#                Esquema = es,
#                Intervalo = esj,
#                R2_original = R2_modelo,
#                Inferior = intervalo[1],
#                Superior = intervalo[2],
#                Longitud = longitud,
#                R2_en_intervalo = R2_intervalo
#              ))
#            }
#          }
#          
#          # Mostrar el data frame
#          rownames(resultados_df) <- NULL  # Eliminar nombres automáticos de filas
#          print(resultados_df)
#          
#          
#          
#          
#          # Procesando ganadores
#          resultados_ganadores <- data.frame(
#            Esquema = integer(),
#            IntervaloGanador = integer(),
#            GanadorPor = integer(),
#            TamanoIntervalo = numeric()
#          )
#          
#          for (es in 1:length(resultadosInter)) {
#            esquema <- resultadosInter[[es]]
#            intervalos_ganadores <- data.frame(
#              Intervalo = integer(),
#              Inferior = numeric(),
#              Superior = numeric(),
#              Longitud = numeric(),
#              R2_en_intervalo = integer()
#            )
#            
#            for (esj in 1:length(esquema)) {
#              intervalo <- esquema[[esj]]
#              longitud <- intervalo[2] - intervalo[1]
#              R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
#              
#              if (R2_intervalo == 1) {
#                intervalos_ganadores <- rbind(intervalos_ganadores, data.frame(
#                  Intervalo = esj,
#                  Inferior = intervalo[1],
#                  Superior = intervalo[2],
#                  Longitud = longitud,
#                  R2_en_intervalo = R2_intervalo
#                ))
#              }
#            }
#            
#            if (nrow(intervalos_ganadores) == 0) {
#              intervalo_ganador <- 0
#              ganador_por <- 0
#              tamano_intervalo <- 0
#            } else if (nrow(intervalos_ganadores) == 1) {
#              intervalo_ganador <- intervalos_ganadores$Intervalo
#              ganador_por <- 1
#              tamano_intervalo <- intervalos_ganadores$Longitud
#            } else if (nrow(intervalos_ganadores) == 2) {
#              intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
#              ganador_por <- 2
#              tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
#            } else if (nrow(intervalos_ganadores) == 3) {
#              intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
#              ganador_por <- 3
#              tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
#            } else {
#              intervalo_ganador <- intervalos_ganadores$Intervalo[which.min(intervalos_ganadores$Longitud)]
#              ganador_por <- 3
#              tamano_intervalo <- intervalos_ganadores$Longitud[which.min(intervalos_ganadores$Longitud)]
#            }
#            
#            resultados_ganadores <- rbind(resultados_ganadores, data.frame(
#              Esquema = es,
#              IntervaloGanador = intervalo_ganador,
#              GanadorPor = ganador_por,
#              TamanoIntervalo = tamano_intervalo
#            ))
#            
#            
#          }
#          
#          # Mostrar el data frame
#          rownames(resultados_ganadores) <- NULL  # Eliminar nombres automáticos de filas
#          print(resultados_ganadores)
#          
#          
#          #Filtrar ganadores
#          resultados_ganadores_filtrados <- subset(resultados_ganadores, IntervaloGanador!=0)
#          
#          #Seleccionar al mejor
#          if(nrow(resultados_ganadores_filtrados)>0){
#            ganador_absoluto <-resultados_ganadores_filtrados[which.min(resultados_ganadores_filtrados$TamanoIntervalo),]
#          }else{
#            ganador_absoluto <- resultados_ganadores[1,0] #fallbac en caso de que sean 0
#          }
#          
#          #datos
#          datos_ganador <- ganador_absoluto[,c("Esquema","IntervaloGanador","GanadorPor")]
#          print(datos_ganador)
#          
#          
#         m <- m + 1
#         if (m == limit_model) {
#           break
#         }
#       } # Fin procesando modelo
#       if (m == limit_model) {
#         break
#       }
#     } # Fin proceso bloques de 1000
#   } # Fin replicas 
# }




titulo <- paste0("Efica", N,".csv")
ruta<-paste("/home/irving/Documentos/tesis/Evaluacion-Presicion-Bootstrap/Resultados", titulo, sep = "/")
write.csv(conteos_totales, ruta, row.names = FALSE)


# Filtrar ganadores que no tienen intervalo ganador igual a 0
resultados_ganadores_filtrados <- subset(resultados_ganadores, IntervaloGanador != 0)

# Seleccionar el mejor
if (nrow(resultados_ganadores_filtrados) > 0) {
  ganador_absoluto <- resultados_ganadores_filtrados[which.min(resultados_ganadores_filtrados$TamanoIntervalo),]
} else {
  ganador_absoluto <- data.frame(Esquema = 0, IntervaloGanador = 0, GanadorPor = 0)# Fallback en caso de que sean 0
}

# Datos del ganador
datos_ganador <- ganador_absoluto[, c("Esquema", "IntervaloGanador", "GanadorPor")]
print(datos_ganador)

# Actualizar los conteos de la réplica actual
conteo_replica$modelosEvaluados <- conteo_replica$modelosEvaluados + 1

#Acumular el esquema ganador
if (datos_ganador$Esquema >= 1 && datos_ganador$Esquema <= 8) {
  conteo_replica[[paste0("esquema", datos_ganador$Esquema, "Ganador")]] <- conteo_replica[[paste0("esquema", datos_ganador$Esquema, "Ganador")]] + 1
} else {
  conteo_replica$sinesquemaGanador <- conteo_replica$sinesquemaGanador + 1
}

#Acumular el intervalo ganador
if (datos_ganador$IntervaloGanador >= 1 && datos_ganador$IntervaloGanador <= 8) {
  conteo_replica[[paste0("intervalo", datos_ganador$IntervaloGanador, "Ganador")]] <- conteo_replica[[paste0("intervalo", datos_ganador$IntervaloGanador, "Ganador")]] + 1
} else {
  conteo_replica$sinintervaloGanador <- conteo_replica$sinintervaloGanador + 1
}

#Acumular el tipo de ganador
if (datos_ganador$GanadorPor == 1) {
  conteo_replica$ganadorPorUnico <- conteo_replica$ganadorPorUnico + 1
} else if (datos_ganador$GanadorPor == 2) {
  conteo_replica$ganadorEmpate <- conteo_replica$ganadorEmpate + 1
} else if (datos_ganador$GanadorPor == 3) {
  conteo_replica$ganadorEmpateTriple <- conteo_replica$ganadorEmpateTriple + 1
}else{
  conteo_replica$ningunGanador<- conteo_replica$ningunGanador + 1
}






#BCa and ABC
set.seed(1) #for reproducibility
require('boot' )
x<-rnorm(15,mean = 20,sd = 2)
fboot <- function(x, i) mean(x[i]) #compute estimate given

bs <- boot(x, fboot, R = 1000) #generate bootstrap estimates
boot.ci(bs, type = 'bca' , conf = 0.95) #BCa method 95% C.I.

#Calculations and Intervals on Original Scale
fabc <- function(x, w) w%*%x #ABC uses weighted average
abc.ci(x, fabc, conf = 0.95)



#x<-rnorm(10)
theta<- function(x,w){sum(w*x)/sum(w)}
abc.ci(x,theta,conf=0.95)







CalPrecMuestras <- function(archivos_encontrados, caso, replicas, nivConfianza, N) {
  archivo_muestra <- archivos_encontrados$muestra
  archivo_R2 <- archivos_encontrados$R2
  library(readxl)
  
  # Leer los archivos .xlsx considerando que tienen cabecera
  data_muestra <- read_excel(archivo_muestra, col_names = TRUE)
  data_R2 <- read_excel(archivo_R2, col_names = TRUE)
  
  block <- 1000
  cols_por_model <- 2
  limit_model <- 3
  esquemas <- 8
  
  # Resultados del análisis
  nombre_cols <- c("Replica", "esquema", "FrecEficIB1", "FrecEficIB2", "FrecEficIB3",
                   "FrecEficIB1Unico", "FrecEficIB2Unico", "FrecEficIB3Unico",
                   "FrecEficIB1Emp2", "FrecEficIB2Emp2", "FrecEficIB3Emp2",
                   "FrecEficIB1Emp3", "FrecEficIB2Emp3", "FrecEficIB3Emp3")
  
  conteos_totales <- data.frame(matrix(0, ncol = length(nombre_cols), nrow = replicas * esquemas))
  colnames(conteos_totales) <- nombre_cols
  
  no_entro_ninguno <- 0
  
  # Procesando las réplicas
  for (replica in 1:replicas) {
    print(paste("Replica #", replica))
    m <- 0
    
    # Vectores que contienen los datos por defecto
    conteo_replica <- data.frame(matrix(0, nrow = 8, ncol = 14))
    replica_vector <- rep(replica, each = esquemas)
    esquema_vector <- rep(1:esquemas)
    matriz_inicial <- data.frame(replica_vector, esquema_vector)
    conteo_replica[, 1:2] <- matriz_inicial
    
    # Extraer el bloque de datos para la réplica actual
    fila_inicio <- (replica - 1) * N + 1
    fila_fin <- replica * N
    replica_data <- data_muestra[fila_inicio:fila_fin, ]
    R2_replica <- data_R2[replica, ]
    
    # Procesar cada bloque de 1000 columnas
    for (i in seq(1, ncol(replica_data), by = block)) {
      block_end <- min(i + block - 1, ncol(replica_data))
      block_caso <- replica_data[, i:block_end]
      R2_block <- R2_replica[i:block_end]
      
      for (j in seq(1, ncol(block_caso), by = cols_por_model)) {
        model_end <- min(j + cols_por_model - 1, ncol(block_caso))
        modeloActual <- block_caso[, j:model_end]
        R2_modelo <- R2_block[ceiling(j / cols_por_model)]
        
        # Aquí procesas los datos del modelo actual
        print(paste("Procesando modelo con columnas:", j, "-", model_end))
        resultadosInter <- CalPrecicion(modeloActual, alpha, nivConfianza, caso)
        print("procesado")
        print(R2_modelo)
        
        # Procesando resultados
        resultados_df <- data.frame(
          Esquema = integer(),
          Intervalo = integer(),
          R2_original = numeric(),
          Inferior = numeric(),
          Superior = numeric(),
          Longitud = numeric(),
          R2_en_intervalo = integer()
        )
        
        # Llenar el data frame con los intervalos y calcular sus longitudes
        for (es in 1:length(resultadosInter)) {
          esquema <- resultadosInter[[es]]
          
          for (esj in 1:length(esquema)) {
            intervalo <- esquema[[esj]]
            longitud <- intervalo[2] - intervalo[1]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            resultados_df <- rbind(resultados_df, data.frame(
              Esquema = es,
              Intervalo = esj,
              R2_original = R2_modelo,
              Inferior = intervalo[1],
              Superior = intervalo[2],
              Longitud = longitud,
              R2_en_intervalo = R2_intervalo
            ))
          }
        }
        rownames(resultados_df) <- NULL
        print(resultados_df)
        
        # Procesando ganadores
        resultados_ganadores <- data.frame(
          Esquema = integer(),
          IntervaloGanador = integer(),
          GanadorPor = integer(),
          TamanoIntervalo = numeric()
        )
        
        for (es in 1:length(resultadosInter)) {
          esquema <- resultadosInter[[es]]
          intervalos_ganadores <- data.frame(
            Intervalo = integer(),
            Inferior = numeric(),
            Superior = numeric(),
            Longitud = numeric()
          )
          
          # Procesando los intervalos por esquema
          for (esj in 1:length(esquema)) {
            intervalo <- esquema[[esj]]
            longitud <- intervalo[2] - intervalo[1]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            
            if (R2_intervalo == 1) {
              intervalos_ganadores <- rbind(intervalos_ganadores, data.frame(
                Intervalo = esj,
                Inferior = intervalo[1],
                Superior = intervalo[2],
                Longitud = longitud
              ))
            }
          }
          
          # Obtener al mejor intervalo por esquema y como ganó
          if (nrow(intervalos_ganadores) == 0) {
            intervalo_ganador <- 0
            ganador_por <- 0
            tamano_intervalo <- 0
          } else {
            # Determinar el mejor intervalo en caso de empate
            mejor_intervalo <- intervalos_ganadores[which.min(intervalos_ganadores$Longitud), ]
            intervalo_ganador <- mejor_intervalo$Intervalo
            ganador_por <- nrow(intervalos_ganadores)
            tamano_intervalo <- mejor_intervalo$Longitud
          }
          
          resultados_ganadores <- rbind(resultados_ganadores, data.frame(
            Esquema = es,
            IntervaloGanador = intervalo_ganador,
            GanadorPor = ganador_por,
            TamanoIntervalo = tamano_intervalo
          ))
        }
        rownames(resultados_ganadores) <- NULL
        print(resultados_ganadores)
        
        # Procesando los intervalos por esquema
        for (numEsquema in 1:length(resultadosInter)) {
          resultados_esquema <- resultadosInter[[numEsquema]]
          
          intervalos_ganadores <- list()
          for (numIntervalo in 1:length(resultados_esquema)) {
            intervalo <- resultados_esquema[[numIntervalo]]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            
            if (R2_intervalo == 1) {
              if (numIntervalo == 1) conteo_replica[numEsquema, 3] <- conteo_replica[numEsquema, 3] + 1
              if (numIntervalo == 2) conteo_replica[numEsquema, 4] <- conteo_replica[numEsquema, 4] + 1
              if (numIntervalo == 3) conteo_replica[numEsquema, 5] <- conteo_replica[numEsquema, 5] + 1
              
              intervalos_ganadores[[length(intervalos_ganadores) + 1]] <- list(Intervalo = numIntervalo, Longitud = intervalo[2] - intervalo[1])
            }
          }
          
          # Determinar los ganadores únicos, empates y triple empate
          if (length(intervalos_ganadores) == 1) {
            if (intervalos_ganadores[[1]]$Intervalo == 1) conteo_replica[numEsquema, 6] <- conteo_replica[numEsquema, 6] + 1
            if (intervalos_ganadores[[1]]$Intervalo == 2) conteo_replica[numEsquema, 7] <- conteo_replica[numEsquema, 7] + 1
            if (intervalos_ganadores[[1]]$Intervalo == 3) conteo_replica[numEsquema, 8] <- conteo_replica[numEsquema, 8] + 1
          } else if (length(intervalos_ganadores) == 2) {
            mejor_intervalo <- intervalos_ganadores[[which.min(sapply(intervalos_ganadores, function(x) x$Longitud))]]
            if (mejor_intervalo$Intervalo == 1) conteo_replica[numEsquema, 9] <- conteo_replica[numEsquema, 9] + 1
            if (mejor_intervalo$Intervalo == 2) conteo_replica[numEsquema, 10] <- conteo_replica[numEsquema, 10] + 1
            if (mejor_intervalo$Intervalo == 3) conteo_replica[numEsquema, 11] <- conteo_replica[numEsquema, 11] + 1
          } else if (length(intervalos_ganadores) == 3) {
            mejor_intervalo <- intervalos_ganadores[[which.min(sapply(intervalos_ganadores, function(x) x$Longitud))]]
            if (mejor_intervalo$Intervalo == 1) conteo_replica[numEsquema, 12] <- conteo_replica[numEsquema, 12] + 1
            if (mejor_intervalo$Intervalo == 2) conteo_replica[numEsquema, 13] <- conteo_replica[numEsquema, 13] + 1
            if (mejor_intervalo$Intervalo == 3) conteo_replica[numEsquema, 14] <- conteo_replica[numEsquema, 14] + 1
          } else {
            no_entro_ninguno <- no_entro_ninguno + 1
          }
        }
        
        m <- m + 1
        if (m == limit_model) {
          break
        }
      } # Fin procesando modelo
      if (m == limit_model) {
        break
      }
    } # Fin proceso bloques de 1000
    
    # Agregar conteo_replica a conteos_totales
    fila_inicio <- (replica - 1) * esquemas + 1
    fila_fin <- replica * esquemas
    conteos_totales[fila_inicio:fila_fin, ] <- conteo_replica
  } # Fin replicas
  print(conteos_totales)
}
