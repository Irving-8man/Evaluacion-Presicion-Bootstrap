#Tesis, evaluación de la precisión de un modelo con la 
#tecnica de Regresion lineal con Bootstrap y estimadores robustos
#Irving Geyler Cupul Uc

#Librerias
library("robustbase")
library("readxl")
#Función de apoyo para ponderar residuales
#Sea residualesRob <- los residuales del modelo con regresión robusta
#CMERob <- cuadrado medio del error del modleo con regresión robusta
PodResidRobu <- function(residualesRob, CMERob){
  consPes <- 3
  n <- length(residualesRob)
  x <- abs(residualesRob)/sqrt(CMERob)
  w <- rep(1,n)
  xx <- which(x > consPes) 
  w[xx] <- (consPes / w[xx])
  resRobuPonder<- w*residualesRob
  return(resRobuPonder)
}

#Función de apoyo para obtener las muestras Bootstrap de R²
#Sea y <- los observador, z <- los estimados,yAjRob <- ajustados 
#residuales <- obtenenidos de la regresion lineal,residualesRob <- residuales robustos de la regresion lineal, 
#residualesRP <- residuales robustos ponderados,hii<- valor de aplacamiento del modelo con 
#regresion simple, n <- tamaño de la muestra de residuos, 
#B <- repeticiones Bootstrap,tipo <- sea el tipo de esquema Boostrap Robusto. 
#1<-Wu 1, 2<-Wu 2, 3<-Wu 3, 4<-Liu 1, 5<-Liu 2, 6<-Wild
# 7<-residuales balanceados, 8<-residuales pareados
CalcularR2Bootstrap <- function(y, z, yAjRob, residuales, residualesRob, residualesRP, hii, n, B, tipo) {
  sqrt_hii <- sqrt(1 - hii)
  RsBoot <- numeric(B)
  
  if((tipo == 7) || (tipo ==8)){
    N     <- rep(1:n,B)
    NPerm <- sample(N)
  }
  
  for (i in 1:B) {
    residualBT <- switch(
      tipo,
      {
        tt <- rnorm(n)
        (tt * residualesRP) / sqrt_hii
      },
      {
        ai <- (residuales - mean(residuales)) / sd(residuales)
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
      },
      {
        posI=(i-1)*n+1
        posF=i*n
        VPosi=NPerm[posI:posF]
        residualesRP[VPosi]
      },
      {
        posI=(i-1)*n+1
        posF=i*n
        NPerm[posI:posF]
      },
      stop("Esquema no válido")
    )
    
    if(tipo != 8){
      yBoots <- yAjRob + residualBT
      modeloBoots <- lm(yBoots ~ z)
    }else{
      YP <- y[residualBT]
      ZP <- z[residualBT]
      modeloBoots <- lm(YP ~ ZP)
    }
    
    RsBoot[i] <- summary(modeloBoots)$r.squared
  }
  
  return(RsBoot)
}


#Función de apoyo para construir el intervalo de confianza  la muestra de R² Bootstraps del modelo
#Sea data <- los valores z e y del modelo, R2 <- la R² estimada del modelo, 
#muestrasR2Boot <- la muestra de R² Bootstrap, B <- repeticiones Bootstrap, nivConfianza<- nivel de confianza,
#tipo<- sea el tipo de intervalo de confiaza que desea construir,opciones:
#1<- percentil, 2<-BCa
ContruirIntervBoot <- function(data, R2, muestrasR2Boot, B=100, nivConfianza=0.95, tipo){
  intervalo <- numeric(2) 
  alpha <- 1-nivConfianza
  z <- as.numeric(data[[1]])
  y <- as.numeric(data[[2]])
  vectorR2Bootstrap <- muestrasR2Boot
  
  intervalo <- switch(
    tipo,
    {
      puntosCriticos <- quantile(vectorR2Bootstrap, c(alpha/2, 1 - alpha/2))
      as.vector(puntosCriticos)
    },
    {
      n <- length(z)
      z0 <- qnorm(mean(vectorR2Bootstrap < R2))
      suma0=0
      suma02=0
      suma03=0
      for(i in 1:n)
      {R2MI=summary(lm(y[-i]~z[-i]))$r.squared
      suma0=R2MI+suma0
      }
      R2PMI=suma0/n
      for(i in 1:n)
      {Dif0=R2PMI-summary(lm(y[-i]~z[-i]))$r.squared
      suma02=Dif0^2+suma02
      suma03=Dif0^3+suma03
      }
      a=suma03/(6*(suma02^1.5))
      
      z_alfa1 <- qnorm(alpha / 2)
      z_alfa2 <- qnorm(1 - alpha / 2)
      alfa1 <- pnorm(z0 + (z0 + z_alfa1) / (1 - a * (z0 + z_alfa1)))
      alfa2 <- pnorm(z0 + (z0 + z_alfa2) / (1 - a * (z0 + z_alfa2)))
      
      ICInfBootBCa <- quantile(muestrasR2Boot, alfa1)
      ICSupBootBCa <- quantile(muestrasR2Boot, alfa2)
      as.vector(c(ICInfBootBCa, ICSupBootBCa))
    },
    stop("Intervalo no válido")
  )
  
  return (intervalo)
}



#Función propuesta para evaluar la precisión de un modelo creando
#intervalos de confianza con la tecnica de regresion lineal
#con estimadores robustos y esquemas Booststrap
#Sea data<- el modelo con las columnas z e y
#caso <- 1 #Normalidad- homocedasticidad, 2 #No normalidad-homocedasticidad,
#3 #Normalidad-heterocidasticidad y 4 #No normalidad-heterocidastecidad
EvalPrecisionModel <- function (data,alpha,nivConfianza,caso){
  z <- as.numeric(data[[1]])
  y <- as.numeric(data[[2]])
  n <- length(z)
  B <- 100 #Remuestras bootstrap
  
  #Regresion lineal simple
  modeloLineal <- lm(y~z)
  residuales <- residuals(modeloLineal)
  R2 <- summary(modeloLineal)$r.squared
  hii <- hatvalues(modeloLineal)
  yAju <- fitted(modeloLineal)
  
  #Casos 
  NVC <- 1 #Normalidad- homocedasticidad
  NNVC <- 2 #No normalidad-homocedasticidad
  NNC <- 3 #Normalidad-heterocidasticidad
  NNNC <- 4 #No normalidad-heterocidastecidad 
  numRemues <- 8 #Numero de tipos de remuestreos implementados
  valoresBootR <- matrix(0,nrow = B, ncol = numRemues)
  
  #Uso de los estimadores dependendiendo el caso
  if(NVC == caso){
    #Minimos cuadrados
    residualesRob <- residuales
    yAjRob <- yAju
    residualesRP <- residuales
    
  }else{
    # Estimador robusto MM  
    modeloLinealRob <- lmrob(y ~ z, method = "MM")
    residualesRob <- modeloLinealRob$residuals
    yAjRob <- modeloLinealRob$fitted.values
    
    #Segundo caso, no cumple normalidad
    if( NNVC==caso){
      # Residuales robustos sin ponderacion
      residualesRP <- residualesRob
    }
    
    #Tercer caso, no hay homocedasticidad
    if((NNNC==caso) || (NNC==caso)){
      # Residuales robustos con ponderacion
      CMERob <- modeloLinealRob$scale**2
      residualesRP <- PodResidRobu (residualesRob,CMERob)
    }
  }
  
  for(i in 1:numRemues){
    BootR <- CalcularR2Bootstrap(y,z, yAjRob,residuales,residualesRob,residualesRP,hii,n,B,tipo=i)
    valoresBootR[, i] <- BootR
  }
  resultadosInter <- vector("list", numRemues)
  
  #Calculo de los intervalos
  for(i in 1:numRemues){
    muestrasR2Boot <- valoresBootR[,i]
    perc <- ContruirIntervBoot(data,R2,muestrasR2Boot,B,nivConfianza=0.95,tipo=1)
    bca <- ContruirIntervBoot(data,R2,muestrasR2Boot,B,nivConfianza=0.95,tipo=2)
    resultadosInter[[i]] <- list(perc, bca)
  }
  
  return(resultadosInter)
}


##################################################################################################
#Procesamiento de datos y cálculo de datos de interes

#Función principal para el procesamiento del archivo y/o archivos
#con los modelos y sus replicas dado los parametros,
#modelo<-1-Precisos, 2-Imprecisos
#sea para el caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD,
#e indicando tamaño de muestra de los modelos.
ProcesarModelsData <- function(modelo,caso,numMuestras = c(10,15,20,25,30,35) ){
  #Modificar dependiendo de la carpeta donde este alojada
  directorio <- "../Evaluacion-Presicion-Bootstrap/Modelos Simulados"

  carp_model <- switch(modelo,
                       "Precisos" = "Precisos",
                       "Imprecisos" = "Imprecisos",
                       stop("Modelo no válido"))
  arch_model <- switch(modelo,
                       "Precisos" = "EP",
                       "Imprecisos" = "EI",
                       stop("Modelo no válido"))
  ruta_model <- file.path(directorio, carp_model)
  
  arch_caso <- c("NVC","NVD","NNVC","NNVD")
  tipo_caso <- switch(caso,{arch_caso[1]},{arch_caso[2]},{arch_caso[3]},{arch_caso[4]},stop())
  model_caso <- paste(arch_model, tipo_caso, sep = "")
  ruta_model_caso <- file.path(ruta_model, model_caso)
  
  # Recuperando los archivos de la ruta
  archivos <- sort(list.files(ruta_model_caso))
  
  # Encontrar los archivos de muestra y R2 correspondientes
  encontrar_archivos <- function(num) {
    patron_muestra <- paste0("^", model_caso, " ", num)
    patron_R2 <- paste0("^R2", arch_model, tipo_caso, " ", num)
    archivo_muestra <- archivos[grepl(patron_muestra, archivos)]
    archivo_R2 <- archivos[grepl(patron_R2, archivos)]
    if (length(archivo_muestra) == 1 && length(archivo_R2) == 1) {
      list(
        muestra = paste(ruta_model_caso, archivo_muestra, sep = "/"),
        R2 =paste(ruta_model_caso, archivo_R2, sep = "/"))
    } else {
      NULL
    }
  }
  
  replicas <- 5 #replicas de estudio
  #Procesamiento de las muestras solicitadas
  for(muestra in numMuestras){
    archivos_encontrados <- encontrar_archivos(muestra)
    if (!is.null(archivos_encontrados)) {
      ProcesarModels(archivos_encontrados,caso, replicas, nivConfianza = 0.95, N = muestra,MODELO=arch_model, CASO=tipo_caso)
    } else {
      cat("No se encontró un archivo de muestra o R2 correspondiente para el tamaño de muestra:", muestra, "\n")
    }
    
  }
  
}




# Función para procesar y capturar los resultados sobre la precisión de las muestras con sus I.C.
# Dado archivos_encontrados <-  list(muestra,R2) con rutas de archivos
# para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
# para replicas sea el número de replicas
# con nivConfinza entre 0 - 1
# y sea para N el tamaño de las muestras
ProcesarModels <- function(archivos_encontrados, caso, replicas, nivConfianza, N, MODELO, CASO) {
  archivo_muestra <- archivos_encontrados$muestra
  archivo_R2 <- archivos_encontrados$R2
  library(readxl)
  data_muestra <- read_excel(archivo_muestra, col_names = TRUE)
  data_R2 <- read_excel(archivo_R2, col_names = TRUE)
  
  block <- 1000
  cols_por_model <- 2
  limit_model <- 3 #modelos a procesar
  esquemas <- 8
  
  # Resultados del analisis
  #Matriz de conteo de resultados
  nombre_cols <- c("Replica","esquema","FrecEficIB1","FrecEficIB2",
                   "FrecEficIB1Unico","FrecEficIB2Unico",
                   "FrecEficIB1Emp2", "FrecEficIB2Emp2")
  conteos_totales <- matrix(0, ncol = length(nombre_cols), nrow = replicas * esquemas)
  colnames(conteos_totales) <- nombre_cols
  
  #Matriz de conteo de ceros en las replicas por esquemas
  nombre_cols_cer <- c("Replicas", "NumMod", "Esq1", "Esq2", "Esq3", "Esq4", "Esq5", "Esq6", "Esq7", "Esq8")
  conteo_ceros <- matrix(ncol=length(nombre_cols_cer),nrow = replicas)
  colnames(conteo_ceros) <- nombre_cols_cer
  no_entro_ninguno <-0
  
  #Nombre archivos
  nombre_archivo_conteos <- paste("resultados_conteos__", MODELO, "__", CASO, "__N", N, ".csv", sep = "")
  nombre_archivo_ceros <- paste("resultados_ceros__", MODELO, "__", CASO, "__N", N, ".csv", sep = "")

  # Procesando las replicas
  for (replica in 1:replicas) {
    print(paste("Replica #", replica))
    m <- 0
    
    # Vectores que contienen los datos por defecto
    conteo_replica <- matrix(0, nrow = 8, ncol = length(nombre_cols))
    replica_vector <- rep(replica, each = esquemas)
    esquema_vector <- rep(1:esquemas)
    matriz_inicial <- cbind(replica_vector, esquema_vector)
    conteo_replica[, 1:2] <- matriz_inicial
    # Vector para almacenar el conteo de ceros de la réplica actual
    conteo_ceros_replica <- numeric(length(nombre_cols_cer))
    conteo_ceros_replica[1] <- replica
    
    # Extraer el bloque de datos para la réplica actual
    fila_inicio <- (replica - 1) * N + 1
    fila_fin <- replica * N
    replica_data <- data_muestra[fila_inicio:fila_fin, ]
    R2_replica <- data_R2[replica, ]
    
    # Procesar cada bloque de 1000 columnas
    for (i in seq(1, ncol(replica_data), by = block)) {
      block_end <- min(i + block - 1, ncol(replica_data))
      block_caso <- replica_data[, i:block_end]
      R2_block <- R2_replica[i:500]
      Rmod <- 0
      
      for (j in seq(1, ncol(block_caso), by = cols_por_model)) {
        model_end <- min(j + cols_por_model - 1, ncol(block_caso))
        modeloActual <- block_caso[, j:model_end]
        R2_modelo <- R2_block[Rmod+1]
        Rmod <- Rmod + 1
        resultadosInter <- EvalPrecisionModel(modeloActual, alpha, nivConfianza, caso)
        
        # Procesando los intervalos por esquema
        for (numEsquema in 1:length(resultadosInter)) {
          resultados_esquema <- resultadosInter[[numEsquema]]
          
          intervalos_ganadores <- list()
          for (numIntervalo in 1:length(resultados_esquema)) {
            intervalo <- resultados_esquema[[numIntervalo]]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            
            if (R2_intervalo == 1) {
              if (numIntervalo == 1) conteo_replica[numEsquema,3] <- conteo_replica[numEsquema,3] + 1
              if (numIntervalo == 2) conteo_replica[numEsquema,4] <- conteo_replica[numEsquema,4] + 1
              
              intervalos_ganadores[[length(intervalos_ganadores) + 1]] <- list(Intervalo = numIntervalo, Longitud = intervalo[2] - intervalo[1])
            }
          }
          
          # Determinar los ganadores
          if (length(intervalos_ganadores) == 1) {
            if (intervalos_ganadores[[1]]$Intervalo == 1) conteo_replica[numEsquema,5] <- conteo_replica[numEsquema,5] + 1
            if (intervalos_ganadores[[1]]$Intervalo == 2) conteo_replica[numEsquema,6] <- conteo_replica[numEsquema,6] + 1
          } else if (length(intervalos_ganadores) == 2) {
            mejor_intervalo <- intervalos_ganadores[[which.min(sapply(intervalos_ganadores, function(x) x$Longitud))]]
            if (mejor_intervalo$Intervalo == 1) conteo_replica[numEsquema, 7] <- conteo_replica[numEsquema, 7] + 1
            if (mejor_intervalo$Intervalo == 2) conteo_replica[numEsquema, 8] <- conteo_replica[numEsquema, 8] + 1
          } else {
            no_entro_ninguno <- no_entro_ninguno + 1
            conteo_ceros_replica[numEsquema + 2] <- conteo_ceros_replica[numEsquema + 2] + 1
          }
        }#Fin de conteo
        
        m <- m + 1
        if (m %% limit_model == 0) {
          print(paste("Procesado", m, "modelos de replicas:", replica))
        }
        if (m == limit_model) {
          break
        }
      } # Fin procesando modelo
    } # Fin proceso bloques de 1000
    
    fila_inicio <- (replica - 1) * esquemas + 1
    fila_fin <- replica * esquemas
    conteos_totales[fila_inicio:fila_fin, ] <- conteo_replica
    conteo_ceros_replica[2] <- m
    conteo_ceros[replica, ] <- conteo_ceros_replica
    # Guardar resultados tras procesar cada réplica
    #write.csv(x = conteos_totales, file = nombre_archivo_conteos, row.names = FALSE)
    #write.csv(x = conteo_ceros, file = nombre_archivo_ceros, row.names = FALSE)
  }
  cat("Fin de cálculos")
}