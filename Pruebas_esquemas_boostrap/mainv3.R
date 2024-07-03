#Version data frame y version nueva de Bca
#Tesis Regresion lineal con esquemas robustos
#Irving Geyler Cupul Uc

#cargando librerias
library("robustbase")
library("boot")
library("bootstrap")
library("readxl")


#Función de apoyo para ponderar residuales
#Sea residualesRob <- los residuales robustos del modelo
#CMERob <- 
PodResidRobu <- function(residualesRob,CMERob){
  consPes <- 3
  n <- length(residualesRob)
  x <- abs(residualesRob)/sqrt(CMERob)
  w <- rep(1,n)
  xx <- which(x > consPes) 
  w[xx] <- (consPes / w[xx])
  resRobuPonder<- w*residualesRob
  return(resRobuPonder)
}

CalcularR2Bootstrap <- function(y,z, yAjRob,residuales,residualesRob,residualesRP, hii,n, B, tipo) {
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
        #agregar residuales balanceados (Mtro Luis Colorado)
        posI=(i-1)*n+1
        posF=i*n
        VPosi=NPerm[posI:posF]
        residualesRP[VPosi]
      },
      {
        #agregar residuales pareados (Mtro Luis Colorado)
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


#Intervalos de confianza para las muestras Bootstrap 
ContruirIntervBoot <- function(data,R2,muestrasR2Boot,B=100,nivConfianza=0.95,tipo){
  intervalo <- numeric(2) 
  alpha <- 1-nivConfianza
  z <- as.numeric(data[[1]])
  y <- as.numeric(data[[2]])
  x <- vectorR2Bootstrap <- muestrasR2Boot
  
  intervalo <- switch(
    tipo,
    {
      #percentil fuente(An-introduction-to-bootstrap_bradley-efron.pdf/pag92)
      puntosCriticos <- quantile(vectorR2Bootstrap, c(alpha/2, 1 - alpha/2)) # Aproximación bootstrap de los puntos críticos
      as.vector(puntosCriticos)
    },
    {
      #fuente(tesis-balam) Bca
      n<-length(z)
      # Paso 2: Calcular el factor de corrección de sesgo (z0)
      z0 <- qnorm(mean(vectorR2Bootstrap < R2))
      
      # Paso 3: Calcular el factor de aceleración (a)
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
      
      # Paso 4: Ajustar percentiles
      z_alfa1 <- qnorm(alpha / 2)
      z_alfa2 <- qnorm(1 - alpha / 2)
      alfa1 <- pnorm(z0 + (z0 + z_alfa1) / (1 - a * (z0 + z_alfa1)))
      alfa2 <- pnorm(z0 + (z0 + z_alfa2) / (1 - a * (z0 + z_alfa2)))
      
      # Paso 5: Determinar el intervalo de confianza
      ICInfBootBCa <- quantile(muestrasR2Boot, alfa1)
      ICSupBootBCa <- quantile(muestrasR2Boot, alfa2)
      as.vector(c(ICInfBootBCa, ICSupBootBCa))
    },
    {
      #ABC fuente(An%20Introduction%20to%20Bootstrap%20Methods%20with%20Applications%20to%20R.pdf/107)
      #fuen(An-introduction-to-bootstrap_bradley-efron.pdf/pag.101)
      #duda aun con esto
      #Agregar funcion de ejemplo
      fabc <-function(x,w) x%*%w #ABC usando sus pesos
      w <- rep( (1/length(x)), length(x)) 
      intervalo<- abc.ci(x, fabc, conf = nivConfianza) #ABC method C.I.
      as.vector(c(intervalo[2],intervalo[3]))#Revisíon
    },
    stop("Intervalo no válido")
  )
  
  return (intervalo)
}



#Calcular la precisión de un módelo
CalPrecicion <- function (data,alpha,nivConfianza,caso){
  z <- as.numeric(data[[1]])
  y <- as.numeric(data[[2]])
  n <- length(z)
  print("Entrando en cálculos")
  
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
  #Num de remuestreos
  numRemues <- 8
  
  valoresBootR <- matrix(0,nrow = B, ncol = numRemues) #todas los R2 generados por muestra 
  #Minimos cuadrados
  if(NVC == caso){
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
  #print(valoresBootR) #Fin de esquemas boot
  
  ###Empieza creacion de intervalos de confianza
  resultadosInter <- vector("list", numRemues)
  
  #Calculo de los intervalos
  for(i in 1:numRemues){
    muestrasR2Boot <- valoresBootR[,i]
    #Percentil
    perc <- ContruirIntervBoot(data,R2,muestrasR2Boot,B,nivConfianza=0.95,tipo=1)
    #bca
    bca <- ContruirIntervBoot(data,R2,muestrasR2Boot,B,nivConfianza=0.95,tipo=2)
    #Percentil-simetrico
    abc <- ContruirIntervBoot(data,R2,muestrasR2Boot,B,nivConfianza=0.95,tipo=3)
    resultadosInter[[i]] <- list(perc, bca,abc)
  }
  
  #print(resultadosInter)
  return(resultadosInter)
}




##################################################################################################
#Procesamiento de datos y cálculo de datos de interes



#Funcion para la recuperacion de muestras y sus R2 estimadas
#Sea para el parametro modelo, 1-Precisos, 2-Imprecisos
#Y sea para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
ObtenerArchivoData <- function(modelo,caso,numMuestras = c(10,15,20,25,30,35) ){
  #Modificar dependiendo de la carpeta donde este alojada
  directorio <- "../Evaluacion-Presicion-Bootstrap/Modelos Simulados"
  #Elige un modelo
  carp_model <- switch(modelo,
                       "Precisos" = "Precisos",
                       "Imprecisos" = "Imprecisos",
                       stop("Modelo no válido"))
  arch_model <- switch(modelo,
                       "Precisos" = "EP",
                       "Imprecisos" = "EI",
                       stop("Modelo no válido"))
  ruta_model <- file.path(directorio, carp_model)
  
  
  #Elige un caso
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
  
  #Procesamiento de la replicas del caso 
  replicas <- 2
  for(muestra in numMuestras){
    archivos_encontrados <- encontrar_archivos(muestra)
    if (!is.null(archivos_encontrados)) {
      CalPrecMuestrasv2(archivos_encontrados,caso, replicas, nivConfianza = 0.95, N = muestra)
    } else {
      cat("No se encontró un archivo de muestra o R2 correspondiente para el tamaño de muestra:", muestra, "\n")
    }
    
  }
  
}






#Función para implementar el calculo de las precision de las muestras
# Dado archivos_encontrados <-  list(muestra,R2) con rutas de archivos
# para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
# para replicas sea el número de replicas
# con nivConfinza entre 0 - 1
# y sea para N el tamaño de las muestras

#Version dataframe
CalPrecMuestrasv2 <- function(archivos_encontrados, caso, replicas, nivConfianza, N) {
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
  
  # Resultados del analisis
  nombre_cols <- c("Replica","esquema","FrecEficIB1","FrecEficIB2","FrecEficIB3",
                   "FrecEficIB1Unico","FrecEficIB2Unico", "FrecEficIB3Unico",
                   "FrecEficIB1Emp2", "FrecEficIB2Emp2", "FrecEficIB3Emp2",
                   "FrecEficIB1Emp3", "FrecEficIB2Emp3","FrecEficIB3Emp3")
  nombre_cols_cer <- c("Replicas", "NumMod", "Esq1", "Esq2", "Esq3", "Esq4", "Esq5", "Esq6", "Esq7", "Esq8")

  conteos_totales <- matrix(0, ncol = length(nombre_cols), nrow = replicas * esquemas)
  conteo_ceros <- matrix(0,ncol=length(nombre_cols_cer),nrow = replicas)
  colnames(conteos_totales) <- nombre_cols
  
  no_entro_ninguno <-0
  
  # Procesando las replicas
  for (replica in 1:replicas) {
    print(paste("Replica #", replica))
    m <- 0
    
    # Vectores que contienen los datos por defecto
    conteo_replica <- matrix(0, nrow = 8,ncol = 14)
    replica_vector <- rep(replica, each = esquemas)
    esquema_vector <- rep(1:esquemas)
    matriz_inicial <- cbind(replica_vector, esquema_vector)
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
      R2_block <- R2_replica[i:500]
      
      Rmod <- 0
      for (j in seq(1, ncol(block_caso), by = cols_por_model)) {
        model_end <- min(j + cols_por_model - 1, ncol(block_caso))
        modeloActual <- block_caso[, j:model_end]
        R2_modelo <- R2_block[Rmod+1]
        Rmod<-Rmod+1
        
        # Aquí procesas los datos del modelo actual
        print(paste("Procesando modelo con columnas:", j, "-", model_end))
        resultadosInter <- CalPrecicion(modeloActual, alpha, nivConfianza, caso)
        print("procesado")
        print(R2_modelo)
        
        
        #Procesando resultados
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
        rownames(resultados_df) <- NULL
        print(resultados_df)
        
        # Procesando ganadores
        resultados_ganadores <- data.frame(
          Esquema = integer(),
          IntervaloGanador = integer(),
          GanadorPor = integer(),
          TamanoIntervalo = numeric()
        )
        #Procesando los equemas
        for (es in 1:length(resultadosInter)) {
          esquema <- resultadosInter[[es]]
          intervalos_ganadores <- data.frame(
            Intervalo = integer(),
            Inferior = numeric(),
            Superior = numeric(),
            Longitud = numeric()
          )
          
          #Procesando los intervalos por esquema
          for (esj in 1:length(esquema)) {
            intervalo <- esquema[[esj]]
            longitud <- intervalo[2] - intervalo[1]
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            
            #Filtrar a los que si tienen al R2
            if (R2_intervalo == 1) {
              intervalos_ganadores <- rbind(intervalos_ganadores, data.frame(
                Intervalo = esj,
                Inferior = intervalo[1],
                Superior = intervalo[2],
                Longitud = longitud
              ))
            }
          }
          
          #Obtener al mejor intervalo por esquema y como gano
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
          } 
          
          #Agregando al mejor por esquema
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
              if (numIntervalo == 1) conteo_replica[numEsquema,3] <- conteo_replica[numEsquema,3] + 1
              if (numIntervalo == 2) conteo_replica[numEsquema,4] <- conteo_replica[numEsquema,4] + 1
              if (numIntervalo == 3) conteo_replica[numEsquema,5] <- conteo_replica[numEsquema,5] + 1
              
              intervalos_ganadores[[length(intervalos_ganadores) + 1]] <- list(Intervalo = numIntervalo, Longitud = intervalo[2] - intervalo[1])
              #print(intervalos_ganadores)  
            }
          }
          
          # Determinar los ganadores únicos, empates y triple empate
          if (length(intervalos_ganadores) == 1) {
            #Definiendo ganador unico
            if (intervalos_ganadores[[1]]$Intervalo == 1) conteo_replica[numEsquema,6] <- conteo_replica[numEsquema,6] + 1
            if (intervalos_ganadores[[1]]$Intervalo == 2) conteo_replica[numEsquema,7] <- conteo_replica[numEsquema,7] + 1
            if (intervalos_ganadores[[1]]$Intervalo == 3) conteo_replica[numEsquema,8] <- conteo_replica[numEsquema,8] + 1
          } else if (length(intervalos_ganadores) == 2) {
            #Definiendo ganador por empate
            mejor_intervalo <- intervalos_ganadores[[which.min(sapply(intervalos_ganadores, function(x) x$Longitud))]]
            if (mejor_intervalo$Intervalo == 1) conteo_replica[numEsquema, 9] <- conteo_replica[numEsquema, 9] + 1
            if (mejor_intervalo$Intervalo == 2) conteo_replica[numEsquema, 10] <- conteo_replica[numEsquema, 10] + 1
            if (mejor_intervalo$Intervalo == 3) conteo_replica[numEsquema, 11] <- conteo_replica[numEsquema, 11] + 1
            
          } else if (length(intervalos_ganadores) == 3) {
            #Definiendo ganador por triple empate
            mejor_intervalo <- intervalos_ganadores[[which.min(sapply(intervalos_ganadores, function(x) x$Longitud))]]
            if (mejor_intervalo$Intervalo == 1) conteo_replica[numEsquema, 12] <- conteo_replica[numEsquema, 12] + 1
            if (mejor_intervalo$Intervalo == 2) conteo_replica[numEsquema, 13] <- conteo_replica[numEsquema, 13] + 1
            if (mejor_intervalo$Intervalo == 3) conteo_replica[numEsquema, 14] <- conteo_replica[numEsquema, 14] + 1
            
          }else{
            #Ningun ganador
            no_entro_ninguno <-no_entro_ninguno+1
          }
          
          
        }#Fin de conteo
        
        
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
  #write.csv(x = conteos_totales, file = "conteos.csv", row.names = FALSE) 
}