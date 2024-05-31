#Tesis Regresion lineal con esquemas robustos
#Irving Geyler Cupul Uc

#cargando librerias
library("robustbase")
library("boot")
library("bootstrap")
library("readxl")


#Funcion de apoyo para ponderar residuales
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


#Obtener los residualesBT para tratar
#Que sea tipoEsq
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


#Funcion para contruir los intervalos de confianza Bootstrap
ContruirIntervBoot <- function(data,R2,muestrasR2Boot,B=100,nivConfianza=0.95,tipo){
  intervalo <- numeric(2) 
  alpha <- 1-nivConfianza
  z <- data[,1]
  y <- data[,2]
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
      n <- length(y)
      p0 <- length(which(vectorR2Bootstrap >= R2))/B #proporcion de e Rˆ2i’s ≥ Rˆ2 duda
      vectorR2i <- numeric(n)
      for (i in 1:n){
        vectorR2i[i] <- summary(lm(y[-i]~z[-i]))$r.squared
      }
      R2iPromedio <- mean(vectorR2i)
      diferenciasR2iPi <- R2iPromedio - vectorR2i
      sumaDiferenciasCuad <- sum(diferenciasR2iPi**2)
      sumaDiferenciasCubo <- sum(diferenciasR2iPi**3)
      a <- sumaDiferenciasCubo/(6*(sumaDiferenciasCuad**1.5))
      
      ZL <- qnorm(1-p0) + (qnorm(1-p0) - qnorm(1-alpha/2) )/(1-a *(qnorm(1-p0) - qnorm(1-alpha/2)))
      ZU <- qnorm(1-p0) + (qnorm(1-p0) + qnorm(1-alpha/2) )/(1-a *(qnorm(1-p0) + qnorm(1-alpha/2))) 
      
      ICInfBootBCa <- quantile(vectorR2Bootstrap,probs= pnorm(ZL))
      ICSupBootBCa <- quantile(vectorR2Bootstrap,probs= pnorm(ZU))
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



#funcion principal
CalPrecicion <- function (data,alpha,nivConfianza,caso){
  z <- data[,1]
  y <- data[,2]
  n <- length(z)
  print("Entrando en calculos")
  print(data)
  
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




###############################################################################################
#Recuperacion de las muestras y la data
#Sea para el parametro modelo, 1-Precisos, 2-Imprecisos
#Y sea para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
ObtenerArchivoData <- function(modelo,caso,numMuestras = c(10,15,20,25,30,35) ){
  directorio <- "../Evaluacion-Presicion-Bootstrap/Modelos Simulados"
  #Elige modelo y preparar la ruta
  carp_model <- switch(modelo,
                       "Precisos" = "Precisos",
                       "Imprecisos" = "Imprecisos",
                       stop("Modelo no válido"))
  arch_model <- switch(modelo,
                       "Precisos" = "EP",
                       "Imprecisos" = "EI",
                       stop("Modelo no válido"))
  ruta_model <- file.path(directorio, carp_model)

  
  #Elige el caso, preparar la carpeta y la ruta
  arch_caso <- c("NVC","NVD","NNVC","NNVD")
  tipo_caso <- switch(caso,{arch_caso[1]},{arch_caso[2]},{arch_caso[3]},{arch_caso[4]},stop())
  model_caso <- paste(arch_model, tipo_caso, sep = "")

  #ruta preparada
  ruta_model_caso <- file.path(ruta_model, model_caso)
  
  # Obteniendo archivos en la carpeta del caso y sus R2
  archivos <- sort(list.files(ruta_model_caso))
  #archivos_R2 <- sort(archivos[grepl("^R2", archivos)])
  #archivos_Muestras <- sort(archivos[!grepl("^R2", archivos)])
  
  # Función para encontrar los archivos de muestra y R2 correspondientes
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
  
  #Comenzando el procesamiento fuerte
  replicas <- 5
  for(muestra in numMuestras){
    archivos_encontrados <- encontrar_archivos(muestra)
    if (!is.null(archivos_encontrados)) {
      CalPrecMuestras(archivos_encontrados,caso, replicas, nivConfianza = 0.95, N = muestra)
    } else {
      cat("No se encontró un archivo de muestra o R2 correspondiente para el tamaño de muestra:", muestra, "\n")
    }
    
  }
  
}


# ######################################3
# #Función para implementar el calculo de las precision de las muestras
# # Dado archivos_encontrados <-  list(muestra,R2) con rutas de archivos
# # para el parametro caso, 1-NVC, 2-NVD, 3-NNVC, 4-NNVD
# # para replicas sea el número de replicas
# # con nivConfinza entre 0 - 1
# # y sea para N el tamaño de las muestras
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
  limit_model <- 10
  
  # Inicializar el data frame de conteos
  conteos_totales <- data.frame(
    modelosEvaluados = integer(),
    esquema1Ganador = integer(),
    esquema2Ganador = integer(),
    esquema3Ganador = integer(),
    esquema4Ganador = integer(),
    esquema5Ganador = integer(),
    esquema6Ganador = integer(),
    esquema7Ganador = integer(),
    esquema8Ganador = integer(),
    sinesquemaGanador = integer(),
    intervalo1Ganador = integer(),
    intervalo2Ganador = integer(),
    intervalo3Ganador = integer(),
    sinintervaloGanador = integer(),
    ganadorPorUnico = integer(),
    ganadorEmpate = integer(),
    ganadorEmpateTriple = integer(),
    ningunGanador = integer()
  )
  
  # Ajustando por replicas
  for (replica in 1:replicas) {
    print(paste("Replica #", replica))
    
    # Inicializar contadores para la réplica actual
    conteo_replica <- list(
      modelosEvaluados = 0,
      esquema1Ganador = 0,
      esquema2Ganador = 0,
      esquema3Ganador = 0,
      esquema4Ganador = 0,
      esquema5Ganador = 0,
      esquema6Ganador = 0,
      esquema7Ganador = 0,
      esquema8Ganador = 0,
      sinesquemaGanador = 0,
      intervalo1Ganador = 0,
      intervalo2Ganador = 0,
      intervalo3Ganador = 0,
      sinintervaloGanador = 0,
      ganadorPorUnico = 0,
      ganadorEmpate = 0,
      ganadorEmpateTriple = 0,
      ningunGanador = 0
    )
    
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
        
        #Mostrar los mejores
        rownames(resultados_ganadores) <- NULL  # Eliminar nombres automáticos de filas
        print(resultados_ganadores)
        
        
        
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
        
        if (conteo_replica$modelosEvaluados == limit_model) {
          break
        }
      } # Fin procesando modelo
      if (conteo_replica$modelosEvaluados == limit_model) {
        break
      }
    } # Fin proceso bloques de 1000
    
    # Añadir los conteos de la réplica actual al data frame total
    conteos_totales <- rbind(conteos_totales, conteo_replica)
  } # Fin replicas
  
  titulo <- paste0("Efica", N,".csv")
  ruta<-paste("/home/irving/Documentos/tesis/Evaluacion-Presicion-Bootstrap/Resultados", titulo, sep = "/")
  write.csv(conteos_totales, ruta, row.names = FALSE)
}

