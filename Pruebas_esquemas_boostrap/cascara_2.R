conteo_replica <- matrix(0, nrow = 8, ncol = length(nombre_cols))
replica_vector <- rep(replica, each = esquemas)
esquema_vector <- rep(1:esquemas)
matriz_inicial <- cbind(replica_vector, esquema_vector)
conteo_replica[, 1:2] <- matriz_inicial

conteo_ceros_replica <- numeric(length(nombre_cols_cer))
conteo_ceros_replica[1] <- replica

fila_inicio <- (replica - 1) * N + 1
fila_fin <- replica * N
replica_data <- data_muestra[fila_inicio:fila_fin, ]
R2_replica <- as.numeric(data_R2[replica, ])

for (i in seq(1, ncol(replica_data), by = block)) {
  block_end <- min(i + block - 1, ncol(replica_data))
  block_caso <- replica_data[, i:block_end]
  R2_block <- as.numeric(R2_replica[i:500])
  
  Rmod <- 0
  
  for (j in seq(1, ncol(block_caso), by = cols_por_model)) {
    model_end <- min(j + cols_por_model - 1, ncol(block_caso))
    modeloActual <- block_caso[, j:model_end]
    R2_modelo <- R2_block[Rmod+1]
    if (is.na(R2_modelo)) {
      stop("El valor de R2_modelo es NA.")
    }
    Rmod <- Rmod + 1
    resultadosInter <- EvalPrecisionModel(modeloActual, alpha, nivConfianza, caso)#Calcular intervalos
    
    #Procesando resultados por esquema
    for (numEsquema in 1:length(resultadosInter)) {
      resultados_esquema <- resultadosInter[[numEsquema]]
      intervalos_ganadores <- list()
      
      #Validando intervalos por esquema
      for (numIntervalo in 1:length(resultados_esquema)) {
        intervalo <- resultados_esquema[[numIntervalo]]
        ############seguir aqui
        if (!any(is.na(intervalo)) && length(intervalo) == 2) {
          
          if (!is.na(R2_modelo)) {
            R2_intervalo <- ifelse(R2_modelo >= intervalo[1] & R2_modelo <= intervalo[2], 1, 0)
            if (R2_intervalo == 1) {
              if (numIntervalo == 1) conteo_replica[numEsquema, 3] <- conteo_replica[numEsquema, 3] + 1
              if (numIntervalo == 2) conteo_replica[numEsquema, 4] <- conteo_replica[numEsquema, 4] + 1
              intervalos_ganadores[[length(intervalos_ganadores) + 1]] <- list(Intervalo = numIntervalo, Longitud = intervalo[2] - intervalo[1])
            }
          } else {
            warning(paste("Valor NA en R2_modelo:", R2_modelo, "en esquema", numEsquema, "intervalo", numIntervalo))
          }
        } else {
          warning(paste("Intervalo inválido en esquema", numEsquema, "intervalo", numIntervalo, 
                        "Intervalo:", paste(intervalo, collapse = ","), 
                        "con R2_modelo:", R2_modelo))
        }
        
      }#Procesando 
      
      
      
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
    
    
    
    
    
    
    
    
    
    
    #Segumiento de modelos
    num_m <- num_m + 1
    if (num_m %% 50 == 0) {
      print(paste("Procesados", m, "modelos de replica", replica, "/",replicas))
    }
    if (num_m == limit_model) {
      break
    }
  } # Fin procesando modelo
} # Fin proceso bloques de 1000

#Guardar resultados obtenidos de replica
fila_inicio <- (replica - 1) * esquemas + 1
fila_fin <- replica * esquemas
conteos_totales[fila_inicio:fila_fin, ] <- conteo_replica
conteo_ceros_replica[2] <- m
conteo_ceros[replica, ] <- conteo_ceros_replica

#Guardar resultados en xlsx
writeData(wb_conteos, "Conteos", conteos_totales, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb_ceros, "Ceros", conteo_ceros, startCol = 1, startRow = 1, rowNames = FALSE)
saveWorkbook(wb_conteos, nombre_archivo_conteos, overwrite = TRUE)
saveWorkbook(wb_ceros, nombre_archivo_ceros, overwrite = TRUE)