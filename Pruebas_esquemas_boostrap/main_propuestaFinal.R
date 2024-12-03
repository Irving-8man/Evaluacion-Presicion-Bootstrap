data1 <- read.csv("C:/Users/geyle/Downloads/Irving/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_NVC.csv")
data2 <- read.csv("C:/Users/geyle/Downloads/Irving/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_NNVC.csv")
data3 <- read.csv("C:/Users/geyle/Downloads/Irving/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_NVD.csv")
data4 <- read.csv("C:/Users/geyle/Downloads/Irving/tesis/Evaluacion-Presicion-Bootstrap/Data_pruebas/Datos_Caso_NNVD.csv")


#Propuesta final
#Sean los parametros: data una matriz,
# con la columna 1 las "z" la columna 2 las "y" 
# y niConfinza un nivel de confinza entre 0-1
PropFCalcPrecModl <- function(data, nivConfianza=0.95){
  library("nortest")
  library("lmtest")
  library("robustbase")
  
  z <<- as.numeric(data[[1]])
  y <<- as.numeric(data[[2]])
  n <<- length(z)
  B <<- 2500 #Remuestras bootstrap necesarias
  caso <<- 0
  IC_proces <<- 1 #Construir IC, metodo percentil como inicial
  alpha <<- 1-nivConfianza
  #Casos posibles
  NVC <- 1 #Normalidad- homocedasticidad
  NVD <- 2 #Normalidad-heterocidasticidad
  NNVC <- 3 #No normalidad-homocedasticidad
  NNVD <- 4 #No normalidad-heterocidastecidad 
  #Verificación de supuestos
  hay_normalidad <- FALSE
  varianza_constante <- FALSE
  media_cero <- FALSE
  hay_independencia <-FALSE
  conceptosClave <<- rep(0,4) #Banderas de apoyo en las 4 pruebas
  limite_precision <<-0.7
  
  #Regresion lineal simple
  modeloLineal <- lm(y~z)
  residuales <- residuals(modeloLineal)
  R2 <- summary(modeloLineal)$r.squared
  hii <- hatvalues(modeloLineal)
  yAju <- fitted(modeloLineal)
  
  # Función para imprimir texto con colores 
  print_colored <- function(text, text_color = 37, bg_color = NULL) {
    if (!is.null(bg_color)) {
      cat(sprintf("\033[%d;%dm%s\033[0m", text_color, bg_color, text))# Texto con fondo
    } else {
      cat(sprintf("\033[%dm%s\033[0m", text_color, text)) # Solo texto coloreado
    }
  }
  
  #Supuesto de normalidad
  {
    pValShap <- shapiro.test(residuales)$p.value  #Shapiro-Wilk
    ValCShap <- shapiro.test(residuales)$statistic
    pValLILLIE <- lillie.test(residuales)$p.value   #Lilliefort
    ValCLILLIE <- lillie.test(residuales)$statistic
    #Resultados
    pValor <- c(pValShap, pValLILLIE)
    ValorCal <- c(ValCShap, ValCLILLIE)
    Estadistica <- c("Shapiro-Wilk","Liliefort")
    tablaNormal <- data.frame(Estadistica,ValorCal,pValor)
    print_colored("\nPRUEBA DE NORMALIDAD\n", 37, 40)
    print(tablaNormal)
    
    pVal_Min <- min(pValor)
    hay_normalidad <-  pVal_Min > alpha
    
    if (hay_normalidad) {
      cat("\nConclusión: Se cumple el supuesto de normalidad con Shapiro y Lilliefort al",alpha*100,"%.\n")
      conceptosClave[1] <- 1
    }else{
      cat("\nConclusión: No se cumple el supuesto de normalidad con Shapiro y Lilliefort al",alpha*100,"%.\n")
      conceptosClave[1] <- 2
    }
  }#Fin prueba normalidad
  
  
  #Media cero
  {
    #Prueba
    pValor <-t.test(residuales)$p.value #test T-student
    ValorCal <- t.test(residuales)$statistic
    #Construccion de la tabla
    Estadistica = c("T-Student")
    tablaVar <- data.frame(Estadistica, ValorCal, pValor)
    print_colored("\nPRUEBA T-STUDENT PARA MEDIA CERO EN LOS RESIDUALES\n", 37, 40)
    print(tablaVar)
    media_cero <- alpha < pValor
    
    if(media_cero){
      cat("\nConclusión: Se cumple el supuesto de media cero en los residuales al ",alpha*100,"%.\n")
      conceptosClave[3] <- 1
    }else{
      cat("\nConclusión: No se cumple el supuesto de media cero en los residuales al",alpha*100,"%.\n")
      conceptosClave[3] <- 2
    }
  }#Fin de prueba de media cero
  
  
  
  #Prueba de varianzas
  {
    pValor <- 0
    prueba_var <- NULL
    aviso_varianza <- NULL
    ValorCal <- 0
    
    if(hay_normalidad){
      #Aplicar Breusch-Pagan Test
      pValor <- bptest(modeloLineal)$p.value
      ValorCal <- bptest(modeloLineal)$statistic
      conceptosClave[2] <- 1
      Estadistica <- c("Breush-Pagan")
    }else{
      #Aplicar White test
      pValor <- bptest(modeloLineal, varformula = ~ I(z^2))$p.value
      ValorCal <- bptest(modeloLineal, varformula = ~ I(z^2))$statistic
      conceptosClave[2] <- 2
      Estadistica <- c("White")
    }
    
    valorTestVa <-c(pValor)
    tablaVarianza <- data.frame(Estadistica,ValorCal, pValor)
    
    aviso_varianza <- switch(
      conceptosClave[2],
      {
        print_colored("\nPRUEBA DE IGUALDAD DE VARIANZAS CON BREUSH-PAGAN\n", 37, 40)
      },
      {
        print_colored("\nPRUEBA DE IGUALDAD DE VARIANZAS CON WHITE\n", 37, 40)
      },
      stop("Prueba invalida")
    )
    
    print(tablaVarianza)
    
    varianza_constante <- alpha < pValor
    
    if(varianza_constante){
      aviso_varianza <- switch(
          conceptosClave[2],
          {
            cat("\nConclusión: Se cumple el supuesto de varianza constante con el estadístico de Breush-Pagan al",alpha*100,"%.\n")   
          },
          {
            cat("\nConclusión: Se cumple el supuesto de varianza constante con el estadístico de White al",alpha*100,"%.\n")
          },
          stop("Prueba invalida")
        )
    }else{
      aviso_varianza <- switch(
        conceptosClave[2],
        {
          cat("\nConclusión: No se cumple el supuesto de varianza constante con el estadístico de Breush-Pagan al",alpha*100,"%.\n")   
        },
        {
          cat("\nConclusión: No se cumple el supuesto de varianza constante con el estadístico de White al",alpha*100,"%.\n")
        },
        stop("Prueba invalida")
      )
    }
  }#Fin prueba varianzas
  
  
  
  #Independencia
  {
    #Prueba
    pValor <- dwtest(modeloLineal)$p.value #Durbin-Watson test
    ValorCal <- dwtest(modeloLineal)$statistic
    Estadistica = c("Durbin-Watson test")
    tablaVar <- data.frame(Estadistica,ValorCal, pValor)
    print_colored("\nPRUEBA DE DURBIN-WATSON PARA INDEPENDENCIA\n", 37, 40)
    print(tablaVar)
    hay_independencia <- alpha < pValor
    
    if(hay_independencia){
      cat("\nConclusión: Se cumple el supuesto de independencia con DURBIN-WATSON al",alpha*100,"%.\n")
      conceptosClave[4] <- 1
    }else{
      cat("\nConclusión: No se cumple el supuesto de independencia con DURBIN-WATSON al",alpha*100,"%.\n")
      conceptosClave[4] <- 2
    }
  }#Fin de prueba independencia
  
  
  

  #Decision de caso que nos encontramos
  if(hay_normalidad && varianza_constante){
    caso <- NVC
    print_colored("\n ****** CASO: NORMALIDAD - HOMOCEDASTICIDAD (NVC) ******\n", 37, 34)
  }else{
    # Uso de estimador robusto MM  
    modeloLinealRob <<- lmrob(y ~ z, method = "MM")
    
    if (hay_normalidad && !varianza_constante) {
      caso <- NVD
      print_colored("\n ****** CASO: NORMALIDAD - HETEROCIDASTICIDAD (NVD) ******\n", 37, 34)
    }
    if(!hay_normalidad && varianza_constante){
      caso <- NNVC
      print_colored("\n ****** CASO: NO NORMALIDAD - HOMOCEDASTICIDAD (NNVC) ******\n", 37, 34)
    }
    if(!hay_normalidad && !varianza_constante){
      caso <- NNVD
      print_colored("\n ****** CASO: NO NORMALIDAD - HETEROCIDASTICIDAD (NNVD) ******\n", 37, 34)
      IC_proces <- 2
    }
  }
  
  RsBoot <<- numeric(B)#Remuestra de R^2
  #Comienza remuestreos Bootstrap para el coeficiente de determinación
  {
    if(caso != NNVD){

      residuales_utilizar <- switch(
        caso,
        {
          residuales #Usar minimos cuadrados -> NVC
        },
        {
          resTemp <- modeloLinealRob$residuals 
          CMERob <- modeloLinealRob$scale**2
          consPes <- 3
          x <- abs(resTemp)/sqrt(CMERob)
          w <- rep(1,n)
          xx <- which(x > consPes) 
          w[xx] <- (consPes / w[xx])
          w*resTemp #Usar residuales robustos ponderados MM-estimador->NVD
        },
        {
          modeloLinealRob$residuals #Usar residuales robustosMM-estimador ->NNVC
        },
        stop("Caso no valido para residuales")
      )
      
      yajus_usar <- if (caso == NVC ) yAju else modeloLinealRob$fitted.values
      
      #Efectuando Esquema Liu 2
      sqrt_hii <- sqrt(1 - hii)
      for (i in 1:B) {
        media1 <- 0.5 * sqrt(17 / 6) + sqrt(1 / 6)
        media2 <- 0.5 * sqrt(17 / 6) - sqrt(1 / 6)
        H <- rnorm(n, media1, sqrt(0.5))
        D <- rnorm(n, media2, sqrt(0.5))
        tt <- H * D - media1 * media2
        resBoot <- (tt * residuales_utilizar) / sqrt_hii
        yBoots <- yajus_usar + resBoot
        modeloBoots <- lm(yBoots ~ z)
        RsBoot[i] <- summary(modeloBoots)$r.squared
      }
    }else{#Caso NNVD
      N     <- rep(1:n,B)
      NPerm <- sample(N)
      for (i in 1:B) {
        posI <- (i-1)*n+1
        posF <- i*n
        pares_boots<- NPerm[posI:posF]
        YP <- y[pares_boots]
        ZP <- z[pares_boots]
        modeloBoots <- lm(YP ~ ZP)
        RsBoot[i] <- summary(modeloBoots)$r.squared
      }
    }
  }#Fin Bootstrap
  
  intervalo <- numeric(2)#Intervalo
  #Construir intervalo de confinza para R^2
  {
    intervalo <- switch(
      IC_proces,
      {
        puntosCriticos <- quantile(RsBoot, c(alpha/2, 1 - alpha/2))
        as.vector(puntosCriticos)
      },
      {
        z0 <- qnorm(mean(RsBoot < R2))
        suma0 <-0
        suma02<-0
        suma03<-0
        for(i in 1:n){
          R2MI<-summary(lm(y[-i]~z[-i]))$r.squared
          suma0<-R2MI+suma0
        }
        R2PMI<- suma0/n
        for(i in 1:n){
          Dif0<-R2PMI-summary(lm(y[-i]~z[-i]))$r.squared
          suma02<-Dif0^2+suma02
          suma03<-Dif0^3+suma03
        }
        a<-suma03/(6*(suma02^1.5))
        z_alfa1 <- qnorm(alpha / 2)
        z_alfa2 <- qnorm(1 - alpha / 2)
        alfa1 <- pnorm(z0 + (z0 + z_alfa1) / (1 - a * (z0 + z_alfa1)))
        alfa2 <- pnorm(z0 + (z0 + z_alfa2) / (1 - a * (z0 + z_alfa2)))
        
        ICInfBCa <- quantile(RsBoot, alfa1)
        ICSupBCa <- quantile(RsBoot, alfa2)
        as.vector(c(ICInfBCa, ICSupBCa))
      },
      stop("Intervalo no válido")
    )
  }
  
  #Resultado final
  {
    resultado <- switch(
      IC_proces,
      {
        print_colored("\nPRECISION (R2) CON EL ICB PERCENTIL-ESQUEMA BOOTSTRAP LIU2   \n", 37, 40)
      },
      {
        print_colored("\nPRECISION (R2) CON EL ICB BCa-ESQUEMA BOOTSTRAP PAREADO BALANCEADO \n", 37, 40)
      },
      stop("Prueba invalida")
    )
    
    Atributos <- c("R2","R2BootMedia","DesvEstR2Boots","LIR2","LSR2")
    Valores <- c(R2,mean(RsBoot),sd(RsBoot),intervalo[1],intervalo[2])
    tabla_IC <-data.frame(Atributos,Valores)
    print(tabla_IC)
    
    R2_preciso <- ifelse( (limite_precision>=intervalo[1] & limite_precision <= intervalo[2]) || (intervalo[1] >= limite_precision) , 1, 0)
    
    conclusion <- switch(
      IC_proces,
      {
        if (R2_preciso == 1) {
          cat("\nConclusión: El modelo es preciso con el ICB Percentil al ",nivConfianza*100,"%.\n")   
        }else{
          cat("\nConclusión: El modelo es impreciso con el ICB Percentil al",nivConfianza*100,"%.\n")
        }
      },{
        if (R2_preciso == 1) {
          cat("\nConclusión: El modelo es preciso con el ICB BCa al",nivConfianza*100,"%.\n")   
        }else{
          cat("\nConclusión: El modelo es impreciso con el ICB BCa al",nivConfianza*100,"%.\n")
        }
      }
    )
    
    cat("\n Criterio: Si el ",limite_precision, " está contenido en el ICB o es mayor igual al limite inferior del ICB.\n")
  }#Fin de bloque final
  
}#Fin propuesta final





#Caso NVC
PropFCalcPrecModl(data = data1, nivConfianza = 0.95)

#Caso NNVD
PropFCalcPrecModl(data = data4, nivConfianza = 0.95)






#Caso NNVC
PropFCalcPrecModl(data = data2, nivConfianza = 0.95)


#Caso NVD
PropFCalcPrecModl(data = data3, nivConfianza = 0.95)





