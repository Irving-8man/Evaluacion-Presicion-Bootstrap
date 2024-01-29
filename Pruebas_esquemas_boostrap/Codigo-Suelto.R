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