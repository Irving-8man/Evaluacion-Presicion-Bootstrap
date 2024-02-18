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