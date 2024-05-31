datos_texto <- "0.360 0.361
0.340 0.296
0.270 0.325
0.289 0.287
0.248 0.259
0.228 0.269
0.298 0.284
0.237 0.206
0.270 0.214
0.257 0.236
0.265 0.241
0.263 0.263
0.283 0.278
0.285 0.294
0.201 0.202
0.220 0.183
0.194 0.196
0.198 0.188
0.200 0.207
0.221 0.207
0.200 0.207
0.189 0.207
0.211 0.207
0.242 0.207
0.276 0.268
0.283 0.268
0.216 0.233
0.218 0.234"

# Leer datos y separar en dos vectores
datos <- read.table(text = datos_texto)
vector1 <- datos$V1
vector2 <- datos$V2

datos <- data.frame(ValObsKgY = vector1, ValPreKgZ = vector2)

# Guardar el data frame como un archivo CSV
write.csv(datos, file = "N-VC.csv", row.names = FALSE)