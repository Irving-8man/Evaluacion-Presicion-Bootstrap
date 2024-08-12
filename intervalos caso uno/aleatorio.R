# Definir los tamaños de muestra
tamaños_muestra <- c(10, 15, 20, 25, 30, 35)

# Definir los tipos de modelos
tipos_modelo <- c("Preciso", "Impreciso")

# Definir los tipos de casos
tipos_caso <- c("NVC" = "Normalidad y Varianza Constante", 
                "NVD" = "Normalidad y Varianza No Constante", 
                "NNVC" = "No Normalidad y Varianza Constante", 
                "NNVD" = "No Normalidad y Varianza No Constante")

# Función para seleccionar un caso al azar
seleccionar_caso <- function(tamaño) {
  if (tamaño %in% tamaños_muestra) {
    modelo_seleccionado <- sample(tipos_modelo, 1)
    caso_seleccionado <- sample(names(tipos_caso), 1)
    return(list(tamaño_muestra = tamaño, 
                modelo_seleccionado = modelo_seleccionado, 
                caso_seleccionado = caso_seleccionado, 
                descripcion_caso = tipos_caso[caso_seleccionado]))
  } else {
    stop("Tamaño de muestra no válido.")
  }
}

# Ejemplo de uso: seleccionar un caso para tamaño 25
resultado <- seleccionar_caso(35)

# Imprimir el resultado
print(resultado)

