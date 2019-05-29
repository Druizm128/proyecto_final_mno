
##### LEER DATOS #####


#muestra_netflix <- readRDS("datos/muestra_pelis_todos_usuarios")

library(Matrix)
library(tidyverse)

functionALS <- function(datos, lambda = 0.03, n_factores = 5, n_iter = 100, delta = 0.0001){
  
  start_time <- Sys.time()
  # Convertimos los datos a una matrix de formato: usarios X peliculas = R
  R <- sparseMatrix(
    as.integer(datos$nuevos_usuarios_ids), # renglón
    as.integer(datos$nuevos_ids),    # columna 
    x = datos$calif)              # rating
  paste("\nSe convirtiron los datos a una matriz rala")
  
  
  
  # Obtenemos las dimensiones de la matriz R
  n_u <- dim(R)[1]
  n_m <- dim(R)[2]
  paste("Tenemos una matriz rala de dimensiones ", n_u, " x ", n_m, sep="")
  
  
  ##### OBTENEMOS LAS DIMENSIONES DE LOS DATOS #####
  
  # Obtenemos matriz de pesos W
  
  # La matriz W tiene 0 (NA) cuando el usuario no ha calificado la película, y tiene 1 para las películas que sí calificó (a esta matriz le llaman I en el paper)
  W <- R
  W@x[!is.na(W@x)] <- 1
  
  ##### DEFINIMOS LA FUNCION DE COSTOS #####
  
  # Función de costo
  fn_costo <- function (R, U, M, W, lambda, n_u_i, n_m_j) {
    sum(W * ((R - (U %*% M)) ^ 2)) + lambda * (sum(n_u_i %*% (U ^ 2)) + sum((M ^2) %*% n_m_j))
  }
  
  
  ##### FIJAR HIPERPARÁMETROS #####
  
  # Semilla
  set.seed(28882)
  
  # Regularización
  #lambda = 0.03
  
  # Número de factores latentes
  #n_factores = 5
  
  # Número de iteraciones
  #n_iter = 1000
  
  
  ##### INICIALIZAMOS EL ALGORITMO #####
  
  # Factor para hacer chicos los números.
  delta <- 0.0001
  
  # Inicializamos la matriz de películas y factores latentes
  M <- delta * matrix(runif(n_m * n_factores),nrow = n_factores, ncol = n_m)
  
  # Como se sugiere en el artículo, se inicializará la primera fila de la matriz M con los promedios de calificaciones de las películas. Si hay películas SIN calificación, se les pone la calificación promedio de todas las peliculas 
  M[1, ] <- colSums(R, na.rm = TRUE) / colSums(W)
  calif_promedio <- mean(M[1, ], na.rm = TRUE)
  M[1, ][is.na(M[1, ])] <- calif_promedio
  
  # Inicializamos la matriz de Us con 1.
  U <- matrix(0, nrow = n_u, ncol = n_factores)
  U[, 1] <- 1
  
  # Número de ratings para cada usuario
  n_u_i <- rowSums(W)
  
  # Número de usuarios para cada pelicula
  n_m_j <- colSums(W)
  
  
  ###### EJECUCION DEL ALGORITMO #######
  
  ## vector de kk
  iteracion <- vector(mode = "numeric", length = n_iter)
  costo1 <- vector(mode = "numeric", length = n_iter)
  costo2 <- vector(mode = "numeric", length = n_iter)
  
  
  # Este for aproxima U%*%M a R, considerando la matriz de pesos W para calcular el error
  for (kk in 1:n_iter) {
    # El primer paso es fijar M y minimizar la fn de costo. Se itera sobre los usuarios 
    
    for (ii in 1:n_u) {
      # Primero hay que desechar las columnas M y filas R que son irrelevantes para 
      #el usuario ii, para poder incrementar la velocidaddel algorítmo
      
      M_selected <- M[, W[ii,] == 1, drop = FALSE]
      R_selected <- R[ii,][W[ii,] == 1, drop = FALSE]
      
      A<- M_selected %*% t(M_selected) + lambda * n_u_i[ii] * diag(n_factores)
      V <-  (M_selected %*% R_selected)
      
      
      # Actualiza U para cada usuario ii
      U[ii,] <- t(solve(A)%*%V)
      
      
    }
    cost <- fn_costo(R, U, M, W, lambda = lambda, n_u_i, n_m_j)
    print(paste0(kk, " iteración, paso 1: Función de costo = ", cost))
    costo1[kk] <- cost
    
    # Minimizar U %*% M para una U fija iterando sobre n usuarios (items)
    for (jj in 1:n_m) {
      if (sum(W[, jj] == 1) > 0) {
        U_selected <- U[W[, jj] == 1, , drop = FALSE]
        R_col_selected <- R[, jj][W[, jj] == 1, drop = FALSE]
        
        A <- t(U_selected) %*% U_selected + lambda * n_m_j[jj] * diag(n_factores)
        V <- t(U_selected) %*% R_col_selected
        # Actualiza M para cada película jj
        M[, jj] <- solve(A)%*%V
      }
    }
    cost <- fn_costo(R, U, M, W, lambda = lambda, n_u_i, n_m_j)
    print(paste0(kk, " iteración, paso 2: Función de costo = ", cost))
    costo2[kk] <- cost
    
    # Guardar iteración
    iteracion[kk] <- kk
    
  }
  end_time <- Sys.time()

  ###### RESULTADO #######
  monitoreo <- data_frame(iteracion, costo1,costo2)
  
  # Aproximamos la matriz de calificaciones R realizando el producto de matrices U * M
  ratings <-as(U %*% M, "dgCMatrix")

  #round(ratings)
  
  # Error de la aproximación RMSE
  
  recm <- function(calif, pred){
    sqrt(mean((calif - pred)^2))
  }
  rmse <- recm(R,ratings)
  
  # Tiempo de ejecuciopn
  tiempo <- end_time-start_time
  
  resultados <- list(monitoreo = monitoreo, 
                     rmse = rmse, 
                     ratings = round(ratings), 
                     tiempo = tiempo, 
                     n_iter = n_iter, 
                     lambda = lambda, 
                     n_factores = n_factores, 
                     delta = delta,
                     matriz_resultado = ratings)
  
  return(resultados)
  
  
}








