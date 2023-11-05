# Librerías
library(fields)
library(MASS)
library(imager)
library(writexl)


######### Funciones para tratar imágenes #########

# Abre una imagen, la convierte a escala de grises y devuelve una matriz.
abrirImagen <- function(path){
  img <- load.image(path)
  img <- grayscale(img)
  img_matrix <- as.array(img)[,,1,1]
  return(img_matrix)
}

# Rotación de una imagen cuadrada.
rotarImagen <- function(img, rot){
  l <- dim(img)[1]
  im <- as.cimg(img)
  im_rot <- imrotate(im, rot)
  img_rot <- as.array(im_rot)[,,1,1]
  l_rot <- dim(img_rot)[1]
  w_rot <- dim(img_rot)[2]
  l_dif <- l_rot - l
  w_dif <- w_rot - l
  img_rot <- img_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
  return(img_rot)
}

# Reescalamiento de una imagen cuadrada a una matriz de lxl
escalarImagen <- function(img, l){
  img <- as.cimg(img)
  img <- resize(img, l, l)
  img_matrix <- as.array(img)[,,1,1]
  return(img_matrix)
}

######### Estimadores del Variograma Angular, Variograma Angular Cruzado y Pseudo Variograma Angular Cruzado #########

# Estimador del Variograma Angular
VarAng <- function(img){
  var_ang <- numeric(180)
  l <- dim(img)[1]
  im <- as.cimg(img)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:180){
    im_rot <- imrotate(im, i)
    img_rot <- as.array(im_rot)[,,1,1]
    im_ones_rot <- imrotate(im_ones, i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot <- img_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang[i] <- sum((img[ones_rot == 1]-img_rot[ones_rot == 1])^2)/n
  }
  return(var_ang)
}

# Estimador del Variograma Angular Cruzado
VarAngCr <- function(img1, img2){
  var_ang_cr <- numeric(180)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:180){
    im_rot1 <- imrotate(im1, i)
    im_rot2 <- imrotate(im2, i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang_cr[i] <- sum((img1[ones_rot == 1]-img_rot1[ones_rot == 1])*(img2[ones_rot == 1]-img_rot2[ones_rot == 1]))/n
  }
  return(var_ang_cr)
}

# Estimador del Pseudo Variograma Angular Cruzado
PseudoVarAngCr <- function(img1, img2){
  pseudo_var_ang_cr <- numeric(180)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:180){
    im_rot1 <- imrotate(im1, i)
    im_rot2 <- imrotate(im2, i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    pseudo_var_ang_cr[i] <- sum((img_rot1[ones_rot == 1]-img2[ones_rot == 1])^2)/n
  }
  return(pseudo_var_ang_cr)
}


######### Estimadores considerando particiones de 2 en 2 #########

VarAng2 <- function(img){
  var_ang <- numeric(90)
  l <- dim(img)[1]
  im <- as.cimg(img)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:90){
    im_rot <- imrotate(im, 2*i)
    img_rot <- as.array(im_rot)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 2*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot <- img_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang[i] <- sum((img[ones_rot == 1]-img_rot[ones_rot == 1])^2)/n
  }
  return(var_ang)
}

VarAngCr2 <- function(img1, img2){
  var_ang_cr <- numeric(90)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:90){
    im_rot1 <- imrotate(im1, 2*i)
    im_rot2 <- imrotate(im2, 2*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 2*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang_cr[i] <- sum((img1[ones_rot == 1]-img_rot1[ones_rot == 1])*(img2[ones_rot == 1]-img_rot2[ones_rot == 1]))/n
  }
  return(var_ang_cr)
}

PseudoVarAngCr2 <- function(img1, img2){
  pseudo_var_ang_cr <- numeric(90)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:90){
    im_rot1 <- imrotate(im1, 2*i)
    im_rot2 <- imrotate(im2, 2*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 2*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    pseudo_var_ang_cr[i] <- sum((img_rot1[ones_rot == 1]-img2[ones_rot == 1])^2)/n
  }
  return(pseudo_var_ang_cr)
}

######### Estimadores considerando particiones de 3 en 3 #########

VarAng3 <- function(img){
  var_ang <- numeric(60)
  l <- dim(img)[1]
  im <- as.cimg(img)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:60){
    im_rot <- imrotate(im, 3*i)
    img_rot <- as.array(im_rot)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 3*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot <- img_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang[i] <- sum((img[ones_rot == 1]-img_rot[ones_rot == 1])^2)/n
  }
  return(var_ang)
}

VarAngCr3 <- function(img1, img2){
  var_ang_cr <- numeric(60)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:60){
    im_rot1 <- imrotate(im1, 3*i)
    im_rot2 <- imrotate(im2, 3*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 3*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang_cr[i] <- sum((img1[ones_rot == 1]-img_rot1[ones_rot == 1])*(img2[ones_rot == 1]-img_rot2[ones_rot == 1]))/n
  }
  return(var_ang_cr)
}

PseudoVarAngCr3 <- function(img1, img2){
  pseudo_var_ang_cr <- numeric(60)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:60){
    im_rot1 <- imrotate(im1, 3*i)
    im_rot2 <- imrotate(im2, 3*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 3*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    pseudo_var_ang_cr[i] <- sum((img_rot1[ones_rot == 1]-img2[ones_rot == 1])^2)/n
  }
  return(pseudo_var_ang_cr)
}

######### Estimadores considerando particiones de 5 en 5 #########

VarAng5 <- function(img){
  var_ang <- numeric(36)
  l <- dim(img)[1]
  im <- as.cimg(img)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:36){
    im_rot <- imrotate(im, 5*i)
    img_rot <- as.array(im_rot)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 5*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot <- img_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang[i] <- sum((img[ones_rot == 1]-img_rot[ones_rot == 1])^2)/n
  }
  return(var_ang)
}

VarAngCr5 <- function(img1, img2){
  var_ang_cr <- numeric(36)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:36){
    im_rot1 <- imrotate(im1, 5*i)
    im_rot2 <- imrotate(im2, 5*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 5*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    var_ang_cr[i] <- sum((img1[ones_rot == 1]-img_rot1[ones_rot == 1])*(img2[ones_rot == 1]-img_rot2[ones_rot == 1]))/n
  }
  return(var_ang_cr)
}

PseudoVarAngCr5 <- function(img1, img2){
  pseudo_var_ang_cr <- numeric(36)
  l <- dim(img1)[1]
  im1 <- as.cimg(img1)
  im2 <- as.cimg(img2)
  ones <- matrix(1, l, l)
  im_ones <- as.cimg(ones)
  for(i in 1:36){
    im_rot1 <- imrotate(im1, 5*i)
    im_rot2 <- imrotate(im2, 5*i)
    img_rot1 <- as.array(im_rot1)[,,1,1]
    img_rot2 <- as.array(im_rot2)[,,1,1]
    im_ones_rot <- imrotate(im_ones, 5*i)
    ones_rot <- as.array(im_ones_rot)[,,1,1]
    l_rot <- dim(ones_rot)[1]
    w_rot <- dim(ones_rot)[2]
    l_dif <- l_rot - l
    w_dif <- w_rot - l
    img_rot1 <- img_rot1[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    img_rot2 <- img_rot2[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    ones_rot <- ones_rot[(l_dif%/%2+1):(l_rot-(l_dif%/%2+l_dif%%2)),(w_dif%/%2+1):(w_rot-(w_dif%/%2+w_dif%%2))]
    n <- sum(ones_rot == 1)
    pseudo_var_ang_cr[i] <- sum((img_rot1[ones_rot == 1]-img2[ones_rot == 1])^2)/n
  }
  return(pseudo_var_ang_cr)
}

######### Funciones de Covarianza #########

cov_function_exp <- function(h, a=1, scale=1) {
  return(scale * exp(-a * h))
}

cov_function_spherical <- function(h, a=1) {
  cov_values <- ifelse(h <= a, 
                       1 - 1.5 * (h / a) + 0.5 * (h / a)^3,
                       0)
  return(cov_values)
}

cov_function_gaussian <- function(h, a=1) {
  exp(- (h^2) / (2 * a^2))
}


######### SIMULACIÓN DE IMÁGENES CON COVARIANZA EXPONENCIAL #########

# Parámetros de la grilla
nx <- 50  # Número de píxeles a lo largo del eje x
ny <- 50  # Número de píxeles a lo largo del eje y

# Crear una malla
x <- seq(0, 1, length = nx)
y <- seq(0, 1, length = ny)
grid <- expand.grid(x = x, y = y)

# Calcular la matriz de covarianza
distances <- as.matrix(dist(grid))
cov_matrix_exp <- cov_function_exp(distances)

# Iniciar un dataframe para guardar todas las imágenes
all_images1 <- c()
startTime <- Sys.time()
for(i in 1:100){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images1 <- cbind(all_images1, sim_img)
}
write.csv(all_images1, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img1.csv", row.names = FALSE)

all_images2 <- c()
for(i in 101:200){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images2 <- cbind(all_images2, sim_img)
}
write.csv(all_images2, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img2.csv", row.names = FALSE)

all_images3 <- c()
for(i in 201:300){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images3 <- cbind(all_images3, sim_img)
}
write.csv(all_images3, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img3.csv", row.names = FALSE)

all_images4 <- c()
for(i in 301:400){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images4 <- cbind(all_images4, sim_img)
}
write.csv(all_images4, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img4.csv", row.names = FALSE)

all_images5 <- c()
for(i in 401:500){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images5 <- cbind(all_images5, sim_img)
}
write.csv(all_images5, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img5.csv", row.names = FALSE)

all_images6 <- c()
for(i in 501:600){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images6 <- cbind(all_images6, sim_img)
}
write.csv(all_images6, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img6.csv", row.names = FALSE)

all_images7 <- c()
for(i in 601:700){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images7 <- cbind(all_images7, sim_img)
}
write.csv(all_images7, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img7.csv", row.names = FALSE)

all_images8 <- c()
for(i in 701:800){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images8 <- cbind(all_images8, sim_img)
}
write.csv(all_images8, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img8.csv", row.names = FALSE)

all_images9 <- c()
for(i in 801:900){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images9 <- cbind(all_images9, sim_img)
}
write.csv(all_images9, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img9.csv", row.names = FALSE)

all_images10 <- c()
for(i in 901:1000){
  cat("Imagen número:",i)
  set.seed(i)
  sim_img <- mvrnorm(n = 1, mu = rep(0, nx * ny), Sigma = cov_matrix_exp)
  all_images10 <- cbind(all_images10, sim_img)
}
write.csv(all_images10, "C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img10.csv", row.names = FALSE)

endTime <- Sys.time()
print(endTime - startTime)


# Lectura de los archivos con las imágenes
for(i in 1:10){
  eval(parse(text = paste0("sim_img_exp_df",i," <- read.csv('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/Geodésica/sim_img",i,".csv', header = TRUE)","")))
}
sim_img_exp_df <- cbind(sim_img_exp_df1, sim_img_exp_df2, sim_img_exp_df3, sim_img_exp_df4, sim_img_exp_df5, sim_img_exp_df6, sim_img_exp_df7, sim_img_exp_df8, sim_img_exp_df9, sim_img_exp_df10)

image_exp_list <- list()

# Convertir cada columna en una matriz 2D y almacenar en la lista
for (i in 1:ncol(sim_img_exp_df)) {
  image_vector <- as.vector(sim_img_exp_df[, i])
  image_matrix <- matrix(image_vector, nrow = nx)
  image_exp_list[[i]] <- image_matrix
}

# Resultados del Variograma Angular Cruzado para los ángulo 3, 30 y 45
ang_min_3 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- VarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 3))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  ang_min_3[i] <- rbind(as.integer(ind))
}
mean(ang_min_3)
sd(ang_min_3)
summary(ang_min_3)

ang_min_30 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- VarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 30))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  ang_min_30[i] <- rbind(as.integer(ind))
}
mean(ang_min_30)
median(ang_min_30)
sd(ang_min_30)

ang_min_45 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- VarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 45))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  ang_min_45[i] <- rbind(as.integer(ind))
}
mean(ang_min_45)
median(ang_min_45)
sd(ang_min_45)


# Resultados del Pseudo Variograma Angular Cruzado para los ángulo 3, 30 y 45
pseudo_ang_min_3 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- PseudoVarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 3))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  pseudo_ang_min_3[i] <- rbind(as.integer(ind))
}
mean(pseudo_ang_min_3)
median(pseudo_ang_min_3)
sd(pseudo_ang_min_3)

pseudo_ang_min_30 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- PseudoVarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 30))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  pseudo_ang_min_30[i] <- rbind(as.integer(ind))
}
mean(pseudo_ang_min_30)
median(pseudo_ang_min_30)
sd(pseudo_ang_min_30)

pseudo_ang_min_45 <- c()
for(i in 1:1000){
  cat("Imagen número:",i,"\n")
  vac <- PseudoVarAngCr(image_exp_list[[i]], rotarImagen(image_exp_list[[i]], 45))
  ind <- which(sapply(vac, function(x) x == min(vac)))
  pseudo_ang_min_45[i] <- rbind(as.integer(ind))
}
mean(pseudo_ang_min_45)
median(pseudo_ang_min_45)
sd(pseudo_ang_min_45)


######### ANÁLISIS DE IMÁGENES CON DISTORSIONES (TID2013) #########

# Se abren las imágenes de referencia y con distorsión, y se reescalan a 300x300.
for(i in 1:25){
  cat("Imagen número: ", i, " \n")
  if(i < 10){
    eval(parse(text = paste0("img", i, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/reference_images/I0", i, ".bmp'), 300)")))
    for(j in 1:24){
      if(j < 10){
        eval(parse(text = paste0("img", i, "_d_", j, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/distorted_images/i0", i, "_0", j, "_5.bmp'), 300)")))
      }
      else{
        eval(parse(text = paste0("img", i, "_d_", j, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/distorted_images/i0", i, "_", j, "_5.bmp'), 300)")))
      }
    }
  }
  else{
    eval(parse(text = paste0("img", i, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/reference_images/I", i, ".bmp'), 300)")))
    for(j in 1:24){
      if(j < 10){
        eval(parse(text = paste0("img", i, "_d_", j, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/distorted_images/i", i, "_0", j, "_5.bmp'), 300)")))
      }
      else{
        eval(parse(text = paste0("img", i, "_d_", j, " <- escalarImagen(abrirImagen('C:/Users/basti/OneDrive/Documentos/Google Drive Local/Mi Unidad/12vo Semestre/Memoria/tid2013/distorted_images/i", i, "_", j, "_5.bmp'), 300)")))
      }
    }
  }
}

# Se obtienen los Variogramas Angulares de las imágenes de referencia
for(i in 1:25){
  cat("Variograma Angular nro: ", i, " \n")
  eval(parse(text = paste0("var_ang_", i, " <- VarAng(img", i, ")")))
}

# Se obtienen los Variogramas Angulares Cruzado para los ángulos 3, 30 y 45
for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("var_ang_cr_", i, "_", j, "_3 <- VarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 3))")))
  }
}
var_ang_cr_3 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("var_ang_cr_3 <- cbind(var_ang_cr_3, var_ang_cr_", i, "_", j, "_3)")))
  }
}

for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("var_ang_cr_", i, "_", j, "_30 <- VarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 30))")))
  }
}
var_ang_cr_30 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("var_ang_cr_30 <- cbind(var_ang_cr_30, var_ang_cr_", i, "_", j, "_30)")))
  }
}

for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("var_ang_cr_", i, "_", j, "_45 <- VarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 45))")))
  }
}
var_ang_cr_45 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("var_ang_cr_45 <- cbind(var_ang_cr_45, var_ang_cr_", i, "_", j, "_45)")))
  }
}

# Se obtienen los Pseudo Variogramas Angulares Cruzado para los ángulos 3, 30 y 45
for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("pseudo_var_ang_cr_", i, "_", j, "_3 <- PseudoVarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 3))")))
  }
}
pseudo_var_ang_cr_3 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("pseudo_var_ang_cr_3 <- cbind(pseudo_var_ang_cr_3, pseudo_var_ang_cr_", i, "_", j, "_3)")))
  }
}

for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("pseudo_var_ang_cr_", i, "_", j, "_30 <- PseudoVarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 30))")))
  }
}
pseudo_var_ang_cr_30 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("pseudo_var_ang_cr_30 <- cbind(pseudo_var_ang_cr_30, pseudo_var_ang_cr_", i, "_", j, "_30)")))
  }
}

for(i in 1:25){
  cat("Imagen nro: ", i, " \n")
  for(j in 1:24){
    cat("   Distorsión nro: ", j, " \n")
    eval(parse(text = paste0("pseudo_var_ang_cr_", i, "_", j, "_45 <- PseudoVarAngCr(img", i, ", rotarImagen(img", i, "_d_", j, ", 45))")))
  }
}
pseudo_var_ang_cr_45 <- c()
for(i in 1:25){
  for(j in 1:24){
    eval(parse(text = paste0("pseudo_var_ang_cr_45 <- cbind(pseudo_var_ang_cr_45, pseudo_var_ang_cr_", i, "_", j, "_45)")))
  }
}

# Se obtienen las estimaciones de los ángulos de rotación 3, 30 y 45 para el Variograma Angular Cruzado
ang_min_3_vac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(var_ang_cr_3[,i], function(x) x == min(var_ang_cr_3[,i])))
  ang_min_3_vac[i] <- rbind(as.integer(ind))
}
mean(ang_min_3_vac)
# Se obtienen las estimaciones de los ángulos de rotación para el Variograma Angular Cruzado por imagen y por distorsión
ang_min_3_vac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_3_vac[(24*i-23):(24*i)]), "\n")
  ang_min_3_vac_x_img <- rbind(ang_min_3_vac_x_img, cbind(mean(ang_min_3_vac[(24*i-23):(24*i)]), sd(ang_min_3_vac[(24*i-23):(24*i)])))
}
ang_min_3_vac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_3_vac[(1:600)%%24==i]), "\n")
  ang_min_3_vac_x_dist <- rbind(ang_min_3_vac_x_dist, cbind(mean(ang_min_3_vac[(1:600)%%24==i]), sd(ang_min_3_vac[(1:600)%%24==i])))
}

ang_min_30_vac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(var_ang_cr_30[,i], function(x) x == min(var_ang_cr_30[,i])))
  ang_min_30_vac[i] <- rbind(as.integer(ind))
}
mean(ang_min_30_vac)
# Se obtienen las estimaciones de los ángulos de rotación para el Variograma Angular Cruzado por imagen y por distorsión
ang_min_30_vac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_30_vac[(24*i-23):(24*i)]), "\n")
  ang_min_30_vac_x_img <- rbind(ang_min_30_vac_x_img, cbind(mean(ang_min_30_vac[(24*i-23):(24*i)]), sd(ang_min_30_vac[(24*i-23):(24*i)])))
}
ang_min_30_vac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_30_vac[(1:600)%%24==i]), "\n")
  ang_min_30_vac_x_dist <- rbind(ang_min_30_vac_x_dist, cbind(mean(ang_min_30_vac[(1:600)%%24==i]), sd(ang_min_30_vac[(1:600)%%24==i])))
}

ang_min_45_vac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(var_ang_cr_45[,i], function(x) x == min(var_ang_cr_45[,i])))
  ang_min_45_vac[i] <- rbind(as.integer(ind))
}
mean(ang_min_45_vac)
# Se obtienen las estimaciones de los ángulos de rotación para el Variograma Angular Cruzado por imagen y por distorsión
ang_min_45_vac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_45_vac[(24*i-23):(24*i)]), "\n")
  ang_min_45_vac_x_img <- rbind(ang_min_45_vac_x_img, cbind(mean(ang_min_45_vac[(24*i-23):(24*i)]), sd(ang_min_45_vac[(24*i-23):(24*i)])))
}
ang_min_45_vac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_45_vac[(1:600)%%24==i]), "\n")
  ang_min_45_vac_x_dist <- rbind(ang_min_45_vac_x_dist, cbind(mean(ang_min_45_vac[(1:600)%%24==i]), sd(ang_min_45_vac[(1:600)%%24==i])))
}

# Se obtienen las estimaciones de los ángulos de rotación 3, 30 y 45 para el Variograma Angular Cruzado
ang_min_3_pvac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(pseudo_var_ang_cr_3[,i], function(x) x == min(pseudo_var_ang_cr_3[,i])))
  ang_min_3_pvac[i] <- rbind(as.integer(ind))
}
mean(ang_min_3_pvac)
# Se obtienen las estimaciones de los ángulos de rotación para el Pseudo Variograma Angular Cruzado por imagen y por distorsión
ang_min_3_pvac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_3_pvac[(24*i-23):(24*i)]), "\n")
  ang_min_3_pvac_x_img <- rbind(ang_min_3_pvac_x_img, cbind(mean(ang_min_3_pvac[(24*i-23):(24*i)]), sd(ang_min_3_pvac[(24*i-23):(24*i)])))
}
ang_min_3_pvac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_3_pvac[(1:600)%%24==i]), "\n")
  ang_min_3_pvac_x_dist <- rbind(ang_min_3_pvac_x_dist, cbind(mean(ang_min_3_pvac[(1:600)%%24==i]), sd(ang_min_3_pvac[(1:600)%%24==i])))
}

ang_min_30_pvac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(pseudo_var_ang_cr_30[,i], function(x) x == min(pseudo_var_ang_cr_30[,i])))
  ang_min_30_pvac[i] <- rbind(as.integer(ind))
}
mean(ang_min_30_pvac)
# Se obtienen las estimaciones de los ángulos de rotación para el Pseudo Variograma Angular Cruzado por imagen y por distorsión
ang_min_30_pvac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_30_pvac[(24*i-23):(24*i)]), "\n")
  ang_min_30_pvac_x_img <- rbind(ang_min_30_pvac_x_img, cbind(mean(ang_min_30_pvac[(24*i-23):(24*i)]), sd(ang_min_30_pvac[(24*i-23):(24*i)])))
}
ang_min_30_pvac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_30_pvac[(1:600)%%24==i]), "\n")
  ang_min_30_pvac_x_dist <- rbind(ang_min_30_pvac_x_dist, cbind(mean(ang_min_30_pvac[(1:600)%%24==i]), sd(ang_min_30_pvac[(1:600)%%24==i])))
}

ang_min_45_pvac <- c()
for(i in 1:600){
  cat("Imagen número:",i,"\n")
  ind <- which(sapply(pseudo_var_ang_cr_45[,i], function(x) x == min(pseudo_var_ang_cr_45[,i])))
  ang_min_45_pvac[i] <- rbind(as.integer(ind))
}
mean(ang_min_45_pvac)
# Se obtienen las estimaciones de los ángulos de rotación para el Pseudo Variograma Angular Cruzado por imagen y por distorsión
ang_min_45_pvac_x_img <- c()
for(i in 1:25){
  cat(length(ang_min_45_pvac[(24*i-23):(24*i)]), "\n")
  ang_min_45_pvac_x_img <- rbind(ang_min_45_pvac_x_img, cbind(mean(ang_min_45_pvac[(24*i-23):(24*i)]), sd(ang_min_45_pvac[(24*i-23):(24*i)])))
}
ang_min_45_pvac_x_dist <- c()
for(i in 0:23){
  cat(length(ang_min_45_pvac[(1:600)%%24==i]), "\n")
  ang_min_45_pvac_x_dist <- rbind(ang_min_45_pvac_x_dist, cbind(mean(ang_min_45_pvac[(1:600)%%24==i]), sd(ang_min_45_pvac[(1:600)%%24==i])))
}

# Se guardan los resultados obtenidos para cada ángulo y estimador.
ang_min_x_img <- cbind(ang_min_3_vac_x_img, ang_min_30_vac_x_img, ang_min_45_vac_x_img, ang_min_3_pvac_x_img, ang_min_30_pvac_x_img, ang_min_45_pvac_x_img)
ang_min_x_dist <- cbind(ang_min_3_vac_x_dist, ang_min_30_vac_x_dist, ang_min_45_vac_x_dist, ang_min_3_pvac_x_dist, ang_min_30_pvac_x_dist, ang_min_45_pvac_x_dist)
