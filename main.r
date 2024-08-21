# Definindo as matrizes M1, M2, MA e MB
M1 <- matrix(c(2, 0, 0, -5, -6, 3, 5, 4, -2), nrow = 3, byrow = TRUE)
M2 <- matrix(c(2, 4, 0, 6, -5, 0, 5, 7, -4), nrow = 3, byrow = TRUE)
A <- matrix(c(2, 6, 6, -5, -6, 3, 5, 4, -2), nrow = 3, byrow = TRUE)
B <- matrix(c(6, 8, 0, -8, -6, 3, 5, 7, -9), nrow = 3, byrow = TRUE)
C <- matrix(c(2, 1, -3, 0, 2, 1, 5, 1, 3), nrow = 3, byrow = TRUE)
# Matrizes auxiliares
m2 <- matrix(nrow = 3, ncol = 3)
m3 <- matrix(nrow = 3, ncol = 3)
ma <- matrix(nrow = 3, ncol = 3)
mb <- matrix(nrow = 3, ncol = 3)

aux <- numeric(length = nrow(M1) * ncol(M1))
aux2 <- numeric(length = nrow(M2) * ncol(M2))
auxa <- numeric(length = nrow(A) * ncol(A))
auxb <- numeric(length = nrow(B) * ncol(B))

index <- 1 

# Expansão em cofatores para M1
for (i in 1:nrow(M1)) {
  for (j in 1:ncol(M1)) {
    m2[i, j] <- M1[i, j]
    aux <- m2[2, 2] * M1[3, 3]
    AUX <- m2[3, 2] * M1[2, 3] 
    c11 <- aux - AUX
    aux <- m2[2, 1] * M1[3, 3]
    AUX <- m2[3, 1] * M1[2, 3]
    c12 <- AUX - aux
    aux <- m2[2, 1] * M1[3, 2]
    AUX <- m2[3, 1] * M1[2, 2]
    c13 <- aux - AUX
  }
}

# Expansão em cofatores para M2
for (i in 1:nrow(M2)) {
  for (j in 1:ncol(M2)) {
    m3[i, j] <- M2[i, j]
    aux2 <- m3[2, 2] * M2[3, 3]
    AUX2 <- m3[3, 2] * M2[2, 3] 
    cl11 <- aux2 - AUX2
    aux2 <- m3[2, 1] * M2[3, 3]
    AUX2 <- m3[3, 1] * M2[2, 3]
    cl12 <- AUX2 - aux2
    aux2 <- m3[2, 1] * M2[3, 2]
    AUX2 <- m3[3, 1] * M2[2, 2]
    cl13 <- aux2 - AUX2
  }
}

# Expansão em cofatores para MA
for (i in 1:nrow(A)) {
  for (j in 1:ncol(A)) {
    ma[i, j] <- A[i, j]
    auxa <- ma[2, 2] * A[3, 3]
    AUXA <- ma[3, 2] * A[2, 3] 
    ca11 <- auxa - AUXA
    auxa <- ma[2, 1] * A[3, 3]
    AUXA <- ma[3, 1] * A[2, 3]
    ca12 <- AUXA - auxa
    auxa <- ma[2, 1] * A[3, 2]
    AUXA <- ma[3, 1] * A[2, 2]
    ca13 <- auxa - AUXA
  }
}

# Expansão em cofatores para MB
for (i in 1:nrow(B)) {
  for (j in 1:ncol(B)) {
    mb[i, j] <- B[i, j]
    auxb <- mb[2, 2] * B[3, 3]
    AUXB <- mb[3, 2] * B[2, 3] 
    cb11 <- auxb - AUXB
    auxb <- mb[2, 1] * B[3, 3]
    AUXB <- mb[3, 1] * B[2, 3]
    cb12 <- AUXB - auxb
    auxb <- mb[2, 1] * B[3, 2]
    AUXB <- mb[3, 1] * B[2, 2]
    cb13 <- auxb - AUXB
  }
}

# Soma de A e B (MA + MB)
Maux <- A + B
det_A <- det(A)
det_B <- det(B)
resultado_a <- det_A + det_B
soma_AB <- A + B
det_soma_AB <- det(soma_AB)
produto_AB <- A %*% B
det_produto_AB <- det(produto_AB)
cofactorMatrix <- function(A) {
  n <- nrow(A)
  cof <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      minor <- A[-i, -j]
      cof[i, j] <- (-1)^(i+j) * det(minor)
    }
  }
  return(cof)
}
adjoint <- function(A) {
  return(t(cofactorMatrix(A)))
}

# d) Adj(A)
adj_A <- adjoint(A)
print("Adj(A):")
print(adj_A)

# e) Adj(B)
adj_B <- adjoint(B)
print("Adj(B):")
print(adj_B)

# f) det(A^-1)
det_inv_A <- 1 / det_A
print(paste("det(A^-1) =", det_inv_A))

# g) det(B^-1)
det_inv_B <- 1 / det_B
print(paste("det(B^-1) =", det_inv_B))

det_C <- det(C)
print(paste("det(C) =", det_C))

# b) Calculando a inversa de C
inv_C <- solve(C)
# Resultados
print('============Quest 1===============')
print('Dada a matriz M1:') 
print(M1)
print('Expansão em cofatores ao longo da primeira linha da matriz M1:')
print(c11)
print(c12)
print(c13) 

print('============Quest 2===============')
print('Dada a matriz M2:') 
print(M2)
print('Expansão em cofatores ao longo da primeira linha da matriz M2:')
print(cl11)
print(cl12)
print(cl13)

print('============Quest 3===============')
print('Dadas as matrizes A e B, calcule A+B:')
print('Matriz A:')
print(A)
print('Matriz B:')
print(B)
print('Resposta:')
print(Maux)
print(paste("det(A) + det(B) =", resultado_a))
print(paste("det(A + B) =", det_soma_AB))
print(paste("det(AB) =", det_produto_AB))
print('============Quest 4===============')
print("inv(C):")
print(inv_C)


