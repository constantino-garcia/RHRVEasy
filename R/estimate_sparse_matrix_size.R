# estimate sparse matrix size per entry
library("Matrix")

matrix_sizes <- round(10 ^ seq(1, 4, by = 0.5))
ps <- c(0.01, 0.1, 0.5, 0.75)
results = sapply(matrix_sizes, \(matrix_size) {
  sapply(ps, \(p) {
    ss <- max(1, round(matrix_size * matrix_size * p))
    ind <- sample(matrix_size * matrix_size, ss, replace = FALSE)
    i <- floor((ind - 1) / matrix_size) + 1
    j <- ((ind - 1) %% matrix_size) + 1
    A <- sparseMatrix(i, j, dims = c(matrix_size, matrix_size))
    object.size(A) / length(ind)
  })
})
colnames(results) <- matrix_sizes
rownames(results) <- ps
print(results)
