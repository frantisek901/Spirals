#### Script for coumputing Hausdorf fractal dimension of opinion space



# create a dataset to calculate de Hausdorff-Besicovitch dimension
mat <- matrix(runif(512*512),nrow = 512,ncol = 512)

mat[mat<=0.5] <- 0
mat[mat>0.5] <- 1

cant <- sum(mat)

fragment <- rep(2,10)**(0:9)
Table <- data.frame(Delta = rep(512,10)/(fragment ), N = fragment**2)
Table$LogDelta <- log(Table$Delta)

for(i in 2:10){
  delta_aux <- Table$Delta[i]

  for(j in 1:fragment [i]){
    row_id <- ((j-1)*delta_aux+1):(j*delta_aux)
    for(k in 1:fragment [i]){
      col_id <- ((k-1)*delta_aux+1):(k*delta_aux)
      if(sum(mat[row_id,col_id]) == 0){
        Table$N[i] <- Table$N[i] - 1
      }
    }
  }
}

Table$LogN <- log(Table$N)
lm_dim <- lm(Table$LogN ~ Table$LogDelta)

plot(Table$LogN ~ Table$LogDelta)
abline(lm_dim)

print('The box-counting dimension is:')
print(-lm_dim$coefficients[2])

# without the borders
Table <- Table[2:nrow(Table),]
lm_dim <- lm(Table$LogN ~ Table$LogDelta)

plot(Table$LogN ~ Table$LogDelta)
abline(lm_dim)

print('The box-counting dimension is:')
print(-lm_dim$coefficients[2])