install.packages("dbscan")
library(dbscan)

data("iris")
x <- as.matrix(iris[, 1:4])
x <- as.matrix(iris[, 3])

db <- dbscan(x, eps = .1, minPts = 2)
db

cbind(x, db$cluster)[order(x), ]

pairs(x, col = db$cluster + 1L)

lof <- lof(x, k = 4)
pairs(x, cex = lof)

opt <- optics(x, eps = 1, minPts = 4)
opt

