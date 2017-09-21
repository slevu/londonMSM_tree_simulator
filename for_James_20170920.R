## Say I have a matrix like below. It represents tips of a tree in row/col, and values = node corresponding to most recent common ancestor (MRCA) for each pair of tips.
m <- structure(c(1, 5, 5, 5, 5, 2, 7, 6, 5, 7, 3, 6, 5, 6, 6, 4), .Dim = c(4L, 4L), .Dimnames = list(c("t2", "t3", "t4", "t1"), c("t2", "t3", "t4", "t1")))
m

## And a vector of heights of each tips/node, in same order as matrix
heights <- c(1.9, 0.5, 0, 1.1, 2.4, 1.5, 1)

## And a threshold of distance between node and tips
k <- 1.2

## I want a new matrix that gives 1 if two daughter tips are both less than the threshold away from their MRCA, and 0 otherwise (and in diagonal). So I do this:
ml <- m; ml[] <-  0
for (i in 1:dim(m)[1]) for (j in 1:dim(m)[2]){
  ml[i,j] <- ifelse( i!=j && (heights[m[i,j]] - heights[i]) < k &&
                       (heights[m[j,i]] - heights[j]) < k, 1, 0)
}
ml

## Which should give me this matrix:
to_be_found <- structure(c(0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(c("t2", "t3", "t4", "t1"), c("t2", "t3", "t4", "t1")))
