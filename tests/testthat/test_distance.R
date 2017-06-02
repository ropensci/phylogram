library(phylogram)
context("distance computation and embedding")

# simulate a sequence set
set.seed(999)
bases <- c("a", "c", "g", "t")
x <- list(sample(bases, replace = TRUE, size = 100))
evolve <- function(a) if(runif(1) > 0.95) sample(bases, 1) else a
for(i in 2:10) x[[i]] <- unname(sapply(x[[i - 1]], evolve))

# generate a k-mer distance matrix
x.dist <- kdistance(x)

# embed sequences
x.mbed <- mbed(x, counts = TRUE)

# build a divisive tree
x.tree <- topdown(x.mbed)

test_that("objects have correct classes", {
  expect_is(x.dist, "dist")
  expect_is(x.mbed, "mbed")
  expect_is(x.tree, "dendrogram")
})

