library(phylogram)
context("distance computation and embedding")

# simulate a sequence set
set.seed(999)
bases <- c("a", "c", "g", "t")
x <- list(sample(bases, replace = TRUE, size = 100))
evolve <- function(a) if(runif(1) > 0.95) sample(bases, 1) else a
for(i in 2:10) x[[i]] <- unname(sapply(x[[i - 1]], evolve))
names(x) <- paste0("S", 1:10)

# create DNAbin object
rawbases <- as.raw(c(136, 40, 72, 24))
xDNA <- lapply(x, function(s) rawbases[match(s, bases)])
class(xDNA) <- "DNAbin"

# generate a k-mer distance matrix
x.dist <- kdistance(x, method = "edgar")
xDNA.dist <- kdistance(xDNA, method = "edgar")

# embed sequences
x.mbed <- mbed(x, counts = TRUE)
xDNA.mbed <- mbed(x, counts = TRUE)

# build a divisive tree
set.seed(999)
x.tree <- topdown(x)
set.seed(999)
xDNA.tree <- topdown(xDNA)

test_that("objects have correct classes", {
  expect_is(x.dist, "dist")
  expect_is(x.mbed, "mbed")
  expect_is(x.tree, "dendrogram")
})

test_that("character and DNAbin formats give same results", {
  expect_equal(x.dist, xDNA.dist)
  expect_equal(x.mbed, xDNA.mbed)
  expect_equal(x.tree, xDNA.tree)
})
