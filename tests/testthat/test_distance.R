library(phylogram)
context("distance computation and tree building")

# simulate a DNA sequence dataset
set.seed(999)
bases <- c("A", "C", "G", "T")
x <- list(sample(bases, replace = TRUE, size = 100))
evolve <- function(a) if(runif(1) > 0.95) sample(bases, 1) else a
for(i in 2:10) x[[i]] <- unname(sapply(x[[i - 1]], evolve))
names(x) <- paste("Sequence", 1:10)
# convert to DNAbin object
rawbases <- as.raw(c(136, 40, 72, 24))
xDNA <- lapply(x, function(s) rawbases[match(s, bases)])
class(xDNA) <- "DNAbin"

# simulate an AA sequence dataset this time an alignment
set.seed(999)
aminos <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
y <- matrix(sample(aminos, replace = TRUE, size = 100), nrow = 1)
evolve <- function(a) if(runif(1) > 0.95) sample(aminos, 1) else a
for(i in 2:10) y <- rbind(y, sapply(y[i - 1,], evolve))
rownames(y) <- paste("Sequence", 1:10)
# convert to AAbin object
rawaminos <- as.raw((65:89)[-c(2, 10, 15, 21, 24, 26)])
yAA <- apply(y, c(1, 2), function(s) rawaminos[match(s, aminos)])
class(yAA) <- "AAbin"

# count k-mers
x.kcounts <- kcount(x)
xDNA.kcounts <- kcount(xDNA)
y.kcounts <- kcount(y, k = 2)
yAA.kcounts <- kcount(yAA, k = 2)

# generate k-mer distance matrices
x.dist <- kdistance(x, method = "edgar")
xDNA.dist <- kdistance(xDNA, method = "edgar")
y.dist <- kdistance(y, method = "edgar", k = 2)
yAA.dist <- kdistance(yAA, method = "edgar", k = 2)

# embed sequences
set.seed(999)
x.mbed <- mbed(x)
set.seed(999)
xDNA.mbed <- mbed(xDNA)
set.seed(999)
y.mbed <- mbed(y, k = 2)
set.seed(999)
yAA.mbed <- mbed(y, k = 2)

# build divisive trees
set.seed(999)
x.tree <- topdown(x, nstart = 20)
set.seed(999)
xDNA.tree <- topdown(xDNA, nstart = 20)
set.seed(999)
y.tree <- topdown(y, nstart = 20, k = 2)
set.seed(999)
yAA.tree <- topdown(yAA, nstart = 20, k = 2)

# prune trees
x.tree <- prune(x.tree, pattern = "^Sequence 10$")
xDNA.tree <- prune(xDNA.tree, pattern = "^Sequence 10$")
y.tree <- prune(y.tree, pattern = "^Sequence 10$")
yAA.tree <- prune(yAA.tree, pattern = "^Sequence 10$")

# arrange trees in ladder form
x.tree <- ladder(x.tree, decreasing = FALSE)
xDNA.tree <- ladder(xDNA.tree, decreasing = FALSE)
y.tree <- ladder(y.tree, decreasing = TRUE)
yAA.tree <- ladder(yAA.tree, decreasing = TRUE)

# reposition trees
x.tree <- reposition(x.tree, shift = 0.01)
xDNA.tree <- reposition(xDNA.tree, shift = 0.01)
y.tree <- reposition(y.tree, shift = 0.01)
yAA.tree <- reposition(yAA.tree, shift = 0.01)

test_that("objects have correct classes", {
  expect_is(x.dist, "dist")
  expect_is(y.dist, "dist")
  expect_is(x.mbed, "mbed")
  expect_is(y.mbed, "mbed")
  expect_is(x.tree, "dendrogram")
  expect_is(y.tree, "dendrogram")
})

test_that("character and DNAbin formats give same results", {
  expect_equal(x.kcounts, xDNA.kcounts)
  expect_equal(y.kcounts, yAA.kcounts)
  expect_equal(x.dist, xDNA.dist)
  expect_equal(y.dist, yAA.dist)
  expect_equal(x.mbed[,], xDNA.mbed[,])
  expect_equal(y.mbed[,], yAA.mbed[,])
  expect_equal(x.tree, xDNA.tree)
  expect_equal(y.tree, yAA.tree)
})

test_that("object dimensions are correct", {
  expect_equal(ncol(x.kcounts), 1024)
  expect_equal(ncol(y.kcounts), 400)
  expect_equal(length(x.dist), length(y.dist))
  expect_equal(ncol(x.mbed), 9)
  expect_equal(ncol(y.mbed), 9)
  expect_equal(length(x.tree), 2)
  expect_equal(length(y.tree), 2)
})
