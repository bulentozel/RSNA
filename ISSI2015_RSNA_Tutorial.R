#########################################
# ISSI 2015 - Introduction to SNA with R
# Bulent Ozel, ozel@uji.es 
# Universitat Jaume I, Castello, Spain         
#########################################


##### PART 1. A BRIEF INTRODUCTION to R ENVIRONMENT
#### 1.0 Initializing and Using R Packages:

## Installing (once only) packages ::
install.packages(c("sna", "network"))

## Using packages:
library(sna)
library(network)
# or:
require(sna)
require(network)


### 1.1. Introduction to basic R syntax
a <- 3
sqrt(a)
b <- a
a == b
a != sqrt(a)
a
b
ls()
help(sqrt)
`?`(sqrt)
help.start()
help.search("cent")
apropos("deg")
rm(a)
example(rgraph)
example(hist)

#### 1.3. Vectors and Matrices in R

# Creating vectors using the concatenation operator

x <- c(1, 3, 5)
y <- c("one", "three", "five")
c(x, y)
# entries must be of the same type, they are forced:
z <- c(T, F)
c(x, y, z)

# Sequences and replication
rep(1, 5)
rep(1:5, 2)
seq(from = 1, to = 5, by = 2)
1:5
# create a vector by concatenation
rep(1:5, each = 2)
rep(1:5, times = 5:1)

# Any and all (with vectors)
x <- 1:5
x > 2
any(x > 2)
all(x > 2)

is.element(4, x)

# From vectors to matrices
vaListues <- 1:24
vaListues
A <- matrix(vaListues, nrow = 6, ncol = 4, byrow = F)
A
A <- matrix(vaListues, nrow = 6, ncol = 4, byrow = T)
A

# Retrieving vaListues within a matrix:
A[1, 2]
A[1, ]
A[, 2]
A[2:3, 3:4]

# Filtering out rows and/or columns:
A[-1, ]
A[-2, -2]


# From vectors to matrices:

B <- cbind(1:5, 1:5)
B

C <- rbind(1:5, 1:5)
C

# Concatinating matrices:
cbind(B, t(C))

dim(B)
dim(C)
NROW(B)
NCOL(B)
cbind(B, B)
t(B)
cbind(t(B), C)


#### 1.4. Elementwise operations

# Most arithmatic operators are applied elementwise

x <- 1:5
x + 1
x * 2
x/3
x - 4
x^5

# if both operands are vector, then they must of same size:
x + x
x * x

A <- rbind(1:5, 2:6)
B <- rbind(3:7, 4:8)

A + B
A/B

# Inner product
x %*% x

A
B
A %*% t(B)


# LogicaList operators
A > 0
A == B
A != B
!(A == B)

(A > 2) | (B > 4)
(A > 2) || (B > 4)

(A > 2) & (B > 4)
(A > 2) && (B > 4)

# ... many other basic transformations
log(A)
exp(B)
sqrt(A + B)

# nested statements!
log((sqrt(A + B) + A) * B)

#### 1.5 Lists, data frames, and arrays

## Lists: Handling collection of an "irregular" set of data:::

# Creating lists:
aList <- list(1:5)
aList
aList <- list(1:5, letters[1:3])
aList


B <- matrix(1:3, nrow = 3, ncol = 2, byrow = T)

aList <- list(1:5, letters[1:3], B)
aList

# Retrieving data from lists:

aList[[1]]
aList[[2]][3]
aList[[3]][1, 2]

# Naming items within a list:

aList <- list(author.ids = 1:4, tot.pub = 7)
names(aList)


# multiple ways to access:
aList[["author.ids"]]
aList[[1]]
aList$author.ids

# elementwise operators on lists:
aList[[1]] + 3
aList[[2]] <- aList[[2]] * 4


# adding new items on the goo:
aList$description <- "Members of a research team."
aList[["gender.attr"]] <- c("F", "M", "M", "X")
aList



## Data Frames: Handling tabulated data sets::::
dF <- data.frame(income.level = 1:5, is.corrupt = c(F, F, F, F, T), name = LETTERS[1:5])
dF

# retrieve values:
dF[1, 2]
dF[2, 1:2]
dF[-1, ]
dF[, -2]
names(dF)
dF[2]
dF[[2]]


# Update entries:

dF$is.corrupt[5] <- FALSE
dF

dF[5, 2] <- "TRUE"
dF

# Add new column:

dF$education <- c("PhD", "BSc", "Pre", "Pre", "Pre")
dF

# converting matrices to data frames:

dF <- as.data.frame(cbind(1:5, rnorm(5)))
dF

is.data.frame(dF)
is.matrix(dF)


# changing column names in a data frame
colnames(dF) <- c("Sample", "Prob.")
dF
colnames(dF)[2] <- "Chance"
dF


## Arrays: When a use of 2 dimensional repres. is not enough::::
a <- array(1:18, dim = c(2, 3, 3))
a
a[1, 2, 3]
a[1, 2, ]
a[1, , ]
a[-1, 2:3, 1:2]
a * 5
a <- array(dim = c(2, 3, 2, 5, 6))
dim(a)

#### 1.6. Finding built-in data sets.
data()
data(flo)
# view the object 
class(flo)
flo

###### 1.7. Elementary visuaListizations
data(USArrests)
class(USArrests)
colnames(USArrests)

plot(USArrests$UrbanPop, USArrests$Murder)

# Adding labels:
plot(USArrests$UrbanPop, USArrests$Murder, xlab = "Population Size", ylab = "Murder Cases", main = "Urban Crime")


plot(USArrests$UrbanPop, USArrests$Murder, xlab = "Population Size", ylab = "Murder Cases", main = "Urban Crime", type = "n")

text(USArrests$UrbanPop, USArrests$Murder, rownames(USArrests))

# Checking frequency distributions
hist(USArrests$Murder)

# Comparing distibutions of variables:
boxplot(USArrests)

# Network visualization:
gplot(flo)

gplot(flo, label = colnames(flo))

gplot(flo, label = colnames(flo))

node.size <- (degree(flo) + 1)^0.5
gplot(flo, label = colnames(flo), vertex.cex = node.size)

### 1.8. Importing and exporting data
`?`(read.table)
`?`(read.csv)
`?`(scan)
apropos("read")
`?`(load)
`?`(save)
`?`(write.table)
apropos("write")


########## PART 2: NETWORK OBJECTS: IMPORT, EXPLORATION, MANIPULATION ####

### 2.1 Using built-in datasets #
library(network)
data(package = "network")
data(flo)
flo
# for more information....
`?`(data)
`?`(flo)

node.size <- (degree(flo) + 1)^0.5
gplot(flo, label = colnames(flo), vertex.cex = node.size)

### 2.2 Importing network data #


# Locate and read an adjacency matrix
getwd()
setwd("/Users/bulent/Desktop/ISSI2015-RSNA")
list.files()

# Importing co-autgorship informaion from an adjaceny matrix

AxA <- read.table("AxA_Adj.csv", header = FALSE, sep = ";")
class(AxA)
AxA.matrix <- as.matrix(AxA, matrix.type = "adjaceny")
Collaboration.Net <- as.network(AxA.matrix)
Collaboration.Net
gplot(Collaboration.Net)


# Importing the same data recorded in an edgelist format.

AxA.Edges <- read.table("AxA_Edg.csv", header = TRUE, sep = ",")
AxA.Edges
Collaboration.Net <- as.network(AxA.Edges[, 1:2], matrix.type = "edgelist", directed = FALSE)
Collaboration.Net
gplot(Collaboration.Net)

# for more information....
`?`(read.paj)
`?`(read.table)
`?`(network)
`?`(as.network.matrix)

### 2.3 Creating network objects from scratch #

# an empty network
g <- network.initialize(5)
g

# accessing vertex attributes
g %v% "vertex.names"

#Adding edges
g[1, 2] <- 1
g[2, c(1, 3, 5)] <- 1
g
gplot(g, label = g %v% "vertex.names")

# deleting edges
g[2, 5] <- 0
g
gplot(g, label = g %v% "vertex.names")
g[, ] <- 0
g

# adding/re-creating from new random edges
m <- matrix(rbinom(25, 1, 0.5), 5)
m
diag(m) <- 0
m
g
g[m > 0] <- 1
g

# setting edge weights
m <- matrix(1:5^2, nrow = 5, ncol = 5)
mg <- as.matrix(g)
m[mg == 0] <- 0
g
g %e% "link.weights" <- m

# retrieving edge values
list.edge.attributes(g)
g %e% "link.weights"

# storing link weights
g.adj.matrix <- as.sociomatrix(g, attrname = "link.weights")
g.adj.matrix

#For more information....
`?`(network.extraction)
`?`(ego.extract)
`?`(add.edge)
`?`(delete.edges)
`?`(delete.vertices)
`?`(get.edges)



### 2.4 Basic network descriptives  #
summary(Collaboration.Net)

print(Collaboration.Net)

network.dyadcount(Collaboration.Net)

network.edgecount(Collaboration.Net)

network.size(Collaboration.Net)

as.sociomatrix(Collaboration.Net)

Collaboration.Net[, ]

plot(Collaboration.Net, displaylabels = T, boxed.labels = F)
plot(Collaboration.Net, displaylabels = T, mode = "circle")

node.size <- (degree(Collaboration.Net) + 1)^0.5
node.label <- Collaboration.Net %v% "vertex.names"
gplot(Collaboration.Net, label = node.label, vertex.cex = node.size)

# for more information
`?`(summary.network)
`?`(network.dyadcount)
`?`(network.edgecount)
`?`(as.sociomatrix)
`?`(as.matrix.network)
`?`(is.directed)
`?`(plot.network)

### 2.5 Adding network and vertex attributes #

# Adding attributes
Collaboration.Net %v% "Degrees" <- degree(Collaboration.Net)
Collaboration.Net %n% "Description" <- "A sample co-autorship network."

# Listing attributes
list.vertex.attributes(Collaboration.Net)
list.network.attributes(Collaboration.Net)
list.edge.attributes(Collaboration.Net)

# Retrieving attributes
Collaboration.Net %v% "Degrees"
Collaboration.Net %n% "Description"

# for more information
`?`(attribute.methods)


########## PART 3: CLASSICAL NETWORK ANALYSIS WITH SNA PACKAGE ####

require(sna)

### 3.1 Network level descriptives #

## Connectivity:::

Connectivity <- gden(Collaboration.Net)

## Centralization:::

Centrality.InDeg <- centralization(Collaboration.Net, degree, cmode = "indegree")
Centrality.Eigen <- centralization(Collaboration.Net, evcent)


## Reciprocity:::

Reciprocity <- grecip(Collaboration.Net)

## Transitivity:::

Transitivity <- gtrans(Collaboration.Net)


## Sub-groups: Connected components:::


N.Groups <- components(Collaboration.Net)
Core <- component.largest(Collaboration.Net, result = "graph")
dev.off()
gplot(Core)
Groups <- component.dist(Collaboration.Net)
Groups
plot(1:length(Groups$cdist), Groups$cdist, xlab = "Size", ylab = "Frequency")

## More on connectivity, distance measurement, and cohesion:::

g <- rgraph(10, tp = 3/19)
g
is.connected(g)
gplot(g)
is.connected(g, connected = "weak")
geodist(g)
reachability(g)
symmetrize(g)
symmetrize(g, rule = "strong")


# for more information and choices...
`?`(centralization)
`?`(gden)
`?`(grecip)
`?`(gtrans)
`?`(hierarchy)
`?`(components)
`?`(component.dist)
`?`(component.largest)
`?`(clique.census)
`?`(components)
`?`(component.dist)
`?`(dyad.census)
`?`(is.isolate)
`?`(isolates)
`?`(kcycle.census)
`?`(kpath.census)
`?`(triad.census)



### 3.2 Node level analysis #


## Degree centralities:::

# Total degree
degree(Collaboration.Net)

# In degrees
ideg <- degree(Collaboration.Net, cmode = "indegree")

# Out degrees
odeg <- degree(Collaboration.Net, cmode = "outdegree")


# Visiualizing in-degree vs out-degrees with basic visualization techniques:
plot(ideg, odeg, type = "n", xlab = "Incoming", ylab = "Outgoing")

# Plotting ideg by odeg
abline(0, 1, lty = 3)
text(jitter(ideg), jitter(odeg), network.vertex.names(Collaboration.Net), cex = 0.75, col = 2)

# Plotting simple histograms of the degree distributions:
par(mfrow = c(1, 3)) # Set up a 1x3 display
hist(ideg, xlab = "Indegree", main = "Indegree Distribution", prob = TRUE)
hist(odeg, xlab = "Outdegree", main = "Outdegree Distribution", prob = TRUE)
hist(ideg + odeg, xlab = "Total Degree", main = "Total Degree Distribution", prob = TRUE)
par(mfrow = c(1, 1))


## Betweenness centralities::::

bet <- betweenness(Collaboration.Net, gmode = "graph")
bet
gplot(Collaboration.Net, vertex.cex = 1 + sqrt(bet)/2, gmode = "graph")


## Closeness centralities::::

clo <- closeness(Collaboration.Net)
clo

# Developing an alternative closeness:
# e.g: total geo-distances to all others
closeness.new <- function(x) {
	geo <- 1/geodist(x)$gdist
	diag(geo) <- 0
	# sum up the rows:
	apply(geo, 1, sum)
}

clo.new <- closeness.new(Collaboration.Net)
hist(clo.new, xlab = "Alt. Closeness", prob = TRUE)

# Running a correlation test on metrics:
cor(clo.new, bet)
plot(clo.new, bet)
cor.test(clo.new, bet)

#For more information....
`?`(betweenness)
`?`(bonpow)
`?`(closeness)
`?`(degree)
`?`(evcent)
`?`(graphcent)
`?`(infocent)
`?`(prestige)
`?`(stresscent)



########## PART 4: STORING R-OBJECTS ####

save.image("session.RData")

load("session.RData")


