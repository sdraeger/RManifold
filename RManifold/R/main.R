library(class)
library(tidyverse)
library(rgl)
library(RDRToolbox)
library(cccd)
library(Rfast)


isomap <- function(X, k, d) {
  # Compute k neighbor graph
  graph = nng(data, k=k)
  # Compute distance matrix
  ds = dist(data)
  dsm = as.matrix(ds)
  # Compute shortest pairwise paths between vertices
  shortest_paths = floyd(dsm)
  # Apply MDS on distance matrix
  fit = cmdscale(shortest_paths, eig=TRUE, k=d)
  fit
}


data = SwissRoll(N = 1000)
df = data.frame(data)

colnames(df) = c('x', 'y', 'z')

ggplot(df, aes(x=x, y=y)) +
  geom_point(aes(colour=z))

k = 5
d = 2
fit = isomap(data, k, d)
emb_df = data.frame(fit$points)
colnames(emb_df) = c('x', 'y')
ggplot(emb_df, aes(x=x, y=y)) + geom_point()
