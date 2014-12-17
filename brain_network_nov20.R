m1 = readLines("data/brain/new_20141120/mdcbc810_ascssc_085.network")
m2 = readLines("data/brain/new_20141120/pfc35_ascssc_07.network")

gene_rows_m1 = min(grep("connect", m1)) - 1
genes1 = m1[1:gene_rows_m1]
conns1 = m1[-1:-gene_rows_m1]


gene_rows_m2 = min(grep("connect", m2)) - 1
genes2 = m2[1:gene_rows_m2]
conns2 = m2[-1:-gene_rows_m2]

length(intersect(genes1, genes2))

gene = intersect(genes1, genes2)

n1 = matrix(0, nrow = length(gene), ncol = length(gene))
n2 = matrix(0, nrow = length(gene), ncol = length(gene))

for(j in seq_along(conns1)) {
  gs = strsplit(conns1[j], split = "\tconnect\t")[[1]]
  ms = match(gs, gene)
  if (!(any(is.na(ms)))) {
    n1[ms[1], ms[2]] = 1
  }
}
#n1 = n1 + t(n1)

for(j in seq_along(conns2)) {
  gs = strsplit(conns2[j], split = "\tconnect\t")[[1]]
  ms = match(gs, gene)
  if (!(any(is.na(ms)))) {
    n2[ms[1], ms[2]] = 1
  }
}
#n2 = n2 + t(n2)

length(conns1)
sum(n1)/2
length(conns2)
sum(n2)/2
table(apply(n1 + n2, 1, sum))

rmat = matrix(0, nrow = 10, ncol = 10)
rmat[2:5, 1] = 1
rmat[1,2:5] = 1
rmat[3:6, 7:8] = 1
rmat[7:8, 3:6] = 1
rmat[9,10] = 1
rmat[10,9] = 1
conComp(rmat)
library(ggm)


library(igraph)
g1 = graph.adjacency(n1)
g2 = graph.adjacency(n2)
V(g1)$name = ""
V(g2)$name = ""
par(mfrow = c(1,2))
plot(g1, edge.arrow.width = 0, vertex.size = 0.2, main = "Network 1")
plot(g2, edge.arrow.width = 0, vertex.size = 0.2, main = "Network 2")

sum(apply(n2, 1, sum) + apply(n1, 1, sum) == 0)
## load sbm data

ncomb = n1 + n2

fit1 = fit_SBM(A = n1, Nobs = 1, q = 5, niter = 25, stop_thres = 0.000001)
#|----##--Reparameterizing this function --Wed Dec 17 15:10:23 2014--
fit2 = fit_SBM(A = n2, Nobs = 1, q = 5, niter = 25, stop_thres = 0.000001)
#|----##--Reparameterizing this function --Wed Dec 17 15:10:23 2014--
fitc = fit_SBM(A = ncomb, Nobs = 1, q = 5, niter = 25, stop_thres = 0.000001)
#|----##--Reparameterizing this function --Wed Dec 17 15:10:23 2014--







