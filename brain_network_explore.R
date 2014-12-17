m1 = readLines("data/brain/mdcbc815.network")
m2 = readLines("data/brain/pfc17.network")

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
n1 = n1 + t(n1)

for(j in seq_along(conns2)) {
  gs = strsplit(conns2[j], split = "\tconnect\t")[[1]]
  ms = match(gs, gene)
  if (!(any(is.na(ms)))) {
    n2[ms[1], ms[2]] = 1
  }
}
n2 = n2 + t(n2)

length(conns1)
sum(n1)/2
length(conns2)
sum(n2)/2

library(igraph)
g1 = graph.adjacency(n1)
g2 = graph.adjacency(n2)
V(g1)$name = ""
V(g2)$name = ""
par(mfrow = c(1,2))
plot(g1, edge.arrow.width = 0, vertex.size = 0.2, main = "Network 1")
plot(g2, edge.arrow.width = 0, vertex.size = 0.2, main = "Network 2")

sum(apply(n2, 1, sum) + apply(n1, 1, sum) == 0)

gen_network = sample_generating_models(Nnodes = 5475, mode = "block", K = 20, is_null = TRUE)
#|----##Edit paramter K. Make it into model_params, a list taking in possible model parameters, passed into specific models. --Sat Nov 29 02:46:46 2014--
new_probmat = matrix(runif(n = 400, min = 0.0000001, max = 0.01), nrow = 20, ncol = 20)

nm = sample_generating_models(Nnodes = 5000, mode = "block", K = 3, is_null = TRUE)
#|----##Edit paramter K. Make it into model_params, a list taking in possible model parameters, passed into specific models. --Sat Nov 29 02:46:46 2014--
ms = sample_network_pair(gen_model = nm)
perform_hyptest(ms[[1]], ms[[2]], mode = "block", n_models = 10, pl = list(cc_adj = 2, thres_ignore = 5, alphas = 0.05), pval_adj_fx = mult_pearson) -> test
nodes = c(100, 500, 1000, 2000, 5000)
timings = c(0.003, 0.02, 0.047, 0.557, 2.433)
plot(nodes, timings)



perform_hyptest(n1, n2, mode = "block", n_models = 100, pl = list(cc_adj = 2, thres_ignore = 5, alphas = 0.05), pval_adj_fx = mult_pearson) -> test

perform_hyptest(n1, n2, mode = "block", n_models = 100, pl = list(cc_adj = 2, thres_ignore = 5, alphas = 0.05), pval_adj_fx = mult_highcrit) -> test
