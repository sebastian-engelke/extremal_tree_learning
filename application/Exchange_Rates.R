# BoE_website: https://www.bankofengland.co.uk/boeapps/database/Rates.asp?Travel=NIxASx&into=GBP

## Fetch data from the Bank of England website
# 
#library(pdfetch)
#BoE_exchange_rates <- as.matrix(pdfetch_BOE(identifiers = c("XUDLADS","XUDLCDS","XUDLBK89", "XUDLBK25", "XUDLDKS","XUDLERS","XUDLHDS","XUDLBK33","XUDLBK97","XUDLBK78","XUDLJYS","XUDLBK83","XUDLNDS","XUDLNKS","XUDLBK47","XUDLBK85","XUDLSRS","XUDLSGS","XUDLZRS","XUDLBK93","XUDLSKS","XUDLSFS","XUDLTWS","XUDLBK87","XUDLBK95","XUDLUSS"), 
#                                            from =  "2005-10-01", to = "2020-09-30"))
#colnames(BoE_exchange_rates) <- ctry_codes

# save(BoE_exchange_rates, file = "~/Dropbox/structurelearning/code/Application/Exchange_Data/BoE_exchange_rates.RData")

load(file = "data/BoE_exchange_rates.RData")
data <- BoE_exchange_rates
ctry_codes <- c("AUS", "CAN", "CHN","CZE","DNK","EUR","HKG","HUN","IND","ISR","JPN","MYS","NZL","NOR","POL","RUS","SAU","SGP","ZAF","KOR","SWE","CHE","TWN","THA","TUR","USA")

library("igraph")
library("graphicalExtremes")
library("timeSeries")
library("fGarch")
library("tidyverse")
library("latex2exp")

coords_tree <- read_rds("data/coords_exchange_rates.rds")

# create negative log ratios
nlr <-  apply(data, 2, FUN = function(x) diff(log(x))) #ndifflog

d <-  dim(nlr)[2]
n <-  dim(nlr)[1]

x <- nlr
fit.fGarch <- list()
residus <-  matrix(NA,n,d)
colnames(residus) <-  colnames(nlr)

for(i in 1:d){
  ga <- NULL
  form <- paste("~arma(0,2)+garch(1,1)")
  ga <- garchFit(formula=as.formula(form), data=x[,i], trace=FALSE, cond.dist="norm")
  residus[,i] <- ga@residuals/ga@sigma.t
  fit.fGarch[[i]] <- ga
  cat("Done for ",i,"/",d,"\n")
}


colnames(residus) <-  colnames(nlr)
abs.residus <- abs(residus)
acf(abs.residus[,10])

data <- abs.residus 
graph.full <- make_full_graph(d)

p <- 0.97
k <- n * (1-p)
B <- 100
edge.mat <- matrix(0, nrow=d, ncol=d)
for(i in 1:B){
  G.est <-  emp_vario(data=data[sample(1:nrow(data), size = nrow(data), replace = TRUE),], p = p)
  MST.est <- igraph::mst(graph=graph.full, weights = 2- Gamma2chi(G.est[ends(graph.full,E(graph.full))]), algorithm = "prim")
  edge.mat[ends(MST.est,E(MST.est))] <- edge.mat[ends(MST.est,E(MST.est))] + 1
}

edge.mat.sym <- edge.mat + t(edge.mat)
MST.graph <- graph_from_adjacency_matrix(edge.mat.sym!=0, mode="undirected")
vertex_attr(MST.graph) <- list(name = colnames(data))
MST.graph <- graphicalExtremes:::set_graph_parameters(MST.graph)
E(MST.graph)$width <- edge.mat.sym[ends(MST.graph,E(MST.graph), names=FALSE)] /10
igraph::V(MST.graph)$size <- 13
igraph::V(MST.graph)$color <- grDevices::adjustcolor(col = "#AAC1D8")


pdf(file = paste("Exchange_rate_Tree_bootstrap.pdf", sep=""), width = 7)
par(mar=c(4,4,1.6,1.5), mgp=c(2,0.6,0), pty="s", cex.lab=1.6, cex.axis = 1.3, cex.main=2, pch=1, cex=0.9, lwd=1)
plot(MST.graph, layout=coords_tree)
dev.off()



G.est = emp_vario(data=data, p = p)
colnames(G.est) <- rownames(G.est) <- ctry_codes
MST.est <- igraph::mst(graph=graph.full, weights = 2- Gamma2chi(G.est[ends(graph.full,E(graph.full))]), algorithm = "prim")
vertex_attr(MST.est) <- list(name = colnames(data))
MST.est <- graphicalExtremes:::set_graph_parameters(MST.est)
E(MST.est)$width <- Gamma2chi(G.est[ends(MST.est,E(MST.est), names=FALSE)] ) * 10
igraph::V(MST.est)$color <- grDevices::adjustcolor(col = "#AAC1D8")

igraph::V(MST.est)$size <- 13

pdf(file = paste("Exchange_rate_Tree.pdf", sep=""), width = 7)
par(mar=c(4,4,1.6,1.5), mgp=c(2,0.6,0), pty="s", cex.lab=1.6, cex.axis = 1.3, cex.main=2, pch=1, cex=.9, lwd=1)
plot(MST.est, layout=coords_tree)
dev.off()


chi.tree <- Gamma2chi(complete_Gamma(Gamma=G.est, graph=MST.est))
chi.est = emp_chi(data=data, p = p)

pdf(file = paste("Exchange_rate_Tree_ECs.pdf", sep=""), width = 7)
par(mar=c(4,4,1.6,1.5), mgp=c(2,0.6,0), pty="s", cex.lab=1.6, cex.axis = 1.3, cex.main=2, pch=1, cex=1.3, lwd=1)
plot(chi.tree, chi.est, ylab="Empirical", main="", las=1, xlim=c(0,1), ylim = c(0,1),
     xlab="Fitted model", pch=19, lwd=2)
abline(0,1)
dev.off()



### Exploratory analysis

library("SpatialADAI")
qseq <- seq(.8, .99, length.out = 100)

i <- "EUR"
j <- "SGP"

pdf(file = paste("chiplot_",i,j,".pdf", sep=""), width = 7)
par = par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,6,2) +.1)
chiu.est(data[,c(i,j)], u = qseq, nbs=100,blocklength = 1, xlim=c(min(qseq),1),xlab="1-q", ylab=expression(hat(chi)(q)))
dev.off()



#### Choice of k


G_tree_err <- function(G, tree){
  G_compl <- complete_Gamma(Gamma = G, graph = tree)
  chi <- Gamma2chi(G)
  chi_compl <- Gamma2chi(G_compl)
  return(mean(abs( (chi[upper.tri(chi)]-chi_compl[upper.tri(chi)]) )^2))  
}



p_vec <- seq(.85, 0.997, length.out = 100)
emst_list <- lapply(p_vec, FUN = function(p) emst(data = data, p=p, method = "vario")$graph)
G_list <- lapply(p_vec, FUN = function(p) emp_vario(data = data, p=p))

ll <- length(emst_list)
tree_err <- numeric(ll)
for(i in 1:ll){
  tree_err[i] <- G_tree_err(G = G_list[[i]], tree = emst_list[[i]]) 
}


pdf(file = paste("chi_error.pdf", sep=""), width = 7)
par(mar=c(4,4,1.6,1.5), mgp=c(2,0.6,0), pty="s", cex.lab=1.6, cex.axis = 1.3, cex.main=2, pch=1, cex=1.3, lwd=1)
plot(p_vec, tree_err*10^2, xlab="1-q", ylab=TeX("$\\hat{\\Delta}[chi](q) \\times 10^2$"), las=1)
dev.off()

