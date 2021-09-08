
library(ape)

eucalypts <- read.tree("bamm/ML1_modified.tre")
plot(myrteae, cex=0.3)
axisPhylo()
myrteae

# Estimando taxas de especiação e extinção independente de caracter #
# Referencias: (Rabosky, 2014; et al., 2017) #

# 1) Baixar o programa em http://bamm-project.org; ¡¡LEIA O TUTORIAL COM CALMA DEPOIS!!
# 2) Colocar o programa na mesma pasta com todas as outras coisas (filogenia e "control file")
# 3) Excluir outgroups da filogenia (tudo que não é Myrteae) - várias formas de se fazer isso
myrteae_noOG<-drop.tip(myrteae, c("Leptospermum_scoparium","Eucalyptus_perriniana","Syzygium_maire",
                                  "Syzygium_gustavioides","Syzygium_buxifolium","Syzygium_jambos",
                                  "Syzygium_paniculatum","Syzygium_amplifolium","Syzygium_muellerii",
                                  "Syzygium_guineense","Metrosideros_stipularis","Medrosideros_perforata",
                                  "Metrosideros_nervulosa","Xanthomyrtus_montivaga","Xanthomyrtus_compacta",
                                  "Metrosideros_perforata"))

write.tree(myrteae_noOG, file = "myrteae_noOG")

# 4) Adicionar o prior ao "control file"
library(BAMMtools)
setBAMMpriors(read.tree("myrteae_noOG"), total.taxa = 2500) # porque estimamos que existam 2500 especies de Myrteae, mesmo que elas não estejam amostradas
# vai salvar um arquivo .txt com os prior empiricos no seu diretorio

# 5) Configurar o "control file"
# Abra o template do control file e ajuste o nome do arquivo da arvore, do sample size por clado e os priors 
# Vamos colocar poucas gerações de MCMC apenas para esse tutorial, mas aconselha-se rodar pelo menos 100.000.000 gerações quando pra valer
# 6) Rodar o BAMM
# Após abrir o programa, vá no terminal e execute "bamm -c 'DIREÇÃO_DO_CONTROL_FILE_NO_SEU_COMPUTADOR'"
# Note que rodamos a análise em si em c++ e não no R (seria muito lento)
# Se rodar direitinho, próximo passo é analisar os resultados no R

# Analisando os resultados #
# (1) Cheque se houve convergencia das MCMC
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts) #tem que ser acima de 200
effectiveSize(postburn$logLik) #aqui tambem

# A analise não convergiu porque a gente rodou muito poucas gerações de MCMC


# Vamos comparar com uma analise que eu já tinha rodado antes, com 100 000 000 gerações
mcmcout <- read.csv("mcmc_out_apB2.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts) #tem que ser acima de 200
effectiveSize(postburn$logLik)

xx <- mcmcout

# (2) Quantas mudanças de dinamica macroevolutiva? (rate shifts)
# 
# "Approaches that identify a single best shift configuration (e.g., stepwise AIC, 
# or other approaches that simply maximize the likelihood) are inherently limited 
# by their assumption that the model with the best information theoretic score 
# (AIC etc) is the model, given the candidate set of models. However, for most 
# real datasets, the best rate shift configuration is merely one of a large 
# number of possible rate shift configurations that have similar probabilities. 
# The BAMM philosophy is largely oriented around addressing this. To understand 
# the following examples, you must understand what we mean by distinct shift 
# configurations, credible sets of rate shift configurations, and marginal shift probabilities." 
#
# [Dessa forma, a filosofia do BAMM é semelhante a do Rpanda (Morlon et al. 2014)]


# Postprocess 1: Quanta variação de taxa? 

# Resumo da distribuição de probabilidade posterior para diferentes "shift configuration"
table(xx$N_shifts ) / nrow(xx)

install.packages("BAMMtools")
library(BAMMtools)
# Como isso se compara a informação a priori? (i.e.0 shifts)
plotPrior(xx)


# Postprocess 2: Lendo o "event data" / basic plotting

ed <- getEventData(myrteae, "event_data.txt", nsamples = 200)

summary(ed)
plot.bammdata(ed)

# Várias opções graficas a explorar
#
#  tau
#  pal

z <- plot.bammdata(ed, tau = 0.002, lwd=2)
addBAMMlegend(z)


# Ver uma "shift configuration" individual
# essa é uma
x <- subsetEventData(ed, index = 10)
plot.bammdata(x, tau = 0.002, lwd=3)
addBAMMshifts(x, cex=2)
# essa é outra
x <- subsetEventData(ed, index = 20)
plot.bammdata(x, tau = 0.002, lwd=3)
addBAMMshifts(x, cex=2)


# Loop  um set de amostras da distribuição posterior e plotar: 

pvec <- 1:9 * 10

quartz.options(height=10, width=10)
par(mar=c(1,1,1,1))

plot.new()
par(mfrow=c(3,3))

# get general color map: 
zcol <- plot.bammdata(ed, tau = 0.002, show=F)

for (ii in pvec ){
  x <- subsetEventData(ed, index = ii)
  plot.bammdata(x, tau = 0.002, method = "polar", lwd=1.3, colorbreaks = zcol$colorbreaks)
  addBAMMshifts(x, cex=2, par.reset=F, method= "polar")
  mtext(text = paste("sample ", ii, sep=""), side= 1)
  
}


# Postprocess 3: Taxa por clado 

# Primeiro vamos visualizar os nós na árvore:
plot.phylo(myrteae_noOG, show.tip.label = FALSE)
nodelabels()

tcols <- rep("black", length(myrteae_noOG$tip.label))
dclade <- extract.clade(myrteae_noOG, node = 110)$tip.label
tcols[myrteae_noOG$tip.label %in% dclade] <- "red"

# Destacando o clado do genero Eugenia
plot.phylo(myrteae_noOG, tip.color = tcols, cex = 0.7) 

# getCladeRates = media de taxas (mu, lambda) por clado em cada amostra da distribuição posterior

myrteaerates <- getCladeRates(ed)
names(myrteaerates)
myrteaerates$lambda # para as 200 amostras que estamos analisando


plot.new()
par(mfrow=c(2,1))
hist(myrteaerates$lambda, breaks=50, xlim = c(0, 0.2), col = "red") # especiação
hist(myrteaerates$mu, breaks=50, xlim=c(0, 0.2), col="blue") # lambda

# Agora vamos pegar a taxa media para tres grupos:
#   1. Apenas "Eugenia"
#   2. Apenas não Eugenias
#   3. tudo


# Vamos focar no nó 110 (ancestral comum das Eugenias)
#      por agora:

eugenia <- getCladeRates(ed, node = 110)
mean(eugenia$lambda)

# não-eugenias, usando nodetype = "exclude"
non_eugenia <- getCladeRates(ed, node = 110, nodetype = "exclude")
mean(non_eugenia$lambda)

# setup 3 panel plot:
plot.new()
par(mfrow=c(3,1))
hist(eugenia$lambda, breaks=50, col="red", xlim=c(0,0.35))
hist(non_eugenia$lambda, breaks=50, col = "blue", xlim=c(0,0.35))
hist(myrteaerates$lambda, breaks = 50, col = "lightgreen", xlim=c(0,0.35))



# Postprocess 4: Visualizando taxas ao longo do tempo

# Para Myrteae como um todo:

plotRateThroughTime(ed, avgCol = "blue")

# essa é a trajetoria total do "rate-through-time"

# Mas o mais interessante é ver que diferentes clados tiveram diferentes dinamicas
# então que tal se dividirmos de novo entre Eugenias e não-Eugenia?

ivec <- seq(0.1, 0.9, length.out=100)

plot.new()
par(mfrow=c(3,1))
plotRateThroughTime(ed, node = 110, intervals = ivec, intervalCol = "red", start.time = 36, ylim=c(0, 0.35))
mtext(side=3, text = "Eugenia", line = -2)

plotRateThroughTime(ed, node = 110, intervals = ivec, avgCol = "blue", start.time = 36, ylim=c(0, 0.35), nodetype = "exclude")
mtext(side=3, text = "Não-Eugenia", line = -2)

# e tudo junto?
plotRateThroughTime(ed, intervals = ivec, intervalCol = "darkgreen", avgCol = "darkgreen", start.time = 36, ylim=c(0, 0.35))
mtext(side=3, text = "Myrteae", line = -2)



# Tudo no mesmo plot:
par(new=T)
plot.new() 
plotRateThroughTime(ed, node = 110, intervals = ivec, intervalCol = "red", start.time = 36, ylim=c(0, 0.35))

plotRateThroughTime(ed, node = 110, intervals = ivec, avgCol = "blue", start.time = 36, ylim=c(0, 0.35), nodetype = "exclude", add=T)


# Postprocess 5
# 
# "shift-configurations" - Demonstrando algumas amostras diferentes de "shift-configuration"
# Vamos plotar o 95% "credible set" de diferentes "shift-configurations"

css <- credibleShiftSet(ed, expectedNumberOfShifts = 1)
plot(css, lwd=1.5)

# qual é o "melhor" total "shift-configuration"?

best <- credibleShiftSet(ed, expectedNumberOfShifts = 1, threshold = 5)

#?getBestShiftConfiguration
  
z <- plot(best, lwd=2)
addBAMMshifts(best, cex=3)
addBAMMlegend(z)

# Postprocess 6
#  "marginal shift probabilities":

mst <- marginalShiftProbsTree(ed)

plot(mst, cex = 0.7)
add.scale.bar(length = 0.5, lcol="red")

# lado a lado com a arvore
par(new=TRUE)
plot.new()
par(mfrow=c(1,2))

plot(myrteae, cex = 0.7)
plot(mst, cex = 0.7)


# Postprocess 7
# cohorts - o quanto diferentes clados compartilham uma dinamica macroevolutiva? 

cmat <- getCohortMatrix(ed)
cohorts(cmat, ed, use.plot.bammdata = TRUE, lwd=1.5)


res <- getTipRates(ed)
res$lambda


# Referências ####
# Paradis, E. (2011). Analysis of Phylogenetics and Evolution with R. Springer Science & Business Media.
# Rabosky, D. et al. (2014) BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees.Methods in Ecology and Evolution 5.7: 701-707.
  # Tutoriais no BAMM por http://bamm-project.org
# Revell, L. J. (2012). phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217-223.
  # Tutoriais do phytools e blog http://blog.phytools.org
# Vasconcelos, T. N. et al.(2017). Myrteae phylogeny, calibration, biogeography and diversification patterns: 
  # Increased understanding in the most species rich tribe of Myrtaceae. Molecular phylogenetics and evolution, 109, 113-137.

