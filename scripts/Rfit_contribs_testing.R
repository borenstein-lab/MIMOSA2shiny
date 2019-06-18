### Jaeckel dispersion functions
getWScores = function(nSamps){
  return(Rfit::getScores(Rfit::wscores, seq_len(nSamps)/(nSamps + 1)))
}

j.disp <- function(beta, x, y, scrs, ...) {
  e <- y - x %*% beta #E is residuals
  return(drop(crossprod(e[order(e)], scrs)))
}
j.disp.fit = function(yhat, y, scrs = NULL){
  if(is.null(scrs)) scrs = getWScores(length(y))
  e = y - yhat
  return(drop(crossprod(e[order(e)], scrs))) ###Note this is sum(e[order(e)]*scrs)
}
j.disp.mean = function(y, scrs = NULL){
  if(is.null(scrs)) scrs = getWScores(length(y))
  e = y - mean(y, na.rm = T)
  return(drop(crossprod(e[order(e)], scrs)))
}
j.disp.rsq = function(beta, x, y){
  sse = j.disp.fit(beta*x, y) #E is residuals)
  tot_e = j.disp.mean(y)
  return(1 - sse/tot_e)
}
j.disp.rsq.adj = function(null_disp, model_disp, tauhat, df2, df1 = 1){
  ## Adjusted R-squared from summary: 
  Fstat = ((null_disp-model_disp))/(tauhat/2)
  R2 = (df1/df2 * Fstat)/(1 + df1/df2 * Fstat)
  return(R2)
}


function (x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], scores = Rfit::wscores, 
          ...) 
{
  scrs <- getScores(scores, seq_len(length(y))/(length(y) + 
                                                  1)) ## Just normalized ranks basically
  j.grad <- function(x, y, beta, scores, ...) {
    x <- as.matrix(x)
    e <- y - x %*% beta
    r <- rank(e, ties.method = "first")/(length(e) + 1) #Rank of the residuals
    crossprod(x, -1 * getScores(scores, r))
  }

  sd.y <- sd(y)
  ystar <- y/sd.y
  fit0 <- optim(beta0/sd.y, j.disp, method = "BFGS", x = x, 
                y = ystar, scrs = scrs, gr = j.grad, scores = scores, 
                ...)
  optim(fit0$par * sd.y, j.disp, method = "BFGS", x = x, y = y, 
        scrs = scrs, gr = j.grad, scores = scores, ...)
}


ten_spec_dat = process_abundances("data/testData/sim_data/allSpeciesEnv3.txt", "data/testData/sim_data/allMetabolitesEnv3.txt", fluxes_file = "data/testData/sim_data/allMetFluxesEnv3.txt", simulated = T)
config_table = data.table(V1 = c("database", "genomeChoices","metType", "kegg_prefix", "data_prefix", "vsearch_path", "revRxns"), 
                          V2 = c("Greengenes 13_5 or 13_8", "AGORA genomes and models", "KEGG Compound IDs", "data/KEGGfiles/", "data/", "bin/vsearch", T))
config_table = rbind(config_table, data.table(V1 = "manualAGORA", V2 = T))
config_table = config_table[!V1 %in% c("rxnEdit", "revRxns")]
network_results = build_metabolic_model(ten_spec_dat[[1]], config_table, manual_agora = T, degree_filt = 0)
network = network_results[[1]]
species = network_results[[2]]
mets_melt = melt(ten_spec_dat[[2]], id.var = "compound", variable.name = "Sample")
indiv_cmps = get_species_cmp_scores(species, network, normalize = T, leave_rxns = F, manual_agora = ifelse("manualAGORA" %in% config_table[,V1], T, F), kos_only = F)
indiv_cmps = indiv_cmps[compound %in% mets_melt[,compound]]



tot_cmps = indiv_cmps[,sum(CMP), by=list(compound, Sample)]
tot_cmps = merge(tot_cmps, mets_melt[,list(compound, Sample, value)], by = c("compound", "Sample"))
all_comps = tot_cmps[,unique(compound)]
model_dat = data.table(compound = all_comps, Intercept = 0, Slope = 0, Rsq = 0, PVal = 0) #Make all other columns numeric
resid_dat = data.table(expand.grid(compound = all_comps, Sample = tot_cmps[,unique(Sample)]))
for(x in 1:length(all_comps)){
  scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], Rfit::rfit(value~V1)], error=function(e){ NA})
}
y = tot_cmps[compound==all_comps[x], value]
x = tot_cmps[compound==all_comps[x], V1]
scrs <- getScores(wscores, seq_len(length(y))/(length(y) + 1))
betas = coef(scaling_mod)[2]
j.disp(betas, as.matrix(x), y, scrs)

jaeckel.rsq = 1-(scaling_mod$D1/scaling_mod$D0) #Also equal to
j.disp.rsq(betas, x, y, scrs)

contribs2 = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table, nperm = 200)
#We can try this as well but my suspicion is that it doesn't matter
contribs2a = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table = config_table, adj_rsq = T, nperm = 200)
#Real question - relationship between dispersion and covariance
contribs_compare = merge(contribs2, contribs2a, by = c("compound", "Species", "NullDisp"))
contribs_compare[abs(value.x) < 10e-10, value.x:=0]
contribs_compare[abs(value.y) < 10e-10, value.y:=0]
contribs_compare[,RelContrib.x:=value.x/TrueRsq.x]
contribs_compare[,RelContrib.y:=value.y/TrueRsq.y]
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point()
contribs_compare[abs(RelContrib.y) > 1]
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point() + xlim(-0.5, 1.5) + ylim(-0.5, 1.5)
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point() + facet_wrap(~compound, scales = "free")

mets_melt = transform_mets(mets_melt, met_transform = "zscore")
contribs3 = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table, nperm = 200)
#We can try this as well but my suspicion is that it doesn't matter
contribs3a = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table = config_table, adj_rsq = T, nperm = 200)
#Real question - relationship between dispersion and covariance
contribs_compare = merge(contribs3, contribs3a, by = c("compound", "Species", "NullDisp"))
contribs_compare[abs(value.x) < 10e-10, value.x:=0]
contribs_compare[abs(value.y) < 10e-10, value.y:=0]
contribs_compare[,RelContrib.x:=value.x/TrueRsq.x]
contribs_compare[,RelContrib.y:=value.y/TrueRsq.y]
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point()
contribs_compare[abs(RelContrib.y) > 1]
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point() + xlim(-0.5, 1.5) + ylim(-0.5, 1.5)
ggplot(contribs_compare, aes(x=RelContrib.x, y = RelContrib.y)) + geom_point() + facet_wrap(~compound, scales = "free")
#I think it has to do with the tauhat - i don't think we want this
contribs_compare[RelContrib.x > 1]

nperms = c(50, 100, 200, 300, 500)
contribs_count_all = data.table()
time_dat = data.table(Nperm = nperms)
for(j in 1:length(nperms)){
  time_start = Sys.time()
  contribs_count = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table, nperm = nperms[j], return_perm = T) 
  time_end = Sys.time()
  contribs_count[,NPerm:=nperms[j]]
  time_dat[j, Time:=time_end-time_start]
  contribs_count_all = rbind(contribs_count_all, contribs_count, fill = T)
  print(time_dat)
}
contribs_count_all[abs(value) < 10e-10, value:=0]
perm_stats = contribs_count_all[,list(sd(value), sd(value)/abs(mean(value)), mean(value), sd(value[NPerm > 100])/mean(value[NPerm > 100])), by=list(compound, Species, TrueRsq)]
perm_stats[V1 !=0, summary(V2)]
perm_stats[V4 > 0.1]
## Do this better - return_perm for a really long run and then just look at subsets
contribs_long = rank_based_rsq_contribs(indiv_cmps, mets_melt, config_table, nperm = 2500, return_perm = T) 
perm_vals = data.table()
for(k in 1:25){
  contribs_sub = contribs_long$PermMat[(1+(k-1)*100):(k*100)]
  sub_dat = melt(contribs_sub[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)], id.var = "compound")
  perm_vals = rbind(perm_vals, sub_dat)
}
perm_vals[,sd(value), by=list(variable, compound)][,hist(V1)]
perm_vals[,sd(value)/abs(mean(value)), by=list(variable, compound)][,hist(V1)]
perm_vals = data.table()
for(k in 1:10){
  contribs_sub = contribs_long$PermMat[(1+(k-1)*250):(k*250)]
  sub_dat = melt(contribs_sub[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)], id.var = "compound")
  perm_vals = rbind(perm_vals, sub_dat)
}
perm_vals[,sd(value), by=list(variable, compound)][,hist(V1)]
perm_vals[,sd(value)/abs(mean(value)), by=list(variable, compound)][V1 != 0][,hist(V1)]

for(k in 1:){
  contribs_sub = contribs_long$PermMat[(1+(k-1)*500):(k*500)]
  sub_dat = melt(contribs_sub[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)], id.var = "compound")
  perm_vals = rbind(perm_vals, sub_dat)
}
perm_vals[,sd(value), by=list(variable, compound)][,hist(V1)]
perm_vals[,sd(value)/abs(mean(value)), by=list(variable, compound)][V1 != 0][,hist(V1)]
perm_vals[,sd(value)/abs(mean(value)), by=list(variable, compound)][V1 > 10]
perm_vals[,Species:=gsub("Marg_", "", variable)]
perm_vals = merge(perm_vals, contribs_long$Contribs, by = c("compound", "Species"))
perm_vals[,hist(abs((value.x-value.y)/value.y))]
## Let's try this on real data

#Use un-normalized mets to avoid confusion
mets_melt = melt(ten_spec_dat[[2]], id.var = "compound", variable.name = "Sample")
comp_test1 = "ac[e]"
test_cmps = indiv_cmps[compound==comp_test1]
test_mets = mets_melt[compound==comp_test1]
test_mets[,CenteredScaled:=scale(value), by=compound]
test_cmps = merge(test_cmps, test_mets, by = c("Sample", "compound"))
tot_test_cmps = test_cmps[,sum(CMP), by=list(Sample, value, CenteredScaled)]

test_mod = tot_test_cmps[,lm(value~V1)]
test_mod_scale = tot_test_cmps[,lm(CenteredScaled~V1)]
test_cmps[,ModPred:=CMP*coef(test_mod)[2]]
test_cmps[,ModPredScaled:=CMP*coef(test_mod_scale)[2]]
test_cmps[,cov(ModPred, value), by=Species]
test_cmps[,cov(ModPred, value)/var(value), by=Species][,sum(V1)] #Rsq 
test_cmps[,cov(ModPredScaled,CenteredScaled)/var(CenteredScaled), by=Species][,sum(V1)] #Ok phew

mod1 = fit_cmp_mods(test_cmps, mets_melt, rank_based = F)
test_cmps_resids = add_residuals(test_cmps, mod1[[1]], mod1[[2]])
test_cmps_resids[,sum(newValue)+ mod1[[1]][,Intercept] , by=Sample]#Just offset by intercept
var_shares_test = calculate_var_shares(test_cmps_resids)
test_cmps[,cov(ModPred, value), by=Species]
var_shares_test[,V1] #yepppp - why is this?
var_shares_compare = merge(test_cmps[,cov(ModPred, value), by=Species], var_shares_test, by = "Species", all = T)
var_shares_compare[,V1.x-V1.y] #0
## These are no longer equal depending on metabolite scaling!! what? 

tot_test_cmps[,cor(V1, value)*cor(V1, value)]

tot_test_cmps[,FittedValues:=test_mod$fitted.values]
tot_test_cmps[,Resid:=value-FittedValues]
test_cmps = merge(test_cmps, tot_test_cmps[,list(Sample, FittedValues, Resid)], by="Sample")

test_cmps[,cov(Resid, CMP), by=Species] #why are these additive? vars are not independnet
tot_test_cmps[,cov(Resid, V1)] #Basically 0
tot_test_cmps[,cov(Resid, value)/(var(value))] #1-Rsq 

#Wait so what is this function in dispersion-land? This is just the model dispersion right?

test_mod2 = tot_test_cmps[,rfit(value~V1)]
test_mod2$D0
test_mod2$D1
(test_mod2$D0-test_mod2$D1)/test_mod2$D0
test_mod2_scaled = tot_test_cmps[,rfit(CenteredScaled~V1)]
test_mod2_scaled$D0
test_mod2_scaled$D1
(test_mod2_scaled$D0-test_mod2_scaled$D1)/test_mod2_scaled$D0


tot_test_cmps[,FittedValues2:=test_mod2$fitted.values]
tot_test_cmps[,Resid2:=value-FittedValues2]

test_cmps = merge(test_cmps, tot_test_cmps[,list(FittedValues, Resid, FittedValues2, Resid2, Sample)], by="Sample")

j.disp.rsq(coef(test_mod2)[2], tot_test_cmps[,V1], tot_test_cmps[,value]) # ok, fixed

tot_test_cmps[,cov(Resid, V1)] #0
tot_test_cmps[,j.disp.fit(Resid2, V1)]
tot_test_cmps[,j.disp.fit(FittedValues2, value)]

test_cmps[,FittedValues:=CMP*coef(test_mod)[2]]
test_cmps[Species=="BcaccaeAGORA", FittedValues:=FittedValues + coef(test_mod)[1]] #Make it add up
test_cmps[,FittedValues2:=CMP*coef(test_mod2)[2]]
test_cmps[Species=="BcaccaeAGORA", FittedValues2:=FittedValues2 + coef(test_mod2)[1]] #Make it add up

setnames(test_cmps, c("Resid", "Resid2"), c("TotResid", "TotResid2"))

tot_test_cmps[,cov(Resid, V1)]#0
tot_test_cmps[,cov(V1, value)]
test_cmps[,cov(TotResid,FittedValues), by=Species] [,sum(V1)]#0
test_cmps[,var(value)-sum((value-FittedValues)^2), by=Species]

test_cmps[,j.disp.fit(FittedValues2, TotResid2), by=Species]
test_cmps[,j.disp.fit(FittedValues2, TotResid2), by=Species][,sum(V1)] #Nope
test_cmps[,j.disp.fit(FittedValues2, value), by=Species]
null_disp = j.disp.mean(tot_test_cmps[,value])

test_cmps[,null_disp-j.disp.fit(FittedValues2, value), by=Species][,sum(V1)] # not equal to explained_disp - not additive

explained_disp = null_disp - j.disp.fit(tot_test_cmps[,FittedValues2], tot_test_cmps[,value])
#Resid_disp = unexplained_disp


test_cmps_wide = dcast(test_cmps, Sample~Species, value.var = "ModPred")
test_cmps_wide[,Resid:=tot_test_cmps[,Resid2]]

all_species = names(test_cmps_wide)[2:ncol(test_cmps_wide)]
pairwise_disps = sapply(all_species, function(x){
  sapply(all_species, function(y){
    test_cmps_wide[,j.disp.fit(get(x), get(y))]
  })
})
pairwise_disps = data.table(pairwise_disps, Species = all_species)
pairwise_disps[]
test_cmps[,j.disp.fit(Resid2, CMP), by=Species] #[,sum(V1)] #Not additive

colSums(pairwise_disps[,1:(ncol(pairwise_disps)-1), with=F])
rowSums(pairwise_disps[,1:(ncol(pairwise_disps)-1), with=F])
#I don't think these are meaningful
rowMeans(pairwise_disps[,1:(ncol(pairwise_disps)-1), with=F])
pairwise_disps[,1:(ncol(pairwise_disps)-1), with=F]/null_disp

x = c(1,3,4)
y = c(1,3,4)
z = x+y
cov(x,z)
cov(y,z)
var(z)
cov(x,z) + cov(y,z)

##
cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt, rank_based = F)
indiv_cmps2 = add_residuals(indiv_cmps, model_dat = cmp_mods[[1]], resid_dat = cmp_mods[[2]])
var_shares1 = calculate_var_shares(indiv_cmps2)
var_shares2 = calculate_var_shares(indiv_cmps2, model_cov = T)
var_shares_comp = merge(var_shares1, var_shares2, by = c("compound", "Species"), all = T)
var_shares_comp[is.na(VarShare.y) & VarShare.x != 0]
var_shares_comp[is.na(VarShare.x) & VarShare.y != 0]
var_shares_comp[VarShare.x != VarShare.y]
var_shares_comp[abs(VarShare.x-VarShare.y) > 10e-12] #Ok then
var_shares_comp[abs(V1.x - V1.y) > 10e-12]

met_fluxes = fread("data/testData/sim_data/allMetFluxesEnv3.txt")
flux_contribs_test1 = getContributions(met_fluxes)
flux_contribs_test2 = getContributions(met_fluxes, alt_calculate = T)
flux_contribs_compare = merge(flux_contribs_test1, flux_contribs_test2, by = c("compound", "Species"), all = T)
flux_contribs_compare[abs(VarShare.x - VarShare.y) > 10e-12]
flux_contribs_compare[abs(V1.x - V1.y) > 10e-14] #Ok i think
