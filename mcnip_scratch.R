e = 0.00026
cn_ex = 27.6
cn_s = 27.6
cn_m = 10
cn_enz = 3
p = 0.5
q = 0.5
a = 0.5
r = 0.008314
r_death = 0.00015
r_ecloss = 0.001
mic_to_soc = 0.5
mic_to_son = 0.5
doc_input = 0.0005
litter_c = 0.0005
don_input = 0.0005/27.6
litter_n = 0.0005/27.6
a_upt_c = 1e8
vmax_0 = 1e8
ea_upt_c = 48
ea_up = 48
b_upt_c = 0.1
m_upt_c = 0.01
km_0 = 500
km_slope = 5
b_cue = 0.63
m_cue = -0.016
mic_c = 1.1957
mic_n = 0.1196
soc = 144.5986
son = 5.4413
doc = 0.00091631
don = 0.00049421
ec = 0.0325

pars <- c(e = 0, cn_ex = 27.6, cn_s = 27.6, cn_m = 10, cn_enz = 3, p = 0.5, q = 0.5, a = 0.5, r = 0.008314, r_death = 0.00015, r_ecloss = 0.001, mic_to_soc = 0.5, mic_to_son = 0.5, doc_input = 0.0005, litter_c = 0.0005, don_input = 0.0005/27.6, litter_n = 0.0005/27.6, a_upt_c = 1e8, vmax_0 = 1e8, ea_upt_c = 48, ea_up = 48, b_upt_c = 0.1, m_upt_c = 0.01, km_0 = 500, km_slope = 5, b_cue = 0.63, m_cue = -0.016, temp = 20)
state <- c(mic_c = 0.5, mic_n = 0.05, soc = 100, son = 3.6232, doc = 0.5, don = 0.0333, ec = 0.01)


derivs <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    #relevant model equations
    #calculate uptake kinetics
    vmax_upt_c <- vmax_upt_n <- a_upt_c * exp(-ea_upt_c / (r * (temp + 273)))
    km_upt_c <- km_upt_n <- b_upt_c + m_upt_c * temp #not in damm-mcnip
    #calculate depolymerization kinetics
    vmax_c <- vmax_n <- vmax_0 * exp(-ea_up / (r * (temp + 273)))
    km_c <- km_n <- km_0 + km_slope * temp #not in damm-mcnip
    cue <- b_cue + m_cue * temp #not in damm-mcnip
    upt_c <- mic_c * vmax_upt_c * (doc / (km_upt_c + doc)) #microbial C uptake
    cmin <- upt_c * (1 - cue) #C mineralization
    upt_n <- mic_n * vmax_upt_n * (don / (km_upt_n + don)) #microbial N uptake
    death_c <- r_death * mic_c #microbial C turnover
    death_n <- death_c / cn_m #r_death * mic_n #microbial N turnover
    enz_c <- p * cue * upt_c #amount of C available for enzyme production after C allocated to mineralization
    enz_n <- q * upt_n #amount of N available for enzyme production
    eprod <- ifelse(enz_c / cn_enz >= enz_n, enz_n, enz_c / cn_enz) #if enzyme production is N-limited, then number of enzyme produced = N cost for enzyme production (1 enzyme made for every 1 N used), else enzyme production C-limited
    growth_c <- (1 - p) * upt_c * cue + enz_c - cn_enz * eprod #available C for biomass growth
    growth_n <- (1 - q) * upt_n + enz_n - eprod #available N for biomass growth
    growth <- ifelse(growth_c / cn_m >= growth_n, growth_n, growth_c / cn_m) #if microbes N-limited, then amount of microbial biomass growth equals cost of N to produce microbes, else growth is C-limited
    overflow_c <- growth_c - cn_m * growth #extra C after microbial growth goes to overflow metabolism
    nmin <- growth_n - growth #N mineralization
    dmic_c <- cn_m * growth - death_c #microbial biomass C pool, growth multiplied by C:N of microbes because growth = N cost to build a microbe, but C cost is greater
    dmic_n <- growth - death_n #microbial biomass N pool
    eloss <- r_ecloss * ec #enzyme turnover
    dec <- eprod - eloss #enzyme pool
    decom_c <- vmax_c * a * ec * soc / (km_c + soc) #depolymerization of soc by enzymes
    decom_n <- vmax_n * (1 - a) * ec * son / (km_n + son) #depolymerization of son by enzymes
    dsoc <- litter_c + death_c * mic_to_soc - decom_c
    dson <- litter_n + death_n * mic_to_son - decom_n
    ddoc <- doc_input + e + decom_c + death_c * (1 - mic_to_soc) + (cn_enz / (1 + cn_enz)) * eloss - upt_c #doc pool, enzyme turnover deposited here.
    ddon <- don_input + e / cn_ex + decom_n + death_n * (1 - mic_to_son) + (1 / cn_enz) * eloss - upt_n
    #don pool, enzyme turnover deposited here. because enzymes are not split up into separate pools within the model, the amount of C and N created by the turnover of an enzyme is calculated using the C:N of enzymes
    list(c(dmic_c, dmic_n, dec, dsoc, dson, ddoc, ddon))
  })
}

#use rootSolve
library(rootSolve)
test  <- runsteady(y = state, times = c(0,100), func = derivs, parms = pars)
