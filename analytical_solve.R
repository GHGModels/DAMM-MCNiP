temp = 20
soilM = 0.229
R = 0.008314
CN_s = 27.6
CN_m = 10
CN_enz =3
p = 0.5
q = 0.5
a = 0.5
r_death = 0.00015
r_ECloss = 0.001
MIC_to_SOC = 0.5
MIC_to_SON = 0.5
DOC_input = 0.0005
DON_input = DOC_input/CN_s
Litter_C = 0.0005
Litter_N = Litter_C/CN_s
A_UPT_C = 100000000
Ea_UPT_C = 48
b_UPT_C = 0.1
m_UPT_C = 0.01
Vmax_0 = 100000000
Ea_up = 48
Km_0 = 0.0015
Km_slope = 0.00005
Km_O2 = 0.121
Dgas = 1.67
O2airfrac = 0.209
BD = 0.8
PD = 2.52
porosity = 1 - BD/PD;
frac = 0.000414
Dliq = 3.17

O2 <- Dgas*O2airfrac*((porosity-soilM)^(4/3))
vmax_UPT_C <- A_UPT_C*exp(-Ea_UPT_C/(R*(temp+273)))
km_UPT_C <- b_UPT_C + m_UPT_C*temp
CUE <- 0.31
vmax_upt_N = vmax_UPT_C
km_upt_N = km_UPT_C
Vmax_C = Vmax_0*exp(-Ea_up/(R*(temp+273)))
Km_C = 0.0025
Vmax_N = Vmax_C
Km_N = Km_C


