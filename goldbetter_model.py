__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import numpy as np

# Constants
k1 = 0.4  # 1/h
k2 = 0.2  # 1/h
k3 = 0.4  # 1/nMh
k4 = 0.2  # 1/h
k5 = 0.4  # 1/h
k6 = 0.2  # 1/h
k7 = 0.5  # 1/nMh
k8 = 0.1  # 1/h
K_AP = 0.7  # nM
K_AC = 0.6  # nM
K_IB = 2.2  # nM
k_dmb = 0.01  # 1/h
k_dmc = 0.01  # 1/h
k_dmp = 0.01  # 1/h
k_dn = 0.01  # 1/h
k_dnc = 0.12  # 1/h
K_d = 0.3  # nM
K_dp = 0.1  # nM
K_p = 0.1  # nM
K_mB = 0.4  # nM
K_mC = 0.4  # nM
K_mP = 0.31  # nM
k_sb = 0.12  # 1/h
k_sC = 1.6  # 1/h
k_sP = 0.6  # 1/h
m = 2  #
n = 4  #
v_1b = 0.5  # nM/h
v_1c = 0.6  # nM/h
v_1p = 0.4  # nM/h
v_1pc = 0.4  # nM/h
v_2b = 0.1  # nM/h
v_2c = 0.1  # nM/h
v_2p = 0.3  # nM/h
v_2pc = 0.1  # nM/h
v_3b = 0.5  # nM/h
v_3pc = 0.4  # nM/h
v_4b = 0.2  # nM/h
v_4pc = 0.1  # nM/h
v_phos = 0.4  # nM/h
v_dbc = 0.5  # nM/h
v_dbn = 0.6  # nM/h
v_dcc = 0.7  # nM/h
v_din = 0.8  # nM/h
v_dpc = 0.7  # nM/h
v_dpcc = 0.7  # nM/h
v_dpcn = 0.7  # nM/h
v_mb = 0.8  # nM/h
v_mC = 1  # nM/h
v_mP = 1.1  # nM/h
v_sB = 1  # nM/h
v_sC = 1.1  # nM/h
v_sP = 1.5  # nM/h

def hill(L, Ka, coeff):
    theta = L**coeff/(Ka**coeff+L**coeff)
    return theta

def mm(substrate, Km):
    rate = substrate/(Km+substrate)
    return rate

# Equations
# mRNAs of Per, Cry, and Bmal1
dMP_dt = v_sP * hill(B_n, K_AP, n) - v_mP * mm(M_P, K_mP) - k_dmp * M_P
dMC_dt = v_sC * hill(B_n, K_AC, n) - v_mC * mm(M_C, K_mC) - k_dmc * M_C
dMB_dt = v_sB * hill(B_n, K_IB, m) - v_mP * mm(M_B, K_mB) - k_dmb * M_B #todo: This equation seems off, maybe typo in supplement

# Phosphorylated and nonphosphorylated proteins PER and CRY in cytosol
dPC_dt = k_sP * M_P - v_1p * mm(P_C, K_p) + v_2p * mm(P_CP, K_dp) + k4 * PC_C - k3 * P_C * C_C - k_dn * P_C
dCC_dt = k_sC * M_C - v_1c * mm(C_C, K_p) + v_2c * mm(C_CP, K_dp) + k4 * PC_C - k3 * P_C * C_C - k_dnc * C_C
dPCP_dt = v_1p * mm(P_C, K_p) - v_2p * mm(P_CP, K_dp) - v_dpc * mm(P_CP, K_d) - k_dn * P_CP
dCCP_dt = v_1c * mm(C_C, K_p) - v_2c * mm(C_CP, K_dp) - v_dcc * mm(C_CP, K_d) - k_dn * C_CP

