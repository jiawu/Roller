__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def hill(L, Ka, coeff):
    theta = (L**coeff)/(Ka**coeff+L**coeff)
    return theta

def mm(substrate, Km):
    rate = substrate/(Km+substrate)
    return rate

def goldbetter_ode(y, time):
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
    v_mB = 0.8  # nM/h
    v_mC = 1  # nM/h
    v_mP = 1.1  # nM/h
    v_sB = 1  # nM/h
    v_sC = 1.1  # nM/h
    v_sP = 1.5  # nM/h

    # Unpack y
    B_N, M_P, M_C, M_B, P_C, P_CP, PC_C, C_C, C_CP, PC_CP, PC_N, PC_NP, I_N, B_C, B_CP, B_NP = y

    # Equations
    # mRNAs of Per, Cry, and Bmal1
    dMP_dt = v_sP * hill(B_N, K_AP, n) - v_mP * mm(M_P, K_mP) - k_dmp * M_P
    dMC_dt = v_sC * hill(B_N, K_AC, n) - v_mC * mm(M_C, K_mC) - k_dmc * M_C

    # It appears as if there is a typo in the supplement, this should be the correct equation though
    dMB_dt = v_sB * hill(B_N, K_IB, m) - v_mB * mm(M_B, K_mB) - k_dmb * M_B
    #dMB_dt = v_sB * (K_IB**m/(K_IB**m+B_N**m)) - v_mB * mm(M_B, K_mB) - k_dmb * M_B

    # Phosphorylated and nonphosphorylated proteins PER and CRY in cytosol
    dPC_dt = k_sP * M_P - v_1p * mm(P_C, K_p) + v_2p * mm(P_CP, K_dp) + k4 * PC_C - k3 * P_C * C_C - k_dn * P_C
    dCC_dt = k_sC * M_C - v_1c * mm(C_C, K_p) + v_2c * mm(C_CP, K_dp) + k4 * PC_C - k3 * P_C * C_C - k_dnc * C_C
    dPCP_dt = v_1p * mm(P_C, K_p) - v_2p * mm(P_CP, K_dp) - v_dpc * mm(P_CP, K_d) - k_dn * P_CP
    dCCP_dt = v_1c * mm(C_C, K_p) - v_2c * mm(C_CP, K_dp) - v_dcc * mm(C_CP, K_d) - k_dn * C_CP

    # Phosphorylated and nonphosphorylated PER-CRY complex in cytosol and nucleus
    dPCC_dt = -v_1pc * mm(PC_C, K_p) + v_2pc * mm(PC_CP,
                                                  K_dp) - k4 * PC_C + k3 * P_C * C_C + k2 * PC_N - k1 * PC_C - k_dn * PC_C
    dPCN_dt = -v_3pc * mm(PC_N, K_p) + v_4pc * mm(PC_NP,
                                                  K_dp) - k2 * PC_N + k1 * PC_C - k7 * B_N * PC_N + k8 * I_N - k_dn * PC_N
    dPCCP_dt = v_1pc * mm(PC_C, K_p) - v_2pc * mm(PC_CP, K_dp) - v_dpcc * mm(PC_CP, K_d) - k_dn * PC_CP
    dPCNP_dt = v_3pc * mm(PC_N, K_p) - v_4pc * mm(PC_NP, K_dp) - v_dpcn * mm(PC_NP, K_d) - k_dn * PC_NP

    # Phosphorylated and nonphosphorylated protein BMAL1 in the cytosol and nucleus:
    dBC_dt = k_sb * M_B - v_1b * mm(B_C, K_p) + v_2b * mm(B_CP, K_dp) - k5 * B_C + k6 * B_N - k_dn * B_C
    dBCP_dt = v_1b * mm(B_C, K_p) - v_2b * mm(B_CP, K_dp) - v_dbc * mm(B_CP, K_d) - k_dn * B_CP
    dBN_dt = -v_3b * mm(B_N, K_p) + v_4b * mm(B_NP, K_dp) + k5 * B_C - k6 * B_N - k7 * B_N * PC_N + k8 * I_N - k_dn * B_N
    dBNP_dt = v_3b * mm(B_N, K_p) - v_4b * mm(B_NP, K_dp) - v_dbn * mm(B_NP, K_d) - k_dn * B_NP

    # Inactive complex between PER-CRY and CLOCK-BMAL1 in nucleus
    dIN_dt = -k8 * I_N + k7 * B_N * PC_N - v_din * mm(I_N, K_d) - k_dn * I_N

    return (dBN_dt, dMP_dt, dMC_dt, dMB_dt, dPC_dt, dPCP_dt, dPCC_dt,
            dCC_dt, dCCP_dt, dPCCP_dt, dPCN_dt, dPCNP_dt, dIN_dt,
            dBC_dt, dBCP_dt, dBNP_dt)

if __name__ == "__main__":
    t = np.linspace(0,100,1000)
    states = 16
    y0 = tuple([4]*states)
    Y = integrate.odeint(goldbetter_ode, y0, t)
    plt.plot(t, Y)
    plt.show()