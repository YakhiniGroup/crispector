from collections import defaultdict
import numpy as np
from math import exp
from numpy import exp, log1p
#from scipy.special import logsumexp
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
from crispector.utils.constants_and_types import TransDf, AlgResultDf, EDIT_PERCENT, TransResultDf, SITE_NAME, R_SITE, L_SITE, \
    FREQ, \
    TX_READ_NUM, MOCK_READ_NUM, SITE_A, SITE_B, TX_TRANS_READ, MOCK_TRANS_READ, TRANS_PVAL, TRANS_FDR, IS_TRANS, \
    TX_TRANS_BACKGROUND_READ, MOCK_TRANS_BACKGROUND_READ

import mpmath 
from mpmath import exp

def translocations_test(result_df, tx_df, mock_df, FDR_thresh: float,
                        edit_thresh: float):
    """

    :param result_df: NHEJ activity results
    :param tx_df: Treatment translocations
    :param mock_df: Mock translocations
    :param FDR_thresh: Threshold to count FDR value as active
    :param edit_thresh: Threshold to count NHEJ activity as active
    :return: TransResultDf - translocation result
    """
    if tx_df.shape[0] == 0:
        return pd.DataFrame(columns=[SITE_A, SITE_B, TX_TRANS_READ, MOCK_TRANS_READ, TRANS_PVAL, TRANS_FDR])

    trans_d = defaultdict(list)

    # Check translocation only for active editing site.
    active_df_t = result_df.loc[result_df[EDIT_PERCENT] >= edit_thresh]
    active_df = active_df_t.reset_index(drop=True)
    # Run HG test for all translocations
    for row_idx, row_a in active_df.iterrows():
        #print(row_a)
        for _, row_b in active_df.iloc[row_idx+1:].iterrows():
            site_a = row_a[SITE_NAME]
            site_b = row_b[SITE_NAME]

            # Get all translocation between site A & B
            tx_a_b_df = tx_df.loc[((tx_df[R_SITE] == site_a) & (tx_df[L_SITE] == site_b)) |
                                  ((tx_df[R_SITE] == site_b) & (tx_df[L_SITE] == site_a))]
            tx_trans = tx_a_b_df[FREQ].sum()

            # Make sure mock_df is empty with all the relevant fields
            if (R_SITE in mock_df.columns) and (L_SITE in mock_df.columns):
                mock_a_b_df = mock_df.loc[((mock_df[R_SITE] == site_a) & (mock_df[L_SITE] == site_b)) |
                                          ((mock_df[R_SITE] == site_b) & (mock_df[L_SITE] == site_a))]
                mock_trans = mock_a_b_df[FREQ].sum()
            else:
                mock_trans = 0

            # Get number for HG test
            tx_a_read_num = row_a[TX_READ_NUM]
            mock_a_read_num = row_a[MOCK_READ_NUM]
            tx_b_read_num = row_b[TX_READ_NUM]
            mock_b_read_num = row_b[MOCK_READ_NUM]

            # Tx & Mock background number are Geometric average of site A & B
            tx_read_num = (tx_a_read_num * tx_b_read_num) ** 0.5
            mock_read_num = (mock_a_read_num * mock_b_read_num) ** 0.5

            # Compute p-value
            b = tx_trans
            B = tx_trans + mock_trans
            n = tx_read_num
            N = tx_read_num + mock_read_num

            # If both sites have no translocations than move to the next couple
            if B == 0:
                continue

            #p_val = hypergeometric_cdf(b - 1, N, B, n)
            p_val=hypergeometric_cdf_log(b - 1, N, B, n) #changed tp mpmath implementation
            # Store values
            trans_d[SITE_A].append(site_a)
            trans_d[SITE_B].append(site_b)
            trans_d[TX_TRANS_READ].append(tx_trans)
            trans_d[MOCK_TRANS_READ].append(mock_trans)
            trans_d[TX_TRANS_BACKGROUND_READ].append(int(tx_read_num))
            trans_d[MOCK_TRANS_BACKGROUND_READ].append(int(mock_read_num))
            trans_d[TRANS_PVAL].append(p_val)

    # Compute FDR correction
    _, trans_d[TRANS_FDR] = fdrcorrection(trans_d[TRANS_PVAL])
    if len(trans_d) != 0:
        trans_df: TransResultDf = pd.DataFrame.from_dict(trans_d)
    else:
        trans_df = pd.DataFrame(columns=[SITE_A, SITE_B, TX_TRANS_READ, MOCK_TRANS_READ, TRANS_PVAL, TRANS_FDR])

    trans_df = trans_df.sort_values(by=[TRANS_FDR]).reset_index(drop=True)

    # Set Is_translocation column in tx_read and mock_read
    mock_df[IS_TRANS] = False
    tx_df[IS_TRANS] = False
    active_trans_df = trans_df.loc[trans_df[TRANS_FDR] < FDR_thresh]

    for _, row in active_trans_df.iterrows():
        site_a = row[SITE_A]
        site_b = row[SITE_B]
        tx_df.loc[((tx_df[R_SITE] == site_a) & (tx_df[L_SITE] == site_b)) |
                  ((tx_df[R_SITE] == site_b) & (tx_df[L_SITE] == site_a)), IS_TRANS] = True

    return trans_df


def hypergeom_pmf(k, M, n, N):
    '''
    implementing the probability mass function using mpmath package
    '''
    
    #k = number of translocations in Tx sample
    #M = total reads (Mock+Tx)
    #n = Total number of translocations (Tx + M)
    #N = total number of Tx reads 

    mock_n = M - n #number of mock reads 
    pmf = (mpmath.beta(n+1, 1) * mpmath.beta(mock_n +1,1) * mpmath.beta(M-N+1, N+1) /
          (mpmath.beta(k+1, n-k+1) * mpmath.beta(N-k+1,mock_n-N+k+1) * mpmath.beta(M+1, 1)))
    
    return pmf

def hypergeometric_cdf_log(k, M, n, N):
    """
      log pmf calculation using mpmath
    """
    max_val=min(N+1, n+1) #hypergemotric test upper limit
    k2 = np.arange(k + 1, max_val) 
    all_pmfs=[mpmath.log(hypergeom_pmf(k_ind,M,n,N)) for k_ind in k2] #calclulate the logpmf
 
    log_res = mpmath.log(mpmath.fsum(mpmath.exp(val) for val in all_pmfs)) #log sum all distribution values
    return exp(log_res)
