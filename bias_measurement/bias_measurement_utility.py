import numpy as np

from scipy.optimize import curve_fit
import scipy.odr
import scipy.stats

def f_lin(x,A,B):
    return A*x+B

def f_1(x,B):
    return x+B

def f_wrapper_for_odr(beta, x): # parameter order for odr
    return f_lin(x, *beta)


'''
p_test takes two required (x,y) and two optional (label, p_th) arguments

label and p_th are optional. Ideally there should be a way to write
this function in a manner in which the order of p_th and label does not
matter, because they are optional variables. However, I do not have the
required python knowledge to write it in such a manner.

'''

def p_test(x,y,label='',p_th=0.05):
    parameters, cov= curve_fit(f_lin, x, y)

    model = scipy.odr.odrpack.Model(f_wrapper_for_odr)
    data = scipy.odr.odrpack.Data(x,y)
    myodr = scipy.odr.odrpack.ODR(data, model, beta0=parameters,  maxit=0)
    myodr.set_job(fit_type=2)
    parameterStatistics = myodr.run()
    df_e = len(x) - len(parameters) # degrees of freedom, error
    cov_beta = parameterStatistics.cov_beta # parameter covariance matrix from ODR
    sd_beta = parameterStatistics.sd_beta * parameterStatistics.sd_beta
    ci = []
    t_df = scipy.stats.t.ppf(0.975, df_e)
    ci = []
    for i in range(len(parameters)):
        ci.append([parameters[i] - t_df * parameterStatistics.sd_beta[i], parameters[i] + t_df * parameterStatistics.sd_beta[i]])

    # print('sd_beta = ', parameterStatistics.sd_beta)

    # Significance test for A == 1

    tstat_beta = (parameters[0]-1) / parameterStatistics.sd_beta[0] # coeff t-statistics
    pstat_beta_A = (1.0 - scipy.stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0    # coef. p-values

    # Significance test for B == 0

    tstat_beta = (parameters[1]) / parameterStatistics.sd_beta[1] # coeff t-statistics
    pstat_beta_B = (1.0 - scipy.stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0    # coef. p-values


    if pstat_beta_A < p_th or pstat_beta_B < p_th:

        if label =='':
            print("Y is a biased estimator of X")
        else:
            print("%s : biased."%label)
    else:
        if label =='':
            print("Y is an unbiased estimator of X")
        else:
            print("%s : unbiased."%label)
