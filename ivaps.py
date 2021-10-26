"""Treatment estimation functions"""
from linearmodels.iv import IV2SLS
from linearmodels.system.model import SUR
from statsmodels.multivariate.multivariate_ols import _MultivariateOLS
from multiprocessing import Pool
import numpy as np
import pandas as pd

def estimate_aps(predict, X, C, S = 100, delta = 0.1, nprocesses = 1, chunksize = None):
    """Estimate APS for given dataset and prediction function

    Parameters
    -----------
    predict: function
        Function taking a 2D design matrix and returning a 1D vector of predictions
    X: array-like
        2D design matrix
    C: array-like
        Integer column indices for continuous variables
    S: int, default: 100
        Number of draws for each APS estimation
    delta: float, default: 0.1
        Radius of sampling ball
    nprocesses: int, default: 1
        Number of processes used to parallelize APS estimation
    chunksize: int, default: None
        Task chunk size used to parallelize APS estimation

    Returns
    -----------
    np.ndarray
        Array of estimated APS for each observation in sample

    Notes
    ------
    Approximate propensity score (APS) is the average predicted class for a given observation over :math:`S` samples.

    .. math::
        p^s(X_i;\\delta) = \\frac{1}{S} \\sum_{s=1}^{S} Predict(X_i^s)

    :math:`X_i^s` resamples continuous features of :math:`X_i`.
    Discrete and categorical features are kept unchanged from :math:`X_i`.

    Let :math:`X_{ci}` denote the vector of continuous features in :math:`X_i`.
    :math:`X_{ci}^s` is drawn uniformly at random from a ball centered at :math:`X_{ci}` with radius :math:`\\delta`.

    Continuous features are normalized to have mean zero and standard deviation one during the resampling step.
    This transformation in undone prior to prediction.
    """

    X = np.array(X)
    X_c = X[:, C].astype(float)
    c_std = np.std(X_c, axis=0)

    with Pool(processes=nprocesses) as pool:
        return sum(pool.starmap(estimate_aps_helper, [(i, delta, X_c, c_std, X, C, predict) for i in range(S)], chunksize=chunksize))/S

def estimate_aps_helper(i, delta, X_c, c_std, X, C, predict):
    # Resample continuous features
    dev = np.random.uniform(-delta, delta, X_c.shape)
    X_c_s = np.copy(X_c) + c_std * dev
    X_s = np.copy(X)
    X_s[:, C] = X_c_s
    return predict(X_s)

def estimate_treatment_effect(aps, Y, Z, D, W = None, saturated_aps = False, cov_type = "robust", weights = None):
    """Main treatment effect estimation function

    Parameters
    -----------
    aps: array-like
        Array of estimated APS values
    Y: array-like
        Array of outcome variables
    Z: array-like
        Array of treatment recommendations
    D: array-like
        Array of treatment assignments
    W: array-like, default: None
        Array of control variables
    saturated_aps: bool, default: False
        Convert APS variable into a full set of dummy variables
    cov_type: str, default: "robust"
        Covariance type of IV2SLS.
    weights: array-like, default: None
        Observation weights used in estimation

    Returns
    -----------
    tuple(IVResults, dict(D, dict(stat_label, value)))
        Tuple containing the fitted IV model and a dictionary containing the results for the treatment effect.

    Notes
    -----
    Treatment effect is estimated using IV estimation.

    .. math::
        D_i = \\gamma_0(1-I) + \\gamma_1 Z_i + \\gamma_2 p^s(X_i;\\delta) + \\gamma_3 W_i + v_i
    .. math::
        Y_i = \\beta_0(1-I) + \\beta_1 D_i + \\beta_2 p^s(X_i;\\delta) + \\beta_3 W_i + \\epsilon_i

    :math:`\\beta_1` is our causal estimate of the treatment effect. :math:`I` is an indicator for APS taking only a single nondegenerate value in the sample.
    """

    aps = np.array(aps)
    Y = np.array(Y)
    Z = np.array(Z)
    D = np.array(D)
    W = np.array(W)
    weights = np.array(weights)

    # Use only observations where aps is nondegenerate
    obs_tokeep = np.nonzero((aps > 0) & (aps < 1))
    print(f"We will fit on {len(obs_tokeep[0])} values out of {len(Y)} from the dataset for which the APS estimation is nondegenerate.")
    assert len(obs_tokeep[0]) > 0

    aps = aps[obs_tokeep[0]]
    Y = Y[obs_tokeep[0]]
    Z = Z[obs_tokeep[0]]
    D = D[obs_tokeep[0]]
    if W.any():
        W = W[obs_tokeep[0]]
    if weights.any():
        weights = weights[obs_tokeep[0]]

    cols = {"aps":aps, "Y":Y, "Z":Z, "D":D}
    exog = []

    # Check for single non-degeneracy
    constant = len(np.unique(aps)) == 1

    if len(W.shape) > 1:
        for i in range(W.shape[1]):
            cols["W"+str(i)] =  W[:,i]
            exog.append("W"+str(i))
            constant = (len(np.unique(W[:,i])) == 1) | constant
    elif len(W.shape) == 1:
        cols["W"] =  W
        exog.append("W")
        constant = (len(np.unique(W)) == 1) | constant

    # Add constant to specification if not provided
    if not constant:
        cols["const"] = np.ones(len(Y))
        exog.append("const")

    df = pd.DataFrame(cols)

    if saturated_aps:
        df["aps"] = df.aps.astype('category')
        dummy_df = pd.get_dummies(df.aps, prefix = "aps")
        aps_cols = list(dummy_df.columns)[1:]
        df = df.join(dummy_df[aps_cols])
        exog.extend(aps_cols)
    else:
        exog.append("aps")

    if weights.any():
        results = IV2SLS(df['Y'], df[exog], df['D'], df['Z'], weights = weights).fit(cov_type=cov_type)
    else:
        results = IV2SLS(df['Y'], df[exog], df['D'], df['Z']).fit(cov_type=cov_type)

    # Compile results
    res_dict = {"D":{}}
    res_dict["D"]['coef'] = results.params["D"]
    res_dict["D"]['stderr'] = results.std_errors["D"]
    res_dict["D"]['t'] = results.tstats["D"]
    res_dict["D"]['p'] = results.pvalues["D"]
    res_dict["D"]['n'] = results.nobs

    return results, res_dict

def covariate_balance_test(aps, X, Z, W = None, saturated_aps = False, cov_type = "robust"):
    """Covariate Balance Test

    Parameters
    -----------
    aps: array-like
        Array of estimated APS values
    X: array-like
        Array of covariates to test
    Z: array-like
        Array of treatment recommendations
    W: array-like, default: None
        Array of control variables
    saturated_aps: bool, default: False
        Convert APS variable into a full set of dummy variables
    cov_type: str, default: "robust"
        Covariance type of SUR.

    Returns
    -----------
    tuple(SystemResults, dict(X, dict(stat_label, value)))
        Tuple containing the fitted SUR model and a dictionary containing results of covariate balance estimation for each covariate as well as the joint hypothesis.

    Notes
    -----
    This function estimates a system of Seemingly Unrelated Regression (SUR) as defined in the linearmodels package.
    """

    aps = np.array(aps)
    X = np.array(X)
    Z = np.array(Z)
    W = np.array(W)

    # Use only observations where aps is nondegenerate
    obs_tokeep = np.nonzero((aps > 0) & (aps < 1))
    print(f"We will fit on {len(obs_tokeep[0])} values out of {len(X)} from the dataset for which the APS estimation is nondegenerate.")
    assert len(obs_tokeep[0]) > 0

    aps = aps[obs_tokeep[0]]
    X = X[obs_tokeep[0]]
    Z = Z[obs_tokeep[0]]
    if W.any():
        W = W[obs_tokeep[0]]

    cols = {"aps":aps, "Z":Z}

    dep = []
    if len(X.shape) > 1:
        for i in range(X.shape[1]):
            cols["X"+str(i)] =  X[:,i]
            dep.append("X"+str(i))
    else:
        cols["X"] =  X
        dep.append("X")

    exog = ["Z"]

    # Check for single non-degeneracy
    constant = len(np.unique(aps)) == 1

    if len(W.shape) > 1:
        for i in range(W.shape[1]):
            cols["W"+str(i)] =  W[:,i]
            exog.append("W"+str(i))
            constant = (len(np.unique(W[:,i])) == 1) | constant
    elif len(W.shape) == 1:
        cols["W"] =  W
        exog.append("W")
        constant = (len(np.unique(W)) == 1) | constant

    # Add constant to specification if not provided
    if not constant:
        cols["const"] = np.ones(len(aps))
        exog.append("const")

    df = pd.DataFrame(cols)

    if saturated_aps:
        df["aps"] = df.aps.astype('category')
        dummy_df = pd.get_dummies(df.aps, prefix = "aps")
        aps_cols = list(dummy_df.columns)[1:]
        df = df.join(dummy_df[aps_cols])
        exog.extend(aps_cols)
    else:
        exog.append("aps")

    # Covariate balance test
    mv_ols_res = SUR.multivariate_ls(df[dep], df[exog]).fit(cov_type=cov_type)

    # Joint hypothesis test: use multivariate_OLS from statsmodels
    # Edge case: single variable then joint test is the same as the original
    if len(dep) > 1:
        mv_ols_joint = _MultivariateOLS(df[dep], df[exog]).fit()
        L = np.zeros((1,len(exog)))
        L[:,0] = 1
        mv_test_res = mv_ols_joint.mv_test([("Z", L)])
    else:
        mv_test_res = None

    # Compile results
    res_dict = {}
    for x_var in dep:
        res_dict[x_var] = {}
        res_dict[x_var]['coef'] = mv_ols_res.params[f"{x_var}_Z"]
        res_dict[x_var]['stderr'] = mv_ols_res.std_errors[f"{x_var}_Z"]
        res_dict[x_var]['t'] = mv_ols_res.tstats[f"{x_var}_Z"]
        res_dict[x_var]['p'] = mv_ols_res.pvalues[f"{x_var}_Z"]
        res_dict[x_var]['n'] = int(mv_ols_res.nobs/len(dep))
    if mv_test_res:
        res_dict['joint'] = {}
        res_dict['joint']['p'] = mv_test_res.results['Z']['stat'].iloc[0, 4]
        res_dict['joint']['f'] = mv_test_res.results['Z']['stat'].iloc[0, 3]
    else:
        res_dict['joint'] = {}
        res_dict['joint']['p'] = mv_ols_res.pvalues[f"{dep[0]}_Z"]
        res_dict['joint']['t'] = mv_ols_res.tstats[f"{dep[0]}_Z"]

    return mv_ols_res, res_dict
