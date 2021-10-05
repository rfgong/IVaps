"""
This example code computes APS for a cutoff rule on three variables.
The import requires ivaps.py to be placed in the same folder as this file.

Edit this file as needed.
"""
from ivaps import *

def predict(X):
    return ((X[:,0] >= 0.202) & (X[:,1] >= 25000) & (X[:,2] <= 0.03)).astype(int)

df = pd.read_stata("data/safety_net_elig.dta")
df["aps"] = estimate_aps(predict = predict, X = df[["sum_pctg_ssi_mdcd_days", "ucc_per_bed", "profit_margin"]], C = [0,1,2], S = 10000, delta = 0.05, nprocesses = 2)
df.to_stata("data/safety_net_elig.dta")
