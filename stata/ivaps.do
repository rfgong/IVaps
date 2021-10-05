// EXAMPLE WORKFLOW


// Generate APS Values
/*
Option 1: Configure & Run gen_aps.py
Option 2: Configure & Run gen_aps.R
*/


// Load Data
use "data/safety_net_elig.dta", clear 


// 2SLS Example ivreg2: http://www.repec.org/bocode/i/ivreg2.html
// Include "aps if aps > 0 & aps <1" to use nondegenerate aps values
ivreg2 tot_con_sus2020_07_31 (safety_dollars_adj = safety_net) aps if aps > 0 & aps <1, first robust


// Covariate Balance Example mvreg: https://www.stata.com/manuals/mvmvreg.pdf
// Include "aps if aps > 0 & aps <1" to use nondegenerate aps values
mvreg ucc_per_bed profit_margin sum_pctg_ssi beds employees_FTEs occupancy operating_margin costs_per_discharge mdcd_net_rev = safety_net aps if aps > 0 & aps < 1
test safety_net
