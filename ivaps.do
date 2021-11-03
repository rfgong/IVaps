net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace
mata mata mlib index
ssc install ivreg2

prog def estimate_aps
	args pred C S delta nprocesses
	/*    
	pred: Do file
		Produce predictions from data
		(1) resampled continuous variables will have "_samp" appended to their names
		(2) must replace the variable "pred" with predictions
	C: str
		String with continuous variable names separated by spaces
	S: int
		Number of draws for each APS estimation
	delta: float
		Radius of sampling ball
	nprocesses: int
		Number of processes used to parallelize APS estimation
	Returns
	-----------
	APS stored in "aps" variable
	*/
	quietly {
		foreach varname in `C' {
			quietly su `varname'
			global sd_`varname' = r(sd)
			gen `varname'_samp = .
		}
		gen pred = .
		gen aps = 0
	}
	parallel initialize `nprocesses', f
	parallel, prog(estimate_aps_helper): estimate_aps_helper "`pred'" "`C'" `S' `delta'
end

prog def estimate_aps_helper
	args pred C S delta
	forval i=1/`S' {
		foreach varname in `C' {
			replace `varname'_samp = `varname' + runiform(-`delta',`delta')*${sd_`varname'}
		}
		do "`pred'.do"
		replace aps = aps + pred
	}
	replace aps = aps/`S'
end
