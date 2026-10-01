********************************************************************************
* RDROBUST STATA PACKAGE -- rdrobust
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Rocio Titiunik
********************************************************************************
*! version 11.1.1 01oct2026

capture program drop rdrobust
program define rdrobust, eclass
	version 16.0
	syntax anything [if] [in] [, c(real 0) fuzzy(string) deriv(real 0) p(string) q(real 0) h(string) b(string) rho(real 0) covs(string) covs_drop(string) kernel(string) weights(string) bwselect(string) vce(string) level(real 95) all scalepar(real 1) scaleregul(real 1) nochecks nowarnings masspoints(string) bwcheck(real 0) bwrestrict(string) stdvars(string) detail vleverage PRECision(string)]
	marksample touse
	capture mata: mata describe rdrobust_bw()
	if _rc quietly mata: mata mlib index

	* Snapshot Mata externals so we can drop only the variables WE create,
	* leaving the user's Mata workspace and rdrobust Mata functions untouched.
	mata: _mtx = direxternal("*"); st_local("_mata_before", rows(_mtx) ? invtokens(_mtx') : "")
	mata: mata drop _mtx

	* M1: do all work in an isolated temp frame instead of preserve+keep.
	* Only the touse=1 subset is copied (in-memory, no disk I/O), and the
	* user's frame is restored even on error or Ctrl-Break via nobreak.
	local _orig_frame `c(frame)'
	tempname _work_frame
	capture frame drop `_work_frame'
	frame put * if `touse', into(`_work_frame')

	nobreak {
	cwf `_work_frame'
	capture noisily {
	tokenize "`anything'"
	local y `1'
	local x `2'

	******************** Set PRECision *********************
	if ("`precision'"=="") local precision = "double"
	else {
		local precision = lower("`precision'")
		if ("`precision'"!="double" & "`precision'"!="single") {
			di as err `"precision(): incorrectly specified: options(single, double)"'
			exit 198
		}
	}
	local storage_type = "double"
	if ("`precision'"=="single") local storage_type = "float"

	local kernel   = lower("`kernel'")
	local bwselect = lower("`bwselect'")

	* Normalize the remaining string options here, before anything branches on
	* them. These used to be compared case-sensitively against lowercase
	* literals, so stdvars(ON) silently behaved as stdvars(off) -- a 43-fold
	* difference in the selected bandwidth on scaled data, with no warning.
	* Empty values are left empty so the DEFAULTS block below still fires.
	local masspoints = lower("`masspoints'")
	local stdvars    = lower("`stdvars'")
	local bwrestrict = lower("`bwrestrict'")
	local covs_drop  = lower("`covs_drop'")

	******************** Set VCE ***************************
	local nnmatch = 3
	local cr_method = ""
	tokenize `vce'
	local w : word count `vce'
	* Normalize the vce TYPE (first token) only. The remaining tokens are a
	* cluster variable name and an nnmatch count, and Stata variable names are
	* case sensitive, so they must not be lowercased. Without this, vce(NN)
	* failed with rc=7 while kernel() and bwselect() accepted any case.
	if `w' >= 1 local 1 = lower(`"`1'"')
	if `w' == 1 {
		local vce_select `"`1'"'
	}
	if `w' == 2 {
		local vce_select `"`1'"'
		if ("`vce_select'"=="nn") local nnmatch `"`2'"'
		if inlist("`vce_select'","cluster","nncluster","cr1","cr2","cr3","hc0","hc1","hc2","hc3") local clustvar `"`2'"'
	}
	if `w' == 3 {
		local vce_select `"`1'"'
		local clustvar   `"`2'"'
		local nnmatch    `"`3'"'
		if !inlist("`vce_select'","cluster","nncluster","cr1","cr2","cr3") {
			di as error  "{err}{cmd:vce()} incorrectly specified"
			exit 125
		}
	}
	if `w' > 3 {
		di as error "{err}{cmd:vce()} incorrectly specified"
		exit 125
	}

	* Disallow vce(nncluster ...): warn and shift to cr1 (default when clusters)
	if ("`vce_select'"=="nncluster") {
		di as text "Warning: vce(nncluster) is not supported. Switching to vce(cr1) (the default when clusters)."
		local vce_select = "cr1"
	}

	* With a cluster variable, map hc0/hc1/hc2/hc3 to cr1/cr1/cr2/cr3.
	* Per cluster_validation design: hc0+cluster is a silent remap to cr1
	* (the default); hc1/2/3 produce a warning so the user knows their
	* requested HC variant is being upgraded to its cluster analogue.
	if ("`clustvar'"!="") {
		if ("`vce_select'"=="hc0") local vce_select = "cr1"
		if ("`vce_select'"=="hc1") {
			di as text "Warning: vce(hc1 `clustvar') is not a cluster option. Switching to vce(cr1 `clustvar')."
			local vce_select = "cr1"
		}
		if ("`vce_select'"=="hc2") {
			di as text "Warning: vce(hc2 `clustvar') is not a cluster option. Switching to vce(cr2 `clustvar')."
			local vce_select = "cr2"
		}
		if ("`vce_select'"=="hc3") {
			di as text "Warning: vce(hc3 `clustvar') is not a cluster option. Switching to vce(cr3 `clustvar')."
			local vce_select = "cr3"
		}
		* bare vce(cluster clustvar) is equivalent to cr1 (default)
		if ("`vce_select'"=="cluster") local vce_select = "cr1"
	}

	* Preserve the raw vce option for e(vce_select) before internal remapping
	local vce_raw = "`vce_select'"
	if ("`vce_raw'"=="") local vce_raw = "nn"

	* Display label
	local vce_type = "NN"
	if ("`vce_select'"=="hc0")        local vce_type = "HC0"
	if ("`vce_select'"=="hc1")        local vce_type = "HC1"
	if ("`vce_select'"=="hc2")        local vce_type = "HC2"
	if ("`vce_select'"=="hc3")        local vce_type = "HC3"
	if ("`vce_select'"=="cr1")        local vce_type = "CR1"
	if ("`vce_select'"=="cr2")        local vce_type = "CR2"
	if ("`vce_select'"=="cr3")        local vce_type = "CR3"

	* Cluster / CR mapping
	if inlist("`vce_select'","cr1","cr2","cr3") {
		if ("`clustvar'"=="") {
			di as error "{err}{cmd:vce(`vce_select' clustervar)} requires a cluster variable"
			exit 125
		}
		local cluster = "cluster"
		if ("`vce_select'"=="cr1") local cr_method = "cr1"
		if ("`vce_select'"=="cr2") local cr_method = "crv2"
		if ("`vce_select'"=="cr3") local cr_method = "crv3"
		local vce_select = "hc0"
	}
	if ("`vce_select'"=="")              local vce_select = "nn"

	******************** Set BW ***************************
	tokenize `h'	
	local w : word count `h'
	if `w' == 1 {
		local h_l `"`1'"'
		local h_r `"`1'"'
	}
	if `w' == 2 {
		local h_l `"`1'"'
		local h_r `"`2'"'
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:h()} only accepts two inputs"  
		exit 125
	}
	
	tokenize `b'	
	local w : word count `b'
	if `w' == 1 {
		local b_l `"`1'"'
		local b_r `"`1'"'
	}
	if `w' == 2 {
		local b_l `"`1'"'
		local b_r `"`2'"'
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:b()} only accepts two inputs"  
		exit 125
	}
	
	* Transport buffers for the coefficient / variance matrices. These used to
	* be created as GLOBAL matrices literally named `b' and `V', which silently
	* destroyed any user matrices of those names (they were dropped again at the
	* end, so the originals were gone for good).
	tempname bmat Vmat

	*** Manual bandwidth
	if ("`h'"!="") {
		local bwselect = "Manual"
		if ("`b'"=="") {
			local b_r = `h_r'
			local b_l = `h_l'
		}		
		
		* Was `scalar rho = ...`, a GLOBAL scalar that clobbered any user scalar
		* named rho and was never dropped. The option value is used directly,
		* without rounding, so that b = h/rho holds for any positive rho.
		if (`rho' > 0)  {
			* rho() silently overrode an explicit b(). Say so.
			if ("`b'" != "") {
				di as text "Note: both b() and rho() were specified; rho() takes precedence and b() is ignored (b = h/rho)."
			}
			local b_l = `h_l'/`rho'
			local b_r = `h_r'/`rho'
		}
	}	
	
	*** Default bandwidth 
	if ("`h'"=="" & "`bwselect'"=="") local bwselect= "mserd"
	
	******************** Set Fuzzy***************************
	tokenize `fuzzy'	
	local w : word count `fuzzy'
	if `w' == 1 {
		local fuzzyvar `"`1'"'
	}
	if `w' == 2 {
		local fuzzyvar `"`1'"'
		local sharpbw  `"`2'"'
		if `"`2'"' != "sharpbw" {
			di as error  "{err}fuzzy() only accepts sharpbw as a second input" 
			exit 125
		}
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:fuzzy()} only accepts two inputs"  
		exit 125
	}
	
	**** DROP MISSINGS **********************************************
	if ("`covs'"!="") {
		qui ds `covs', alpha
		local covs_list = r(varlist)
		local ncovs: word count `covs_list'
	}
	* Combine all missing-value drops into a single scan.
	local drop_cond "mi(`y') | mi(`x')"
	if ("`fuzzy'"!="")   local drop_cond "`drop_cond' | mi(`fuzzyvar')"
	if ("`cluster'"!="") local drop_cond "`drop_cond' | mi(`clustvar')"
	foreach z of local covs_list {
		local drop_cond "`drop_cond' | mi(`z')"
	}
	if ("`weights'"!="") local drop_cond "`drop_cond' | mi(`weights') | `weights'<=0"
	qui drop if `drop_cond'

	**** Convert string clustvar to numeric ************************
	* vce(cluster var) traditionally accepts string/categorical clusters in
	* Stata. rdrobust's Mata path uses st_data(), which is numeric-only, so
	* map a string clustvar to an integer factor via egen group() in a
	* tempvar used only by the Mata call. The original `clustvar' is kept
	* unchanged so display and ereturn report the user's variable name.
	local clustvar_num "`clustvar'"
	if ("`cluster'"!="") {
		cap confirm numeric variable `clustvar'
		if (_rc) {
			tempvar _clustvar_num
			qui egen `storage_type' `_clustvar_num' = group(`clustvar')
			local clustvar_num "`_clustvar_num'"
		}
	}

	**** CHECK colinearity ******************************************
	local covs_drop_coll = 0	
	if ("`covs_drop'"=="") local covs_drop = "pinv"	
	if ("`covs'"!="") {	
		
	if ("`covs_drop'"=="invsym")  local covs_drop_coll = 1
	if ("`covs_drop'"=="pinv")    local covs_drop_coll = 2
	
	if ("`covs_drop'"!="off") {	
		
		qui _rmcoll `covs_list'
		local nocoll_controls_cat `r(varlist)'
		local nocoll_controls ""
		foreach myString of local nocoll_controls_cat {
			if ~strpos("`myString'", "o."){
				if ~strpos("`myString'", "MYRUNVAR"){
					local nocoll_controls "`nocoll_controls' `myString'"
				}
				}
			}			
		local covs_new `nocoll_controls'
		qui ds `covs_new', alpha
		local covs_list_new = r(varlist)
		local ncovs_new: word count `covs_list_new'
		
		if (`ncovs_new'<`ncovs') {
				local ncovs = "`ncovs_new'"
				local covs_list = "`covs_list_new'"
				di as error  "{err}Multicollinearity issue detected in {cmd:covs}. Redundant covariates were removed." 				
			}	
		}	
	}
	
				
	**** DEFAULTS ***************************************
	if ("`masspoints'"=="") local masspoints = "adjust"
	if ("`stdvars'"=="")    local stdvars    = "on"
	if ("`bwrestrict'"=="") local bwrestrict = "on"

	* Validate the on/off-style options against their whitelists. Previously an
	* unrecognized value fell through to the "not on" branch and was silently
	* treated as off, which quietly defeated the scale-robustness default.
	if !inlist("`masspoints'","adjust","check","off") {
		di as error "{err}{cmd:masspoints()} incorrectly specified (received '`masspoints''); allowed: adjust, check, off."
		exit 198
	}
	if !inlist("`stdvars'","on","off") {
		di as error "{err}{cmd:stdvars()} incorrectly specified (received '`stdvars''); allowed: on, off."
		exit 198
	}
	if !inlist("`bwrestrict'","on","off") {
		di as error "{err}{cmd:bwrestrict()} incorrectly specified (received '`bwrestrict''); allowed: on, off."
		exit 198
	}
	if !inlist("`covs_drop'","off","invsym","pinv") {
		di as error "{err}{cmd:covs_drop()} incorrectly specified (received '`covs_drop''); allowed: off, invsym, pinv."
		exit 198
	}
	*****************************************************************
	
	* Only compute what we need. `su x, d` also computes skewness/kurtosis/etc.
	qui su `x'
	local N     = r(N)
	local x_min = r(min)
	local x_max = r(max)
	local x_sd  = r(sd)
	qui _pctile `x', p(25 75)
	local x_iq = r(r2) - r(r1)
	local range_l = abs(`c'-`x_min')
	local range_r = abs(`x_max'-`c')
	
	if (`deriv'>0 & "`p'"=="" & `q'==0) local p = `deriv'+1
	if ("`p'"=="")  local p = 1
	if ("`q'"=="0") local q = `p'+1

	**************************** BEGIN ERROR CHECKING ************************************************
	if ("`checks'"=="") {
			if (`c'<=`x_min' | `c'>=`x_max'){
			 di as error  "{err}{cmd:c()} should be set within the range of `x'"  
			 exit 125
			}
						
			if (`N'<20 & "`h'"=="") {
			 di as error "{err}Not enough observations to perform bandwidth calculations. Using the maximum distance from the cutoff for h."
			 local bwselect = "Manual"
			 local bw_range = max(`range_l',`range_r')
			 local h   = `bw_range'
			 local h_l = `bw_range'
			 local h_r = `bw_range'
			 if (`rho'>0) {
			   local b_l = `bw_range'/`rho'
			   local b_r = `bw_range'/`rho'
			 }
			 else if ("`b'"=="") {
			   local b_l = `bw_range'
			   local b_r = `bw_range'
			 }
			}
			
			if (!inlist("`kernel'","uni","uniform","tri","triangular","epa","epanechnikov","")){
			 di as error  "{err}{cmd:kernel()} incorrectly specified"
			 exit 7
			}

			if (inlist("`bwselect'","CCT","IK","CV","cct","ik","cv")){
				di as error  "{err}{cmd:bwselect()} options IK, CCT and CV have been deprecated. Please see help for new options"
				exit 7
			}

			if (!inlist("`bwselect'","mserd","msetwo","msesum","msecomb1","msecomb2") & !inlist("`bwselect'","cerrd","certwo","cersum","cercomb1","cercomb2") & "`bwselect'"!="Manual"){
				di as error  "{err}{cmd:bwselect()} incorrectly specified"
				exit 7
			}

			if (!inlist("`vce_select'","nn","","cluster","nncluster") & !inlist("`vce_select'","cr1","cr2","cr3") & !inlist("`vce_select'","hc0","hc1","hc2","hc3")){
			 di as error  "{err}{cmd:vce()} incorrectly specified"
			 exit 7
			}

			if (`p'<0 | `q'<=0 | `deriv'<0 | `nnmatch'<=0){
			 di as error  "{err}{cmd:p()}, {cmd:q()}, {cmd:deriv()}, {cmd:nnmatch()} should be positive"
			 exit 411
			}

			if (`p'>=`q' & `q'>0){
			 di as error  "{err}{cmd:q()} should be higher than {cmd:p()}"
			 exit 125
			}

			if (`deriv'>`p' & `deriv'>0){
			 di as error  "{err}{cmd:deriv()} can not be higher than {cmd:p()}"
			 exit 125
			}

			if (`p'>0) {
				local p_round = round(`p')/`p'
				local q_round = round(`q')/`q'
				local d_round = round(`deriv'+1)/(`deriv'+1)
				local m_round = round(`nnmatch')/`nnmatch'

				if (`p_round'!=1 | `q_round'!=1 |`d_round'!=1 |`m_round'!=1 ){
				 di as error  "{err}{cmd:p()}, {cmd:q()}, {cmd:deriv()} and {cmd:nnmatch()} should be integers"  
				 exit 126
				}
			}
			if (`level'>=100 | `level'<=0){
			 di as error  "{err}{cmd:level()} should be a number in (0, 100)"
			 exit 125
			}
			if (`rho' < 0) {
			 di as error  "{err}{cmd:rho()} should be strictly positive"
			 exit 125
			}
			if (`bwcheck' < 0 | (`bwcheck' != round(`bwcheck'))) {
			 di as error  "{err}{cmd:bwcheck()} must be a non-negative integer (0 = unset)"
			 exit 125
			}
			if ("`masspoints'" != "" & ///
			    !inlist("`masspoints'", "check", "adjust", "off")) {
			 di as error  "{err}{cmd:masspoints()} must be one of check, adjust, off"
			 exit 125
			}
	}

	* Fail early, and say why, when one side cannot support the polynomial fits.
	* Done here, in the validation block, because an exit from inside the Mata
	* work blocks does not surface its own return code (see the note below).
	mata: _x0 = st_data(., ("`y' `x'"), 0)[,2]; st_local("_M0_l", strofreal(rows(uniqrows(select(_x0, _x0:<`c'))))); st_local("_M0_r", strofreal(rows(uniqrows(select(_x0, _x0:>=`c')))))
	mata: mata drop _x0
	if (min(`_M0_l', `_M0_r') < `q'+1) {
		local _side = cond(`_M0_l' < `q'+1, "left", "right")
		di as error "{err}Not enough distinct running-variable values on the `_side' side of the cutoff (" min(`_M0_l', `_M0_r') ") to fit a polynomial of order q = `q'."
		exit 2001
	}
	*********************** END ERROR CHECKING ************************************************************
	}
	* End of validation-only capture block. Splitting validation from the
	* Mata-work block sidesteps a Stata quirk where `exit N` inside an `if{}`
	* inside `capture noisily { ... }` does NOT halt the noisily block when
	* later `mata { ... }` blocks exist at the same scope; control jumps past
	* the first mata block and lands in the second, leaking transient Mata
	* externals and surfacing a misleading rc=3499 "X_l not found" instead of
	* the real validation error.
	local _rc = _rc
	if `_rc' {
		cwf `_orig_frame'
		frame drop `_work_frame'
		exit `_rc'
	}
	capture noisily {
	if ("`vce_select'"=="nn" | "`masspoints'"=="check" | "`masspoints'"=="adjust") {
		sort `x', stable
		if ("`vce_select'"=="nn") {
			tempvar dups dupsid
			* Use the TEMPVAR macros, not the literal names: `gen dups = _N`
			* created a permanent variable called `dups`, so a user variable of
			* that name broke every default vce(nn) run (and the rc=110 was
			* masked into a misleading rc=3499 by the capture-noisily+mata path).
			by `x': gen `storage_type' `dups' = _N
			by `x': gen `storage_type' `dupsid' = _n
		}
	}

	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
		local kernel_type = "Epanechnikov"
		local C_c = 2.34
	}
	else if ("`kernel'"=="uniform" | "`kernel'"=="uni") {
		local kernel_type = "Uniform"
		local C_c = 1.843
	}
	else {
		local kernel_type = "Triangular"
		local C_c = 2.576
	}
	
	*** Start MATA ********************************************************

	mata{
	
	*** Preparing data
		YX = st_data(., ("`y' `x'"), 0)
		Y = YX[,1]; X = YX[,2]
		ind_l = selectindex(X:<`c'); ind_r = selectindex(X:>=`c')
		X_l = X[ind_l];	X_r = X[ind_r]
		Y_l = Y[ind_l];	Y_r = Y[ind_r]
		dZ=dT=dC=Z_l=Z_r=T_l=T_r=C_l=C_r=fw_l=fw_r=g_l=g_r=dups_l=dups_r=dupsid_l=dupsid_r=g_l=g_r=eT_l=eT_r=eZ_l=eZ_r=indC_l=indC_r=eC_l=eC_r=0
		
		N   = length(X);	N_l = length(X_l);	N_r = length(X_r)

		// Unique-value counts per side, computed HERE because this block is
		// common to the manual-h and auto-bandwidth paths. They used to be set
		// only inside the auto-bandwidth branch, so h() + masspoints(check)
		// crashed rc=111 at the display -- or, after an earlier run in the same
		// session, silently displayed STALE counts from that run.
		st_numscalar("M_l", length(uniqrows(X_l)))
		st_numscalar("M_r", length(uniqrows(X_r)))
				
		if ("`covs'"!="") {
			Z   = st_data(.,tokens("`covs_list'"), 0); dZ  = cols(Z)
			Z_l = Z[ind_l,];	Z_r = Z[ind_r,]
		}
	
		if ("`fuzzy'"!="") {
			T = st_data(.,("`fuzzyvar'"), 0);	T_l = T[ind_l];	T_r = T[ind_r]; dT = 1
			// Reject fully degenerate first stage (no variation, no jump).
			// One-sided non-compliance falls through to the perf_comp branch.
			if (variance(T_l)==0 & variance(T_r)==0 & abs(mean(T_l) - mean(T_r)) < sqrt(epsilon(1))) {
				_error("Fuzzy RD: first-stage variable has no variation and no jump at the cutoff. The fuzzy estimator is not identified.")
			}
			if (variance(T_l)==0 | variance(T_r)==0){
				T_l = T_r = 0
				st_local("perf_comp","perf_comp")
			}
			if ("`sharpbw'"!=""){
				T_l = T_r = 0
				st_local("sharpbw","sharpbw")
			}
		}
	
		if ("`cluster'"!="") {
			C  = st_data(.,("`clustvar_num'"), 0)
			C_l  = C[ind_l]; C_r  = C[ind_r]
			indC_l = order(C_l,1);  indC_r = order(C_r,1) 
			g_l = rows(panelsetup(C_l[indC_l],1));	g_r = rows(panelsetup(C_r[indC_r],1))
			st_numscalar("g_l",  g_l);     st_numscalar("g_r",   g_r)
		}	
	
		if ("`weights'"!="") {
			fw = st_data(.,("`weights'"), 0)
			fw_l = fw[ind_l];	fw_r = fw[ind_r]
		}
		
		if ("`vce_select'"=="nn") {
			dups      = st_data(.,("`dups'"), 0); dupsid    = st_data(.,("`dupsid'"), 0)
			dups_l    = dups[ind_l];    dups_r    = dups[ind_r]
			dupsid_l  = dupsid[ind_l];  dupsid_r  = dupsid[ind_r]
		}
		
		
		h_l = `h_l'
		h_r = `h_r'
		b_l = `b_l'
		b_r = `b_r'

	***********************************************************************
	******** Computing bandwidth selector *********************************
	***********************************************************************		
masspoints_found = 0
	
	if ("`h'"=="") {	
	
		BWp = min((`x_sd',`x_iq'/1.349))
		x_sd = y_sd = 1
		c = `c'
		*** Starndardized ******************
		if  ("`stdvars'"=="on")  {	
			y_sd = sqrt(variance(Y))
			x_sd = sqrt(variance(X))
			X_l = X_l/x_sd;	X_r = X_r/x_sd
			Y_l = Y_l/y_sd;	Y_r = Y_r/y_sd
			c = `c'/x_sd
			BWp = min((1, (`x_iq'/x_sd)/1.349))
		}		
		x_l_min = min(X_l);	x_l_max = max(X_l)
		x_r_min = min(X_r);	x_r_max = max(X_r)
	
		range_l = c - x_l_min
		range_r = x_r_max - c
		************************************		
		
		mN = `N'
		bwcheck = `bwcheck'
		covs_drop_coll = `covs_drop_coll'

		// Always compute unique-value vectors so the bwcheck and de-standardize
		// paths below have them available (even when masspoints=="off").
		X_uniq_l = sort(uniqrows(X_l),-1)
		X_uniq_r = uniqrows(X_r)
		M_l = length(X_uniq_l)
		M_r = length(X_uniq_r)
		M = M_l + M_r
		if ("`masspoints'"=="check" | "`masspoints'"=="adjust") {
			st_numscalar("M_l", M_l); st_numscalar("M_r", M_r)
			mass_l = 1-M_l/N_l
			mass_r = 1-M_r/N_r
			if (mass_l>=0.2 | mass_r>=0.2){
				masspoints_found = 1
				display("{err}Mass points detected in the running variable.")
				if ("`masspoints'"=="adjust" & "`bwcheck'"=="0") bwcheck = 10
				if ("`masspoints'"=="check") display("{err}Try using option {cmd:masspoints(adjust)}.")
			}
		}
				
		c_bw = `C_c'*BWp*mN^(-1/5)
		if ("`masspoints'"=="adjust") c_bw = `C_c'*BWp*M^(-1/5)
		if  ("`bwrestrict'"=="on") {
			bw_max = max((range_l,range_r))
			c_bw = min((c_bw, bw_max))
		}
		if (bwcheck > 0) {
			bwcheck_l = min((bwcheck, M_l))
			bwcheck_r = min((bwcheck, M_r))
			bw_min_l = abs(X_uniq_l:-c)[bwcheck_l]*(1+sqrt(epsilon(1)))
			bw_min_r = abs(X_uniq_r:-c)[bwcheck_r]*(1+sqrt(epsilon(1)))
			c_bw = max((c_bw, bw_min_l, bw_min_r))
		}		
		
			
		// T1: per-side V-fit caches reused across all pilot calls.
		vcache_l = asarray_create("string")
		vcache_r = asarray_create("string")
		// Set when any pilot fit is not identified (rdrobust_bw returns missing).
		bw_undef = 0

		*** Step 1: d_bw
		C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`q'+1, nu=`q'+1, o_B=`q'+2, h_V=c_bw, h_B=range_l*(1+sqrt(epsilon(1))), 0, "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
		bw_undef = max((bw_undef, hasmissing(C_d_l[1..3])))
		C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`q'+1, nu=`q'+1, o_B=`q'+2, h_V=c_bw, h_B=range_r*(1+sqrt(epsilon(1))), 0, "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
		bw_undef = max((bw_undef, hasmissing(C_d_r[1..3])))
		if (C_d_l[1]==0 | C_d_l[2]==0 | C_d_r[1]==0 | C_d_r[2]==0 |C_d_l[1]==. | C_d_l[2]==. | C_d_l[3]==. |C_d_r[1]==. | C_d_r[2]==. | C_d_r[3]==.) printf("{err}Not enough variability to compute the preliminary bandwidth. Consider using option {cmd:stdvars(on)} (now the default) to standardize the running variable before bandwidth selection; or check for mass points with {cmd:masspoints(check)}.\n")
	
		*** BW-TWO
		if  ("`bwselect'"=="msetwo" |  "`bwselect'"=="certwo" | "`bwselect'"=="msecomb2" | "`bwselect'"=="cercomb2" )  {		
			* Preliminar bw
			d_bw_l = (  (C_d_l[1]              /   C_d_l[2]^2))^C_d_l[4]
			d_bw_r = (  (C_d_r[1]              /   C_d_r[2]^2))^C_d_r[4]
			if  ("`bwrestrict'"=="on") {
			d_bw_l = min((d_bw_l, range_l))
			d_bw_r = min((d_bw_r, range_r))
			}
			if (bwcheck > 0) {
				d_bw_l = max((d_bw_l, bw_min_l))
				d_bw_r = max((d_bw_r, bw_min_r))
			}
			* Bias bw
			C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_l, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_b_l[1..3])))
			b_bw_l = (  (C_b_l[1]              /   (C_b_l[2]^2 + `scaleregul'*C_b_l[3])))^C_b_l[4]
			C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_r, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_b_r[1..3])))
			b_bw_r = (  (C_b_r[1]              /   (C_b_r[2]^2 + `scaleregul'*C_b_r[3])))^C_b_r[4]
			if  ("`bwrestrict'"=="on") {
			b_bw_l = min((b_bw_l, range_l))
			b_bw_r = min((b_bw_r, range_r))
			}
			* Main bw
			C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_l, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_h_l[1..3])))
			h_bw_l = (  (C_h_l[1]              /   (C_h_l[2]^2 + `scaleregul'*C_h_l[3])))^C_h_l[4]
			C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_r, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_h_r[1..3])))
			h_bw_r = (  (C_h_r[1]              /   (C_h_r[2]^2 + `scaleregul'*C_h_r[3])))^C_h_r[4]
			if  ("`bwrestrict'"=="on") {
			h_bw_l = min((h_bw_l, range_l))
			h_bw_r = min((h_bw_r, range_r))
			}
		}
		
		*** BW-SUM
		if  ("`bwselect'"=="msesum" | "`bwselect'"=="cersum" |  "`bwselect'"=="msecomb1" | "`bwselect'"=="msecomb2" |  "`bwselect'"=="cercomb1" | "`bwselect'"=="cercomb2")  {
			* Preliminar bw
			d_bw_s = ( ((C_d_l[1] + C_d_r[1])  /  (C_d_r[2] + C_d_l[2])^2))^C_d_l[4]
			if  ("`bwrestrict'"=="on")  d_bw_s = min((d_bw_s, bw_max))
			if (bwcheck > 0) d_bw_s = max((d_bw_s, bw_min_l, bw_min_r))		
			* Bias bw
			C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_s, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_b_l[1..3])))
			C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_s, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_b_r[1..3])))
			b_bw_s = ( ((C_b_l[1] + C_b_r[1])  /  ((C_b_r[2] + C_b_l[2])^2 + `scaleregul'*(C_b_r[3]+C_b_l[3]))))^C_b_l[4]
			if  ("`bwrestrict'"=="on") b_bw_s = min((b_bw_s, bw_max))
			* Main bw
			C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_s, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_h_l[1..3])))
			C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_s, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_h_r[1..3])))
			h_bw_s = ( ((C_h_l[1] + C_h_r[1])  /  ((C_h_r[2] + C_h_l[2])^2 + `scaleregul'*(C_h_r[3] + C_h_l[3]))))^C_h_l[4]
			if  ("`bwrestrict'"=="on") h_bw_s = min((h_bw_s, bw_max))
		}
		
		*** RD
		if  ("`bwselect'"=="mserd" | "`bwselect'"=="cerrd" | "`bwselect'"=="msecomb1" | "`bwselect'"=="msecomb2" | "`bwselect'"=="cercomb1" | "`bwselect'"=="cercomb2" | "`bwselect'"=="") {
			* Preliminar bw
			d_bw_d = ( ((C_d_l[1] + C_d_r[1])  /  (C_d_r[2] - C_d_l[2])^2))^C_d_l[4]
			if  ("`bwrestrict'"=="on") d_bw_d = min((d_bw_d, bw_max))
			
			if (bwcheck > 0) d_bw_d = max((d_bw_d, bw_min_l, bw_min_r))		
			* Bias bw
			C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_d, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_b_l[1..3])))
			C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`q', nu=`p'+1, o_B=`q'+1, h_V=c_bw, h_B=d_bw_d, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_b_r[1..3])))
			b_bw_d = ( ((C_b_l[1] + C_b_r[1])  /  ((C_b_r[2] - C_b_l[2])^2 + `scaleregul'*(C_b_r[3] + C_b_l[3]))))^C_b_l[4]
			if  ("`bwrestrict'"=="on") b_bw_d = min((b_bw_d, bw_max))
			
			* Main bw
			C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_d, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_l, dupsid_l, covs_drop_coll, "`cr_method'", vcache_l)
			bw_undef = max((bw_undef, hasmissing(C_h_l[1..3])))
			C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=`p', nu=`deriv', o_B=`q', h_V=c_bw, h_B=b_bw_d, `scaleregul', "`vce_select'", `nnmatch', "`kernel'", dups_r, dupsid_r, covs_drop_coll, "`cr_method'", vcache_r)
			bw_undef = max((bw_undef, hasmissing(C_h_r[1..3])))
			h_bw_d = ( ((C_h_l[1] + C_h_r[1])  /  ((C_h_r[2] - C_h_l[2])^2 + `scaleregul'*(C_h_r[3] + C_h_l[3]))))^C_h_l[4]
			if  ("`bwrestrict'"=="on") h_bw_d = min((h_bw_d, bw_max))
			
		}	
		


		if (C_b_l[1]==0 | C_b_l[2]==0 | C_b_r[1]==0 | C_b_r[2]==0 |C_b_l[1]==. | C_b_l[2]==. | C_b_l[3]==. | C_b_r[1]==. | C_b_r[2]==. | C_b_r[3]==.) printf("{err}Not enough variability to compute the bias bandwidth (b). Consider using option {cmd:stdvars(on)} (now the default) to standardize the running variable; or check for mass points with {cmd:masspoints(check)}.\n")
		if (C_h_l[1]==0 | C_h_l[2]==0 | C_h_r[1]==0 | C_h_r[2]==0 |C_h_l[1]==. | C_h_l[2]==. | C_h_l[3]==. | C_h_r[1]==. | C_h_r[2]==. | C_h_r[3]==.) printf("{err}Not enough variability to compute the loc. poly. bandwidth (h). Consider using option {cmd:stdvars(on)} (now the default) to standardize the running variable; or check for mass points with {cmd:masspoints(check)}.\n")
		// Stopped at the top of the estimation block below.
		if (bw_undef) st_local("_bw_undef", "1")
	
		cer_h = mN^(-(`p'/((3+`p')*(3+2*`p'))))
		if ("`cluster'"!="") cer_h = (g_l+g_r)^(-(`p'/((3+`p')*(3+2*`p'))))
		cer_b = 1	
			
		if  ("`bwselect'"=="mserd" | "`bwselect'"=="cerrd" | "`bwselect'"=="msecomb1" | "`bwselect'"=="msecomb2" | "`bwselect'"=="cercomb1" | "`bwselect'"=="cercomb2") {
			h_l = h_r = h_mserd = x_sd*h_bw_d
			b_l = b_r = b_mserd = x_sd*b_bw_d
		}	
		if  ("`bwselect'"=="msesum" | "`bwselect'"=="cersum" |  "`bwselect'"=="msecomb1" | "`bwselect'"=="msecomb2" |  "`bwselect'"=="cercomb1" | "`bwselect'"=="cercomb2")  {
			h_l = h_r = h_msesum = x_sd*h_bw_s
			b_l = b_r = b_msesum = x_sd*b_bw_s
		}
		if  ("`bwselect'"=="msetwo" |  "`bwselect'"=="certwo" | "`bwselect'"=="msecomb2" | "`bwselect'"=="cercomb2")  {		
			h_l = h_msetwo_l = x_sd*h_bw_l
			h_r = h_msetwo_r = x_sd*h_bw_r
			b_l = b_msetwo_l = x_sd*b_bw_l
			b_r = b_msetwo_r = x_sd*b_bw_r
		}
		if  ("`bwselect'"=="msecomb1" | "`bwselect'"=="cercomb1") {
			h_l = h_r = h_msecomb1 = min((h_mserd,h_msesum))
			b_l = b_r = b_msecomb1 = min((b_mserd,b_msesum))
		}
		if  ("`bwselect'"=="msecomb2" | "`bwselect'"=="cercomb2") {
			h_l = (sort((h_mserd,h_msesum,h_msetwo_l)',1))[2]
			h_r = (sort((h_mserd,h_msesum,h_msetwo_r)',1))[2]
			b_l = (sort((b_mserd,b_msesum,b_msetwo_l)',1))[2]
			b_r = (sort((b_mserd,b_msesum,b_msetwo_r)',1))[2]
		}		
		if  ("`bwselect'"=="cerrd" | "`bwselect'"=="cersum" | "`bwselect'"=="certwo" | "`bwselect'"=="cercomb1" | "`bwselect'"=="cercomb2"){
			h_l = h_l*cer_h
			h_r = h_r*cer_h
			b_l = b_l*cer_b
			b_r = b_r*cer_b
		}
		
		rho = `rho'
		if (rho>0)  {
			b_l = h_l/rho
			b_r = h_r/rho
		}
					
		*** De-Starndardized *********************************
		c = `c'*x_sd
		X_uniq_l = X_uniq_l*x_sd
		X_uniq_r = X_uniq_r*x_sd
		X_l = X_l*x_sd;	X_r = X_r*x_sd
		Y_l = Y_l*y_sd;	Y_r = Y_r*y_sd
		range_l = range_l*x_sd
		range_r = range_r*x_sd  
		*****************************************************

		
		} /* close if for bw selector */
	
	}
	

	mata{
	
		if (st_local("_bw_undef")=="1") {
			display("{err}Not enough variability in the running variable to compute the bandwidth. Check for mass points with option {cmd:masspoints(check)}; if the running variable is discrete, an RD design may not be identified at this sample size.")
			exit(2001)
		}

		*** Estimation and Inference
		
		c = strtoreal("`c'")
	
		w_h_l = rdrobust_kweight(X_l,`c',h_l,"`kernel'");	w_h_r = rdrobust_kweight(X_r,`c',h_r,"`kernel'")
		w_b_l = rdrobust_kweight(X_l,`c',b_l,"`kernel'");	w_b_r = rdrobust_kweight(X_r,`c',b_r,"`kernel'")

		// Cluster-robust variances need many clusters; with a handful per side
		// they are unreliable, and with p+1 or fewer the variance is not
		// identified (the standard error collapses to zero up to rounding).
		if ("`cluster'"!="" & st_local("warnings")=="") {
			gh_l = rows(uniqrows(select(C_l, w_h_l:>0))); gh_r = rows(uniqrows(select(C_r, w_h_r:>0)))
			if (min((gh_l, gh_r)) <= `p'+1) printf("{txt}Warning: only %g (left) and %g (right) clusters within the bandwidth. With p+1 = %g or fewer clusters on a side the cluster-robust variance is not identified and the standard error can be zero.\n", gh_l, gh_r, `p'+1)
			else if (min((gh_l, gh_r)) < 10) printf("{txt}Warning: only %g (left) and %g (right) clusters within the bandwidth. Cluster-robust standard errors are unreliable with fewer than 10 clusters on a side.\n", gh_l, gh_r)
		}
		
		if ("`weights'"!="") {
			w_h_l = fw_l:*w_h_l;	w_h_r = fw_r:*w_h_r
			w_b_l = fw_l:*w_b_l;	w_b_r = fw_r:*w_b_r			
		}
		
		ind_h_l = selectindex(w_h_l:> 0);		ind_h_r = selectindex(w_h_r:> 0)
		ind_b_l = selectindex(w_b_l:> 0);		ind_b_r = selectindex(w_b_r:> 0)
		N_h_l = length(ind_h_l);	N_b_l = length(ind_b_l)
		N_h_r = length(ind_h_r);	N_b_r = length(ind_b_r)

		if (rows(uniqrows(sort(X_l[ind_h_l],1))) < `p'+1 |
		    rows(uniqrows(sort(X_r[ind_h_r],1))) < `p'+1 |
		    rows(uniqrows(sort(X_l[ind_b_l],1))) < `q'+1 |
		    rows(uniqrows(sort(X_r[ind_b_r],1))) < `q'+1) {
			display("{err}Not enough distinct running-variable values with positive weight to fit the requested polynomials on each side of the cutoff.")
			exit(2001)
		}
		
		if (N_h_l<10 | N_h_r<10 | N_b_l<10 | N_b_r<10){
		 display("{err}Estimates might be unreliable due to low number of effective observations.")
		}
		
		ind_l = ind_b_l; ind_r = ind_b_r
		if (h_l>b_l) ind_l = ind_h_l   
		if (h_r>b_r) ind_r = ind_h_r   
		eN_l = length(ind_l); eN_r = length(ind_r)
		eY_l  = Y_l[ind_l];	eY_r  = Y_r[ind_r]
		eX_l  = X_l[ind_l];	eX_r  = X_r[ind_r]
		W_h_l = w_h_l[ind_l];	W_h_r = w_h_r[ind_r]
		W_b_l = w_b_l[ind_l];	W_b_r = w_b_r[ind_r]
		
		edups_l = edups_r = edupsid_l= edupsid_r = 0	
		if ("`vce_select'"=="nn") {
			edups_l   = dups_l[ind_l];	  	    edups_r   = dups_r[ind_r]
			edupsid_l = dupsid_l[ind_l];	    edupsid_r = dupsid_r[ind_r]
		}
		
		u_l = (eX_l:-`c')/h_l;	u_r = (eX_r:-`c')/h_r;
		// Q1: Vandermonde via successive multiplication.
		R_q_l = J(eN_l,(`q'+1),1); R_q_r = J(eN_r,(`q'+1),1)
		if (`q' >= 1) {
			uq_l = eX_l :- `c'; uq_r = eX_r :- `c'
			for (j=2; j<=(`q'+1); j++)  {
				R_q_l[.,j] = R_q_l[.,j-1] :* uq_l
				R_q_r[.,j] = R_q_r[.,j-1] :* uq_r
			}
		}
		R_p_l = R_q_l[,1::(`p'+1)]; R_p_r = R_q_r[,1::(`p'+1)]
	
		********************************************************************************
		************ Computing RD estimates ********************************************
		********************************************************************************
		L_l = quadcross(R_p_l:*W_h_l,u_l:^(`p'+1)); L_r = quadcross(R_p_r:*W_h_r,u_r:^(`p'+1)) 
			invG_q_l  = cholinv(quadcross(R_q_l,W_b_l,R_q_l));	invG_q_r  = cholinv(quadcross(R_q_r,W_b_r,R_q_r))
			invG_p_l  = cholinv(quadcross(R_p_l,W_h_l,R_p_l));	invG_p_r  = cholinv(quadcross(R_p_r,W_h_r,R_p_r)) 
		
		if (rank(invG_p_l)==. | rank(invG_p_r)==. | rank(invG_q_l)==. | rank(invG_q_r)==. ){
		display("{err}Invertibility problem: check variability of running variable around cutoff. Try checking for mass points with option {cmd:masspoints(check)}.")
			* rc=1 is Stata's --Break--, so under `qui` the user saw ONLY
			* "--Break--" with no hint of the real cause. 506 is the standard
			* "matrix not positive definite" code.
			exit(506)
		}
		
		e_p1 = J((`q'+1),1,0); e_p1[`p'+2]=1
		e_v  = J((`p'+1),1,0); e_v[`deriv'+1]=1
		Q_q_l = ((R_p_l:*W_h_l)' - h_l^(`p'+1)*(L_l*e_p1')*((invG_q_l*R_q_l')':*W_b_l)')'
		Q_q_r = ((R_p_r:*W_h_r)' - h_r^(`p'+1)*(L_r*e_p1')*((invG_q_r*R_q_r')':*W_b_r)')'
		D_l = eY_l; D_r = eY_r		
		
		if ("`fuzzy'"!="") {
			T    = st_data(.,("`fuzzyvar'"), 0);	dT = 1
			T_l  = select(T,X:<`c');  eT_l  = T_l[ind_l]
			T_r  = select(T,X:>=`c'); eT_r  = T_r[ind_r]
			D_l  = D_l,eT_l; D_r = D_r,eT_r
		}
		
		if ("`covs'"!="") {
			eZ_l = Z_l[ind_l,]; eZ_r = Z_r[ind_r,]
			D_l  = D_l,eZ_l; D_r = D_r,eZ_r
			U_p_l = quadcross(R_p_l:*W_h_l,D_l); U_p_r = quadcross(R_p_r:*W_h_r,D_r)
		}
		
		if ("`cluster'"!="") {
			eC_l  = C_l[ind_l];	     eC_r  = C_r[ind_r]
			indC_l = order(eC_l,1);  indC_r = order(eC_r,1) 
			g_l = rows(panelsetup(eC_l[indC_l],1));	g_r = rows(panelsetup(eC_r[indC_r],1))
		}
		
		beta_p_l = invG_p_l*quadcross(R_p_l:*W_h_l,D_l); beta_q_l = invG_q_l*quadcross(R_q_l:*W_b_l,D_l); beta_bc_l = invG_p_l*quadcross(Q_q_l,D_l) 
		beta_p_r = invG_p_r*quadcross(R_p_r:*W_h_r,D_r); beta_q_r = invG_q_r*quadcross(R_q_r:*W_b_r,D_r); beta_bc_r = invG_p_r*quadcross(Q_q_r,D_r)
		beta_p  = beta_p_r  - beta_p_l
		beta_q  = beta_q_r  - beta_q_l
		beta_bc = beta_bc_r - beta_bc_l
		
		************************************ No Covariates **********************************
		
		if (dZ==0) {		
				tau_cl = tau_Y_cl = `scalepar'*factorial(`deriv')*beta_p[(`deriv'+1),1]
				tau_bc = tau_Y_bc = `scalepar'*factorial(`deriv')*beta_bc[(`deriv'+1),1]
				s_Y = 1				
				tau_Y_cl_l = `scalepar'*factorial(`deriv')*beta_p_l[(`deriv'+1),1]
				tau_Y_cl_r = `scalepar'*factorial(`deriv')*beta_p_r[(`deriv'+1),1]
				tau_Y_bc_l = `scalepar'*factorial(`deriv')*beta_bc_l[(`deriv'+1),1]
				tau_Y_bc_r = `scalepar'*factorial(`deriv')*beta_bc_r[(`deriv'+1),1]				
				bias_l = tau_Y_cl_l - tau_Y_bc_l
				bias_r = tau_Y_cl_r - tau_Y_bc_r 		
	
				beta_Y_p_l = `scalepar'*factorial(`deriv')*beta_p_l[,1]
				beta_Y_p_r = `scalepar'*factorial(`deriv')*beta_p_r[,1]

				*** Fuzzy RD ********************
				if (dT>0) {
					tau_T_cl =  factorial(`deriv')*beta_p[(`deriv'+1),2]
					tau_T_bc = 	factorial(`deriv')*beta_bc[(`deriv'+1),2]
					
					beta_T_p_l = factorial(`deriv')*beta_p_l[,2]
					beta_T_p_r = factorial(`deriv')*beta_p_r[,2]
					
					s_Y = (1/tau_T_cl \ -(tau_Y_cl/tau_T_cl^2))
					B_F = tau_Y_cl-tau_Y_bc \ tau_T_cl-tau_T_bc
					tau_cl = tau_Y_cl/tau_T_cl
					tau_bc = tau_cl - s_Y'*B_F
					sV_T = 0 \ 1										
					tau_T_cl_l =  factorial(`deriv')*beta_p_l[(`deriv'+1),2]
					tau_T_cl_r =  factorial(`deriv')*beta_p_r[(`deriv'+1),2]
					tau_T_bc_l =  factorial(`deriv')*beta_bc_l[(`deriv'+1),2]
					tau_T_bc_r =  factorial(`deriv')*beta_bc_r[(`deriv'+1),2]					
					B_F_l = tau_Y_cl_l-tau_Y_bc_l \ tau_T_cl_l-tau_T_bc_l
					B_F_r = tau_Y_cl_r-tau_Y_bc_r \ tau_T_cl_r-tau_T_bc_r					
					bias_l = s_Y'*B_F_l
					bias_r = s_Y'*B_F_r		
		
				
				}	
				
		}
		
		*********************************** Covariates **********************************
				
		if (dZ>0) {				
			ZWD_p_l  = quadcross(eZ_l,W_h_l,D_l)
			ZWD_p_r  = quadcross(eZ_r,W_h_r,D_r)
			colsZ = (2+dT)::(2+dT+dZ-1)
			UiGU_p_l =  quadcross(U_p_l[,colsZ],invG_p_l*U_p_l) 
			UiGU_p_r =  quadcross(U_p_r[,colsZ],invG_p_r*U_p_r) 
			ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
			ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
			ZWY_p_l = ZWD_p_l[,1::1+dT] - UiGU_p_l[,1::1+dT] 
			ZWY_p_r = ZWD_p_r[,1::1+dT] - UiGU_p_r[,1::1+dT]     
			ZWZ_p = ZWZ_p_r + ZWZ_p_l
			ZWY_p = ZWY_p_r + ZWY_p_l
			if ("`covs_drop_coll'"=="0") gamma_p = cholinv(ZWZ_p)*ZWY_p
			if ("`covs_drop_coll'"=="1") gamma_p =  invsym(ZWZ_p)*ZWY_p
			if ("`covs_drop_coll'"=="2") gamma_p =    pinv(ZWZ_p)*ZWY_p
	
			s_Y = (1 \  -gamma_p[,1])
			
			*** Sharp RD ********************
			if (dT==0) {
				* factorial(deriv) converts the local-polynomial coefficient
				* into the derivative estimate. It was missing on these tau
				* lines while beta_Y_p_l/r below and V both carry it, so for
				* deriv >= 2 tau was 1/deriv! of the correct value and the
				* returns were internally inconsistent.
				tau_cl = `scalepar'*factorial(`deriv')*s_Y'*beta_p[(`deriv'+1),]'
				tau_bc = `scalepar'*factorial(`deriv')*s_Y'*beta_bc[(`deriv'+1),]'				
				tau_Y_cl_l = `scalepar'*factorial(`deriv')*s_Y'*beta_p_l[(`deriv'+1),]'
				tau_Y_cl_r = `scalepar'*factorial(`deriv')*s_Y'*beta_p_r[(`deriv'+1),]'
				tau_Y_bc_l = `scalepar'*factorial(`deriv')*s_Y'*beta_bc_l[(`deriv'+1),]'
				tau_Y_bc_r = `scalepar'*factorial(`deriv')*s_Y'*beta_bc_r[(`deriv'+1),]'				
				bias_l = tau_Y_cl_l-tau_Y_bc_l
				bias_r = tau_Y_cl_r-tau_Y_bc_r 		
				
				beta_Y_p_l = `scalepar'*factorial(`deriv')*(s_Y'*beta_p_l')'
				beta_Y_p_r = `scalepar'*factorial(`deriv')*(s_Y'*beta_p_r')	'				
			}
			
			*** Fuzzy RD ********************
			if (dT>0) {
					s_T  = 1 \ -gamma_p[,2]
					sV_T = (0 \ 1 \ -gamma_p[,2] )
					tau_Y_cl   = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_p[  (`deriv'+1),1], beta_p[  (`deriv'+1),colsZ]))
					tau_Y_cl_l = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_p_l[(`deriv'+1),1], beta_p_l[(`deriv'+1),colsZ]))
					tau_Y_cl_r = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_p_r[(`deriv'+1),1], beta_p_r[(`deriv'+1),colsZ]))
			
					tau_Y_bc   = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_bc[  (`deriv'+1),1], beta_bc[  (`deriv'+1),colsZ]))
					tau_Y_bc_l = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_bc_l[(`deriv'+1),1], beta_bc_l[(`deriv'+1),colsZ]))
					tau_Y_bc_r = `scalepar'*factorial(`deriv')*s_Y'*vec((beta_bc_r[(`deriv'+1),1], beta_bc_r[(`deriv'+1),colsZ]))	
					
					tau_T_cl   = factorial(`deriv')*s_T'*vec((beta_p[  (`deriv'+1),2], beta_p[  (`deriv'+1),colsZ]))
					tau_T_cl_l = factorial(`deriv')*s_T'*vec((beta_p_l[(`deriv'+1),2], beta_p_l[(`deriv'+1),colsZ]))
					tau_T_cl_r = factorial(`deriv')*s_T'*vec((beta_p_r[(`deriv'+1),2], beta_p_r[(`deriv'+1),colsZ]))
					
					tau_T_bc =   factorial(`deriv')*s_T'*vec((beta_bc[  (`deriv'+1),2], beta_bc[  (`deriv'+1),colsZ]))
					tau_T_bc_l = factorial(`deriv')*s_T'*vec((beta_bc_l[(`deriv'+1),2], beta_bc_l[(`deriv'+1),colsZ]))
					tau_T_bc_r = factorial(`deriv')*s_T'*vec((beta_bc_r[(`deriv'+1),2], beta_bc_r[(`deriv'+1),colsZ]))
					
					beta_Y_p_l = `scalepar'*factorial(`deriv')*(s_Y'*(beta_p_l[,1], beta_p_l[,colsZ])')'
					beta_Y_p_r = `scalepar'*factorial(`deriv')*(s_Y'*(beta_p_r[,1], beta_p_r[,colsZ])')'
					beta_T_p_l =            factorial(`deriv')*(s_T'*(beta_p_l[,2], beta_p_l[,colsZ])')'
					beta_T_p_r =            factorial(`deriv')*(s_T'*(beta_p_r[,2], beta_p_r[,colsZ])')'
					
					
					B_F = tau_Y_cl-tau_Y_bc \ tau_T_cl-tau_T_bc
					s_Y = 1/tau_T_cl \ -(tau_Y_cl/tau_T_cl^2)
					tau_cl = tau_Y_cl/tau_T_cl
					tau_bc = tau_cl - s_Y'*B_F
					
					B_F_l = tau_Y_cl_l-tau_Y_bc_l \ tau_T_cl_l-tau_T_bc_l
					B_F_r = tau_Y_cl_r-tau_Y_bc_r \ tau_T_cl_r-tau_T_bc_r
					
					bias_l = s_Y'*B_F_l
					bias_r = s_Y'*B_F_r
					
					s_Y = (1/tau_T_cl \ -(tau_Y_cl/tau_T_cl^2) \ -(1/tau_T_cl)*gamma_p[,1] + (tau_Y_cl/tau_T_cl^2)*gamma_p[,2])
			}
		}
			
		**************************************************************************
		************ Computing variance-covariance matrix ************************
		**************************************************************************
		hii_p_l=hii_p_r=hii_q_l=hii_q_r=predicts_p_l=predicts_p_r=predicts_q_l=predicts_q_r=0
		if ("`vce_select'"=="hc0" | "`vce_select'"=="hc1" | "`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
			predicts_p_l=R_p_l*beta_p_l
			predicts_p_r=R_p_r*beta_p_r
			predicts_q_l=R_q_l*beta_q_l
			predicts_q_r=R_q_r*beta_q_r
			
			if ("`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
				hii_p_l = rowsum((R_p_l*invG_p_l):*(R_p_l:*W_h_l))
				hii_p_r = rowsum((R_p_r*invG_p_r):*(R_p_r:*W_h_r))
				
				
				if ("`vleverage'"=="") {
				
					hii_q_l = rowsum((R_q_l*invG_q_l):*(R_q_l:*W_b_l))
					hii_q_r = rowsum((R_q_r*invG_q_r):*(R_q_r:*W_b_r))
				
				}
				else {
				
					hii_q_l = hii_p_l
					hii_q_r = hii_p_r
				
				}
				
			}
			
		}
			
		res_h_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_p_l, hii_p_l, "`vce_select'", `nnmatch', edups_l, edupsid_l, `p'+1)
		res_h_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_p_r, hii_p_r, "`vce_select'", `nnmatch', edups_r, edupsid_r, `p'+1)
		
		if ("`vce_select'"=="nn") {
				res_b_l = res_h_l;	res_b_r = res_h_r
		}
		else {
		
				res_b_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_q_l, hii_q_l, "`vce_select'", `nnmatch', edups_l, edupsid_l, `q'+1)
				res_b_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_q_r, hii_q_r, "`vce_select'", `nnmatch', edups_r, edupsid_r, `q'+1)
		
		}

		// V_cl (main variance): pass invG and sqrtRX for CR2/CR3 hat-matrix
		// adjustment. For CR1 / non-cluster these are ignored.
		sqrtRX_p_l = ("`cr_method'"=="crv2" | "`cr_method'"=="crv3") ? R_p_l:*sqrt(W_h_l) : J(0,0,.)
		sqrtRX_p_r = ("`cr_method'"=="crv2" | "`cr_method'"=="crv3") ? R_p_r:*sqrt(W_h_r) : J(0,0,.)
		invG_p_l_c = ("`cr_method'"=="crv2" | "`cr_method'"=="crv3") ? invG_p_l : J(0,0,.)
		invG_p_r_c = ("`cr_method'"=="crv2" | "`cr_method'"=="crv3") ? invG_p_r : J(0,0,.)
		V_Y_cl_l = invG_p_l*rdrobust_vce(dT+dZ, s_Y, R_p_l:*W_h_l, res_h_l, eC_l, indC_l, invG_p_l_c, sqrtRX_p_l, "`cr_method'", 0)*invG_p_l
		V_Y_cl_r = invG_p_r*rdrobust_vce(dT+dZ, s_Y, R_p_r:*W_h_r, res_h_r, eC_r, indC_r, invG_p_r_c, sqrtRX_p_r, "`cr_method'", 0)*invG_p_r
		// V_bc (robust-bias variance): for CRV2/CRV3 with cluster, use decoupled
		// helper (Q_q sandwich + q-regression cluster leverage). Otherwise CR1
		// with k_override = q+1 to align df correction with the q-regression
		// path used at h=b (so SE is continuous across the h=b boundary).
		if (("`cr_method'"=="crv2" | "`cr_method'"=="crv3") & "`cluster'"!="") {
			V_Y_bc_l = invG_p_l*rdrobust_vce_qq_cluster(Q_q_l, R_q_l, W_b_l, invG_q_l, res_b_l, eC_l, indC_l, dT+dZ, s_Y, "`cr_method'")*invG_p_l
			V_Y_bc_r = invG_p_r*rdrobust_vce_qq_cluster(Q_q_r, R_q_r, W_b_r, invG_q_r, res_b_r, eC_r, indC_r, dT+dZ, s_Y, "`cr_method'")*invG_p_r
		}
		else {
			V_Y_bc_l = invG_p_l*rdrobust_vce(dT+dZ, s_Y, Q_q_l, res_b_l, eC_l, indC_l, J(0,0,.), J(0,0,.), "cr1", `q'+1)*invG_p_l
			V_Y_bc_r = invG_p_r*rdrobust_vce(dT+dZ, s_Y, Q_q_r, res_b_r, eC_r, indC_r, J(0,0,.), J(0,0,.), "cr1", `q'+1)*invG_p_r
		}
		V_tau_cl = (`scalepar')^2*factorial(`deriv')^2*(V_Y_cl_l+V_Y_cl_r)[`deriv'+1,`deriv'+1]
		V_tau_rb = (`scalepar')^2*factorial(`deriv')^2*(V_Y_bc_l+V_Y_bc_r)[`deriv'+1,`deriv'+1]
		se_tau_cl = sqrt(V_tau_cl);	se_tau_rb = sqrt(V_tau_rb)

		if ("`fuzzy'"!="") {
			V_T_cl_l = invG_p_l*rdrobust_vce(dT+dZ, sV_T, R_p_l:*W_h_l, res_h_l, eC_l, indC_l, invG_p_l_c, sqrtRX_p_l, "`cr_method'", 0)*invG_p_l
			V_T_cl_r = invG_p_r*rdrobust_vce(dT+dZ, sV_T, R_p_r:*W_h_r, res_h_r, eC_r, indC_r, invG_p_r_c, sqrtRX_p_r, "`cr_method'", 0)*invG_p_r
			if (("`cr_method'"=="crv2" | "`cr_method'"=="crv3") & "`cluster'"!="") {
				V_T_bc_l = invG_p_l*rdrobust_vce_qq_cluster(Q_q_l, R_q_l, W_b_l, invG_q_l, res_b_l, eC_l, indC_l, dT+dZ, sV_T, "`cr_method'")*invG_p_l
				V_T_bc_r = invG_p_r*rdrobust_vce_qq_cluster(Q_q_r, R_q_r, W_b_r, invG_q_r, res_b_r, eC_r, indC_r, dT+dZ, sV_T, "`cr_method'")*invG_p_r
			}
			else {
				// k_override = q+1: see V_Y_bc comment above (continuity at h=b).
				V_T_bc_l = invG_p_l*rdrobust_vce(dT+dZ, sV_T, Q_q_l, res_b_l, eC_l, indC_l, J(0,0,.), J(0,0,.), "cr1", `q'+1)*invG_p_l
				V_T_bc_r = invG_p_r*rdrobust_vce(dT+dZ, sV_T, Q_q_r, res_b_r, eC_r, indC_r, J(0,0,.), J(0,0,.), "cr1", `q'+1)*invG_p_r
			}
			V_T_cl = factorial(`deriv')^2*(V_T_cl_l+V_T_cl_r)[`deriv'+1,`deriv'+1]
			V_T_rb = factorial(`deriv')^2*(V_T_bc_l+V_T_bc_r)[`deriv'+1,`deriv'+1]
			se_tau_T_cl = sqrt(V_T_cl);	se_tau_T_rb = sqrt(V_T_rb)
		}
		
	
		
		**** Stored results
		st_numscalar("N", N)
		st_numscalar("N_l", N_l)
		st_numscalar("N_r", N_r)
		st_numscalar("x_l_min", x_l_min)
		st_numscalar("x_l_max", x_l_max)
		st_numscalar("x_r_min", x_r_min)
		st_numscalar("x_r_max", x_r_max)
	
		st_numscalar("h_l", h_l)
		st_numscalar("h_r", h_r)
		st_numscalar("b_l", b_l)
		st_numscalar("b_r", b_r)
	
		st_numscalar("quant", -invnormal(abs((1-(`level'/100))/2)))
		st_numscalar("N_h_l", N_h_l);	st_numscalar("N_b_l", N_b_l)
		st_numscalar("N_h_r", N_h_r);	st_numscalar("N_b_r", N_b_r)
		
		st_numscalar("tau_cl", tau_cl); st_numscalar("se_tau_cl", se_tau_cl)
		st_numscalar("tau_bc", tau_bc);	st_numscalar("se_tau_rb", se_tau_rb)
		
		st_numscalar("tau_Y_cl_r", tau_Y_cl_r); st_numscalar("tau_Y_cl_l", tau_Y_cl_l)
		st_numscalar("tau_Y_bc_r", tau_Y_bc_r);	st_numscalar("tau_Y_bc_l", tau_Y_bc_l)
		
		st_numscalar("bias_l", bias_l);  st_numscalar("bias_r", bias_r)
		
		st_matrix("beta_Y_p_r", beta_Y_p_r); st_matrix("beta_Y_p_l", beta_Y_p_l)

		st_numscalar("g_l",  g_l);       st_numscalar("g_r",   g_r)
		* e(b) / e(V): three named coefficients mirroring the printed output.
		*   Conventional:   point = tau_cl,  SE = sqrt(V_tau_cl)
		*   Bias-corrected: point = tau_bc,  SE = sqrt(V_tau_cl)  (not a recommended
		*                                                          inferential object
		*                                                          per CCT 2014; kept
		*                                                          for parity with the
		*                                                          displayed table)
		*   Robust:         point = tau_bc,  SE = sqrt(V_tau_rb)  (CCT-recommended RBC)
		* e(V) is block-diagonal (each row is its own estimand).
		st_matrix("`bmat'", (tau_cl, tau_bc, tau_bc))
		st_matrix("`Vmat'", (V_tau_cl, 0, 0 \ 0, V_tau_cl, 0 \ 0, 0, V_tau_rb))
		st_matrix("V_Y_cl_r", V_Y_cl_r); st_matrix("V_Y_cl_l", V_Y_cl_l)
		st_matrix("V_Y_bc_r", V_Y_bc_r); st_matrix("V_Y_bc_l", V_Y_bc_l)
		st_numscalar("masspoints_found", masspoints_found)
		
		if ("`covs'"!="") {
			st_matrix("gamma_p", gamma_p)
		}
					
					
		if ("`fuzzy'"!="") {
			st_numscalar("tau_T_cl", tau_T_cl); st_numscalar("se_tau_T_cl", se_tau_T_cl)
			st_numscalar("tau_T_bc", tau_T_bc);	st_numscalar("se_tau_T_rb", se_tau_T_rb)	
			
			st_numscalar("tau_T_cl_r", tau_T_cl_r); st_numscalar("tau_T_cl_l", tau_T_cl_l)
			st_numscalar("tau_T_bc_r", tau_T_bc_r);	st_numscalar("tau_T_bc_l", tau_T_bc_l)
			
			st_matrix("beta_T_p_r", beta_T_p_r); st_matrix("beta_T_p_l", beta_T_p_l)

		}
	}
	
	************************************************
	********* OUTPUT TABLE *************************
	************************************************
	local rho_l = scalar(h_l)/scalar(b_l)
	local rho_r = scalar(h_r)/scalar(b_r)
	
	disp ""
	if "`fuzzy'"=="" {
		if ("`covs'"=="") {
			if      ("`deriv'"=="0") disp "Sharp RD estimates using local polynomial regression." 
			else if ("`deriv'"=="1") disp "Sharp Kink RD estimates using local polynomial regression."	
			else                     disp "Sharp RD estimates using local polynomial regression. Derivative of order " `deriv' "."	
		}
		else {
			if      ("`deriv'"=="0") disp "Covariate-adjusted Sharp RD estimates using local polynomial regression." 
			else if ("`deriv'"=="1") disp "Covariate-adjusted Sharp Kink RD estimates using local polynomial regression."	
			else                     disp "Covariate-adjusted Sharp RD estimates using local polynomial regression. Derivative of order " `deriv' "."	
		}
	}
	else {
		if ("`covs'"=="") {
			if      ("`deriv'"=="0") disp "Fuzzy RD estimates using local polynomial regression." 
			else if ("`deriv'"=="1") disp "Fuzzy Kink RD estimates using local polynomial regression."	
			else                     disp "Fuzzy RD estimates using local polynomial regression. Derivative of order " `deriv' "."	
		}
		else {
			if      ("`deriv'"=="0") disp "Covariate-adjusted Fuzzy RD estimates using local polynomial regression." 
			else if ("`deriv'"=="1") disp "Covariate-adjusted Fuzzy Kink RD estimates using local polynomial regression."	
			else                     disp "Covariate-adjusted Fuzzy RD estimates using local polynomial regression. Derivative of order " `deriv' "."			
		}
	}

	disp ""
	disp in smcl in gr "{ralign 18: Cutoff c = `c'}"        _col(19) " {c |} " _col(21) in gr "Left of " in yellow "c"  _col(33) in gr "Right of " in yellow "c"         _col(55) in gr "Number of obs = "  in yellow %10.0f scalar(N)
	disp in smcl in gr "{hline 19}{c +}{hline 22}"                                                                                                                       _col(55) in gr "BW type       = "  in yellow "{ralign 10:`bwselect'}" 
	disp in smcl in gr "{ralign 18:Number of obs}"          _col(19) " {c |} " _col(21) as result %9.0f scalar(N_l)             _col(34) %9.0f  scalar(N_r)                              _col(55) in gr "Kernel        = "  in yellow "{ralign 10:`kernel_type'}" 
	disp in smcl in gr "{ralign 18:Eff. Number of obs}"     _col(19) " {c |} " _col(21) as result %9.0f scalar(N_h_l)           _col(34) %9.0f  scalar(N_h_r)                            _col(55) in gr "VCE method    = "  in yellow "{ralign 10:`vce_type'}" 
	disp in smcl in gr "{ralign 18:Order est. (p)}"         _col(19) " {c |} " _col(21) as result %9.0f `p'             _col(34) %9.0f  `p'         
	disp in smcl in gr "{ralign 18:Order bias (q)}"         _col(19) " {c |} " _col(21) as result %9.0f `q'             _col(34) %9.0f  `q'                              
	disp in smcl in gr "{ralign 18:BW est. (h)}"            _col(19) " {c |} " _col(21) as result %9.3f scalar(h_l)           _col(34) %9.3f  scalar(h_r)                                   
	disp in smcl in gr "{ralign 18:BW bias (b)}"            _col(19) " {c |} " _col(21) as result %9.3f scalar(b_l)           _col(34) %9.3f  scalar(b_r)
	disp in smcl in gr "{ralign 18:rho (h/b)}"              _col(19) " {c |} " _col(21) as result %9.3f `rho_l'         _col(34) %9.3f  `rho_r'
	if ("`masspoints'"=="check" | masspoints_found==1) disp in smcl in gr "{ralign 18:Unique obs}"         _col(19) " {c |} " _col(21) as result %9.0f scalar(M_l)           _col(34) %9.0f  scalar(M_r)                    
	if ("`cluster'"!="")                               disp in smcl in gr "{ralign 18:Number of clusters}" _col(19) " {c |} " _col(21) as result %9.0f scalar(g_l)           _col(34) %9.0f  scalar(g_r)                         
	disp ""
			
	if ("`fuzzy'"!="") {		
		disp in yellow "First-stage estimates. Outcome: `fuzzyvar'. Running variable: `x'."
		disp in smcl in gr "{hline 19}{c TT}{hline 60}"
		
		
		if ("`all'"!="") {
			disp in smcl in gr "{ralign 18:Method}"  _col(19) " {c |} " _col(24) "Coef."  _col(33) `"Std. Err."'   _col(46) "z"    _col(52) "P>|z|"   _col(61) `"[`level'% Conf. Interval]"'
			disp in smcl in gr "{hline 19}{c +}{hline 60}"
			disp in smcl in gr "{ralign 18:Conventional}"      _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_T_cl) _col(33) %7.0g scalar(se_tau_T_cl) _col(43) %5.4f scalar(tau_T_cl/se_tau_T_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_T_cl/se_tau_T_cl))) _col(60) %8.0g  scalar(tau_T_cl - quant*se_tau_T_cl) _col(73) %8.0g scalar(tau_T_cl + quant*se_tau_T_cl)  
			disp in smcl in gr "{ralign 18:Bias-corrected}"    _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_T_bc) _col(33) %7.0g scalar(se_tau_T_cl) _col(43) %5.4f scalar(tau_T_bc/se_tau_T_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_T_bc/se_tau_T_cl))) _col(60) %8.0g  scalar(tau_T_bc - quant*se_tau_T_cl) _col(73) %8.0g scalar(tau_T_bc + quant*se_tau_T_cl) 
			disp in smcl in gr "{ralign 18:Robust}"            _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_T_bc) _col(33) %7.0g scalar(se_tau_T_rb) _col(43) %5.4f scalar(tau_T_bc/se_tau_T_rb) _col(52) %5.3f  scalar(2*normal(-abs(tau_T_bc/se_tau_T_rb))) _col(60) %8.0g  scalar(tau_T_bc - quant*se_tau_T_rb) _col(73) %8.0g scalar(tau_T_bc + quant*se_tau_T_rb) 
		}
		else if ("`detail'"!="") {
			disp in smcl in gr "{ralign 18:Method}"  _col(19) " {c |} " _col(24) "Coef."  _col(33) `"Std. Err."'   _col(46) "z"    _col(52) "P>|z|"   _col(61) `"[`level'% Conf. Interval]"'
			disp in smcl in gr "{hline 19}{c +}{hline 60}"
			disp in smcl in gr "{ralign 18:Conventional}"      _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_T_cl) _col(33) %7.0g scalar(se_tau_T_cl) _col(43) %5.4f scalar(tau_T_cl/se_tau_T_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_T_cl/se_tau_T_cl)))  _col(60) %8.0g  scalar(tau_T_cl) - scalar(quant*se_tau_T_cl) _col(73) %8.0g scalar(tau_T_cl + quant*se_tau_T_cl) 
			disp in smcl in gr "{ralign 18:Robust}"            _col(19) " {c |} " _col(22) in ye %7.0g "    -"  _col(33) %7.0g "    -"     _col(43) %5.4f scalar(tau_T_bc/se_tau_T_rb) _col(52) %5.3f  scalar(2*normal(-abs(tau_T_bc/se_tau_T_rb)))  _col(60) %8.0g  scalar(tau_T_bc - quant*se_tau_T_rb) _col(73) %8.0g scalar(tau_T_bc + quant*se_tau_T_rb) 
		}
		else {
			disp in smcl in gr "{ralign 18:}"                   _col(19) " {c |} " _col(22) "Point"    _col(35) " {c |} "    "Robust Inference" 
			disp in smcl in gr "{ralign 18:}"                   _col(19) " {c |} " _col(22) "Estimate" _col(35) " {c |} "    "z-stat"       _col(52) "P>|z|"    _col(61) `"[`level'% Conf. Interval]"'
			disp in smcl in gr "{hline 19}{c +}{hline 60}"			 
			disp in smcl in gr "{ralign 18:RD Effect}"          _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_T_cl)  _col(35) " {c |} "  %5.4f scalar(tau_T_bc/se_tau_T_rb)  _col(52) %5.3f  scalar(2*normal(-abs(tau_T_bc/se_tau_T_rb))) _col(61) %8.0g scalar(tau_T_cl - quant*se_tau_T_cl)  _col(73)  %8.0g scalar(tau_T_cl + quant*se_tau_T_cl)			

		}
		
			disp in smcl in gr "{hline 19}{c BT}{hline 60}"
			disp ""
	}
	
	if ("`fuzzy'"=="") disp           "Outcome: `y'. Running variable: `x'."
	else               disp in yellow "Treatment effect estimates. Outcome: `y'. Running variable: `x'. Treatment Status: `fuzzyvar'."
		
	disp in smcl in gr "{hline 19}{c TT}{hline 60}"
		
	if ("`all'"!="") {
		disp in smcl in gr "{ralign 18:Method}"         _col(19) " {c |} " _col(24) "Coef."               _col(33) `"Std. Err."'    _col(46) "z"                    _col(52) "P>|z|"                                  _col(61) `"[`level'% Conf. Interval]"'
		disp in smcl in gr "{hline 19}{c +}{hline 60}"
		disp in smcl in gr "{ralign 18:Conventional}"   _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_cl)    _col(33) %7.0g scalar(se_tau_cl) _col(43) %5.4f scalar(tau_cl/se_tau_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_cl/se_tau_cl))) _col(60) %8.0g  scalar(tau_cl - quant*se_tau_cl) _col(73) %8.0g scalar(tau_cl + quant*se_tau_cl)  
		disp in smcl in gr "{ralign 18:Bias-corrected}" _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_bc)    _col(33) %7.0g scalar(se_tau_cl) _col(43) %5.4f scalar(tau_bc/se_tau_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_bc/se_tau_cl))) _col(60) %8.0g  scalar(tau_bc - quant*se_tau_cl) _col(73) %8.0g scalar(tau_bc + quant*se_tau_cl)  
		disp in smcl in gr "{ralign 18:Robust}"         _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_bc)    _col(33) %7.0g scalar(se_tau_rb) _col(43) %5.4f scalar(tau_bc/se_tau_rb) _col(52) %5.3f  scalar(2*normal(-abs(tau_bc/se_tau_rb))) _col(60) %8.0g  scalar(tau_bc - quant*se_tau_rb) _col(73) %8.0g scalar(tau_bc + quant*se_tau_rb)  
	}
	else if ("`detail'"!="") {
		disp in smcl in gr "{ralign 18:Method}"         _col(19) " {c |} " _col(24) "Coef."               _col(33) `"Std. Err."'    _col(46) "z"                    _col(52) "P>|z|"                                  _col(61) `"[`level'% Conf. Interval]"'
		disp in smcl in gr "{hline 19}{c +}{hline 60}"
		disp in smcl in gr "{ralign 18:Conventional}"   _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_cl)    _col(33) %7.0g scalar(se_tau_cl)  _col(43) %5.4f scalar(tau_cl/se_tau_cl) _col(52) %5.3f  scalar(2*normal(-abs(tau_cl/se_tau_cl)))  _col(60) %8.0g scalar(tau_cl - quant*se_tau_cl) _col(73) %8.0g scalar(tau_cl + quant*se_tau_cl) 
		disp in smcl in gr "{ralign 18:Robust}"         _col(19) " {c |} " _col(22) in ye %7.0g "    -"           _col(33) %7.0g "    -"            _col(43) %5.4f scalar(tau_bc/se_tau_rb) _col(52) %5.3f  scalar(2*normal(-abs(tau_bc/se_tau_rb)))  _col(60) %8.0g scalar(tau_bc - quant*se_tau_rb) _col(73) %8.0g scalar(tau_bc + quant*se_tau_rb) 
	} 
	else {
*		disp in smcl in gr "{ralign 18:}"                   _col(19) " {c |} " _col(22) "Estimate"     _col(35) "P>|z|"     _col(47)   `"[`level'% Robust CI]"'
*		disp in smcl in gr "{hline 19}{c +}{hline 60}"
*		disp in smcl in gr "{ralign 18:RD Effect}"   _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_cl)    _col(35)   %5.3f  scalar(2*normal(-abs(tau_bc/se_tau_rb)))   _col(45) %8.0g scalar(tau_bc - quant*se_tau_rb)   _col(55)  %8.0g scalar(tau_bc + quant*se_tau_rb) 
*		disp in smcl in gr "{hline 19}{c BT}{hline 60}"
*		disp ""
*		disp ""
		*disp in smcl in gr "{hline 19}{c TT}{hline 60}"
		disp in smcl in gr "{ralign 18:}"                   _col(19) " {c |} " _col(22) "Point"    _col(35) " {c |} "    "Robust Inference" 
		disp in smcl in gr "{ralign 18:}"                   _col(19) " {c |} " _col(22) "Estimate" _col(35) " {c |} "    "z-stat"       _col(52) "P>|z|"    _col(61) `"[`level'% Conf. Interval]"'
		disp in smcl in gr "{hline 19}{c +}{hline 60}"
		disp in smcl in gr "{ralign 18:RD Effect}"          _col(19) " {c |} " _col(22) in ye %7.0g scalar(tau_cl) _col(35) " {c |} "   %5.4f scalar(tau_bc/se_tau_rb)   _col(52) %5.3f scalar(2*normal(-abs(tau_bc/se_tau_rb))) _col(61) %8.0g scalar(tau_bc - quant*se_tau_rb)  _col(73) %8.0g scalar(tau_bc + quant*se_tau_rb) 
		
	}
		disp in smcl in gr "{hline 19}{c BT}{hline 60}"

	if ("`covs'"!="")        di "Covariate-adjusted estimates. Additional covariates included: `ncovs'"
	if ("`cluster'"!="")     di "Std. Err. adjusted for clusters in " "`clustvar'"
	if ("`scalepar'"!="1")   di "Scale parameter: " `scalepar' 
	if ("`scaleregul'"!="1") di "Scale regularization: " `scaleregul'
	if ("`masspoints'"=="check")  di "Running variable checked for mass points."  
	if ("`masspoints'"=="adjust" & masspoints_found==1) di "Estimates adjusted for mass points in the running variable."  	
	
	if ("`warnings'"=="") {
		if (scalar(h_l)>=`range_l' | scalar(h_r)>=`range_r') disp in red "WARNING: bandwidth {it:h} greater than the range of the data."
		if (scalar(b_l)>=`range_l' | scalar(b_r)>=`range_r') disp in red "WARNING: bandwidth {it:b} greater than the range of the data."
		if (scalar(N_h_l)<20 | scalar(N_h_r)<20)             disp in red "WARNING: bandwidth {it:h} too low."
		if (scalar(N_b_l)<20 | scalar(N_b_r)<20)             disp in red "WARNING: bandwidth {it:b} too low."
		if ("`sharpbw'"!="")                                 disp in red "WARNING: bandwidths automatically computed for sharp RD estimation."
		if ("`perf_comp'"!="")                               disp in red "WARNING: bandwidths automatically computed for sharp RD estimation because perfect compliance was detected on at least one side of the threshold."
	}
	
	local ci_l_rb = round(scalar(tau_bc - quant*se_tau_rb),0.001)
	local ci_r_rb = round(scalar(tau_bc + quant*se_tau_rb),0.001)

	matrix colnames `bmat' = Conventional Bias-corrected Robust
	matrix rownames `Vmat' = Conventional Bias-corrected Robust
	matrix colnames `Vmat' = Conventional Bias-corrected Robust

	}
	local _rc = _rc
	cwf `_orig_frame'
	frame drop `_work_frame'
	if `_rc' {
		* Error path: best-effort per-name Mata cleanup before reraising.
		* Note: in some configurations `exit N` from deeply nested if-blocks
		* inside `capture noisily { ... }` short-circuits past this branch,
		* leaving 1-2 transient Mata externals behind. The leak is bounded
		* (doesn't accumulate — next call's _mata_before captures them so
		* they aren't re-dropped) and harmless. Keep the per-name capture
		* form so that when the branch IS reached, missing names don't
		* halt cleanup of the rest.
		capture mata: _mtx = direxternal("*"); st_local("_mata_after", rows(_mtx) ? invtokens(_mtx') : "")
		capture mata: mata drop _mtx
		local _mata_new : list _mata_after - _mata_before
		foreach _mname of local _mata_new {
			capture mata mata drop `_mname'
		}
		error `_rc'
	}
	}

	ereturn clear

	* ST-6: `touse' came from `marksample' and so marked the [if]/[in] sample
	* ONLY. Every missing-value drop happens inside the work frame, so
	* e(sample) claimed rows that never entered the estimation (verified 1390
	* marked vs e(N)=1244). `drop_cond' is written over variables that exist
	* in this frame too, and locals survive the frame switch, so re-apply the
	* very same condition here. Same fix as rdhte.ado ST-8.
	qui replace `touse' = 0 if `drop_cond'

	ereturn post `bmat' `Vmat', esample(`touse')
	
	ereturn scalar N = `N'
	ereturn scalar N_l = scalar(N_l)
	ereturn scalar N_r = scalar(N_r)
	ereturn scalar N_h_l = scalar(N_h_l)
	ereturn scalar N_h_r = scalar(N_h_r)
	ereturn scalar N_b_l = scalar(N_b_l)
	ereturn scalar N_b_r = scalar(N_b_r)
	
	ereturn scalar c = `c'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	
	ereturn scalar h_l = scalar(h_l)
	ereturn scalar h_r = scalar(h_r)
	ereturn scalar b_l = scalar(b_l)
	ereturn scalar b_r = scalar(b_r)

	* 2x2 bws matrix -- rows {h, b}, cols {left, right}. Mirrors R fit$bws / Py fit.bws.
	tempname _bws_mat
	matrix `_bws_mat' = ( scalar(h_l), scalar(h_r) \ scalar(b_l), scalar(b_r) )
	matrix rownames `_bws_mat' = h b
	matrix colnames `_bws_mat' = left right
	ereturn matrix bws = `_bws_mat'

	ereturn scalar tau_cl   = scalar(tau_cl)
	ereturn scalar tau_cl_l = scalar(tau_Y_cl_l)
	ereturn scalar tau_cl_r = scalar(tau_Y_cl_r)
	ereturn scalar tau_bc   = scalar(tau_bc)
	ereturn scalar tau_bc_l = scalar(tau_Y_bc_l)
	ereturn scalar tau_bc_r = scalar(tau_Y_bc_r)

	ereturn scalar bias_l    = scalar(bias_l)
	ereturn scalar bias_r    = scalar(bias_r)
	ereturn scalar se_tau_cl = scalar(se_tau_cl)
	ereturn scalar se_tau_rb = scalar(se_tau_rb)
	
	ereturn scalar level   = `level'
	ereturn scalar ci_l_cl = scalar(tau_cl - quant*se_tau_cl)
	ereturn scalar ci_r_cl = scalar(tau_cl + quant*se_tau_cl)
	ereturn scalar ci_l_rb = scalar(tau_bc - quant*se_tau_rb)
	ereturn scalar ci_r_rb = scalar(tau_bc + quant*se_tau_rb)
	ereturn scalar pv_cl   = scalar(2*normal(-abs(tau_cl/se_tau_cl)))
	ereturn scalar pv_bc   = scalar(2*normal(-abs(tau_bc/se_tau_cl)))
	ereturn scalar pv_rb   = scalar(2*normal(-abs(tau_bc/se_tau_rb)))
	if ("`cluster'"!="") ereturn scalar n_clust = scalar(g_l) + scalar(g_r)
	
	if ("`fuzzy'"!="") {
		ereturn scalar tau_T_cl    = scalar(tau_T_cl)
		ereturn scalar tau_T_bc    = scalar(tau_T_bc)
		ereturn scalar se_tau_T_cl = scalar(se_tau_T_cl)
		ereturn scalar se_tau_T_rb = scalar(se_tau_T_rb)
		ereturn scalar tau_T_cl_l  = scalar(tau_T_cl_l)
		ereturn scalar tau_T_cl_r  = scalar(tau_T_cl_r)
		ereturn scalar tau_T_bc_l  = scalar(tau_T_bc_l)
		ereturn scalar tau_T_bc_r  = scalar(tau_T_bc_r)
		
		ereturn matrix beta_T_p_r = beta_T_p_r
		ereturn matrix beta_T_p_l = beta_T_p_l
	
	}
	
	ereturn matrix beta_Y_p_r = beta_Y_p_r
	ereturn matrix beta_Y_p_l = beta_Y_p_l
	
	if ("`covs'"!="") {
		ereturn matrix coef_covs = gamma_p
	}
	
	ereturn matrix V_cl_l = V_Y_cl_l
	ereturn matrix V_cl_r = V_Y_cl_r
	ereturn matrix V_rb_l = V_Y_bc_l
	ereturn matrix V_rb_r = V_Y_bc_r

	* 3x2 CI matrix -- rows Conventional/Bias-corrected/Robust, cols ll/ul.
	* Matches e(b) / e(V) and Stata's r(table) matrix-based convention for CIs.
	tempname _ci_mat
	matrix `_ci_mat' = ( scalar(tau_cl - quant*se_tau_cl), scalar(tau_cl + quant*se_tau_cl) \ ///
	                     scalar(tau_bc - quant*se_tau_cl), scalar(tau_bc + quant*se_tau_cl) \ ///
	                     scalar(tau_bc - quant*se_tau_rb), scalar(tau_bc + quant*se_tau_rb) )
	matrix rownames `_ci_mat' = Conventional Bias-corrected Robust
	matrix colnames `_ci_mat' = ll ul
	ereturn matrix ci = `_ci_mat'

	ereturn local ci_rb  [`ci_l_rb' ; `ci_r_rb']
	ereturn local kernel     = "`kernel_type'"
	ereturn local bwselect   = "`bwselect'"
	ereturn local vce_select = "`vce_raw'"
	ereturn local vce_type   = "`vce_type'"
	if ("`covs'"!="")    ereturn local covs "`covs_list'"
	if ("`cluster'"!="") ereturn local clustvar "`clustvar'"
	ereturn local runningvar "`x'"
	ereturn local depvar "`y'"
	* Title for estimates replay / estimates table
	if ("`fuzzy'"=="") local _rd_title "Sharp RD estimates"
	else                local _rd_title "Fuzzy RD estimates"
	if ("`covs'"!="")   local _rd_title "`_rd_title' (covariate-adjusted)"
	ereturn local title    "`_rd_title'"
	ereturn local cmdline  "rdrobust `0'"
	ereturn local cmd      "rdrobust"
	ereturn local precision "`precision'"

	* Drop transient matrices/scalars used as Mata-to-Stata transport buffers
	* so they don't leak into the caller's namespace.
	cap scalar drop h_l h_r b_l b_r quant
	cap scalar drop N_h_l N_h_r N_b_l N_b_r
	cap scalar drop tau_cl tau_bc se_tau_cl se_tau_rb
	cap scalar drop tau_Y_cl_l tau_Y_cl_r tau_Y_bc_l tau_Y_bc_r
	cap scalar drop bias_l bias_r g_l g_r masspoints_found
	cap scalar drop M_l M_r
	cap scalar drop tau_T_cl tau_T_bc se_tau_T_cl se_tau_T_rb
	cap scalar drop tau_T_cl_l tau_T_cl_r tau_T_bc_l tau_T_bc_r

	* Drop only the Mata externals we created (set difference vs the
	* entry-time snapshot). Library functions (rdrobust_*) are not
	* externals and stay loaded; the user's own Mata variables are
	* preserved because they appear in both snapshots.
	mata: _mtx = direxternal("*"); st_local("_mata_after", rows(_mtx) ? invtokens(_mtx') : "")
	mata: mata drop _mtx
	local _mata_new : list _mata_after - _mata_before
	if `"`_mata_new'"' != "" mata mata drop `_mata_new'

	* Normalize _rc on success: the `cap scalar drop tau_T_*` above leaks
	* _rc=111 on sharp RD (those scalars only exist for fuzzy designs),
	* and subsequent mata: statements do not update _rc on success.
	capture local _rc_ok = 0

end
