*!version 8.0.3  06-04-2020

capture program drop rdplot
program define rdplot, eclass
	syntax anything [if] [, c(real 0) p(integer 4) nbins(string) covs(string) covs_eval(string) covs_drop(string) binselect(string) scale(string) kernel(string) weights(string) h(string) k(integer 4) support(string) genvars hide ci(real 0) shade graph_options(string)  nochecks *]

	marksample touse
	tokenize "`anything'"
	local y `1'
	local x `2'
	
	******************** Set BW ***************************
	tokenize `h'	
	local w : word count `h'
	if `w' == 1 {
		local h_r = `"`1'"'
		local h_l = `"`1'"'
	}
	if `w' == 2 {
		local h_l `"`1'"'
		local h_r `"`2'"'
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:h()} accepts at most two inputs"  
		exit 125
	}
	******************** Set scale ***************************
	tokenize `scale'	
	local w : word count `scale'
	if `w' == 1 {
		local scale_r = `"`1'"'
		local scale_l = `"`1'"'
	}
	if `w' == 2 {
		local scale_l `"`1'"'
		local scale_r `"`2'"'
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:scale()} accepts at most two inputs"  
		exit 125
	}
	******************** Set nbins ***************************
	tokenize `nbins'	
	local w : word count `nbins'
	if `w' == 1 {
		local nbins_r = `"`1'"'
		local nbins_l = `"`1'"'
	}
	if `w' == 2 {
		local nbins_l `"`1'"'
		local nbins_r `"`2'"'
	}
	if `w' >= 3 {
		di as error  "{err}{cmd:nbins()} accepts at most two inputs"  
		exit 125
	}	
	******************** Set support ***************************
	tokenize `support'	
	local w : word count `support'
	if `w' == 2 {
		local support_l = `"`1'"'
		local support_r = `"`2'"'
	}
	if (`w' != 2 & "`support'"!="") {
		di as error  "{err}{cmd:support()} only accepts two inputs"  
		exit 125
	}	

	*****************************************
	preserve
	sort `x', stable
	qui keep if `touse'
	
	*****************************************************************
	**** DROP MISSINGS ******************************************
	*****************************************************************
	qui drop if `y'==. | `x'==.	
	if ("`covs'"~="") {
		qui ds `covs'
		local covs_list = r(varlist)
		local ncovs: word count `covs_list'	
		foreach z in `covs_list' {
			qui drop if `z'==.
		}
	}
		
	**** CHECK colinearity ******************************************
	local covs_drop_coll = 0
	if ("`covs_drop'"=="") local covs_drop = "on"	
	if ("`covs'"~="") {	
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
		qui ds `covs_new'
		local covs_list_new = r(varlist)
		local ncovs_new: word count `covs_list_new'
		if (`ncovs_new'<`ncovs') {
			if ("`covs_drop'"=="off") {	
				di as error  "{err}Multicollinearity issue detected in {cmd:covs}. Please rescale and/or remove redundant covariates, or add {cmd:covs_drop} option." 
				exit 125
			} 
			else {
				local ncovs = "`ncovs_new'"
				local covs_list = "`covs_list_new'"
				local covs_drop_coll = 1
			}	
		}
	}
	
	******************************************

	qui su `x'	
	local N = r(N)
	local x_min = r(min)
	local x_max = r(max)
	if ("`support'"!="") {	
		if (`support_l'<`x_min') {
			local x_min = `support_l'
		}
		if (`support_r'>`x_max') {
			local x_max = `support_r'
		}
	}
	local range_l = abs(`c'-`x_min')
	local range_r = abs(`x_max'-`c')
	
	qui su `x' if `x'<`c', d
	local n_l = r(N)
	
	qui su `x' if `x'>=`c', d
	local n_r = r(N)
	local n = `n_r' + `n_l'
	
	qui su `y' if `x'<`c'
	local var_l = r(sd)
	
	qui su `y' if `x'>=`c'
	local var_r = r(sd)
	
	if ("`h_l'"=="" & "`h_r'"=="") {
		local h_l = `range_l'
		local h_r = `range_r'
	}
	if "`kernel'"=="" local kernel = "uni"
	
	qui count if `x'<`c'  & `x'>=`c'-`h_l'
	local n_h_l = r(N)
	qui count if `x'>=`c' & `x'<=`c'+`h_r'
	local n_h_r = r(N)
	
	**************************** ERRORS
	if ("`scale_l'"=="" & "`scale_r'"=="") {
		local scale_r = 1
		local scale_l = 1
	}
	if ("`nbins_l'"=="" & "`nbins_r'"=="") {
		local nbins_r = 0
		local nbins_l = 0
	}
	
	if ("`binselect'"=="") {
		local binselect = "esmv" 
	}
	
	if ("`nochecks'"=="") {
		if (`c'<=`x_min' | `c'>=`x_max'){
			di as error  "{err}{cmd:c()} should be set within the range of `x'"  
		exit 125
		}

		if ("`p'"<"0" | "`nbins_l'"<"0" | "`nbins_r'"<"0"){
			di as error  "{err}{cmd:p()} and {cmd:nbins()} should be a positive integers"  
			exit 411
		}
			
		if ("`k'"<="0"){
			di as error  "{err}{cmd:k()} should be a positive integer"  
			exit 411
		}
		
		if (`n'<20){
			 di as error "{err}Not enough observations to perform bin calculations"  
			 exit 2001
		}
	}

	*******************************
	****** Start MATA *************
	*******************************
	mata{
		n_l=`n_l'
		n_r=`n_r'
		p=`p'
		k=`k'
		n=`n'
		c=`c'
		x_min = `x_min'
		x_max = `x_max'
		h_l     = strtoreal("`h_l'");     h_r     = strtoreal("`h_r'")
		nbins_l = strtoreal("`nbins_l'"); nbins_r = strtoreal("`nbins_r'")
		scale_l = strtoreal("`scale_l'"); scale_r = strtoreal("`scale_r'")
		
		y = st_data(.,("`y'"), 0);	x = st_data(.,("`x'"), 0)
		x_l = select(x,x:<c);	x_r = select(x,x:>=c)
		y_l = select(y,x:<c);	y_r = select(y,x:>=c)
			
		*if ("`hide'"=="" | "`genvars'"!="" ){
		
		************************************************************	
		************ Polynomial curve (order = p) ******************
		************************************************************
		
		if ("`covs'"=="") {
		
		rp_l = J(n_l,(p+1),.)
		rp_r = J(n_r,(p+1),.)
		for (j=1; j<=(p+1); j++) {
			rp_l[.,j] = (x_l:-c):^(j-1)
			rp_r[.,j] = (x_r:-c):^(j-1)
		}

		wh_l = rdrobust_kweight(x_l,c,h_l,"`kernel'")
		wh_r = rdrobust_kweight(x_r,c,h_r,"`kernel'")
		
		if ("`weights'"~="") {
			fw = st_data(.,("`weights'"), 0)
			fw_l = select(fw,x:<c);	fw_r = select(fw,x:>=c)
			wh_l = fw_l:*wh_l;	wh_r = fw_r:*wh_r
		}			

		
		if ("`covs_drop'"=="") {
			gamma_p1_l = cholinv(cross(rp_l,wh_l,rp_l))*cross(rp_l, wh_l, y_l)	
			gamma_p1_r = cholinv(cross(rp_r,wh_r,rp_r))*cross(rp_r, wh_r, y_r)
		} else {
			gamma_p1_l = invsym(cross(rp_l,wh_l,rp_l))*cross(rp_l, wh_l, y_l)	
			gamma_p1_r = invsym(cross(rp_r,wh_r,rp_r))*cross(rp_r, wh_r, y_r)
		}
		

		} else {
		
		Y = st_data(.,("`y'"), 0);	X = st_data(.,("`x'"), 0)
		X_l = select(X,X:<`c');	X_r = select(X,X:>=`c')
		Y_l = select(Y,X:<`c');	Y_r = select(Y,X:>=`c')
		h_l = strtoreal("`h_l'"); h_r = strtoreal("`h_r'")
		w_h_l = rdrobust_kweight(X_l,`c',h_l,"`kernel'");	w_h_r = rdrobust_kweight(X_r,`c',h_r,"`kernel'")
		ind_l = selectindex(w_h_l:> 0);	ind_r = selectindex(w_h_r:> 0)
	
		eY_l  = Y_l[ind_l];	eY_r  = Y_r[ind_r]
		eX_l  = X_l[ind_l];	eX_r  = X_r[ind_r]
		W_h_l = w_h_l[ind_l];	W_h_r = w_h_r[ind_r]
	
		u_l = (eX_l:-`c')/h_l;	u_r = (eX_r:-`c')/h_r;
		R_p_l = J(length(ind_l),(`p'+1),.); R_p_r = J(length(ind_r),(`p'+1),.)
		for (j=1; j<=(`p'+1); j++)  {
			R_p_l[.,j] = (eX_l:-`c'):^(j-1);  R_p_r[.,j] = (eX_r:-`c'):^(j-1)
		}
		
		L_l = quadcross(R_p_l:*W_h_l,u_l:^(`p'+1)); L_r = quadcross(R_p_r:*W_h_r,u_r:^(`p'+1)) 
		
		if ("`covs_drop'"=="") {
			invG_p_l  = cholinv(quadcross(R_p_l,W_h_l,R_p_l));	
			invG_p_r  = cholinv(quadcross(R_p_r,W_h_r,R_p_r)) 
		} else {		
			invG_p_l  = invsym(quadcross(R_p_l,W_h_l,R_p_l));	
			invG_p_r  = invsym(quadcross(R_p_r,W_h_r,R_p_r)) 
		}
		
		Z    = st_data(.,tokens("`covs'"), 0); dZ = cols(Z)
		Z_l  = select(Z,X:<`c');	eZ_l = Z_l[ind_l,]
		Z_r  = select(Z,X:>=`c');	eZ_r = Z_r[ind_r,]
		D_l  = eY_l,eZ_l; D_r = eY_r,eZ_r
		U_p_l = quadcross(R_p_l:*W_h_l,D_l); U_p_r = quadcross(R_p_r:*W_h_r,D_r)
		
		beta_p_l = invG_p_l*quadcross(R_p_l:*W_h_l,D_l) 
		beta_p_r = invG_p_r*quadcross(R_p_r:*W_h_r,D_r)
		
		ZWD_p_l  = quadcross(eZ_l,W_h_l,D_l)
		ZWD_p_r  = quadcross(eZ_r,W_h_r,D_r)
		colsZ = (2)::(2+dZ-1)
			
		UiGU_p_l =  quadcross(U_p_l[,colsZ],invG_p_l*U_p_l) 
		UiGU_p_r =  quadcross(U_p_r[,colsZ],invG_p_r*U_p_r) 
		ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
		ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
		ZWY_p_l = ZWD_p_l[,1] - UiGU_p_l[,1] 
		ZWY_p_r = ZWD_p_r[,1] - UiGU_p_r[,1]     
		ZWZ_p = ZWZ_p_r + ZWZ_p_l
		ZWY_p = ZWY_p_r + ZWY_p_l
		if ("`covs_drop'"=="") {
			gamma_p = cholinv(ZWZ_p)*ZWY_p
		} else {
			gamma_p = invsym(ZWZ_p)*ZWY_p
		}
		
		s_Y = (1 \  -gamma_p[,1])
		gamma_p1_l  = (s_Y'*beta_p_l')'
		gamma_p1_r  = (s_Y'*beta_p_r')'		
		}
		
		st_matrix("gamma_p1_l", gamma_p1_l)
		st_matrix("gamma_p1_r", gamma_p1_r)
		
		*********** Preparte data for polynomial curve plot *****
		nplot = 500
		x_plot_l = rangen(c-h_l,c,nplot)
		x_plot_r = rangen(c,c+h_r,nplot)
		rplot_l = J(nplot,(p+1),.)
		rplot_r = J(nplot,(p+1),.)
		for (j=1; j<=(p+1); j++) {
			rplot_l[.,j] = (x_plot_l:-c):^(j-1)
			rplot_r[.,j] = (x_plot_r:-c):^(j-1)
		}	
		
		gammaZ = 0		
		if ("`covs_eval'"=="mean") gammaZ = mean(Z)*gamma_p
				
		*yhat_x = (R_p_l*gamma_p1_l  \ R_p_r*gamma_p1_r ) :+ gammaZ
		*resid_yz = y-Z*gamma_p
					
		y_plot_l = rplot_l*gamma_p1_l :+ gammaZ
		y_plot_r = rplot_r*gamma_p1_r :+ gammaZ
		
		*}
		
		*******************************************************
		**** Optimal Bins (using polynomial order k) **********
		*******************************************************
		rk_l = J(n_l,(k+1),.)
		rk_r = J(n_r,(k+1),.)
		for (j=1; j<=(k+1); j++) {
			rk_l[.,j] = x_l:^(j-1)
			rk_r[.,j] = x_r:^(j-1)
		}
		gamma_k1_l = invsym(cross(rk_l,rk_l))*cross(rk_l,y_l)		
		gamma_k2_l = invsym(cross(rk_l,rk_l))*cross(rk_l,y_l:^2)	
		gamma_k1_r = invsym(cross(rk_r,rk_r))*cross(rk_r,y_r)		
		gamma_k2_r = invsym(cross(rk_r,rk_r))*cross(rk_r,y_r:^2)
			
		*** Bias w/sample
		mu0_k1_l = rk_l*gamma_k1_l
		mu0_k1_r = rk_r*gamma_k1_r
		mu0_k2_l = rk_l*gamma_k2_l
		mu0_k2_r = rk_r*gamma_k2_r
		drk_l = J(n_l,k,.)
		drk_r = J(n_r,k,.)
		for (j=1; j<=k; j++) {
			drk_l[.,j] = j*x_l:^(j-1)
			drk_r[.,j] = j*x_r:^(j-1)
		}
		
	dxi_l=(x_l[2::length(x_l)]-x_l[1::(length(x_l)-1)])
	dxi_r=(x_r[2::length(x_r)]-x_r[1::(length(x_r)-1)])
	dyi_l=(y_l[2::length(y_l)]-y_l[1::(length(y_l)-1)])
	dyi_r=(y_r[2::length(y_r)]-y_r[1::(length(y_r)-1)])
		
	x_bar_i_l = (x_l[2::length(x_l)]+x_l[1::(length(x_l)-1)])/2
	x_bar_i_r = (x_r[2::length(x_r)]+x_r[1::(length(x_r)-1)])/2
		
	drk_i_l = J(n_l-1,k,.);	rk_i_l  = J(n_l-1,(k+1),.)
	drk_i_r = J(n_r-1,k,.);	rk_i_r  = J(n_r-1,(k+1),.)
				   
	for (j=1; j<=(k+1); j++) {
		rk_i_l[.,j] = x_bar_i_l:^(j-1)
		rk_i_r[.,j] = x_bar_i_r:^(j-1)
	}
	  
	  for (j=1; j<=k; j++) {
		drk_i_l[.,j] = j*x_bar_i_l:^(j-1)
		drk_i_r[.,j] = j*x_bar_i_r:^(j-1)
	  }
	  mu1_i_hat_l = drk_i_l*(gamma_k1_l[2::(k+1)])
	  mu1_i_hat_r = drk_i_r*(gamma_k1_r[2::(k+1)])
	   
	  mu0_i_hat_l = rk_i_l*gamma_k1_l
	  mu0_i_hat_r = rk_i_r*gamma_k1_r
	  mu2_i_hat_l = rk_i_l*gamma_k2_l
	  mu2_i_hat_r = rk_i_r*gamma_k2_r

	  mu0_hat_l = rk_l*gamma_k1_l
	  mu0_hat_r = rk_r*gamma_k1_r
	  mu2_hat_l = rk_l*gamma_k2_l
	  mu2_hat_r = rk_r*gamma_k2_r
	  
	  mu1_hat_l = drk_l*(gamma_k1_l[2::(k+1)])
	  mu1_hat_r = drk_r*(gamma_k1_r[2::(k+1)])
		
	  mu1_i_hat_l = drk_i_l*(gamma_k1_l[2::(k+1)])
	  mu1_i_hat_r = drk_i_r*(gamma_k1_r[2::(k+1)])
	  
	  sigma2_hat_l_bar = mu2_i_hat_l - mu0_i_hat_l:^2
	  sigma2_hat_r_bar = mu2_i_hat_r - mu0_i_hat_r:^2

	  sigma2_hat_l = mu2_hat_l - mu0_hat_l:^2
	  sigma2_hat_r = mu2_hat_r - mu0_hat_r:^2
	  
	  var_y_l = variance(y_l)
	  var_y_r = variance(y_r)

	  B_es_hat_dw = (((c-x_min)^2/(12*n))*sum(mu1_hat_l:^2),((x_max-c)^2/(12*n))*sum(mu1_hat_r:^2))
	  V_es_hat_dw = ((0.5/(c-x_min))*sum(dxi_l:*dyi_l:^2),(0.5/(x_max-c))*sum(dxi_r:*dyi_r:^2))
	  V_es_chk_dw = ((1/(c-x_min))*sum(dxi_l:*sigma2_hat_l_bar),(1/(x_max-c))*sum(dxi_r:*sigma2_hat_r_bar))
	  J_es_hat_dw = ceil((((2*B_es_hat_dw):/V_es_hat_dw)*n):^(1/3))
	  J_es_chk_dw = ceil((((2*B_es_hat_dw):/V_es_chk_dw)*n):^(1/3))
		
	  B_qs_hat_dw = ((n_l^2/(24*n))*sum(dxi_l:^2:*mu1_i_hat_l:^2), (n_r^2/(24*n))*sum(dxi_r:^2:*mu1_i_hat_r:^2))
	  V_qs_hat_dw = ((1/(2*n_l))*sum(dyi_l:^2),(1/(2*n_r))*sum(dyi_r:^2))
	  V_qs_chk_dw = ((1/n_l)*sum(sigma2_hat_l), (1/n_r)*sum(sigma2_hat_r))
	  J_qs_hat_dw = ceil((((2*B_qs_hat_dw):/V_qs_hat_dw)*n):^(1/3))
	  J_qs_chk_dw = ceil((((2*B_qs_hat_dw):/V_qs_chk_dw)*n):^(1/3))
	  
	  J_es_hat_mv  = (ceil((var_y_l/V_es_hat_dw[1])*(n/log(n)^2)), ceil((var_y_r/V_es_hat_dw[2])*(n/log(n)^2)))
	  J_es_chk_mv  = (ceil((var_y_l/V_es_chk_dw[1])*(n/log(n)^2)), ceil((var_y_r/V_es_chk_dw[2])*(n/log(n)^2)))
	  J_qs_hat_mv  = (ceil((var_y_l/V_qs_hat_dw[1])*(n/log(n)^2)), ceil((var_y_r/V_qs_hat_dw[2])*(n/log(n)^2)))
	  J_qs_chk_mv  = (ceil((var_y_l/V_qs_chk_dw[1])*(n/log(n)^2)), ceil((var_y_r/V_qs_chk_dw[2])*(n/log(n)^2)))

	if ("`binselect'"=="es" ) {
		J_star_l_orig = J_es_hat_dw[1]
		J_star_r_orig = J_es_hat_dw[2]
	}
	
	if ("`binselect'"=="esmv" | "`binselect'"=="") {
		J_star_l_orig = J_es_hat_mv[1]
		J_star_r_orig = J_es_hat_mv[2]
	}
	
	if ("`binselect'"=="espr" ) {
		J_star_l_orig = J_es_chk_dw[1]
		J_star_r_orig = J_es_chk_dw[2]
	}
	
	if ("`binselect'"=="esmvpr" ) {
		J_star_l_orig = J_es_chk_mv[1]
		J_star_r_orig = J_es_chk_mv[2]
	}
	
	if ("`binselect'"=="qs" ) {
		J_star_l_orig = J_qs_hat_dw[1]
		J_star_r_orig = J_qs_hat_dw[2]
	}
	
	if ("`binselect'"=="qsmv" ) {
		J_star_l_orig = J_qs_hat_mv[1]
		J_star_r_orig = J_qs_hat_mv[2]
	}
	
	if ("`binselect'"=="qspr" ) {
		J_star_l_orig = J_qs_chk_dw[1]
		J_star_r_orig = J_qs_chk_dw[2]
	}
	
	if ("`binselect'"=="qsmvpr" ) {
		J_star_l_orig = J_qs_chk_mv[1]
		J_star_r_orig = J_qs_chk_mv[2]
	}

	if (nbins_l!=0 & nbins_r!=0) {
		J_star_l_orig = nbins_l
		J_star_r_orig = nbins_r
	}
	
	if (`var_l'==0) {
		J_star_l = 1
		J_star_l_orig = 1
		display("{err}Warning: not enough variability in the outcome variable below the threshold")
	}
	if (`var_r'==0) {
		J_star_r = 1
		J_star_r_orig = 1
		display("{err}Warning: not enough variability in the outcome variable above the threshold")
	}

	J_star_l = round(`scale_l'*J_star_l_orig)
	J_star_r = round(`scale_r'*J_star_r_orig)
	
	st_numscalar("nbins_l", nbins_l)
	st_numscalar("nbins_r", nbins_r)
	st_numscalar("J_star_l", J_star_l)
	st_numscalar("J_star_r", J_star_r)
	st_numscalar("J_star_l_orig", J_star_l_orig)
	st_numscalar("J_star_r_orig", J_star_r_orig)
		
	st_matrix("J_es_hat_dw", J_es_hat_dw)
	st_matrix("J_qs_hat_dw", J_qs_hat_dw)
	st_matrix("J_es_chk_dw", J_es_chk_dw)
	st_matrix("J_qs_chk_dw", J_qs_chk_dw)
	st_matrix("J_es_hat_mv", J_es_hat_mv)
	st_matrix("J_qs_hat_mv", J_qs_hat_mv)
	st_matrix("J_es_chk_mv", J_es_chk_mv)
	st_matrix("J_qs_chk_mv", J_qs_chk_mv)
	}
	
	********************************************************
	**** Generate id and rdplot vars ***********************
	********************************************************	
	local J_star_l = J_star_l
	local J_star_r = J_star_r
	
	qui gen rdplot_id = .
	qui gen rdplot_min_bin = .
	qui gen rdplot_max_bin = .
	qui gen rdplot_mean_bin = .
	
	
	if ("`binselect'"=="qs" | "`binselect'"=="qspr"  | "`binselect'"=="qsmv" | "`binselect'"=="qsmvpr") {
		pctile binsL = `x' if `x'<`c',  nq(`J_star_l')
		pctile binsR = `x' if `x'>=`c', nq(`J_star_r')
	}

mata {
	x_min = `x_min'
	x_max = `x_max'

	if ("`binselect'"=="es" | "`binselect'"=="espr"  | "`binselect'"=="esmv" | "`binselect'"=="esmvpr" | "`binselect'"=="") {
		binsL = rangen(x_min,c    , `J_star_l'+1)
		binsR = rangen(c    ,x_max, `J_star_r'+1)
		bins = binsL[1..length(binsL)-1]\binsR		
	}
	
	if ("`binselect'"=="qs" | "`binselect'"=="qspr"  | "`binselect'"=="qsmv" | "`binselect'"=="qsmvpr") {
		bins = (x_min \ st_data(.,"binsL",0) \ c \ st_data(.,"binsR",0) \ x_max )
	}
	
	st_view(ZZ=.,., "`x' rdplot_id rdplot_min_bin rdplot_max_bin rdplot_mean_bin", "`touse'")
	bin_i = 2
	for(i=1; i<=rows(ZZ); i++) {
		while(ZZ[i,1] >= bins[bin_i] & bin_i < length(bins)) bin_i++
		/* PUT rdplot_id */
		ZZ[i,2] = bin_i - `J_star_l' - 2
		if (ZZ[i,2] >= 0) ZZ[i,2] = ZZ[i,2] + 1
		/* PUT rdplot_min_bin rdplot_max_bin rdplot_mean_bin */
		ZZ[i,3] = bins[bin_i-1] 
		ZZ[i,4] = bins[bin_i]
		ZZ[i,5] = (bins[bin_i]+bins[bin_i-1])/2
	}
	
}

** STATA: Generate inputs for RDPLOT (and possibly for reporting back to user)
if  ("`covs_eval'"=="" | "`covs_eval'"=="0") {
collapse (count) rdplot_N=`x' (mean) rdplot_min_bin rdplot_max_bin rdplot_mean_bin    ///
         (mean) rdplot_mean_x=`x'  rdplot_mean_y=`y'    ///
		 (semean) rdplot_se_y=`y', by(rdplot_id) fast
}

**************************************************************************
**** covs_eval **********************************************************
**************************************************************************
if  ("`covs_eval'"=="mean") {
	tempvar  rdplot_id2  yhat_tmp yhatZ
	qui gen `rdplot_id2' = rdplot_id + `J_star_l'
	qui reg `y' `covs_list' i.`rdplot_id2'
	qui predict `yhatZ'
	
	collapse (count) rdplot_N=`x' (mean) rdplot_min_bin rdplot_max_bin rdplot_mean_bin    ///
         (mean) rdplot_mean_x=`x'   rdplot_mean_y=`yhatZ'     ///
		 (semean) rdplot_se_y=`y', by(rdplot_id) fast
}		 

		qui replace rdplot_N=rdplot_N-1
		qui gen quant = -invt(rdplot_N, abs((1-(`ci'/100))/2))
		qui gen rdplot_ci_l = rdplot_mean_y - quant*rdplot_se_y
		qui gen rdplot_ci_r = rdplot_mean_y + quant*rdplot_se_y
		qui drop quant
	
	mata{
	if ("`genvars'"!="") {
	** MATA: Save rdplot inputs to return to user in original dataset
	rdplot = st_data(.,.)
	}
	}
	
	qui gen bin_length = rdplot_max_bin-rdplot_min_bin
	qui su bin_length if rdplot_id<0, d
	local bin_avg_l = r(mean)
	local bin_med_l = r(p50)
	qui su bin_length if rdplot_id>0, d
	local bin_avg_r = r(mean)
	local bin_med_r = r(p50)
		 
	if ("`binselect'"=="es"){
		local binselect_type="evenly spaced number of bins using spacings estimators."
		scalar J_star_l_IMSE = J_es_hat_dw[1,1]
		scalar J_star_r_IMSE = J_es_hat_dw[1,2]
		scalar J_star_l_MV   = J_es_hat_mv[1,1]
		scalar J_star_r_MV   = J_es_hat_mv[1,2]
	}
	if ("`binselect'"=="espr"){
		local binselect_type="evenly spaced number of bins using polynomial regression."
		scalar J_star_l_IMSE = J_es_chk_dw[1,1]
		scalar J_star_r_IMSE = J_es_chk_dw[1,2]
		scalar J_star_l_MV   = J_es_chk_mv[1,1]
		scalar J_star_r_MV   = J_es_chk_mv[1,2]
	}
    if ("`binselect'"=="esmv" | "`binselect'"==""){
		local binselect_type="evenly spaced mimicking variance number of bins using spacings estimators."
		scalar J_star_l_IMSE = J_es_hat_dw[1,1]
		scalar J_star_r_IMSE = J_es_hat_dw[1,2]
		scalar J_star_l_MV   = J_es_hat_mv[1,1]
		scalar J_star_r_MV   = J_es_hat_mv[1,2]
	}
    if ("`binselect'"=="esmvpr"){
		local binselect_type="evenly spaced mimicking variance number of bins using polynomial regression."
		scalar J_star_l_IMSE = J_es_chk_dw[1,1]
		scalar J_star_r_IMSE = J_es_chk_dw[1,2]
		scalar J_star_l_MV   = J_es_chk_mv[1,1]
		scalar J_star_r_MV   = J_es_chk_mv[1,2]
	}
    if ("`binselect'"=="qs"){
		local binselect_type="quantile spaced number of bins using spacings estimators."
		scalar J_star_l_IMSE = J_qs_hat_dw[1,1]
		scalar J_star_r_IMSE = J_qs_hat_dw[1,2]
		scalar J_star_l_MV   = J_qs_hat_mv[1,1]
		scalar J_star_r_MV   = J_qs_hat_mv[1,2]
	}
    if ("`binselect'"=="qspr"){
		local binselect_type="quantile spaced number of bins using polynomial regression."
		scalar J_star_l_IMSE = J_qs_chk_dw[1,1]
		scalar J_star_r_IMSE = J_qs_chk_dw[1,2]
		scalar J_star_l_MV   = J_qs_chk_mv[1,1]
		scalar J_star_r_MV   = J_qs_chk_mv[1,2]
	}
    if ("`binselect'"=="qsmv"){
		local binselect_type="quantile spaced mimicking variance quantile spaced using spacings estimators."  
		scalar J_star_l_IMSE = J_qs_hat_dw[1,1]
		scalar J_star_r_IMSE = J_qs_hat_dw[1,2]
		scalar J_star_l_MV   = J_qs_hat_mv[1,1]
		scalar J_star_r_MV   = J_qs_hat_mv[1,2]
	}
    if ("`binselect'"=="qsmvpr"){
		local binselect_type="quantile spaced mimicking variance number of bins using polynomial regression."
		scalar J_star_l_IMSE = J_qs_chk_dw[1,1]
		scalar J_star_r_IMSE = J_qs_chk_dw[1,2]
		scalar J_star_l_MV   = J_qs_chk_mv[1,1]
		scalar J_star_r_MV   = J_qs_chk_mv[1,2]
	}
	if (nbins_l!=0 | nbins_r!=0 ) local binselect_type= "RD plot with manually set number of bins."
	
	scalar scale_l = J_star_l / J_star_l_IMSE
	scalar scale_r = J_star_r / J_star_r_IMSE
		
	qui getmata x_plot_l x_plot_r y_plot_l y_plot_r, force
		
	ereturn clear
	ereturn scalar N_l = `n_l'
	ereturn scalar N_r = `n_r'
	ereturn scalar c = `c'
	ereturn scalar J_star_l = J_star_l
	ereturn scalar J_star_r = J_star_r
	ereturn matrix coef_l = gamma_p1_l
	ereturn matrix coef_r = gamma_p1_r
	ereturn local binselect = "`binselect'"

	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") local kernel_type = "Epanechnikov"
	else if ("`kernel'"=="uniform" | "`kernel'"=="uni") local kernel_type = "Uniform"
	else  local kernel_type = "Triangular"
			
	disp ""
	disp in smcl in yellow "RD Plot with " "`binselect_type'" 
	disp ""
	
	disp in smcl in gr "{ralign 21: Cutoff c = `c'}"      _col(22) " {c |} " _col(23) in gr "Left of " in yellow "c"  _col(36) in gr "Right of " in yellow "c" _col(54) in gr "Number of obs  = "  in yellow %10.0f `n'
	disp in smcl in gr "{hline 22}{c +}{hline 22}"                                                                                                             _col(54) in gr "Kernel         = "  in yellow "{ralign 10:`kernel_type'}" 
	disp in smcl in gr "{ralign 21:Number of obs}"        _col(22) " {c |} " _col(23) as result %9.0f `n_l'      _col(37) %9.0f  `n_r'                        
	disp in smcl in gr "{ralign 21:Eff. Number of obs}"   _col(22) " {c |} " _col(23) as result %9.0f `n_h_l'      _col(37) %9.0f  `n_h_r'                        
	disp in smcl in gr "{ralign 21:Order poly. fit (p)}"  _col(22) " {c |} " _col(23) as result %9.0f `p'        _col(37) %9.0f  `p'                              
	disp in smcl in gr "{ralign 21:BW poly. fit (h)}"     _col(22) " {c |} " _col(23) as result %9.3f `h_l'      _col(37) %9.3f  `h_r' 
	disp in smcl in gr "{ralign 21:Number of bins scale}" _col(22) " {c |} " _col(23) as result %9.3f `scale_l'  _col(37) %9.3f  `scale_r' 
	disp ""
	disp "Outcome: `y'. Running variable: `x'."
	disp in smcl in gr "{hline 22}{c TT}{hline 22}"
	disp in smcl in gr   _col(22) " {c |} " _col(23) in gr "Left of " in yellow "c"  _col(36) in gr "Right of " in yellow "c" 
	disp in smcl in gr "{hline 22}{c +}{hline 22}"                                                                                
	disp in smcl in gr "{ralign 21:Bins selected}"          _col(22) " {c |} " _col(23) as result %9.0f e(J_star_l)      _col(37) %9.0f  e(J_star_r)
	disp in smcl in gr "{ralign 21:Average bin length}"    _col(22) " {c |} " _col(23) as result %9.3f `bin_avg_l'   _col(37) %9.3f  `bin_avg_r'
	disp in smcl in gr "{ralign 21:Median bin length}"     _col(22) " {c |} " _col(23) as result %9.3f `bin_med_l'   _col(37) %9.3f  `bin_med_r'
	disp in smcl in gr "{hline 22}{c +}{hline 22}"  
	disp in smcl in gr "{ralign 21:IMSE-optimal bins}"      _col(22) " {c |} " _col(23) as result %9.0f J_star_l_IMSE    _col(37) %9.0f  J_star_r_IMSE
	disp in smcl in gr "{ralign 21:Mimicking Var. bins}"    _col(22) " {c |} " _col(23) as result %9.0f J_star_l_MV      _col(37) %9.0f  J_star_r_MV
	disp in smcl in gr "{hline 22}{c +}{hline 22}"    
	disp in smcl in gr "{lalign 1:Rel. to IMSE-optimal:}"   _col(22) " {c |} " 
	disp in smcl in gr "{ralign 21:Implied scale}"          _col(22) " {c |} " _col(23) as result %9.3f scale_l                 _col(37) %9.3f  scale_r
	disp in smcl in gr "{ralign 21:WIMSE var. weight}"      _col(22) " {c |} " _col(23) as result %9.3f 1/(1+scale_l^3)         _col(37) %9.3f  1/(1+scale_r^3)
   	disp in smcl in gr "{ralign 21:WIMSE bias weight}"      _col(22) " {c |} " _col(23) as result %9.3f scale_l^3/(1+scale_l^3) _col(37) %9.3f  scale_r^3/(1+scale_r^3)
	disp in smcl in gr "{hline 22}{c BT}{hline 22}"    	
	disp ""
		if ("`covs'"!="")    disp "Covariate-adjusted estimates. Additional covariates included: `ncovs'"
		if (`covs_drop_coll'==1) di as error  "{err}Variables dropped due to multicollinearity."


	if ("`hide'"==""){
		if (`"`graph_options'"'=="" ) local graph_options = `"title("Regression function fit", color(gs0)) "'
		
		if (`ci'==0) {
				twoway (scatter rdplot_mean_y rdplot_mean_bin, sort msize(small)  mcolor(gs10)) ///
				(line y_plot_l x_plot_l, lcolor(black) sort lwidth(medthin) lpattern(solid) ) ///
				(line y_plot_r x_plot_r, lcolor(black) sort lwidth(medthin) lpattern(solid) ),  ///
				xline(`c', lcolor(black) lwidth(medthin)) xscale(r(`x_min' `x_max'))  legend(cols(2) order(1 "Sample average within bin" 2 "Polynomial fit of order `p'" )) `graph_options'
		}
		else {
			if ("`shade'"==""){
				twoway (rcap rdplot_ci_l rdplot_ci_r rdplot_mean_bin, color(gs11)) ///
				(scatter rdplot_mean_y rdplot_mean_bin, sort msize(small)  mcolor(gs10)) ///
				(line y_plot_l x_plot_l, lcolor(black) sort lwidth(medthin) lpattern(solid))   ///
				(line y_plot_r x_plot_r, lcolor(black) sort lwidth(medthin) lpattern(solid)),  ///
				xline(`c', lcolor(black) lwidth(medthin)) xscale(r(`x_min' `x_max')) legend(cols(2) order(2 "Sample average within bin" 3 "Polynomial fit of order `p'" )) `graph_options'
			}
			else {
				twoway (rarea rdplot_ci_l rdplot_ci_r rdplot_mean_bin if rdplot_id<0, sort color(gs11)) ///
				       (rarea rdplot_ci_l rdplot_ci_r rdplot_mean_bin if rdplot_id>0, sort color(gs11)) ///
				(scatter rdplot_mean_y rdplot_mean_bin, sort msize(small)  mcolor(gs10)) ///
				(line y_plot_l x_plot_l, lcolor(black) sort lwidth(medthin) lpattern(solid)) ///
				(line y_plot_r x_plot_r, lcolor(black) sort lwidth(medthin) lpattern(solid)) ,  ///
				xline(`c', lcolor(black) lwidth(medthin)) xscale(r(`x_min' `x_max')) legend(cols(2) order(2 "Sample average within bin" 3 "Polynomial fit of order `p'" )) `graph_options'
			}						
		}
	}
	
restore

****************************
** PART 2: genvars=TRUE
****************************
if ("`genvars'"!="") {
	qui for any id N min_bin max_bin mean_bin mean_x mean_y se_y ci_l ci_r hat_y: qui gen rdplot_X = .
}
	
mata{
	if ("`genvars'"~="") {
		st_view(ZZ=.,., "`x' rdplot_id rdplot_N rdplot_min_bin rdplot_max_bin rdplot_mean_bin rdplot_mean_x rdplot_mean_y rdplot_se_y rdplot_ci_l rdplot_ci_r rdplot_hat_y", "`touse'")
		for (i=1; i<=rows(ZZ); i++) {
		if (ZZ[i,1]!=.) {
			bin_i = 2; while(ZZ[i,1] >= bins[bin_i] & bin_i < length(bins)) bin_i++
		 	rdplot_i = bin_i - `J_star_l' - 2
			if (rdplot_i >= 0) rdplot_i = rdplot_i + 1
			ZZ[i,2..11] = select(rdplot, rdplot[.,1]:==rdplot_i)
			ZZ[i,12] = 0; for (j=0; j<=p; j++) {
			if (ZZ[i,2] <0) ZZ[i,12] = ZZ[i,12] + ((ZZ[i,1]-c)^j)*gamma_p1_l[j+1] 
			else           ZZ[i,12] = ZZ[i,12] + ((ZZ[i,1]-c)^j)*gamma_p1_r[j+1]
			}		
		}
		}
	}
}

mata mata clear
end


 
