rdrobust = function(y, x, c = NULL, fuzzy = NULL, deriv = NULL,  
                    p = NULL, q = NULL, h = NULL, b = NULL, rho = NULL, 
                    covs = NULL, covs_drop = TRUE, ginv.tol = 1e-20,
                    kernel = "tri", weights = NULL, bwselect = "mserd",
                    vce = "nn", cluster = NULL, nnmatch = 3, level = 95, 
                    scalepar = 1, scaleregul = 1, sharpbw = FALSE, 
                    all = NULL, subset = NULL, masspoints = "adjust",
                    bwcheck = NULL, bwrestrict=TRUE, stdvars=FALSE) {
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  if (is.null(c)) c <- 0
  
  # p
  if (is.null(p) & !is.null(deriv)) {p = deriv+1}
  
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 1
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }
  
  # q
  if (length(q) == 0) {
    flag_no_q <- TRUE
    q <- p + 1
  } else if ((length(q) > 1) | !(q[1]%in%c(0:20)) | (q[1]<p)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  } else {
    flag_no_q <- FALSE
  }
  
  # deriv
  if (length(deriv) == 0) {
    flag_no_deriv <- TRUE
    deriv <- 0
  } else if ((length(deriv) > 1) | !(deriv[1]%in%c(0:20)) | (deriv[1]>p)) {
    stop("Derivative order incorrectly specified.\n")
  } else {
    flag_no_deriv <- FALSE
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- subset(covs,subset)
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  if (!is.null(fuzzy)){
    if (!is.null(subset)) fuzzy <- fuzzy[subset]
    na.ok <- na.ok & complete.cases(fuzzy)
  } 

  if (!is.null(weights)){
    if (!is.null(subset)) weights <- weights[subset]
    na.ok <- na.ok & complete.cases(weights) & weights>=0
  } 
  
  x = as.matrix(x[na.ok])
  y = as.matrix(y[na.ok])
  
  if (!is.null(covs))    covs    = as.matrix(covs)[na.ok, , drop = FALSE]
  if (!is.null(fuzzy))   fuzzy   = as.matrix(  fuzzy[na.ok])
  if (!is.null(cluster)) cluster = as.matrix(cluster[na.ok])
  if (!is.null(weights)) weights = as.matrix(weights[na.ok])
  
  if (is.null(masspoints)) masspoints=FALSE

  if (vce=="nn" | masspoints=="check" |masspoints=="adjust") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    if (!is.null(covs))    covs    =  as.matrix(covs)[order_x,,drop=FALSE]
    if (!is.null(fuzzy))   fuzzy   =   fuzzy[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
    if (!is.null(weights)) weights = weights[order_x,,drop=FALSE]
  }

  kernel   = tolower(kernel)
  bwselect = tolower(bwselect)
  vce      = tolower(vce)
  
  X_l = x[x<c,,drop=FALSE];  X_r = x[x>=c,,drop=FALSE]
  Y_l = y[x<c,,drop=FALSE];  Y_r = y[x>=c,,drop=FALSE]
  x_min = min(x);  x_max = max(x)
  range_l = abs(max(X_l)-min(X_l));   range_r = abs(max(X_r)-min(X_r))
  N_l = length(X_l);   N_r = length(X_r)
  N = N_r + N_l
  quant = -qnorm(abs((1-(level/100))/2))
  
  vce_type = "NN"
  if (vce=="hc0")         vce_type = "HC0"
  if (vce=="hc1")         vce_type = "HC1"
  if (vce=="hc2")         vce_type = "HC2"
  if (vce=="hc3")      	  vce_type = "HC3"
  if (!is.null(cluster))	vce_type = "Cluster"

  
  ############## COLLINEARITY
  covs_drop_coll=dZ=0
  if (covs_drop == TRUE) covs_drop_coll = 1 
  if (!is.null(covs)) {
    covs.names = colnames(covs)
    if (is.null(covs.names)) {
      covs.names = paste("z",1:ncol(covs),sep="")
      colnames(covs) = covs.names
    }
    covs = covs[,order(nchar(covs.names))]
    covs = as.matrix(covs)
    dZ = length(covs.names)
    covs.check = covs_drop_fun(covs)
    if (covs.check$ncovs < dZ & covs_drop==FALSE) {
      print("Multicollinearity issue detected in covs. Please rescale and/or remove redundant covariates, or use covs_drop option.")  
    }
    if (covs.check$ncovs < dZ & isTRUE(covs_drop)) {
      covs  <- as.matrix(covs.check$covs)
      dZ    <- covs.check$ncovs
      #covs_drop_coll <-1
      #print("Multicollinearity issue detected in covs. Redundant covariates dropped.")  
    }
  }

  #####################################################   CHECK ERRORS
  exit=0
  if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit = 1
  }
  
  if  (bwselect!="mserd" & bwselect!="msetwo" & bwselect!="msesum" & bwselect!="msecomb1" & bwselect!="msecomb2" & bwselect!="cerrd" & bwselect!="certwo" & bwselect!="cersum" & bwselect!="cercomb1" & bwselect!="cercomb2" & bwselect!=""){
    print("bwselect incorrectly specified")  
    exit = 1
  }
  
  if (bwselect=="cct" | bwselect=="ik" | bwselect=="cv"){
    print("bwselect options IK, CCT and CV have been depricated. Please see help for new options")  
    exit = 1
  }
  
  if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
    print("vce incorrectly specified")
    exit = 1
  }
    
    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit = 1
    }
    
    if (!is.null(rho)){  
       if (rho<0){
          print("rho should be greater than 0")
          exit = 1
        }
    }
  
    if (exit>0) stop()
    if (!is.null(h)) bwselect = "Manual"
    if (!is.null(h) & is.null(rho) & is.null(b)) {
      rho = 1
      b = h
    }
    if (!is.null(h) & !is.null(rho) ) b = h/rho
    
  
  if (N<20){
			print("Not enough observations to perform bandwidth calculations. Estimates computed using entire sample")
      h = b = max(range_l,range_r)
			bwselect = "Manual"
			}
  
  
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
  }   else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
  }   else  {
    kernel_type = "Triangular"
  }

  
  
  
  
  
  vce_type = "NN"
  if (vce=="hc0")     		vce_type = "HC0"
  if (vce=="hc1")      	  vce_type = "HC1"
  if (vce=="hc2")      	  vce_type = "HC2"
  if (vce=="hc3")      	  vce_type = "HC3"
  if (vce=="cluster")  	  vce_type = "Cluster"
  if (vce=="nncluster") 	vce_type = "NNcluster"
  
  
  mN = N
  M_l = N_l
  M_r = N_r
  if (masspoints=="check" | masspoints=="adjust") {
    X_uniq_l = sort(unique(X_l), decreasing=TRUE)
    X_uniq_r = unique(X_r)
    M_l = length(X_uniq_l)
    M_r = length(X_uniq_r)
    M = M_l + M_r
    mass_l = 1-M_l/N_l
    mass_r = 1-M_r/N_r				
    if (mass_l>=0.1 | mass_r>=0.1){
      print("Mass points detected in the running variable.")
      if (masspoints=="check") print("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & masspoints=="adjust") bwcheck <- 10
    }				
  }
		
		
    ############################################################################################
    #print("Preparing data.") 
    
    if (is.null(h)) {
      invisible(capture.output( rdbws<- rdbwselect(y=y, x=x, c=c, fuzzy=fuzzy,  deriv=deriv, p=p, q=q, 
                                                   covs=covs, covs_drop=covs_drop,
                       kernel=kernel,  weights=weights, bwselect=bwselect,  bwcheck = bwcheck, bwrestrict=bwrestrict,
                       vce=vce, cluster=cluster,  nnmatch=nnmatch,  scaleregul=scaleregul,
                       sharpbw = sharpbw, subset=subset, all=FALSE, masspoints=masspoints, stdvars=stdvars, prchk=FALSE)))
      
      h_l = c(rdbws$bws[1]); b_l = c(rdbws$bws[3])
      h_r = c(rdbws$bws[2]); b_r = c(rdbws$bws[4])
      
      if (!is.null(rho)) {
        b_l = h_l/rho
    		b_r = h_r/rho
      }
      
    } else{
      if (length(h)==1) h_l = h_r = h
      if (length(h)==2) {
        h_l = h[1]
        h_r = h[2]
      }
      if (is.null(b)) {
        b_l = h_l
        b_r = h_r
      } else {
        if (length(b)==1) b_l = b_r = b
        if (length(b)==2) {
          b_l = b[1]
          b_r = b[2]
        }  
      }  
    }

  w_h_l <- rdrobust_kweight(X_l,c,h_l,kernel);	w_h_r <- rdrobust_kweight(X_r,c,h_r,kernel)
  w_b_l <- rdrobust_kweight(X_l,c,b_l,kernel);	w_b_r <- rdrobust_kweight(X_r,c,b_r,kernel)
  
  if (!is.null(weights)) {
    fw_l <- weights[x<c,, drop=FALSE]  
    fw_r <- weights[x>=c,,drop=FALSE]
    w_h_l <- fw_l*w_h_l;	w_h_r <- fw_r*w_h_r
    w_b_l <- fw_l*w_b_l;	w_b_r <- fw_r*w_b_r			
  }

  ind_h_l <- w_h_l> 0;		ind_h_r <- w_h_r> 0
  ind_b_l <- w_b_l> 0;		ind_b_r <- w_b_r> 0
  N_h_l <- sum(ind_h_l); N_b_l <- sum(ind_b_l)
  N_h_r <- sum(ind_h_r); N_b_r <- sum(ind_b_r)
  
  #if (N_h_l<5 | N_h_r<5 | N_b_l<5 | N_b_r<5){
  #  stop("Not enough observations to perform calculations")
  #  exit(1)
  #}
  
  ind_l = ind_b_l; ind_r = ind_b_r
  if (h_l>b_l) ind_l = ind_h_l   
  if (h_r>b_r) ind_r = ind_h_r   
  
  eN_l = sum(ind_l); eN_r = sum(ind_r)
  eY_l  = Y_l[ind_l,,drop=FALSE];	eY_r  = Y_r[ind_r,,drop=FALSE]
  eX_l  = X_l[ind_l,,drop=FALSE];	eX_r  = X_r[ind_r,,drop=FALSE]
  W_h_l = w_h_l[ind_l];	W_h_r = w_h_r[ind_r]
  W_b_l = w_b_l[ind_l];	W_b_r = w_b_r[ind_r]
  
  edups_l = edupsid_l = matrix(0,eN_l,1)
  edups_r = edupsid_r = matrix(0,eN_r,1)
 
  if (vce=="nn") {
    for (i in 1:eN_l) edups_l[i]=sum(eX_l==eX_l[i])
    for (i in 1:eN_r) edups_r[i]=sum(eX_r==eX_r[i])
    i=1
    while (i<=eN_l) {
      edupsid_l[i:(i+edups_l[i]-1)] = 1:edups_l[i]
      i = i+edups_l[i]
    }
    i=1
    while (i<=eN_r) {
      edupsid_r[i:(i+edups_r[i]-1)]=1:edups_r[i]
      i=i+edups_r[i]
    }
  }          
          
  u_l <- (eX_l-c)/h_l;	u_r <-(eX_r-c)/h_r
  R_q_l = matrix(NA,eN_l,(q+1)); R_q_r = matrix(NA,eN_r,(q+1))
  for (j in 1:(q+1))  {
    R_q_l[,j] = (eX_l-c)^(j-1);  R_q_r[,j] = (eX_r-c)^(j-1)
  }
  R_p_l = R_q_l[,1:(p+1)]; R_p_r = R_q_r[,1:(p+1)]

  #display("Computing RD estimates.")
  L_l = crossprod(R_p_l*W_h_l,u_l^(p+1)); L_r = crossprod(R_p_r*W_h_r,u_r^(p+1)) 
  invG_q_l  = qrXXinv((sqrt(W_b_l)*R_q_l));	invG_q_r  = qrXXinv((sqrt(W_b_r)*R_q_r))
  invG_p_l  = qrXXinv((sqrt(W_h_l)*R_p_l));	invG_p_r  = qrXXinv((sqrt(W_h_r)*R_p_r))
  e_p1 = matrix(0,(q+1),1); e_p1[p+2]=1
  e_v  = matrix(0,(p+1),1); e_v[deriv+1]=1
  Q_q_l = t(t(R_p_l*W_h_l) - h_l^(p+1)*(L_l%*%t(e_p1))%*%t(t(invG_q_l%*%t(R_q_l))*W_b_l))
  Q_q_r = t(t(R_p_r*W_h_r) - h_r^(p+1)*(L_r%*%t(e_p1))%*%t(t(invG_q_r%*%t(R_q_r))*W_b_r))
  D_l = eY_l; D_r = eY_r
  eC_l=eC_r=eT_l=eT_r=eZ_l=eZ_r=NULL
  dT =g_l=g_r= 0
  
  if (!is.null(fuzzy)) {
    dT=1
    T_l  = fuzzy[x<c,,drop=FALSE];  eT_l  = T_l[ind_l,,drop=FALSE]
    T_r  = fuzzy[x>=c,,drop=FALSE]; eT_r  = T_r[ind_r,,drop=FALSE]
    D_l  = cbind(D_l,eT_l); D_r = cbind(D_r,eT_r)
  }
  
  if (!is.null(covs)) {
    Z_l  = covs[x<c,,drop=FALSE];	  eZ_l = Z_l[ind_l,,drop=FALSE]
    Z_r  = covs[x>=c,,drop=FALSE];	eZ_r = Z_r[ind_r,,drop=FALSE]
    D_l  = cbind(D_l,eZ_l); D_r = cbind(D_r,eZ_r)
    U_p_l = crossprod(R_p_l*W_h_l,D_l); U_p_r = crossprod(R_p_r*W_h_r,D_r)
  }
              
  if (!is.null(cluster)) {
    C_l  = cluster[x<c,,drop=FALSE]; C_r= cluster[x>=c,,drop=FALSE]
    eC_l  = C_l[ind_l];	     eC_r  = C_r[ind_r]
    g_l = length(unique(eC_l));	g_r = length(unique(eC_r))
  }
                                                     
  beta_p_l = invG_p_l%*%crossprod(R_p_l*W_h_l,D_l); beta_q_l = invG_q_l%*%crossprod(R_q_l*W_b_l,D_l); beta_bc_l = invG_p_l%*%crossprod(Q_q_l,D_l) 
  beta_p_r = invG_p_r%*%crossprod(R_p_r*W_h_r,D_r); beta_q_r = invG_q_r%*%crossprod(R_q_r*W_b_r,D_r); beta_bc_r = invG_p_r%*%crossprod(Q_q_r,D_r)
  beta_p  = beta_p_r  - beta_p_l
  beta_q  = beta_q_r  - beta_q_l
  beta_bc = beta_bc_r - beta_bc_l

  if (is.null(covs)) {	
  tau_cl = tau_Y_cl = scalepar*factorial(deriv)*beta_p[(deriv+1),1]
  tau_bc = tau_Y_bc = scalepar*factorial(deriv)*beta_bc[(deriv+1),1]
  s_Y = 1
  
  tau_Y_cl_l = scalepar*factorial(deriv)*beta_p_l[(deriv+1),1]
  tau_Y_cl_r = scalepar*factorial(deriv)*beta_p_r[(deriv+1),1]
  tau_Y_bc_l = scalepar*factorial(deriv)*beta_bc_l[(deriv+1),1]
  tau_Y_bc_r = scalepar*factorial(deriv)*beta_bc_r[(deriv+1),1]
  bias_l = tau_Y_cl_l-tau_Y_bc_l
  bias_r = tau_Y_cl_r-tau_Y_bc_r 
  
  if (!is.null(fuzzy)) {
     tau_T_cl = factorial(deriv)*beta_p[(deriv+1),2]
     tau_T_bc = factorial(deriv)*beta_bc[(deriv+1),2]
     tau_cl   = tau_Y_cl/tau_T_cl
     s_Y      = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
     B_F      = c(tau_Y_cl-tau_Y_bc , tau_T_cl-tau_T_bc)
     tau_bc   = tau_cl - t(s_Y)%*%B_F
     sV_T     = c(0 , 1)
     
     tau_T_cl_l = factorial(deriv)*beta_p_l[(deriv+1),2]
     tau_T_cl_r = factorial(deriv)*beta_p_l[(deriv+1),2]
     tau_T_bc_l = factorial(deriv)*beta_bc_l[(deriv+1),2]
     tau_T_bc_r = factorial(deriv)*beta_bc_r[(deriv+1),2]
     B_F_l = c(tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l)
     B_F_r = c(tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r)
     bias_l = t(s_Y)%*%B_F_l
		 bias_r = t(s_Y)%*%B_F_r
  }	
  } else {	
    ZWD_p_l  = crossprod(eZ_l*W_h_l,D_l)
    ZWD_p_r  = crossprod(eZ_r*W_h_r,D_r)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU_p_l =  crossprod(matrix(U_p_l[,colsZ],nrow=p+1),invG_p_l%*%U_p_l) 
    UiGU_p_r =  crossprod(matrix(U_p_r[,colsZ],nrow=p+1),invG_p_r%*%U_p_r) 
    ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
    ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
    ZWY_p_l = ZWD_p_l[,1:(1+dT)] - UiGU_p_l[,1:(1+dT)] 
    ZWY_p_r = ZWD_p_r[,1:(1+dT)] - UiGU_p_r[,1:(1+dT)]     
    ZWZ_p = ZWZ_p_r + ZWZ_p_l
    ZWY_p = ZWY_p_r + ZWY_p_l
    if (covs_drop_coll == 0) gamma_p = chol2inv(chol(ZWZ_p))%*%ZWY_p
    if (covs_drop_coll == 1) gamma_p = ginv(ZWZ_p, tol = ginv.tol)%*%ZWY_p
    s_Y = c(1 ,  -gamma_p[,1])
    
    if (is.null(fuzzy)) {
        tau_cl = scalepar*t(s_Y)%*%beta_p[(deriv+1),]
        tau_bc = scalepar*t(s_Y)%*%beta_bc[(deriv+1),]
        
        tau_Y_cl_l = scalepar*t(s_Y)%*%beta_p_l[(deriv+1),]
        tau_Y_cl_r = scalepar*t(s_Y)%*%beta_p_r[(deriv+1),]
        tau_Y_bc_l = scalepar*t(s_Y)%*%beta_bc_l[(deriv+1),]
        tau_Y_bc_r = scalepar*t(s_Y)%*%beta_bc_r[(deriv+1),]
        bias_l = tau_Y_cl_l-tau_Y_bc_l
        bias_r = tau_Y_cl_r-tau_Y_bc_r 

    } else {
      s_T  = c(1,    -gamma_p[,2])
      sV_T = c(0, 1, -gamma_p[,2])
      tau_Y_cl = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p[(deriv+1),1], beta_p[(deriv+1),colsZ])
      tau_T_cl = factorial(deriv)*t(s_T)%*%c(beta_p[(deriv+1),2], beta_p[(deriv+1),colsZ])
      tau_Y_bc = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc[(deriv+1),1],beta_bc[(deriv+1),colsZ])
      tau_T_bc = factorial(deriv)*t(s_T)%*%c(beta_bc[(deriv+1),2],beta_bc[(deriv+1),colsZ])

      tau_Y_cl_l = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p_l[(deriv+1),1], beta_p_l[(deriv+1),colsZ])
      tau_Y_cl_r = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p_r[(deriv+1),2], beta_p_r[(deriv+1),colsZ])
      tau_Y_bc_l = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc_l[(deriv+1),1],beta_bc_l[(deriv+1),colsZ])
      tau_Y_bc_r = scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc_r[(deriv+1),2],beta_bc_r[(deriv+1),colsZ])

      tau_T_cl_l = factorial(deriv)*t(s_T)%*%c(beta_p_l[(deriv+1),1], beta_p_l[(deriv+1),colsZ])
      tau_T_cl_r = factorial(deriv)*t(s_T)%*%c(beta_p_r[(deriv+1),2], beta_p_r[(deriv+1),colsZ])
      tau_T_bc_l = factorial(deriv)*t(s_T)%*%c(beta_bc_l[(deriv+1),1],beta_bc_l[(deriv+1),colsZ])
      tau_T_bc_r = factorial(deriv)*t(s_T)%*%c(beta_bc_r[(deriv+1),2],beta_bc_r[(deriv+1),colsZ])

      tau_cl = tau_Y_cl/tau_T_cl
      B_F   = c(tau_Y_cl-tau_Y_bc,     tau_T_cl-tau_T_bc)
      B_F_l = c(tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l)
      B_F_r = c(tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r)
      
      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
      tau_bc = tau_cl - t(s_Y)%*%B_F
      
      bias_l = t(s_Y)%*%B_F_l
      bias_r = t(s_Y)%*%B_F_r

      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2) , -(1/tau_T_cl)*gamma_p[,1] + (tau_Y_cl/tau_T_cl^2)*gamma_p[,2])
    }
  }

#*display("Computing variance-covariance matrix.")
  
  hii_l=hii_r=predicts_p_l=predicts_p_r=predicts_q_l=predicts_q_r=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_p_l=R_p_l%*%beta_p_l
    predicts_p_r=R_p_r%*%beta_p_r
    predicts_q_l=R_q_l%*%beta_q_l
    predicts_q_r=R_q_r%*%beta_q_r
    if (vce=="hc2" | vce=="hc3") {
      hii_l=matrix(NA,eN_l,1)	
      for (i in 1:eN_l) hii_l[i] = R_p_l[i,]%*%invG_p_l%*%(R_p_l*W_h_l)[i,]
      hii_r=matrix(NA,eN_r,1)	
      for (i in 1:eN_r) hii_r[i] = R_p_r[i,]%*%invG_p_r%*%(R_p_r*W_h_r)[i,]
    }
  }
  						
	res_h_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_p_l, hii_l, vce, nnmatch, edups_l, edupsid_l, p+1)
	res_h_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_p_r, hii_r, vce, nnmatch, edups_r, edupsid_r, p+1)
	if (vce=="nn") {
			res_b_l = res_h_l;	res_b_r = res_h_r
	} 	else {
			res_b_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_q_l, hii_l, vce, nnmatch, edups_l, edupsid_l, q+1)
			res_b_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_q_r, hii_r, vce, nnmatch, edups_r, edupsid_r, q+1)
  }
			                       
	V_Y_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(R_p_l*W_h_l), res_h_l, eC_l)%*%invG_p_l
	V_Y_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(R_p_r*W_h_r), res_h_r, eC_r)%*%invG_p_r
	V_Y_rb_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(Q_q_l),       res_b_l, eC_l)%*%invG_p_l
	V_Y_rb_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(Q_q_r),       res_b_r, eC_r)%*%invG_p_r
	V_tau_cl = scalepar^2*factorial(deriv)^2*(V_Y_cl_l+V_Y_cl_r)[deriv+1,deriv+1]
	V_tau_rb = scalepar^2*factorial(deriv)^2*(V_Y_rb_l+V_Y_rb_r)[deriv+1,deriv+1]
	se_tau_cl = sqrt(V_tau_cl);	se_tau_rb = sqrt(V_tau_rb)

	if (!is.null(fuzzy)) {
		V_T_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(R_p_l*W_h_l), res_h_l, eC_l)%*%invG_p_l
		V_T_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(R_p_r*W_h_r), res_h_r, eC_r)%*%invG_p_r
		V_T_rb_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(Q_q_l), res_b_l, eC_l)%*%invG_p_l
		V_T_rb_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(Q_q_r), res_b_r, eC_r)%*%invG_p_r
		V_T_cl = factorial(deriv)^2*(V_T_cl_l+V_T_cl_r)[deriv+1,deriv+1]
		V_T_rb = factorial(deriv)^2*(V_T_rb_l+V_T_rb_r)[deriv+1,deriv+1]
		se_tau_T_cl = sqrt(V_T_cl);	se_tau_T_rb = sqrt(V_T_rb)
	}
  
  tau = c(tau_cl, tau_bc, tau_bc)
  se  = c(se_tau_cl,se_tau_cl,se_tau_rb)
  t   =  tau/se
  pv  = 2*pnorm(-abs(t))
  ci  = matrix(NA,nrow=3,ncol=2)
  rownames(ci) = c("Conventional","Bias-Corrected","Robust")
  colnames(ci) = c("Lower","Upper")
  ci[1,] = c(tau[1] - quant*se[1], tau[1] + quant*se[1])
  ci[2,] = c(tau[2] - quant*se[2], tau[2] + quant*se[2])
  ci[3,] = c(tau[3] - quant*se[3], tau[3] + quant*se[3])
    
  if (!is.null(fuzzy)) {  
      tau_T = c(tau_T_cl, tau_T_bc, tau_T_bc)
      se_T  = c(se_tau_T_cl, se_tau_T_cl, se_tau_T_rb)
      t_T   = tau_T/se_T
      pv_T  = 2*pnorm(-abs(t_T))
      ci_T  = matrix(NA,nrow=3,ncol=2)
      ci_T[1,] = c(tau_T[1] - quant*se_T[1], tau_T[1] + quant*se_T[1])
      ci_T[2,] = c(tau_T[2] - quant*se_T[2], tau_T[2] + quant*se_T[2])
      ci_T[3,] = c(tau_T[3] - quant*se_T[3], tau_T[3] + quant*se_T[3])
  }

    coef = matrix(tau,3,1)
    se   = matrix(se, 3,1)
    z    = matrix(t,  3,1)
    pv   = matrix(pv, 3,1)
    ci   = ci

  bws=matrix(c(h_l,b_l,h_r,b_r),2,2)
  colnames(bws)=c("left","right")
  rownames(bws)=c("h","b")
  
  rownames(coef)=rownames(se)=rownames(se)=rownames(z)=rownames(pv)=c("Conventional","Bias-Corrected","Robust")
  colnames(coef)="Coeff"
  colnames(se)="Std. Err."
  colnames(z)="z"
  colnames(pv)="P>|z|"
  colnames(bws)=c("left","right")
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("CI Lower","CI Upper")
    
  Estimate=matrix(NA,1,4)
  colnames(Estimate)=c("tau.us","tau.bc","se.us","se.rb")
  Estimate[1,] <- c(tau_cl,tau_bc, se_tau_cl, se_tau_rb) 

  out=list(Estimate=Estimate, bws=bws, coef=coef, se=se, z=z, pv=pv, ci=ci,
           beta_p_l=beta_p_l[,1], beta_p_r=beta_p_r[,1],
           V_cl_l=V_Y_cl_l, V_cl_r=V_Y_cl_r, V_rb_l=V_Y_rb_l, V_rb_r=V_Y_rb_r,
           N=c(N_l,N_r), N_h=c(N_h_l,N_h_r), N_b=c(N_b_l,N_b_r), M=c(M_l,M_r),
           tau_cl=c(tau_Y_cl_l,tau_Y_cl_r), tau_bc=c(tau_Y_bc_l,tau_Y_bc_r),
           c=c, p=p, q=q, bias=c(bias_l,bias_r), kernel=kernel_type, all=all,
           vce=vce_type, bwselect=bwselect, level=level, masspoints=masspoints)
  out$call <- match.call()
  class(out) <- "rdrobust"
  return(out)
}


print.rdrobust <- function(x,...){
  cat("Call: rdrobust\n\n")
  cat(paste("Number of Obs.           ",  format(x$N[1]+x$N[2], width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")
  cat(paste("Number of Obs.           ",  format(x$N[1],  width=10, justify="right"),  "   ", format(x$N[2],   width=10, justify="right"),       "\n", sep=""))
  cat(paste("Eff. Number of Obs.      ",  format(x$N_h[1],width=10, justify="right"),  "   ", format(x$N_h[2], width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order est. (p)           ",  format(x$p,     width=10, justify="right"),  "   ", format(x$p,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q,     width=10, justify="right"),  "   ", format(x$q,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("BW est. (h)              ",  format(sprintf("%10.3f",x$bws[1,1])),         "   ", format(sprintf("%10.3f",x$bws[1,2])),       "\n", sep=""))
  cat(paste("BW bias (b)              ",  format(sprintf("%10.3f",x$bws[2,1])),         "   ", format(sprintf("%10.3f",x$bws[2,2])),       "\n", sep=""))
  cat(paste("rho (h/b)                ",  format(sprintf("%10.3f",x$bws[1,1]/x$bws[2,1])),  "   ", format(sprintf("%10.3f",x$bws[1,2]/x$bws[2,2])),       "\n", sep=""))
  if (x$masspoints=="adjust" | x$masspoints=="check") cat(paste("Unique Obs.              ",  format(x$M[1], width=10, justify="right"), "   ", format(x$M[2],width=10, justify="right"),        "\n", sep=""))
  cat("\n")
}

summary.rdrobust <- function(object,...) {
  x    <- object
  args <- list(...)
  #if (is.null(args[['level']])) { level <- 0.05 } else { level <- args[['level']] }

  cat("Call: rdrobust\n\n")
  cat(paste("Number of Obs.           ",  format(x$N[1]+x$N[2], width=10, justify="right"),"\n", sep=""))
  cat(paste("BW type                  ",  format(x$bwselect, width=10, justify="right"),"\n", sep=""))
  cat(paste("Kernel                   ",  format(x$kernel,   width=10, justify="right"),"\n", sep=""))
  cat(paste("VCE method               ",  format(x$vce,      width=10, justify="right"),"\n", sep=""))
  cat("\n")
  cat(paste("Number of Obs.           ",  format(x$N[1],   width=10, justify="right"),  "   ", format(x$N[2],   width=10, justify="right"),       "\n", sep=""))
  cat(paste("Eff. Number of Obs.      ",  format(x$N_h[1], width=10, justify="right"),  "   ", format(x$N_h[2], width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order est. (p)           ",  format(x$p,      width=10, justify="right"),  "   ", format(x$p,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("Order bias  (q)          ",  format(x$q,      width=10, justify="right"),  "   ", format(x$q,      width=10, justify="right"),       "\n", sep=""))
  cat(paste("BW est. (h)              ",  format(sprintf("%10.3f",x$bws[1,1])),  "   ", format(sprintf("%10.3f",x$bws[1,2])),      "\n", sep=""))
  cat(paste("BW bias (b)              ",  format(sprintf("%10.3f",x$bws[2,1])),  "   ", format(sprintf("%10.3f",x$bws[2,2])),      "\n", sep=""))
  cat(paste("rho (h/b)                ",  format(sprintf("%10.3f",x$bws[1,1]/x$bws[2,1])),  "   ", format(sprintf("%10.3f",x$bws[1,2]/x$bws[2,2])),       "\n", sep=""))
  if (x$masspoints=="adjust" | x$masspoints=="check") cat(paste("Unique Obs.              ",  format(x$M[1], width=10, justify="right"), "   ", format(x$M[2],width=10, justify="right"),        "\n", sep=""))
  cat("\n")

  ### compute CI
  #z    <- qnorm(100 - level / 2)
  z <- -qnorm(abs((1-(x$level/100))/2))
  
  CI_us_l <- x$Estimate[, "tau.us"] - x$Estimate[, "se.us"] * z;
  CI_us_r <- x$Estimate[, "tau.us"] + x$Estimate[, "se.us"] * z;
  CI_bc_l <- x$Estimate[, "tau.bc"] - x$Estimate[, "se.us"] * z;
  CI_bc_r <- x$Estimate[, "tau.bc"] + x$Estimate[, "se.us"] * z;
  CI_rb_l <- x$Estimate[, "tau.bc"] - x$Estimate[, "se.rb"] * z;
  CI_rb_r <- x$Estimate[, "tau.bc"] + x$Estimate[, "se.rb"] * z;
  
  t_us =x$Estimate[, "tau.us"]/x$Estimate[, "se.us"]
  t_bc =x$Estimate[, "tau.bc"]/x$Estimate[, "se.us"]
  t_rb =x$Estimate[, "tau.bc"]/x$Estimate[, "se.rb"]
  
  pv_us = 2*pnorm(-abs(t_us))
  pv_bc = 2*pnorm(-abs(t_bc))
  pv_rb = 2*pnorm(-abs(t_rb))
  
  ### print output
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
  
  cat(format("Method"          , width=14, justify="right"))
  cat(format("Coef."           , width=10, justify="right"))
  cat(format("Std. Err."       , width=10 , justify="right"))
  cat(format("z"               , width=10, justify="right"))
  cat(format("P>|z|"           , width=10, justify="right"))
  cat(format(paste("[ ", x$level, "%", " C.I. ]", sep=""), width=25, justify="centre"))
  cat("\n")
  
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
  
    cat(format("Conventional", width=14, justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[1, "tau.us"]) , width=10, justify="right"))
    cat(format(paste(sprintf("%3.3f", x$Estimate[1, "se.us"]), sep=""), width=10, justify="right"))
    cat(format(sprintf("%3.3f", t_us) , width=10, justify="right"))
    cat(format(sprintf("%3.3f", pv_us), width=10, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", CI_us_l[1]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_us_r[1]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    
    if (is.null(x$all)) {
      cat(format("Robust", width=14, justify="right"))
      cat(format("-", width=10, justify="right"))
      cat(format("-", width=10, justify="right"))
      #cat(format(sprintf("%3.3f", x$Estimate[1, "tau.bc"]) , width=10, justify="right"))
      #cat(format(paste(sprintf("%3.3f", x$Estimate[1, "se.rb"]), sep=""), width=10, justify="right"))
      cat(format(sprintf("%3.3f", t_rb) , width=10, justify="right"))
      cat(format(sprintf("%3.3f", pv_rb), width=10, justify="right"))
      cat(format(paste("[", sprintf("%3.3f", CI_rb_l[1]), " , ", sep="")  , width=14, justify="right"))
      cat(format(paste(sprintf("%3.3f", CI_rb_r[1]), "]", sep=""), width=11, justify="left"))
      cat("\n") 
    } else {
      cat(format("Bias-Corrected", width=14, justify="right"))
      cat(format(sprintf("%3.3f", x$Estimate[1, "tau.bc"]) , width=10, justify="right"))
      cat(format(paste(sprintf("%3.3f", x$Estimate[1, "se.us"]), sep=""), width=10, justify="right"))
      cat(format(sprintf("%3.3f", t_bc) , width=10, justify="right"))
      cat(format(sprintf("%3.3f", pv_bc), width=10, justify="right"))
      cat(format(paste("[", sprintf("%3.3f", CI_bc_l[1]), " , ", sep="")  , width=14, justify="right"))
      cat(format(paste(sprintf("%3.3f", CI_bc_r[1]), "]", sep=""), width=11, justify="left"))
      cat("\n")
  
      cat(format("Robust", width=14, justify="right"))
      cat(format(sprintf("%3.3f", x$Estimate[1, "tau.bc"]) , width=10, justify="right"))
      cat(format(paste(sprintf("%3.3f", x$Estimate[1, "se.rb"]), sep=""), width=10, justify="right"))
      cat(format(sprintf("%3.3f", t_rb) , width=10, justify="right"))
      cat(format(sprintf("%3.3f", pv_rb), width=10, justify="right"))
      cat(format(paste("[", sprintf("%3.3f", CI_rb_l[1]), " , ", sep="")  , width=14, justify="right"))
      cat(format(paste(sprintf("%3.3f", CI_rb_r[1]), "]", sep=""), width=11, justify="left"))
      cat("\n")
    }
    
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
}

