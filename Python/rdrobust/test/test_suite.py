#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 21:56:03 2021

@author: rmasini
"""

from rdrobust import rdbwselect, rdrobust, rdplot
import pandas  as pd
import os
#=========================================================================
# Test Suite


# ISSUES/QUESTIONS
# Multicollinearity issue. Solved with RRQR
# cov_drop = TRUE always in rdbwselect call? Fix in Python
# Strange behaviour with the subset option (Discrepancy in the third decimal point ?!?!)
# The R version does not report in the final table Coef. and S.E. of Robust estimate. Why? 
# Test rdplot with covariates    


# Load Dataset

# # Create the 'RDSenate_aug.csv' for testing the package features
# df = pd.read_csv('RDSenate.csv')
# n = len(df.margin)
# df['fuzzy'] = ((df.margin + np.random.normal(size = n , scale = 2))>=0).astype(int)
# df['cov_abc'] =  np.random.normal(size = n , scale = 2)
# df['some_cov'] =  np.random.normal(size = n , scale = 2)
# df['otherZ'] =  np.random.normal(size = n , scale = 2)
# df['weights'] = np.random.uniform(size = n)
# df['subset'] = np.random.choice([0,1], n)
# df['cluster'] = np.random.choice(np.arange(5), n)
# df.to_csv('RDSenate_aug.csv', index = False)    
    
df = pd.read_csv('RDSenate_aug.csv')

x = df.margin
y = df.vote
fuzzy = df.fuzzy
covs = df.iloc[:,3:6]
weights = df.weights
subset = df.subset.astype(bool)
cluster = df.cluster

##############################################################################
##### Test rdbwselect Function ###############################################
##############################################################################
y = y
x = x
c = None
fuzzy = None
deriv = None
p = None
q = None
covs = None
covs_drop = True
kernel = "tri"
weights = None
bwselect = "mserd"
vce = "nn"
cluster = None
nnmatch = 3
scaleregul = 1
sharpbw = False
all = True
subset = None
masspoints = "adjust"
bwcheck = None
bwrestrict = True
stdvars = False
prchk = True

out1 = rdbwselect(y, x, c, fuzzy, deriv, p, q, covs, covs_drop, kernel,
                  weights, bwselect, vce, cluster, nnmatch, scaleregul,
                  sharpbw, all, subset, masspoints, bwcheck, bwrestrict,
                  stdvars, prchk)

print(out1)

##############################################################################
##### Test Rdrobust Function #################################################
##############################################################################
y = y
x = x
c = None
fuzzy = None
deriv = None
p = None
q = None 
h = None 
b = None 
rho = None
covs = None
covs_drop = True
kernel = "tri"          # tri, epa, uni
weights = None
bwselect = "mserd"      # mserd, msetwo, msesum, msecomb1, msecomb2, cerrd, certwo, cersum, cercomb1, cercomb2
vce = "nn"              # nn, hc0, hc1, hc2, hc3
cluster = None
nnmatch = 3
level = 95
scalepar = 1
scaleregul = 1
sharpbw = False
all = None
subset = None
masspoints = "adjust"   # adjust, check, off
bwcheck = None
bwrestrict = True
stdvars = False        

out2 =  rdrobust(y, x, c, fuzzy, deriv, p, q, h, b, rho, covs, covs_drop,
                kernel, weights, bwselect, vce, cluster, nnmatch, level,
                scalepar, scaleregul, sharpbw, all, subset, masspoints,
                bwcheck, bwrestrict, stdvars)
print(out2)

##############################################################################
##### Test Rplot Function ####################################################
##############################################################################
y = y 
x = x
c = 0
p = 4
nbins = None
binselect = "esmv"
scale = None 
kernel = "uni"
weights = None
h = None
covs = covs
covs_eval = "mean"
covs_drop = True
support = None
subset = None
masspoints = "adjust"
hide = False
ci = None
shade = False
title = None
x_label = None
y_label = None
x_lim = None
y_lim = None
col_dots = None
col_lines = None

out3 =  rdplot(y, x, c, p, nbins, binselect, scale, kernel, weights, h, 
              covs,  covs_eval, covs_drop, support, subset, masspoints,
              hide, ci, shade, title, x_label, y_label, x_lim, y_lim,
              col_dots, col_lines)

print(out3)