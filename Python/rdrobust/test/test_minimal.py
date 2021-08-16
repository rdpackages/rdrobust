#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 23:57:26 2021

@author: rmasini
"""


import numpy as np
from rdrobust import rdrobust, rdbwselect,rdplot

x = np.random.uniform(size = 100)
y = np.random.uniform(size = 100)
print(rdbwselect(x,y, c =0.5))
print(rdrobust(x,y, c =0.5))
print(rdplot(x,y, c =0.5))
