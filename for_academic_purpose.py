7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:34:40 2019

@author: tug81087
"""
import numpy as np
a,b = np.loadtxt('test_data.txt',skiprows=1,max_rows=10,unpack=True)
print('a= ',a)
print('b= ',b)