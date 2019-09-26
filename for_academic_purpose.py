7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:34:40 2019

@author: tug81087
"""
import numpy as np
a = np.arange(4)
a = a[:, np.newaxis]
#b = np.tile(a,(4,1))
#print(b.transpose())
print(np.tile(a,(1,4)))