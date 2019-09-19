#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:22:25 2019

@author: Tanvir Chowdhury
"""
from numpy import sqrt, arctan, log

def beta_model(u,w_p_bar,wg,y,delta,P):
    Iplus=sqrt((1+y ** 2) * (1+u ** 2/wg ** 2)) + u*y/wg
    Iminus=sqrt((1+y ** 2)*(1+u ** 2/wg ** 2)) - u*y/wg
    part1=(wg ** 2-(wg ** 2 + u ** 2) * delta ** 2)/(2*u*sqrt(wg ** 2+u ** 2))
    part2=(2*w_p_bar**2*delta/u**2)*((wg/u)*(arctan(wg*P/u) - arctan(wg/u))+1/P-1)
    e1=1+(w_p_bar ** 2/u ** 2)*((1-delta ** 2)*y/P-log(Iplus/Iminus)*part1)+part2
    return (e1-1)/(e1+1)

def polarize(u,w1,w2,a_0_1,a_0_2):
    a_u_1 = a_0_1 * w1 **2/(w1 ** 2 + u **2)
    a_u_2 = a_0_2 * w2 ** 2/(w2 ** 2 + u **2)
    return a_u_1, a_u_2
   
def diffuse(u,g0,w_sub,w1,d0_perp):
    lamb = sqrt((u**2+ w_sub ** 2*(1-g0))/(2 * w_sub ** 2*g0))
    du_perp = d0_perp*(1- lamb * arctan(1/lamb));
    fw = (u ** 2+(1+g0)* w_sub ** 2)/((u ** 2+w1 ** 2)*(u ** 2+ w_sub ** 2) ** 2)
    return lamb, du_perp, fw