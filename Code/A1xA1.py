# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 17:26:33 2021

@author: Raj
"""


from sympy import Symbol, poly, div, degree_list
from sympy.abc import x, y

#The xa and xb are the canonical generators of the formal group algebra
xa = Symbol('x_a')
xb = Symbol('x_b')
#Xa and Xb are the dual basis elements in the dual of the formal Demazure algebra
Xa = Symbol('X_a')
Xb = Symbol('X_b')
Xab = Symbol('X_ab')
#a11 is a coefficient in the formal group law
a11 = Symbol('a_11')
#c2 in the coefficient in the formal inverse
c2 = -a11
#These are the formal inverses of xa and xb truncated to degree 2
x_nega = poly(-xa-c2*(xa**2),xa,xb)
x_negb = poly(-xb-c2*(xb**2),xa,xb)

def divide_xa(f):
   g = poly(0,xa,xb)
   for i in [1,2]:
       for j in [0,1,2]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**(i-1))*(xb**j),xa,xb)
   return g      
           
def divide_xb(f):
   g = poly(0,xa,xb)
   for i in [0,1,2]:
       for j in [1,2]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**i)*(xb**(j-1)),xa,xb)
   return g  


Da_xa = divide_xa(poly(xa-x_nega,xa,xb))
Da_xb = divide_xa(poly(0,xa,xb))
Da_xaxa = Da_xa*xa+x_nega*Da_xa
Da_xaxb = Da_xa*xb+x_nega*Da_xb
Da_xbxb = Da_xb*xb+xb*Da_xb
Da = [0,Da_xb,Da_xbxb,Da_xa,Da_xaxb,Da_xaxa]

Db_xa = divide_xb(poly(0,xa,xb))
Db_xb = divide_xb(poly(xb-x_negb,xa,xb))
Db_xaxa = Db_xa*xa+xa*Db_xa
Db_xaxb = Db_xb*xa+x_negb*Db_xa
Db_xbxb = Db_xb*xb+x_negb*Db_xb
Db = [0,Db_xb,Db_xbxb,Db_xa,Db_xaxb,Db_xaxa]


def compose_Da(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2]:
       for j in [0,1,2]:
           if ((i+j)<=2):
               g = poly(g + Da[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

def compose_Db(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2]:
       for j in [0,1,2]:
           if ((i+j)<=2):
               g = poly(g + Db[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g
           
DaDb_xa = compose_Da(Db_xa)

DaDb_xb = compose_Da(Db_xb)

da_xa = Da_xa.coeff_monomial(1)
db_xa = Db_xa.coeff_monomial(1)
dadb_xa = DaDb_xa.coeff_monomial(1)


da_xb = Da_xb.coeff_monomial(1)
db_xb = Db_xb.coeff_monomial(1)
dadb_xb = DaDb_xb.coeff_monomial(1)


ev_xa = da_xa*Xa+db_xa*Xb+dadb_xa*Xab
ev_xb = da_xb*Xa+db_xb*Xb+dadb_xb*Xab

print("ev_xa = ", ev_xa)
print("ev_xb =", ev_xb)

#Coefficients q_{w,w'}^v in the coproduct

#Xa
#Xb
#XaXb
q_a_a_ab = 0
q_a_b_ab = 1
q_b_b_ab = -da_xb


#Coefficients of relation M
print("")
print("Second Relations (Products):")
print("")
print("XaXa=", q_a_a_ab*Xab)
print("")
print("XaXb=", q_a_b_ab*Xab)
print("")
print("XaXab=", 0)
print("")
print("XaXba=", 0)
print("")
print("XbXb=", q_b_b_ab*Xab)
print("")
print("XbXab=", 0)
print("")
print("XabXab=", 0)

    
