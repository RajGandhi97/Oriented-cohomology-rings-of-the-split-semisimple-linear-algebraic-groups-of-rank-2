# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 18:40:22 2021

@author: Raj
"""


from sympy import Symbol, poly, div, degree_list, simplify
from sympy.abc import x, y


#The xa and xb are the canonical generators of the formal group algebra
xa = Symbol('x_a')
xb = Symbol('x_b')
#Xa and Xb are the dual basis elements in the dual of the formal Demazure algebra
Xa = Symbol('X_a')
Xb = Symbol('X_b')
Xab = Symbol('X_ab')
Xba = Symbol('X_ba')
Xaba = Symbol('X_aba')
#aij are the coefficients in the formal group law
a11 = Symbol('a_11')
a12 = Symbol('a_12')
a13 = Symbol('a_13')
a22 = Symbol('a_22')
a14 = Symbol('a_14')
a23 = Symbol('a_23')
#ci are the coefficient in the formal inverse
c2 = -a11
c3 = -a11*c2
#These are the formal inverses of xa and xb truncated to degree 3
x_nega = poly(-xa-c2*(xa**2)-c3*(xa**3),xa,xb)
x_negb = poly(-xb-c2*(xb**2)-c3*(xb**3),xa,xb)
#This is x_{a+b} expaneded as a series xa+_Fx_b truncated to degree 3
xab = poly(xa+xb+a11*xa*xb+a12*((xa**2)*xb+xa*(xb**2)),xa,xb)
#This is sa(xb) and sb(xa) expanded as a  polynomial of xa and xb up to degree 3
saxb = xab
sbxa = xab
#These are the constant terms of ka and kb in Xa^2=ka*X and Xb^2=kb*Xb
ka = -a11
kb = -a11

def divide_xa(f):
   g = poly(0,xa,xb)
   for i in [1,2,3]:
       for j in [0,1,2,3]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**(i-1))*(xb**j),xa,xb)
   return g      
           
def divide_xb(f):
   g = poly(0,xa,xb)
   for i in [0,1,2,3]:
       for j in [1,2,3]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**i)*(xb**(j-1)),xa,xb)
   return g  


Da_xa = divide_xa(poly(xa-x_nega,xa,xb))
Da_xb = divide_xa(poly(xb-xab,xa,xb))
Da_xaxa = Da_xa*xa+x_nega*Da_xa
Da_xaxb = Da_xa*xb+x_nega*Da_xb
Da_xbxb = Da_xb*xb+xab*Da_xb
Da_xaxaxa = Da_xa*xa*xa+x_nega*Da_xaxa
Da_xaxaxb = Da_xa*xa*xb+x_nega*Da_xaxb
Da_xaxbxb = Da_xa*xb*xb+x_nega*Da_xbxb
Da_xbxbxb = Da_xb*xb*xb+xab*Da_xbxb
Da = [0,Da_xb,Da_xbxb,Da_xbxbxb,Da_xa,Da_xaxb,Da_xaxbxb,Da_xaxa,Da_xaxaxb,Da_xaxaxa]

Db_xa = divide_xb(poly(xa-xab,xa,xb))
Db_xb = divide_xb(poly(xb-x_negb,xa,xb))
Db_xaxa = Db_xa*xa+xab*Db_xa
Db_xaxb = Db_xb*xa+x_negb*Db_xa
Db_xbxb = Db_xb*xb+x_negb*Db_xb
Db_xaxaxa = Db_xa*xa*xa+xab*Db_xaxa
Db_xaxaxb = Db_xb*xa*xa+x_negb*Db_xaxa
Db_xaxbxb = Db_xb*xa*xb+x_negb*Db_xaxb
Db_xbxbxb = Db_xb*xb*xb+x_negb*Db_xbxb
Db = [0,Db_xb,Db_xbxb,Db_xbxbxb,Db_xa,Db_xaxb,Db_xaxbxb,Db_xaxa,Db_xaxaxb,Db_xaxaxa]


def compose_Da(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3]:
       for j in [0,1,2,3]:
           if ((i+j)<=3):
               g = poly(g + Da[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

def compose_Db(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3]:
       for j in [0,1,2,3]:
           if ((i+j)<=3):
               g = poly(g + Db[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g
           
DbDa_xa = compose_Db(Da_xa)
DaDb_xa = compose_Da(Db_xa)
DaDbDa_xa = compose_Da(DbDa_xa)

DbDa_xb = compose_Db(Da_xb)
DaDb_xb = compose_Da(Db_xb)
DaDbDa_xb = compose_Da(DbDa_xb)

da_xa = Da_xa.coeff_monomial(1)
db_xa = Db_xa.coeff_monomial(1)
dadb_xa = DaDb_xa.coeff_monomial(1)
dbda_xa = DbDa_xa.coeff_monomial(1)
dadbda_xa = DaDbDa_xa.coeff_monomial(1)

da_xb = Da_xb.coeff_monomial(1)
db_xb = Db_xb.coeff_monomial(1)
dadb_xb = DaDb_xb.coeff_monomial(1)
dbda_xb = DbDa_xb.coeff_monomial(1)
dadbda_xb = DaDbDa_xb.coeff_monomial(1)


ev_xa = da_xa*Xa+db_xa*Xb+dadb_xa*Xab+dbda_xa*Xba+dadbda_xa*Xaba
ev_xb = da_xb*Xa+db_xb*Xb+dadb_xb*Xab+dbda_xb*Xba+dadbda_xb*Xaba

print("ev_xa = ", ev_xa)
print("ev_xb =", ev_xb)



Da_sbxa = compose_Da(sbxa)
da_sbxa = Da_sbxa.coeff_monomial(1)


#Coefficients q_{w,w'}^v in the coproduct
#Xa
#Xb
#XaXb
q_a_a_ab = 0
q_a_b_ab = 1
q_b_b_ab = -da_xb
#XbXa 
q_a_a_ba = -db_xa
q_a_b_ba = 1
q_b_b_ba = 0
#XaXbXa
q_a_a_aba = -dadb_xa-ka*db_xa-ka*db_xa
q_a_b_aba = ka
q_a_ab_aba = 1
q_a_ba_aba = 1 - da_sbxa
q_b_b_aba = 0
q_b_ab_aba = 0
q_b_ba_aba = -da_xb
q_ab_ab_aba = 0
q_ab_ba_aba = 0
q_ba_ba_aba = 0

#Coefficients of relation M
print("")
print("Second Relations (Products):")
print("")
print("XaXa=", q_a_a_ab*Xab+q_a_a_ba*Xba+q_a_a_aba*Xaba)
print("")
print("XaXb=", q_a_b_ab*Xab+q_a_b_ba*Xba+q_a_b_aba*Xaba)
print("")
print("XaXab=", q_a_ab_aba*Xaba)
print("")
print("XaXba=", q_a_ba_aba*Xaba)
print("")
print("XaXaba=", 0)
print("")
print("XbXb=", q_b_b_ab*Xab+q_b_b_ba*Xba+q_b_b_aba*Xaba)
print("")
print("XbXab=", q_b_ab_aba*Xaba)
print("")
print("XbXba=", q_b_ba_aba*Xaba)
print("")
print("XbXaba=", 0)
print("")
print("XabXab=", q_ab_ab_aba*Xaba)
print("")
print("XabXba=", q_ab_ba_aba*Xaba)
print("")
print("XabXaba=", 0)
print("")
print("XbaXba=", q_ba_ba_aba*Xaba)
print("")
print("XbaXaba=", 0)
print("")
print("XabaXaba=", 0)

