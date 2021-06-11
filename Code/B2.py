# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:17:53 2020

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
Xba = Symbol('X_ba')
Xaba = Symbol('X_aba')
Xbab = Symbol('Xbab')
Xabab = Symbol('X_abab')
#aij are the coefficients in the formal group law
a11 = Symbol('a_11')
a12 = Symbol('a_12')
a13 = Symbol('a_13')
a22 = Symbol('a_22')
a14 = Symbol('a_14')
a23 = Symbol('a_23')
a15 = Symbol('a_15')
a24 = Symbol('a_24')
a33 = Symbol('a_33')
#ci are the coefficient in the formal inverse
c2 = -a11
c3 = -a11*c2
c4 = -a11*c3+a12*c2-2*a13+a22
#These are the formal inverses of xa and xb truncated to degree 4
x_nega = poly(-xa-c2*(xa**2)-c3*(xa**3)-c4*(xa**4),xa,xb)
x_negb = poly(-xb-c2*(xb**2)-c3*(xb**3)-c4*(xb**4),xa,xb)
#This is x_{b+b} expaneded as a series xb+_Fx_b truncated to degree 4
xaa = poly(2*xa+a11*(xa**2)+2*a12*(xa**3)+2*a13*(xa**4)+a22*(xa**4),xa,xb)
#This is x_{a+b} expaneded as a series xa+_Fx_b truncated to degree 4
xab = poly(xa+xb+a11*xa*xb+a12*((xa**2)*xb+xa*(xb**2))+a13*((xa**3)*xb+xa*(xb**3))+a22*(xa**2)*(xb**2),xa,xb)
#This is x_{a+2b} expaneded as a series xa+_Fx_b+_Fxb, correct up to degree 4
x2ab = poly(xaa+xb+a11*xaa*xb+a12*((xaa**2)*xb+xaa*(xb**2))+a13*((xaa**3)*xb+xaa*(xb**3))+a22*(xaa**2)*(xb**2),xa,xb)
#This is sa(xb), sb(xa), sasb(xa), sbsa(xb) expanded as a  polynomial of xa and xb up to degree 4
sbxa = poly(xab,xa,xb)
saxb = poly(x2ab,xa,xb)
sasbxa = poly(xab,xa,xb)
sbsaxb = poly(x2ab,xa,xb)
#The constant terms of ka and kb in Xa^2=kaXa and Xb^2=kbXb
ka = -a11
kb = -a11

def divide_xa(f):
   g = poly(0,xa,xb)
   for i in [1,2,3,4]:
       for j in [0,1,2,3,4]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**(i-1))*(xb**j),xa,xb)
   return g      
           
def divide_xb(f):
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4]:
       for j in [1,2,3,4]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**i)*(xb**(j-1)),xa,xb)
   return g  

Da_xa = divide_xa(poly(xa-x_nega,xa,xb))
Da_xb = divide_xa(poly(xb-saxb,xa,xb))
Da_xaxa = Da_xa*xa+x_nega*Da_xa
Da_xaxb = Da_xa*xb+x_nega*Da_xb
Da_xbxb = Da_xb*xb+saxb*Da_xb
Da_xaxaxa = Da_xa*xa*xa+x_nega*Da_xaxa
Da_xaxaxb = Da_xa*xa*xb+x_nega*Da_xaxb
Da_xaxbxb = Da_xa*xb*xb+x_nega*Da_xbxb
Da_xbxbxb = Da_xb*xb*xb+saxb*Da_xbxb
Da_xaxaxaxa = Da_xa*xa*xa*xa+x_nega*Da_xaxaxa
Da_xaxaxaxb = Da_xa*xa*xa*xb+x_nega*Da_xaxaxb
Da_xaxaxbxb = Da_xa*xa*xb*xb+x_nega*Da_xaxbxb
Da_xaxbxbxb = Da_xa*xb*xb*xb+x_nega*Da_xbxbxb
Da_xbxbxbxb = Da_xb*xb*xb*xb+saxb*Da_xbxbxb
Da = [0,Da_xb,Da_xbxb,Da_xbxbxb,Da_xbxbxbxb,Da_xa,Da_xaxb,Da_xaxbxb,Da_xaxbxbxb,Da_xaxa,Da_xaxaxb,Da_xaxaxbxb,Da_xaxaxa,Da_xaxaxaxb,Da_xaxaxaxa]

Db_xa = divide_xb(poly(xa-sbxa,xa,xb))
Db_xb = divide_xb(poly(xb-x_negb,xa,xb))
Db_xaxa = Db_xa*xa+sbxa*Db_xa
Db_xaxb = Db_xb*xa+x_negb*Db_xa
Db_xbxb = Db_xb*xb+x_negb*Db_xb
Db_xaxaxa = Db_xa*xa*xa+sbxa*Db_xaxa
Db_xaxaxb = Db_xb*xa*xa+x_negb*Db_xaxa
Db_xaxbxb = Db_xb*xa*xb+x_negb*Db_xaxb
Db_xbxbxb = Db_xb*xb*xb+x_negb*Db_xbxb
Db_xaxaxaxa = Db_xa*xa*xa*xa+sbxa*Db_xaxaxa
Db_xaxaxaxb = Db_xb*xa*xa*xa+x_negb*Db_xaxaxa
Db_xaxaxbxb = Db_xb*xa*xa*xb+x_negb*Db_xaxaxb
Db_xaxbxbxb = Db_xb*xa*xb*xb+x_negb*Db_xaxbxb
Db_xbxbxbxb = Db_xb*xb*xb*xb+x_negb*Db_xbxbxb
Db = [0,Db_xb,Db_xbxb,Db_xbxbxb,Db_xbxbxbxb,Db_xa,Db_xaxb,Db_xaxbxb,Db_xaxbxbxb,Db_xaxa,Db_xaxaxb,Db_xaxaxbxb,Db_xaxaxa,Db_xaxaxaxb,Db_xaxaxaxa]

def compose_Da(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4]:
       for j in [0,1,2,3,4]:
           if ((i+j)<=4):
               g = poly(g + Da[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

def compose_Db(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4]:
       for j in [0,1,2,3,4]:
           if ((i+j)<=4):
               g = poly(g + Db[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

DbDa_xa = compose_Db(Da_xa)
DaDb_xa = compose_Da(Db_xa)
DbDaDb_xa = compose_Db(DaDb_xa)
DaDbDa_xa = compose_Da(DbDa_xa)
DaDbDaDb_xa = compose_Da(DbDaDb_xa)

DbDa_xb = compose_Db(Da_xb)
DaDb_xb = compose_Da(Db_xb)
DbDaDb_xb = compose_Db(DaDb_xb)
DaDbDa_xb = compose_Da(DbDa_xb)
DaDbDaDb_xb = compose_Da(DbDaDb_xb)

da_xa = Da_xa.coeff_monomial(1)
db_xa = Db_xa.coeff_monomial(1)
dadb_xa = DaDb_xa.coeff_monomial(1)
dbda_xa = DbDa_xa.coeff_monomial(1)
dadbda_xa = DaDbDa_xa.coeff_monomial(1)
dbdadb_xa = DbDaDb_xa.coeff_monomial(1)
dadbdadb_xa = DaDbDaDb_xa.coeff_monomial(1)

da_xb = Da_xb.coeff_monomial(1)
db_xb = Db_xb.coeff_monomial(1)
dadb_xb = DaDb_xb.coeff_monomial(1)
dbda_xb = DbDa_xb.coeff_monomial(1)
dadbda_xb = DaDbDa_xb.coeff_monomial(1)
dbdadb_xb = DbDaDb_xb.coeff_monomial(1)
dadbdadb_xb = DaDbDaDb_xb.coeff_monomial(1)

ev_xa = da_xa*Xa+db_xa*Xb+dadb_xa*Xab+dbda_xa*Xba+dadbda_xa*Xaba+dbdadb_xa*Xbab+dadbdadb_xa*Xabab
ev_xb = da_xb*Xa+db_xb*Xb+dadb_xb*Xab+dbda_xb*Xba+dadbda_xb*Xaba+dbdadb_xb*Xbab+dadbdadb_xb*Xabab

print("ev_xa = ", ev_xa)
print("ev_xb =", ev_xb)


Da_sbxa = compose_Da(sbxa)
Db_saxb = compose_Db(saxb)
Da_sbsaxb = compose_Da(sbsaxb)
DaDb_saxb = compose_Da(compose_Db(saxb))
DaDb_xasaxb = compose_Da(compose_Db(xa*saxb))
db_saxb = Db_saxb.coeff_monomial(1)
da_sbxa = Da_sbxa.coeff_monomial(1)
da_sbsaxb = Da_sbsaxb.coeff_monomial(1)
dadb_saxb = DaDb_saxb.coeff_monomial(1)
dadb_xasaxb = DaDb_xasaxb.coeff_monomial(1)

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
#XbXaXb
q_a_a_bab = 0
q_a_b_bab = kb
q_a_ab_bab = -db_xa
q_a_ba_bab = 0
q_b_b_bab = -dbda_xb-kb*da_xb-kb*da_xb
q_b_ab_bab = 1-db_saxb
q_b_ba_bab = 1
q_ab_ab_bab = 0
q_ab_ba_bab = 0
q_ba_ba_bab = 0
#XaXbXaXb
q_a_a_abab = Symbol('q_a_a_abab')
q_a_b_abab = Symbol('q_a_b_abab')
q_a_ab_abab = Symbol('q_a_ab_abab')
q_a_ba_abab = Symbol('q_a_ba_abab')
q_a_aba_abab = 0
q_a_bab_abab = 1-da_sbxa
q_b_b_abab = Symbol('q_b_b_abab')
q_b_ab_abab = Symbol('q_b_ab_abab')
q_b_ba_abab = Symbol('q_b_ba_abab')
q_b_aba_abab = 1
q_b_bab_abab = -da_xb-da_sbsaxb
q_ab_ab_abab = Symbol('q_ab_ab_abab')
q_ab_ba_abab = Symbol('q_ab_ba_abab')
q_ab_aba_abab = 0
q_ab_bab_abab = 0
q_ba_ba_abab = Symbol('q_ba_ba_abab')
q_ba_aba_abab = 0
q_ba_bab_abab = 0
q_aba_aba_abab = 0
q_aba_bab_abab = 0
q_bab_bab_abab = 0

#Coefficients of relation M
print("")
print("Second Relations (Products):")
print("")
print("XaXa=", q_a_a_ab*Xab+q_a_a_ba*Xba+q_a_a_aba*Xaba+q_a_a_bab*Xbab+q_a_a_abab*Xabab)
print("")
print("XaXb=", q_a_b_ab*Xab+q_a_b_ba*Xba+q_a_b_aba*Xaba+q_a_b_bab*Xbab+q_a_b_abab*Xabab)
print("")
print("XaXab=", q_a_ab_aba*Xaba+q_a_ab_bab*Xbab+q_a_ab_abab*Xabab)
print("")
print("XaXba=", q_a_ba_aba*Xaba+q_a_ba_bab*Xbab+q_a_ba_abab*Xabab)
print("")
print("XaXaba=", q_a_aba_abab*Xabab)
print("")
print("XaXbab=", q_a_bab_abab*Xabab)
print("")
print("XaXabab=", 0)
print("")
print("XbXb=", q_b_b_ab*Xab+q_b_b_ba*Xba+q_b_b_aba*Xaba+q_b_b_bab*Xbab+q_b_b_abab*Xabab)
print("")
print("XbXab=", q_b_ab_aba*Xaba+q_b_ab_bab*Xbab+q_b_ab_abab*Xabab)
print("")
print("XbXba=", q_b_ba_aba*Xaba+q_b_ba_bab*Xbab+q_b_ba_abab*Xabab)
print("")
print("XbXaba=", q_b_aba_abab*Xabab)
print("")
print("XbXbab=", q_b_bab_abab*Xabab)
print("")
print("XbXabab=", 0)
print("")
print("XabXab=", q_ab_ab_aba*Xaba+q_ab_ab_bab*Xbab+q_ab_ab_abab*Xabab)
print("")
print("XabXba=", q_ab_ba_aba*Xaba+q_ab_ba_bab*Xbab+q_ab_ba_abab*Xabab)
print("")
print("XabXaba=", q_ab_aba_abab*Xabab)
print("")
print("XabXbab=", q_ab_bab_abab*Xabab)
print("")
print("XabXabab=", 0)
print("")
print("XbaXba=", q_ba_ba_aba*Xaba+q_ba_ba_bab*Xbab+q_ba_ba_abab*Xabab)
print("")
print("XbaXaba=", q_ba_aba_abab*Xabab)
print("")
print("XbaXbab=", q_ba_bab_abab*Xabab)
print("")
print("XbaXabab=", 0)
print("")
print("XabaXaba=", q_aba_aba_abab*Xabab)
print("")
print("XabaXbab=", q_aba_bab_abab*Xabab)
print("")
print("XabaXabab=", 0)
print("")
print("XbabXbab=", q_bab_bab_abab*Xabab)
print("")
print("XbabXabab=", 0)
print("")
print("XababXabab=", 0)


