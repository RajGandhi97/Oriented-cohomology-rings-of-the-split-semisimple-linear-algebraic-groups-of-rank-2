# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 13:43:39 2021

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
Xbaba = Symbol('X_baba')
Xababa = Symbol('X_ababa')
Xbabab = Symbol('X_babab')
Xababab = Symbol('X_ababab')
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
c5 = -a11*c4+a12*c3+a12*c2**2-4*a13*c2+2*a22*c2
c6 = -a11*c5+a12*c4+2*a12*c2*c3-4*a13*c3-3*a13*c2**2+3*a14*c2-2*a15+2*a22*c3+a22*c2**2-a23*c2+2*a24-a33

#The truncate() function takes as input a polynomial in xa and xb and truncates 
#it to degree 6
def truncate(f):
    g = poly(0,xa,xb)
    for i in [0,1,2,3,4,5,6]:
        for j in [0,1,2,3,4,5,6]:
            if ((i+j) <= 6):
                C = f.coeff_monomial((xa**i)*(xb**j))
                g = poly(g+C*(xa**i)*(xb**j),xa,xb)
    return g
         

#These elements of the formal group algebra truncated to degree 6
x_nega = poly(-xa-c2*(xa**2)-c3*(xa**3)-c4*(xa**4)-c5*(xa**5)-c6*(xa**6),xa,xb)
x_negb = poly(-xb-c2*(xb**2)-c3*(xb**3)-c4*(xb**4)-c5*(xb**5)-c6*(xb**6),xa,xb)            
xbb = poly(2*xb+a11*(xb**2)+2*a12*(xb**3)+2*a13*(xb**4)+a22*(xb**4)+2*a14*(xb**5)+2*a23*(xb**5)+2*a15*(xb**6)+2*a24*(xb**6)+a33*(xb**6),xa,xb)
xaa = poly(2*xa+a11*(xa**2)+2*a12*(xa**3)+2*a13*(xa**4)+a22*(xa**4)+2*a14*(xa**5)+2*a23*(xa**5)+2*a15*(xa**6)+2*a24*(xa**6)+a33*(xa**6),xa,xb)
xab = poly(xa+xb+a11*xa*xb+a12*((xa**2)*xb+xa*(xb**2))+a13*((xa**3)*xb+xa*(xb**3))+a22*(xa**2)*(xb**2)+a14*((xa**4)*xb+xa*(xb**4))+a23*((xa**2)*(xb**3)+(xa**3)*(xb**2))+a15*((xa**5)*xb+xa*(xb**5))+a24*((xa**2)*(xb**4)+(xa**4)*(xb**2))+a33*(xa**3)*(xb**3),xa,xb)
xa3b_long = poly(xbb+xab+a11*xbb*xab+a12*((xbb**2)*xab+xbb*(xab**2))+a13*((xbb**3)*xab+xbb*(xab**3))+a22*(xbb**2)*(xab**2)+a14*((xab)*(xbb**4)+(xab**4)*xbb)+a23*((xab**2)*(xbb**3)+(xab**3)*(xbb**2))+a15*((xbb**5)*xab+xbb*(xab**5))+a24*((xbb**4)*(xab**2)+(xbb**2)*(xab**4))+a33*(xbb**3)*(xab**3),xa,xb)
xa3b = truncate(xa3b_long)
xa2b_long = poly(xbb+xa+a11*xbb*xa+a12*((xbb**2)*xa+xbb*(xa**2))+a13*((xbb**3)*xa+xbb*(xa**3))+a22*(xbb**2)*(xa**2)+a14*((xa)*(xbb**4)+(xa**4)*xbb)+a23*((xa**2)*(xbb**3)+(xa**3)*(xbb**2))+a15*((xbb**5)*xa+xbb*(xa**5))+a24*((xbb**4)*(xa**2)+(xbb**2)*(xa**4))+a33*(xbb**3)*(xa**3),xa,xb)
xa2b = truncate(xa2b_long)
x2a3b_long =  poly(xa+xa3b+a11*xa*xa3b+a12*((xa**2)*xa3b+xa*(xa3b**2))+a13*((xa**3)*xa3b+xa*(xa3b**3))+a22*(xa**2)*(xa3b**2)+a14*(xa*(xa3b**4)+(xa**4)*xa3b)+a23*((xa**2)*(xa3b**3)+(xa**3)*(xa3b**2))+a15*((xa**5)*xa3b+xa*(xa3b**5))+a24*((xa**4)*(xa3b**2)+(xa**2)*(xa3b**4))+a33*(xa**3)*(xa3b**3),xa,xb)
x2a3b = truncate(x2a3b_long)
#These are the reflections applied to the generators truncated to degree 6 
saxb = xab
sbsaxb = xa2b
sasbsaxb = xa2b
sbsasbsaxb = xab
sasbsasbsaxb = xb
sbxa = xa3b
sasbxa = x2a3b
sbsasbxa = x2a3b
sasbsasbxa = xa3b
sbsasbsasb = xa
#The ka and kb in Xa^2=kaXa and Xb^2=kbXb
ka = -a11
kb = -a11

def divide_xa(f):
   g = poly(0,xa,xb)
   for i in [1,2,3,4,5,6]:
       for j in [0,1,2,3,4,5,6]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**(i-1))*(xb**j),xa,xb)
   return g      
           
def divide_xb(f):
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4,5,6]:
       for j in [1,2,3,4,5,6]:
           C = f.coeff_monomial((xa**i)*(xb**j))
           g = poly(g+C*(xa**i)*(xb**(j-1)),xa,xb)
   return g  

Da_xa = truncate(divide_xa(poly(xa-x_nega,xa,xb)))
Da_xb = truncate(divide_xa(poly(xb-saxb,xa,xb)))
Da_xaxa = truncate(Da_xa*xa+x_nega*Da_xa)
Da_xaxb = truncate(Da_xa*xb+x_nega*Da_xb)
Da_xbxb = truncate(Da_xb*xb+saxb*Da_xb)
Da_xaxaxa = truncate(Da_xa*xa*xa+x_nega*Da_xaxa)
Da_xaxaxb = truncate(Da_xa*xa*xb+x_nega*Da_xaxb)
Da_xaxbxb = truncate(Da_xa*xb*xb+x_nega*Da_xbxb)
Da_xbxbxb = truncate(Da_xb*xb*xb+saxb*Da_xbxb)
Da_xaxaxaxa = truncate(Da_xa*xa*xa*xa+x_nega*Da_xaxaxa)
Da_xaxaxaxb = truncate(Da_xa*xa*xa*xb+x_nega*Da_xaxaxb)
Da_xaxaxbxb = truncate(Da_xa*xa*xb*xb+x_nega*Da_xaxbxb)
Da_xaxbxbxb = truncate(Da_xa*xb*xb*xb+x_nega*Da_xbxbxb)
Da_xbxbxbxb = truncate(Da_xb*xb*xb*xb+saxb*Da_xbxbxb)
Da_xaxaxaxaxa = truncate(Da_xa*xa*xa*xa*xa+x_nega*Da_xaxaxaxa)
Da_xaxaxaxaxb = truncate(Da_xa*xa*xa*xa*xb+x_nega*Da_xaxaxaxb)
Da_xaxaxaxbxb = truncate(Da_xa*xa*xa*xb*xb+x_nega*Da_xaxaxbxb)
Da_xaxaxbxbxb = truncate(Da_xa*xa*xb*xb*xb+x_nega*Da_xaxbxbxb)
Da_xaxbxbxbxb = truncate(Da_xa*xb*xb*xb*xb+x_nega*Da_xbxbxbxb)
Da_xbxbxbxbxb = truncate(Da_xb*xb*xb*xb*xb+saxb*Da_xbxbxbxb)
Da_xaxaxaxaxaxa = truncate(Da_xa*xa*xa*xa*xa*xa+x_nega*Da_xaxaxaxaxa)
Da_xaxaxaxaxaxb = truncate(Da_xa*xa*xa*xa*xa*xb+x_nega*Da_xaxaxaxaxb)
Da_xaxaxaxaxbxb = truncate(Da_xa*xa*xa*xa*xb*xb+x_nega*Da_xaxaxaxbxb)
Da_xaxaxaxbxbxb = truncate(Da_xa*xa*xa*xb*xb*xb+x_nega*Da_xaxaxbxbxb)
Da_xaxaxbxbxbxb = truncate(Da_xa*xa*xb*xb*xb*xb+x_nega*Da_xaxbxbxbxb)
Da_xaxbxbxbxbxb = truncate(Da_xa*xb*xb*xb*xb*xb+x_nega*Da_xbxbxbxbxb)
Da_xbxbxbxbxbxb = truncate(Da_xb*xb*xb*xb*xb*xb+saxb*Da_xbxbxbxbxb)
Da = [0,Da_xb,Da_xbxb,Da_xbxbxb,Da_xbxbxbxb,Da_xbxbxbxbxb,Da_xbxbxbxbxbxb,Da_xa,Da_xaxb,Da_xaxbxb,Da_xaxbxbxb,Da_xaxbxbxbxb,Da_xaxbxbxbxbxb,Da_xaxa,Da_xaxaxb,Da_xaxaxbxb,Da_xaxaxbxbxb,Da_xaxaxbxbxbxb,Da_xaxaxa,Da_xaxaxaxb,Da_xaxaxaxbxb,Da_xaxaxaxbxbxb,Da_xaxaxaxa,Da_xaxaxaxaxb,Da_xaxaxaxaxbxb,Da_xaxaxaxaxa,Da_xaxaxaxaxaxb,Da_xaxaxaxaxaxa]
 

Db_xa = truncate(divide_xb(poly(xa-sbxa,xa,xb)))
Db_xb = truncate(divide_xb(poly(xb-x_negb,xa,xb)))
Db_xaxa = truncate(Db_xa*xa+sbxa*Db_xa)
Db_xaxb = truncate(Db_xb*xa+x_negb*Db_xa)
Db_xbxb = truncate(Db_xb*xb+x_negb*Db_xb)
Db_xaxaxa = truncate(Db_xa*xa*xa+sbxa*Db_xaxa)
Db_xaxaxb = truncate(Db_xb*xa*xa+x_negb*Db_xaxa)
Db_xaxbxb = truncate(Db_xb*xa*xb+x_negb*Db_xaxb)
Db_xbxbxb = truncate(Db_xb*xb*xb+x_negb*Db_xbxb)
Db_xaxaxaxa = truncate(Db_xa*xa*xa*xa+sbxa*Db_xaxaxa)
Db_xaxaxaxb = truncate(Db_xb*xa*xa*xa+x_negb*Db_xaxaxa)
Db_xaxaxbxb = truncate(Db_xb*xa*xa*xb+x_negb*Db_xaxaxb)
Db_xaxbxbxb = truncate(Db_xb*xa*xb*xb+x_negb*Db_xaxbxb)
Db_xbxbxbxb = truncate(Db_xb*xb*xb*xb+x_negb*Db_xbxbxb)
Db_xaxaxaxaxa = truncate(Db_xa*xa*xa*xa*xa+sbxa*Db_xaxaxaxa)
Db_xaxaxaxaxb = truncate(Db_xb*xa*xa*xa*xa+x_negb*Db_xaxaxaxa)
Db_xaxaxaxbxb = truncate(Db_xb*xa*xa*xa*xb+x_negb*Db_xaxaxaxb)
Db_xaxaxbxbxb = truncate(Db_xb*xa*xa*xb*xb+x_negb*Db_xaxaxbxb)
Db_xaxbxbxbxb = truncate(Db_xb*xa*xb*xb*xb+x_negb*Db_xaxbxbxb)
Db_xbxbxbxbxb = truncate(Db_xb*xb*xb*xb*xb+x_negb*Db_xbxbxbxb)
Db_xaxaxaxaxaxa = truncate(Db_xa*xa*xa*xa*xa*xa+sbxa*Db_xaxaxaxaxa)
Db_xaxaxaxaxaxb = truncate(Db_xb*xa*xa*xa*xa*xa+x_negb*Db_xaxaxaxaxa)
Db_xaxaxaxaxbxb = truncate(Db_xb*xa*xa*xa*xa*xb+x_negb*Db_xaxaxaxaxb)
Db_xaxaxaxbxbxb = truncate(Db_xb*xa*xa*xa*xb*xb+x_negb*Db_xaxaxaxbxb)
Db_xaxaxbxbxbxb = truncate(Db_xb*xa*xa*xb*xb*xb+x_negb*Db_xaxaxbxbxb)
Db_xaxbxbxbxbxb = truncate(Db_xb*xa*xb*xb*xb*xb+x_negb*Db_xaxbxbxbxb)
Db_xbxbxbxbxbxb = truncate(Db_xb*xb*xb*xb*xb*xb+x_negb*Db_xbxbxbxbxb)
Db = [0,Db_xb,Db_xbxb,Db_xbxbxb,Db_xbxbxbxb,Db_xbxbxbxbxb,Db_xbxbxbxbxbxb,Db_xa,Db_xaxb,Db_xaxbxb,Db_xaxbxbxb,Db_xaxbxbxbxb,Db_xaxbxbxbxbxb,Db_xaxa,Db_xaxaxb,Db_xaxaxbxb,Db_xaxaxbxbxb,Db_xaxaxbxbxbxb,Db_xaxaxa,Db_xaxaxaxb,Db_xaxaxaxbxb,Db_xaxaxaxbxbxb,Db_xaxaxaxa,Db_xaxaxaxaxb,Db_xaxaxaxaxbxb,Db_xaxaxaxaxa,Db_xaxaxaxaxaxb,Db_xaxaxaxaxaxa]

def compose_Da(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4,5,6]:
       for j in [0,1,2,3,4,5,6]:
           if ((i+j)<=6):
               g = poly(g + Da[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

def compose_Db(f):
   p = 0
   g = poly(0,xa,xb)
   for i in [0,1,2,3,4,5,6]:
       for j in [0,1,2,3,4,5,6]:
           if ((i+j)<=6):
               g = poly(g + Db[p]*(f.coeff_monomial((xa**i)*(xb**j))),xa,xb)
               p = p+1
   return g

DbDa_xa = compose_Db(Da_xa)
DaDb_xa = compose_Da(Db_xa)
DbDaDb_xa = compose_Db(DaDb_xa)
DaDbDa_xa = compose_Da(DbDa_xa)
DaDbDaDb_xa = compose_Da(DbDaDb_xa)
DbDaDbDa_xa = compose_Db(DaDbDa_xa)
DaDbDaDbDa_xa = compose_Da(DbDaDbDa_xa)
DbDaDbDaDb_xa = compose_Db(DaDbDaDb_xa)
DaDbDaDbDaDb_xa = compose_Da(DbDaDbDaDb_xa)

DbDa_xb = compose_Db(Da_xb)
DaDb_xb = compose_Da(Db_xb)
DbDaDb_xb = compose_Db(DaDb_xb)
DaDbDa_xb = compose_Da(DbDa_xb)
DaDbDaDb_xb = compose_Da(DbDaDb_xb)
DbDaDbDa_xb = compose_Db(DaDbDa_xb)
DaDbDaDbDa_xb = compose_Da(DbDaDbDa_xb)
DbDaDbDaDb_xb = compose_Db(DaDbDaDb_xb)
DaDbDaDbDaDb_xb = compose_Da(DbDaDbDaDb_xb)

da_xa = Da_xa.coeff_monomial(1)
db_xa = Db_xa.coeff_monomial(1)
dadb_xa = DaDb_xa.coeff_monomial(1)
dbda_xa = DbDa_xa.coeff_monomial(1)
dadbda_xa = DaDbDa_xa.coeff_monomial(1)
dbdadb_xa = DbDaDb_xa.coeff_monomial(1)
dadbdadb_xa = DaDbDaDb_xa.coeff_monomial(1)
dbdadbda_xa = DbDaDbDa_xa.coeff_monomial(1)
dadbdadbda_xa = DaDbDaDbDa_xa.coeff_monomial(1)
dbdadbdadb_xa = DbDaDbDaDb_xa.coeff_monomial(1)
dadbdadbdadb_xa = DaDbDaDbDaDb_xa.coeff_monomial(1)

da_xb = Da_xb.coeff_monomial(1)
db_xb = Db_xb.coeff_monomial(1)
dadb_xb = DaDb_xb.coeff_monomial(1)
dbda_xb = DbDa_xb.coeff_monomial(1)
dadbda_xb = DaDbDa_xb.coeff_monomial(1)
dbdadb_xb = DbDaDb_xb.coeff_monomial(1)
dadbdadb_xb = DaDbDaDb_xb.coeff_monomial(1)
dbdadbda_xb = DbDaDbDa_xb.coeff_monomial(1)
dadbdadbda_xb = DaDbDaDbDa_xb.coeff_monomial(1)
dbdadbdadb_xb = DbDaDbDaDb_xb.coeff_monomial(1)
dadbdadbdadb_xb = DaDbDaDbDaDb_xb.coeff_monomial(1)

ev_xa = da_xa*Xa+db_xa*Xb+dadb_xa*Xab+dbda_xa*Xba+dadbda_xa*Xaba+dbdadb_xa*Xbab+dadbdadb_xa*Xabab+dbdadbda_xa*Xbaba+dadbdadbda_xa*Xababa+dbdadbdadb_xa*Xbabab+dadbdadbdadb_xa*Xababab
ev_xb = da_xb*Xa+db_xb*Xb+dadb_xb*Xab+dbda_xb*Xba+dadbda_xb*Xaba+dbdadb_xb*Xbab+dadbdadb_xb*Xabab+dbdadbda_xb*Xbaba+dadbdadbda_xb*Xababa+dbdadbdadb_xb*Xbabab+dadbdadbdadb_xb*Xababab

print("ev_xa = ", ev_xa)
print("ev_xb =", ev_xb)

Da_sbxa = compose_Da(sbxa)
Db_saxb = compose_Db(saxb)
Db_sasbxa = compose_Db(sasbxa)
Da_sbsaxb = compose_Da(sbsaxb)
Da_sbsasbxa = compose_Da(sbsasbxa)
Da_sbsasbsaxb = compose_Da(sbsasbsaxb)
Db_sasbsaxb = compose_Db(sasbsaxb)
DaDb_saxb = compose_Da(compose_Db(saxb))
DaDb_xasaxb = compose_Da(compose_Db(xa*saxb))
DbDa_xbsaxb = compose_Db(compose_Da(xb*saxb))
DbDa_xbsbxa = compose_Db(compose_Da(xb*sbxa))
DbDa_sbxa = compose_Db(compose_Da(sbxa))

db_saxb = Db_saxb.coeff_monomial(1)
da_sbxa = Da_sbxa.coeff_monomial(1)
db_sasbxa = Db_sasbxa.coeff_monomial(1)
da_sbsaxb = Da_sbsaxb.coeff_monomial(1)
da_sbsasbxa = Da_sbsasbxa.coeff_monomial(1)
da_sbsasbsaxb = Da_sbsasbsaxb.coeff_monomial(1)
db_sasbsaxb = Db_sasbsaxb.coeff_monomial(1)
dadb_saxb = DaDb_saxb.coeff_monomial(1)
dadb_xasaxb = DaDb_xasaxb.coeff_monomial(1)
dbda_xbsaxb = DbDa_xbsaxb.coeff_monomial(1)
dbda_xbsbxa = DbDa_xbsbxa.coeff_monomial(1)
dbda_sbxa = DbDa_sbxa.coeff_monomial(1)



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
#XbXaXbXa
q_a_a_baba = Symbol('q_a_a_baba')
q_a_b_baba = Symbol('q_a_b_baba')
q_a_ab_baba = Symbol('q_a_ab_baba')
q_a_ba_baba = Symbol('q_a_ba_baba')
q_a_aba_baba = -db_xa-db_sasbxa
q_a_bab_baba = 1
q_b_b_baba = Symbol('q_b_b_baba')
q_b_ab_baba = Symbol('q_b_ab_baba')
q_b_ba_baba = Symbol('q_b_ba_baba')
q_b_aba_baba = 1-db_saxb
q_b_bab_baba = 0
q_ab_ab_baba = Symbol('q_ab_ab_baba')
q_ab_ba_baba = Symbol('q_ab_ba_baba')
q_ab_aba_baba = 0
q_ab_bab_baba = 0
q_ba_ba_baba = Symbol('q_ba_ba_baba')
q_ba_aba_baba = 0
q_ba_bab_baba = 0
q_aba_aba_baba = 0
q_aba_bab_baba = 0
q_bab_bab_baba = 0
#XaXbXaXbXa
q_a_a_ababa = Symbol('q_a_a_ababa')
q_a_b_ababa = Symbol('q_a_b_ababa')
q_a_ab_ababa = Symbol('q_a_ab_ababa')
q_a_ba_ababa = Symbol('q_a_ba_ababa')
q_a_aba_ababa = Symbol('q_a_aba_ababa')
q_a_bab_ababa = Symbol('q_a_bab_ababa')
q_a_abab_ababa = 1
q_a_baba_ababa = 1-da_sbxa-da_sbsasbxa
q_b_b_ababa = Symbol('q_b_b_ababa')
q_b_ab_ababa = Symbol('q_b_ab_ababa')
q_b_ba_ababa = Symbol('q_b_ba_ababa')
q_b_aba_ababa = Symbol('q_b_aba_ababa')
q_b_bab_ababa = Symbol('q_b_bab_ababa')
q_b_abab_ababa = 0
q_b_baba_ababa = -da_xb-da_sbsaxb
q_ab_ab_ababa = Symbol('q_ab_ab_ababa')
q_ab_ba_ababa = Symbol('q_ab_ba_ababa')
q_ab_aba_ababa = Symbol('q_ab_aba_ababa')
q_ab_bab_ababa = Symbol('q_ab_bab_ababa')
q_ab_abab_ababa = 0
q_ab_baba_ababa = 0
q_ba_ba_ababa = Symbol('q_ba_ba_ababa')
q_ba_aba_ababa = Symbol('q_ba_aba_ababa')
q_ba_bab_ababa = Symbol('q_ba_bab_ababa')
q_ba_abab_ababa = 0
q_ba_baba_ababa = 0
q_aba_aba_ababa = Symbol('q_aba_aba_ababa')
q_aba_bab_ababa = Symbol('q_aba_bab_ababa')
q_aba_abab_ababa = 0
q_aba_baba_ababa = 0
q_bab_bab_ababa = Symbol('q_bab_bab_ababa')
q_bab_abab_ababa = 0
q_bab_baba_ababa = 0
q_abab_abab_ababa = 0
q_abab_baba_ababa = 0
q_baba_baba_ababa = 0
#XbXaXbXaXb
q_a_a_babab = Symbol('q_a_a_babab')
q_a_b_babab = Symbol('q_a_b_babab')
q_a_ab_babab = Symbol('q_a_ab_babab')
q_a_ba_babab = Symbol('q_a_ba_babab')
q_a_aba_babab = Symbol('q_a_aba_babab')
q_a_bab_babab = Symbol('q_a_bab_babab')
q_a_abab_babab = -db_xa-db_sasbxa
q_a_baba_babab = 0
q_b_b_babab = Symbol('q_b_b_babab')
q_b_ab_babab = Symbol('q_b_ab_babab')
q_b_ba_babab = Symbol('q_b_ba_babab')
q_b_aba_babab = Symbol('q_b_aba_babab')
q_b_bab_babab = Symbol('q_b_bab_babab')
q_b_abab_babab = 1-db_saxb-db_sasbsaxb
q_b_baba_babab = 1
q_ab_ab_babab = Symbol('q_ab_ab_babab')
q_ab_ba_babab = Symbol('q_ab_ba_babab')
q_ab_aba_babab = Symbol('q_ab_aba_babab')
q_ab_bab_babab = Symbol('q_ab_bab_babab')
q_ab_abab_babab = 0
q_ab_baba_babab = 0
q_ba_ba_babab = Symbol('q_ba_ba_babab')
q_ba_aba_babab = Symbol('q_ba_aba_babab')
q_ba_bab_babab = Symbol('q_ba_bab_babab')
q_ba_abab_babab = 0
q_ba_baba_babab = 0
q_aba_aba_babab = Symbol('q_aba_aba_babab')
q_aba_bab_babab = Symbol('q_aba_bab_babab')
q_aba_abab_babab = 0
q_aba_baba_babab = 0
q_bab_bab_babab = Symbol('q_bab_bab_babab')
q_bab_abab_babab = 0
q_bab_baba_babab = 0
q_abab_abab_babab = 0
q_abab_baba_babab = 0
q_baba_baba_babab = 0
#XaXbXaXbXaXb
q_a_a_ababab = Symbol('q_a_a_ababab')
q_a_b_ababab = Symbol('q_a_b_ababab')
q_a_ab_ababab = Symbol('q_a_ab_ababab')
q_a_ba_ababab = Symbol('q_a_ba_ababab')
q_a_aba_ababab = Symbol('q_a_aba_ababab')
q_a_bab_ababab = Symbol('q_a_bab_ababab')
q_a_abab_ababab = Symbol('q_a_abab_ababab')
q_a_baba_ababab = Symbol('q_a_baba_ababab')
q_a_ababa_ababab = 0
q_a_babab_ababab = 1-da_sbxa-da_sbsasbxa
q_b_b_ababab = Symbol('q_b_b_ababab')
q_b_ab_ababab = Symbol('q_b_ab_ababab')
q_b_ba_ababab = Symbol('q_b_ba_ababab')
q_b_aba_ababab = Symbol('q_b_aba_ababab')
q_b_bab_ababab = Symbol('q_b_bab_ababab')
q_b_abab_ababab = Symbol('q_b_abab_ababab')
q_b_baba_ababab = Symbol('q_b_baba_ababab')
q_b_ababa_ababab = 1
q_b_babab_ababab = -da_xb-da_sbsaxb-da_sbsasbsaxb
q_ab_ab_ababab = Symbol('q_ab_ab_ababab')
q_ab_ba_ababab = Symbol('q_ab_ba_ababab')
q_ab_aba_ababab = Symbol('q_ab_aba_ababab')
q_ab_bab_ababab= Symbol('q_ab_bab_ababab')
q_ab_abab_ababab = Symbol('q_ab_abab_ababab')
q_ab_baba_ababab = Symbol('q_ab_baba_ababab')
q_ab_ababa_ababab = 0
q_ab_babab_ababab = 0
q_ba_ba_ababab = Symbol('q_ba_ba_ababab')
q_ba_aba_ababab = Symbol('q_ba_aba_ababab')
q_ba_bab_ababab = Symbol('q_ba_bab_ababab')
q_ba_abab_ababab = Symbol('q_ba_abab')
q_ba_baba_ababab = Symbol('q_ba_baba_ababab')
q_ba_ababa_ababab = 0
q_ba_babab_ababab = 0
q_aba_aba_ababab = Symbol('q_aba_aba_ababab')
q_aba_bab_ababab = Symbol('q_aba_bab_ababab')
q_aba_abab_ababab = Symbol('q_aba_abab_ababab')
q_aba_baba_ababab = Symbol('q_aba_baba_ababab')
q_aba_ababa_ababab = 0
q_aba_babab_ababab = 0
q_bab_bab_ababab = Symbol('q_bab_bab_ababab')
q_bab_abab_ababab = Symbol('q_bab_abab_ababab')
q_bab_baba_ababab = Symbol('q_bab_baba_ababab')
q_bab_ababa_ababab = 0
q_bab_babab_ababab = 0
q_abab_abab_ababab = Symbol('q_abab_abab_ababab')
q_abab_baba_ababab = Symbol('q_abab_baba_ababab')
q_abab_ababa_ababab = 0
q_abab_babab_ababab = 0
q_baba_baba_ababab = Symbol('q_baba_baba_ababab')
q_baba_ababa_ababab = 0
q_baba_babab_ababab = 0
q_ababa_ababa_ababab = 0
q_ababa_babab_ababab = 0
q_babab_babab_ababab = 0

#Coefficients of relation M
print("")
print("Second Relations (Products):")
print("")
print("XaXa=", Symbol.simplify(q_a_a_ab*Xab+q_a_a_ba*Xba+q_a_a_aba*Xaba+q_a_a_bab*Xbab+q_a_a_abab*Xabab+q_a_a_baba*Xbaba+q_a_a_ababa*Xababa+q_a_a_babab*Xbabab+q_a_a_ababab*Xababab))
print("") 
print("XaXb=", Symbol.simplify(q_a_b_ab*Xab+q_a_b_ba*Xba+q_a_b_aba*Xaba+q_a_b_bab*Xbab+q_a_b_abab*Xabab+q_a_b_baba*Xbaba+q_a_b_ababa*Xababa+q_a_b_babab*Xbabab+q_a_b_ababab*Xababab))
print("")   
print("XaXab=", Symbol.simplify(q_a_ab_aba*Xaba+q_a_ab_bab*Xbab+q_a_ab_abab*Xabab+q_a_ab_baba*Xbaba+q_a_ab_ababa*Xababa+q_a_ab_babab*Xbabab+q_a_ab_ababab*Xababab))
print("")   
print("XaXba=", Symbol.simplify(q_a_ba_aba*Xaba+q_a_ba_bab*Xbab+q_a_ba_abab*Xabab+q_a_ba_baba*Xbaba+q_a_ba_ababa*Xababa+q_a_ba_babab*Xbabab+q_a_ba_ababab*Xababab))
print("")   
print("XaXaba=", Symbol.simplify(q_a_aba_abab*Xabab+q_a_aba_baba*Xbaba+q_a_aba_ababa*Xababa+q_a_aba_babab*Xbabab+q_a_aba_ababab*Xababab))
print("")   
print("XaXbab=", Symbol.simplify(q_a_bab_abab*Xabab+q_a_bab_baba*Xbaba+q_a_bab_ababa*Xababa+q_a_bab_babab*Xbabab+q_a_bab_ababab*Xababab))
print("")   
print("XaXabab=", Symbol.simplify(q_a_abab_ababa*Xababa+q_a_abab_babab*Xbabab+q_a_abab_ababab*Xababab))
print("")   
print("XaXbaba=", Symbol.simplify(q_a_baba_ababa*Xababa+q_a_baba_babab*Xbabab+q_a_baba_ababab*Xababab))
print("")   
print("XaXababa=", Symbol.simplify(q_a_ababa_ababab*Xababab))
print("")   
print("XaXbabab=", Symbol.simplify(q_a_babab_ababab*Xababab))
print("")   
print("XaXababab=", 0)
print("")   
print("XbXb=", Symbol.simplify(q_b_b_ab*Xab+q_b_b_ba*Xba+q_b_b_aba*Xaba+q_b_b_bab*Xbab+q_b_b_abab*Xabab+q_b_b_baba*Xbaba+q_b_b_ababa*Xababa+q_b_b_babab*Xbabab+q_b_b_ababab*Xababab))
print("")   
print("XbXab=", Symbol.simplify(q_b_ab_aba*Xaba+q_b_ab_bab*Xbab+q_b_ab_abab*Xabab+q_b_ab_baba*Xbaba+q_b_ab_ababa*Xababa+q_b_ab_babab*Xbabab+q_b_ab_ababab*Xababab))
print("")   
print("XbXba=", Symbol.simplify(q_b_ba_aba*Xaba+q_b_ba_bab*Xbab+q_b_ba_abab*Xabab+q_b_ba_baba*Xbaba+q_b_ba_ababa*Xababa+q_b_ba_babab*Xbabab+q_b_ba_ababab*Xababab))
print("")   
print("XbXaba=", Symbol.simplify(q_b_aba_abab*Xabab+q_b_aba_baba*Xbaba+q_b_aba_ababa*Xababa+q_b_aba_babab*Xbabab+q_b_aba_ababab*Xababab))
print("")   
print("XbXbab=", Symbol.simplify(q_b_bab_abab*Xabab+q_b_bab_baba*Xbaba+q_b_bab_ababa*Xababa+q_b_bab_babab*Xbabab+q_b_bab_ababab*Xababab))
print("")   
print("XbXabab=", Symbol.simplify(q_b_abab_ababa*Xababa+q_b_abab_babab*Xbabab+q_b_abab_ababab*Xababab))
print("")
print("XbXbaba=", Symbol.simplify(q_b_baba_ababa*Xababa+q_b_baba_babab*Xbabab+q_b_baba_ababab*Xababab))
print("")   
print("XbXababa=", Symbol.simplify(q_b_ababa_ababab*Xababab))
print("") 
print("XbXbabab=", Symbol.simplify(q_a_babab_ababab*Xababab))
print("")       
print("XbXababab=", 0)
print("") 
print("XabXab=", Symbol.simplify(q_ab_ab_aba*Xaba+q_ab_ab_bab*Xbab+q_ab_ab_abab*Xabab+q_ab_ab_baba*Xbaba+q_ab_ab_ababa*Xababa+q_ab_ab_babab*Xbabab+q_ab_ab_ababab*Xababab))
print("")   
print("XabXba=", Symbol.simplify(q_ab_ba_aba*Xaba+q_ab_ba_bab*Xbab+q_ab_ba_abab*Xabab+q_ab_ba_baba*Xbaba+q_ab_ba_ababa*Xababa+q_ab_ba_babab*Xbabab+q_ab_ba_ababab*Xababab))
print("")   
print("XabXaba=", Symbol.simplify(q_ab_aba_abab*Xabab+q_ab_aba_baba*Xbaba+q_ab_aba_ababa*Xababa+q_ab_aba_babab*Xbabab+q_ab_aba_ababab*Xababab))
print("")   
print("XabXbab=", Symbol.simplify(q_ab_bab_abab*Xabab+q_ab_bab_baba*Xbaba+q_ab_bab_ababa*Xababa+q_ab_bab_babab*Xbabab+q_ab_bab_ababab*Xababab))
print("")   
print("XabXabab=", Symbol.simplify(q_ab_abab_ababa*Xababa+q_ab_abab_babab*Xbabab+q_ab_abab_ababab*Xababab))
print("")   
print("XabXbaba=", Symbol.simplify(q_ab_baba_ababa*Xababa+q_ab_baba_babab*Xbabab+q_ab_baba_ababab*Xababab))
print("")   
print("XabXababa=", Symbol.simplify(q_ab_ababa_ababab*Xababab))
print("")   
print("XabXbabab=", Symbol.simplify(q_ab_babab_ababab*Xababab))
print("")   
print("XabXababab=", 0)
print("")   
print("XbaXba=", Symbol.simplify(q_ba_ba_aba*Xaba+q_ba_ba_bab*Xbab+q_ba_ba_abab*Xabab+q_ba_ba_baba*Xbaba+q_ba_ba_ababa*Xababa+q_ba_ba_babab*Xbabab+q_ba_ba_ababab*Xababab))
print("")   
print("XbaXaba=", Symbol.simplify(q_ba_aba_abab*Xabab+q_ba_aba_baba*Xbaba+q_ba_aba_ababa*Xababa+q_ba_aba_babab*Xbabab+q_ba_aba_ababab*Xababab))
print("")   
print("XbaXbab=", Symbol.simplify(q_ba_bab_abab*Xabab+q_ba_bab_baba*Xbaba+q_ba_bab_ababa*Xababa+q_ba_bab_babab*Xbabab+q_ba_bab_ababab*Xababab))
print("")   
print("XbaXabab=", Symbol.simplify(q_ba_abab_ababa*Xababa+q_ba_abab_babab*Xbabab+q_ba_abab_ababab*Xababab))
print("")   
print("XbaXbaba=", Symbol.simplify(q_ba_baba_ababa*Xababa+q_ba_baba_babab*Xbabab+q_ba_baba_ababab*Xababab))
print("")   
print("XbaXababa=", Symbol.simplify(q_ba_ababa_ababab*Xababab))
print("")   
print("XbaXbabab=", Symbol.simplify(q_ba_babab_ababab*Xababab))
print("")   
print("XbaXababab=", 0)
print("")   
print("XabaXaba=", Symbol.simplify(q_aba_aba_abab*Xabab+q_aba_aba_baba*Xbaba+q_aba_aba_ababa*Xababa+q_aba_aba_babab*Xbabab+q_aba_aba_ababab*Xababab))
print("")   
print("XabaXbab=", Symbol.simplify(q_aba_bab_abab*Xabab+q_aba_bab_baba*Xbaba+q_aba_bab_ababa*Xababa+q_aba_bab_babab*Xbabab+q_aba_bab_ababab*Xababab))
print("")   
print("XabaXabab=", Symbol.simplify(q_aba_abab_ababa*Xababa+q_aba_abab_babab*Xbabab+q_aba_abab_ababab*Xababab))
print("")   
print("XabaXbaba=", Symbol.simplify(q_aba_baba_ababa*Xababa+q_aba_baba_babab*Xbabab+q_aba_baba_ababab*Xababab))
print("")   
print("XabaXababa=", Symbol.simplify(q_aba_ababa_ababab*Xababab))
print("")   
print("XabaXbabab=", Symbol.simplify(q_aba_babab_ababab*Xababab))
print("")   
print("XabaXababab=", 0)
print("")   
print("XbabXbab=", Symbol.simplify(q_bab_bab_abab*Xabab+q_bab_bab_baba*Xbaba+q_bab_bab_ababa*Xababa+q_bab_bab_babab*Xbabab+q_bab_bab_ababab*Xababab))
print("")   
print("XbabXabab=", Symbol.simplify(q_bab_abab_ababa*Xababa+q_bab_abab_babab*Xbabab+q_bab_abab_ababab*Xababab))
print("")   
print("XbabXbaba=", Symbol.simplify(q_bab_baba_ababa*Xababa+q_bab_baba_babab*Xbabab+q_bab_baba_ababab*Xababab))
print("")   
print("XbabXababa=", Symbol.simplify(q_bab_ababa_ababab*Xababab))
print("")   
print("XbabXbabab=", Symbol.simplify(q_bab_babab_ababab*Xababab))
print("")   
print("XbabXababab=", 0)
print("")   
print("XababXabab=", Symbol.simplify(q_abab_abab_ababa*Xababa+q_abab_abab_babab*Xbabab+q_abab_abab_ababab*Xababab))
print("")   
print("XababXbaba=", Symbol.simplify(q_abab_baba_ababa*Xababa+q_abab_baba_babab*Xbabab+q_abab_baba_ababab*Xababab))
print("")   
print("XababXababa=", Symbol.simplify(q_abab_ababa_ababab*Xababab))
print("")   
print("XababXbabab=", Symbol.simplify(q_abab_babab_ababab*Xababab))
print("")   
print("XababXababab=", 0)
print("")   
print("XbabaXbaba=", Symbol.simplify(q_baba_baba_ababa*Xababa+q_baba_baba_babab*Xbabab+q_baba_baba_ababab*Xababab))
print("")   
print("XbabaXababa=", Symbol.simplify(q_baba_ababa_ababab*Xababab))
print("") 
print("XbabaXbabab=", Symbol.simplify(q_baba_babab_ababab*Xababab))
print("") 
print("XbabaXababab=", 0)
print("") 
print("XababaXababa=", Symbol.simplify(q_ababa_ababa_ababab*Xababab))
print("")   
print("XababaXbabab=", Symbol.simplify(q_ababa_babab_ababab*Xababab))
print("")   
print("XababaXababab=", 0)
print("")   
print("XbababXbabab=", Symbol.simplify(q_babab_babab_ababab*Xababab))
print("")   
print("XbababXababab=", 0)
print("")   
print("XabababXababab=", 0)
print("")   

