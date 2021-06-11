# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:00:08 2020

@author: Raj
"""

#Let $F(x,y)$ be one-dimensional commutative formal group law over a 
#commutative, unital ring $R$. The purpose of this program is to compute the
#coefficients $c_i$ of the formal inverse 
#$g(x)=-c_1*x-c_2*x^2-c_3*x^3-c_4*x^4-c_5*x^5-c_6*x^6$ of $x$ under $F(x,y)$  
#up to degree 6, in terms of the coefficients $a_{i,j}$ of the formal group
#law $F(x,y)=\sum_{i,j=1}^\infty a_{i,j}x^iy^j$.

from sympy import Symbol, poly, degree_list
from sympy.abc import x


#The $c_i$ are the coefficients of the formal inverse 
#$g(x)=-c_1*x-c_2*x^2-c_3*x^3-c_4*x^4-c_5*x^5-c_6*x^6$.
c1 = Symbol('c_1')
c2 = Symbol('c_2')
c3 = Symbol('c_3')
c4 = Symbol('c_4')
c5 = Symbol('c_5')
c6 = Symbol('c_6')

#The $a_{i,j}$ are the coefficients of the formal group law 
#$F(x,y)=\sum_{i,j=1}^\infty a_{i,j}x^iy^j$ up to degree 6 terms.
a11 = Symbol('a_11')
a12 = Symbol('a_12')
a13 = Symbol('a_13')
a22 = Symbol('a_22')
a14 = Symbol('a_14')
a23 = Symbol('a_23')
a15 = Symbol('a_15')
a24 = Symbol('a_24')
a33 = Symbol('a_33')

#$g$ is the formal inverse of $x$ under $F(x,y)$.
g = poly(-c1*x-c2*(x**2)-c3*(x**3)-c4*(x**4)-c5*(x**5)-c6*(x**6),x)


#$hij$ is defined to be $x^i*g(x)^j$.
h11 = poly(x*g,x)
h12 = poly(x*(g**2),x)
h21 = poly((x**2)*g,x)
h13 = poly(x*(g**3),x)
h22 = poly((x**2)*(g**2),x)
h31 = poly((x**3)*g,x)
h14 = poly(x*(g**4),x)
h23 = poly((x**2)*(g**3),x)
h32 = poly((x**3)*(g**2),x)
h41 = poly((x**4)*g,x)
h15 = poly(x*(g**5),x)
h24 = poly((x**2)*(g**4),x)
h33 = poly((x**3)*(g**3),x)
h42 = poly((x**4)*(g**2),x)
h51 = poly((x**5)*g,x)


#$f_i(x)$ is defined by the formula $F(x,g(x))=\sum_{i=1}^\infty f_i(x).
f1 = poly(x+g,x)
f2 = poly(a11*h11)
f3 = poly(a12*(h12+h21),x)
f4 = poly(a13*(h13+h31)+a22*(h22),x)
f5 = poly(a14*(h14+h41)+a23*(h23+h32),x)
f6 = poly(a15*(h15+h51)+a24*(h24+h42)+a33*(h33),x)


#$f$ is the sum $F(x,g(x))$ up to degree $6$ terms.
f = poly(f1+f2+f3+f4+f5+f6,x)

#We collect the coefficient of $x^i$ in $f(x)$ for $i<=6$. 
#Since $F(x,g(x))=0$ by definition, these coefficients must all be $0$.
#This allows us to express the $c_i$ in terms of the $a_{i,j}$.
for i in [1,2,3,4,5,6]:
    print("")
    print("Degree", i)
    print([f.coeff_monomial(x**i)])
