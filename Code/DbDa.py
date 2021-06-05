# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 21:31:02 2021

@author: Raj
"""


from sympy import Symbol, cos
#The Da and Db are the Demazure operators
Da = Symbol('D_a')
Db = Symbol('D_b')
#The xa and xb are the generators of the formal group algebra
xa = Symbol('x_a')
xb = Symbol('x_b')
#sa is a reflection along the root alpha and sb is defined similarly
xasa = Symbol('x_a*s_a')
xbsb = Symbol('x_b*s_b')
sa = Symbol('s_a')
sb = Symbol('s_b')
Aa = Symbol('a')
Bb = Symbol('b')


#The powerset function
def powerset(s):
    x = len(s)
    masks = [1 << i for i in range(x)]
    for i in range(1 << x):
        yield [ss for mask, ss in zip(masks, s) if i & mask]


#X stores the coefficients p_{E1,E2}^I of the coproduct
X = [] 

#Puts the coefficients in X
for E1 in powerset([1,2]):
    for E2 in powerset([1,2]):
        if (1 in E1) & (1 in E2): 
            B1 = -xbsb
        elif (1 not in E1) & (1 not in E2):
            B1 = Db
        else:
            B1 = sb
        if (2 in E1) & (2 in E2): 
            B2 = -xasa
        elif (2 not in E1) & (2 not in E2):
            B2 = Da
        else:
            B2 = sa
        S1 = [];
        S2 = [];
        for x in E1:
            if (x == 1):
                S1.append(Bb)
            elif (x == 2):
                S1.append(Aa)
        for y in E2:
            if (y == 1):
                S2.append(Bb)
            elif (y == 2):
                S2.append(Aa)
        X.append([[S1,S2],[B1,B2]])



#X0 is X minus the coefficients that are clearly 0 in the augmented Demazure algebra                
X0 = []   
sab = {xasa,xbsb}
Dab = {Da,Db}
b = 0
 
#Removes most coefficients that are 0 from X and stores the remaining coefficients in X0                             
for x in X:
    if ((x[1][1] == Da) | (x[1][1] == Db)):
        b = 1
    elif ((x[1][0] == Da) | (x[1][0] == Db)) & (-x[1][1] not in sab):
        b = 1
    elif ((x[1][0] not in Dab) & (x[1][1] not in Dab)) & ((-x[1][0] in sab) | (-x[1][1] in sab)):
        b = 1
    if b == 0:
        X0.append(x)
    b = 0      
    
    
#Y stores the subsets of [1,2,3]; elements not necessarily unique                        
Y = []

for E1 in powerset([Bb,Aa]):
    Y.append(E1)


#Z stores the elements of Y minus any repeated elements
Z = []
for y in Y:
    if y not in Z:
        Z.append(y)


#ZZ1 stores pairs of elements in Z
ZZ1 = [] 
for z1 in Z:
    for z2 in Z:
        ZZ1.append([z1,z2])


#ZZ stores the elements of ZZ1 minus repeated elements
ZZ = []
for zz in ZZ1:
    if zz not in ZZ:
        ZZ.append(zz)

#R contains elements in ZZ that contain ...XaXa... or ...XbXb... in one of 
#the tensor factors in at least one of the sequences.

R = []

for zz in ZZ:
    if zz[0] != []:  
        l = len(zz[0])
        for i in range(l-1):
            if zz[0][i] == zz[0][i+1]:
                R.append(zz)  
                
for zz in ZZ:
    if zz[1] != []:  
        l = len(zz[1])
        for i in range(l-1):
            if zz[1][i] == zz[1][i+1]:
                R.append(zz)
                
#RR contains elements of RR without repetition
RR = []
for r in R:
    if r not in RR:
        RR.append(r)
        
#SS contains elements in ZZ that do NOT contain ...XaXa... or ...XbXb... in
# one of the tensor factors in at least one of the sequences.       
SS = []
for zz in ZZ:
    if zz not in RR:
        SS.append(zz)

            
            
#L1 stores the elements q_{w,w'}^v (the dual product coefficients) that do
#NOT contain ...XaXa... or ...XbXb... in one of the tensor factors 

print("Elements q_{w,w'}^v (the dual product coefficients) that do NOT contain ...XaXa... or ...XbXb... in one of the tensor factors ")                                 
L1 = []
t1 = []

for ss in SS:
    for x in X0:
        if (x[0] == ss):
            t1.append(x[1])
    L1.append([ss,t1])
    t1 = []
print(*L1, sep = "\n")


#L2 stores the elements q_{w,w'}^v (the dual product coefficients) that DO
#contain ...XaXa... or ...XbXb... in one of the tensor factors     
print("")
print("Elements q_{w,w'}^v (the dual product coefficients) that DO contain ...XaXa... or ...XbXb... in one of the tensor factors ")                                 
L2 = []
t2 = []
for rr in RR:
    for x in X0:
        if (x[0] == rr):
            t2.append(x[1])
    L2.append([rr,t2])
    t2 = []
print(*L2, sep = "\n")

