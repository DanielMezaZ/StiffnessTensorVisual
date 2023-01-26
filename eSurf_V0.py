# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:45:44 2023

@author: jesus
"""
import numpy as np
from matplotlib import pyplot as plt
E=1
nu=0.5
G=E/(2*(1+nu))
[  E_1,  E_2,  E_3]=[  E,  E,  E]
[ G_23, G_13, G_12]=[  G,  G,  G]
[nu_23,nu_13,nu_12]=[ nu, nu, nu]
S=np.zeros((3,3,3,3))
S[0][0][0][0]=1/E_1
S[1][1][1][1]=1/E_2
S[2][2][2][2]=1/E_3
S[1][2][1][2]=1/G_23
S[0][2][0][2]=1/G_13
S[0][1][0][1]=1/G_12

S[1][2][2][1]=S[1][2][1][2]
S[0][2][2][0]=S[0][2][0][2]
S[0][1][1][0]=S[0][1][0][1]
S[2][1][1][2]=S[1][2][1][2]
S[2][0][0][2]=S[0][2][0][2]
S[1][0][0][1]=S[0][1][0][1]

S[0][0][1][1]=-nu_12/E_2
S[0][0][2][2]=-nu_13/E_3
S[1][1][2][2]=-nu_23/E_3

S[2][2][0][0]=S[0][0][2][2]
S[1][1][0][0]=S[0][0][1][1]
S[2][2][1][1]=S[1][1][2][2]

S2=np.zeros((6,6))
S2[0][0]=S[0][0][0][0]
S2[1][1]=S[1][1][1][1]
S2[2][2]=S[2][2][2][2]
S2[3][3]=2*S[1][2][1][2]
S2[4][4]=2*S[0][2][0][2]
S2[5][5]=2*S[0][1][0][1]
S2[0][1]=S[0][0][1][1]
S2[0][2]=S[0][0][2][2]
S2[1][2]=S[1][1][2][2]
S2[1][0]=S2[0][1]
S2[2][0]=S2[0][2]
S2[2][1]=S2[1][2]
    
nElev=200
nAzim=200
elev_prev=np.linspace(0,1,nElev,endpoint=False)
elev=np.pi*elev_prev
azim=np.linspace(-np.pi,np.pi,nAzim,endpoint=False)

angles=np.meshgrid(azim,elev)
d_1=np.cos(angles[0])*np.sin(angles[1])
d_2=np.sin(angles[0])*np.sin(angles[1])
d_3=                  np.cos(angles[1])

d=[d_1,d_2,d_3]
#d=[1,0.2,0.2]
E_inv=0

for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                E_inv=E_inv+S[i][j][k][l]*d[i]*d[j]*d[j]*d[l]
E=1/E_inv
E_mag=np.zeros((nElev*nAzim,1))
E_vec=np.zeros((nElev*nAzim,3))
counter=0
for i in range(nElev):
    for j in range(nElev):
        E_vec[counter][0]=E[i][j]*d[0][i][j]
        E_vec[counter][1]=E[i][j]*d[1][i][j]
        E_vec[counter][2]=E[i][j]*d[2][i][j]
        E_mag[counter]= E_vec[counter][0]* E_vec[counter][0]+ \
            E_vec[counter][1]* E_vec[counter][1]+ E_vec[counter][2]* E_vec[counter][2]
        counter=counter+1
   
#E_test=E_vec[:,1]
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
values=E_mag
#surf = ax.scatter(d_1,d_2,d_3, c=values,cmap='coolwarm')
surf = ax.scatter(E_vec[:,0],E_vec[:,1],E_vec[:,2], c=values,cmap='coolwarm')