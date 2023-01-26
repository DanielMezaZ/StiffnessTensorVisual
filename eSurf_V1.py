# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:45:44 2023

@author: jesus
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.tri as mtri
import matplotlib
from  matplotlib.colors import LinearSegmentedColormap

cmap="coolwarm" #LinearSegmentedColormap.from_list('rg',["g", "w", "r"], N=256) 
fig, ax = plt.subplots(ncols=1,nrows=1,subplot_kw={"projection": "3d"})

###############################################################################

# #Play with numbers
# [  E_1,  E_2,  E_3]=[     1.3,    1.3,   1.3]
# [ G_23, G_13, G_12]=[   0.7,  0.7, 0.7]
# [nu_23,nu_13,nu_12]=[   0.25,  0.3, 0.3]

# #Play with numbers
# [  E_1,  E_2,  E_3]=[     2.3,    1.5,   1.0]
# [ G_23, G_13, G_12]=[   0.7,  0.9, 1.2]
# [nu_23,nu_13,nu_12]=[   0.25,  0.3, 0.1]

# #Isotropic
# E=2
# nu=0.2
# G=E/(2*(1+nu))
# K=E/(3*(1-2*nu))
# [  E_1,  E_2,  E_3]=[    E,   E,   E]
# [ G_23, G_13, G_12]=[    G,   G,   G]
# [nu_23,nu_13,nu_12]=[   nu,  nu,  nu]

# Transverse Isotropic
E_p=1.5
E_z=5
nu_p=0.2
nu_pz=0.4
G_p=E_p/(2*(1+nu_p))
G_pz=2
[  E_1,  E_2,  E_3]=[    E_p,   E_p,   E_z]
[ G_23, G_13, G_12]=[    G_pz,   G_pz,   G_p]
[nu_23,nu_13,nu_12]=[   nu_pz,  nu_pz,  nu_p] 

# #Honeycomb
# E_s=4
# t=0.25
# h=1         # "hori"
# l=1       # "vert"
# f=0.3       # Scale independent E- and G-modulus 
# s=0.3         # Scale independent nu
# a=np.pi/180*20
# E_1=E_s*(t/l)**3*np.cos(a)/(h/l+np.sin(a))/(np.sin(a)**2)
# E_2=E_s*(t/l)**3*(h/l+np.sin(a))/(np.cos(a)**3)
# E_3=E_s*f
# nu_12=np.cos(a)**2/((h/l+np.sin(a))*np.sin(a))
# nu_23=0.3*s
# nu_13=0.3*s
# G_12=E_s*(t/l)**3*(h/l+np.sin(a))/((h/l)**2*(1+2*h/l)*np.cos(a))
# G_23=E_s/(2*(1+0.3))*f
# G_13=E_s/(2*(1+0.3))*f

###############################################################################

def componentsC():
    ## Silver
    # [C11,C12,C13,C14,C15,C16]=[1.22,0.92,0.92,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.22,0.92,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.22,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.45,0.00,0.00]
    # [                C55,C56]=[                    0.45,0.00]
    # [                    C66]=[                         0.45]
    
    # # Copper
    # [C11,C12,C13,C14,C15,C16]=[1.68,1.21,1.21,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.68,1.21,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.68,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               1.51,0.00,0.00]
    # [                C55,C56]=[                    1.51,0.00]
    # [                    C66]=[                         1.51]
    
    # # Quartz
    # [C11,C12,C13,C14,C15,C16]=[0.87,0.07,0.12,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     0.87,0.12,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.06,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               1.16,0.00,0.00]
    # [                C55,C56]=[                    1.16,-.36]
    # [                    C66]=[                         0.80]

    # # Gold
    # [C11,C12,C13,C14,C15,C16]=[1.85,1.58,1.58,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.85,1.58,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.85,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.79,0.00,0.00]
    # [                C55,C56]=[                    0.79,0.00]
    # [                    C66]=[                         0.79]

    # # Aluminum phosphate
    # [C11,C12,C13,C14,C15,C16]=[1.05,0.29,0.69,-.18,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.05,0.69,0.18,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.34,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.46,0.00,0.00]
    # [                C55,C56]=[                    0.46,-.25]
    # [                    C66]=[                         0.76]
    
    # # Uranium    
    # [C11,C12,C13,C14,C15,C16]=[2.15,0.46,0.22,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.99,1.08,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          2.67,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               2.49,0.00,0.00]
    # [                C55,C56]=[                    1.47,0.00]
    # [                    C66]=[                         0.89]
    
    # # Magnesium
    # [C11,C12,C13,C14,C15,C16]=[5.65,2.32,1.81,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     5.65,1.81,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          5.87,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               3.36,0.00,0.00]
    # [                C55,C56]=[                    3.36,0.00]
    # [                    C66]=[                         3.36]
    
    # # Titanium
    # [C11,C12,C13,C14,C15,C16]=[1.23,1.00,0.69,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.23,0.69,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.53,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.61,0.00,0.00]
    # [                C55,C56]=[                    0.61,0.00]
    # [                    C66]=[                         0.24]
    
    # # # Tin
    # [C11,C12,C13,C14,C15,C16]=[0.75,0.62,0.44,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     0.75,0.44,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.53,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.44,0.00,0.00]
    # [                C55,C56]=[                    0.44,0.00]
    # [                    C66]=[                         0.47]
    
    
    # # Sulfur
    # [C11,C12,C13,C14,C15,C16]=[2.40,1.33,1.71,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     2.05,1.59,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          4.83,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.86,0.00,0.00]
    # [                C55,C56]=[                    1.74,0.00]
    # [                    C66]=[                         1.52]
 
    # # Feldspar
    # [C11,C12,C13,C14,C15,C16]=[0.62,0.43,0.37,0.00,-.14,0.00]
    # [    C22,C23,C24,C25,C26]=[     1.58,0.22,0.00,-.03,0.00]
    # [        C33,C34,C35,C36]=[          1.00,0.00,-.17,0.00]
    # [            C44,C45,C46]=[               0.28,0.00,-.05]
    # [                C55,C56]=[                    0.41,0.00]
    # [                    C66]=[                         0.72]
 
    # # Seignette salt
    # [C11,C12,C13,C14,C15,C16]=[2.55,1.41,1.16,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     3.81,1.46,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          3.71,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               2.68,0.00,0.00]
    # [                C55,C56]=[                    0.64,0.00]
    # [                    C66]=[                         1.96]
 
    # # Birkenholz (Triplay)
    # [C11,C12,C13,C14,C15,C16]=[0.36,0.11,0.17,0.00,0.00,0.00]
    # [    C22,C23,C24,C25,C26]=[     0.12,0.08,0.00,0.00,0.00]
    # [        C33,C34,C35,C36]=[          1.81,0.00,0.00,0.00]
    # [            C44,C45,C46]=[               0.30,0.00,0.00]
    # [                C55,C56]=[                    0.22,0.00]
    # [                    C66]=[                         0.13]
  

    # Wood moisture
    [E1,E2,E3]=[11.18,2.31,0.56]
    [G12,G13,G23]=[1.37,1.01,0.43]
    [n12,n21,n13,n31,n23,n32]=[0.01,0.04,0.11,2.21,0.26,1.09]
    S=np.array([[1/E1,-n12/E2,-n13/E3,0,0,0],[-n21/E1,1/E2,-n23/E3,0,0,0],\
                [-n31/E1,-n32/E2,1/E3,0,0,0],[0,0,0,1/G23,0,0],\
                [0,0,0,0,1/G13,0],[0,0,0,0,0,1/G12]])
    return S
    # # # Invent
    # [C11,C12,C13,C14,C15,C16]=[0.72,0.43,0.37,0.02,-.14,0.03]
    # [    C22,C23,C24,C25,C26]=[     1.08,0.22,0.03,-.03,0.01]
    # [        C33,C34,C35,C36]=[          1.20,0.01,-.17,0.00]
    # [            C44,C45,C46]=[               0.28,0.02,-.05]
    # [                C55,C56]=[                    0.41,0.02]
    # [                    C66]=[                         0.72]
    
    # C=np.array([[C11,C12,C13,C14,C15,C16],[C12,C22,C23,C24,C25,C26],\
    #             [C13,C12,C33,C34,C35,C36],[C14,C24,C34,C44,C45,C46],\
    #             [C15,C25,C35,C45,C55,C56],[C16,C26,C36,C46,C56,C66]])
    

    #return np.linalg.inv(C)


def components(nu_23,nu_13,nu_12,G_23, G_13, G_12, E_1,  E_2,  E_3):
    S=np.zeros((6,6))
    S[0][0]=1/E_1
    S[1][1]=1/E_2
    S[2][2]=1/E_3
    S[3][3]=1/G_23/2 #4
    S[4][4]=1/G_13/2 #4
    S[5][5]=1/G_12/2 #4
    S[0][1]=-nu_12/E_2 #2
    S[0][2]=-nu_13/E_3 #2
    S[1][2]=-nu_23/E_3 #2
    S[1][0]=S[0][1]
    S[2][0]=S[0][2]
    S[2][1]=S[1][2]
    return S

def generateMesh(nElev,nAzim):
    elev=np.linspace(0,np.pi,nElev,endpoint=True)
    # For 2D
    #elev=np.linspace(np.pi/2-0.01,np.pi/2+0.01,nElev,endpoint=True)
    azim=np.linspace(0,2*np.pi,nAzim,endpoint=True)
    angles=np.meshgrid(azim,elev)

    d_1=np.cos(angles[0])*np.sin(angles[1])
    d_2=np.sin(angles[0])*np.sin(angles[1])
    d_3=                  np.cos(angles[1])
    d1=[d_1,d_2,d_3]

    d_1=-np.sin(angles[0])
    d_2=np.cos(angles[0])
    d_3=                  0
    d2=[d_1,d_2,d_3]
    
    d_1=np.cos(angles[0])*np.sin(angles[1]+np.pi/2)
    d_2=np.sin(angles[0])*np.sin(angles[1]+np.pi/2)
    d_3=                  np.cos(angles[1]+np.pi/2)
    d3=[d_1,d_2,d_3]

    return d1,d2,d3

def generateFlatTensor(d,dt_1,dt_2):
    d0=d[0]*d[0]
    d1=d[1]*d[1]
    d2=d[2]*d[2]
    d3=d[1]*d[2]/np.sqrt(2)*2
    d4=d[2]*d[0]/np.sqrt(2)*2
    d5=d[0]*d[1]/np.sqrt(2)*2
    d_i=[d0,d1,d2,d3,d4,d5]
    
    dt1_0=dt_1[0]*dt_1[0]
    dt1_1=dt_1[1]*dt_1[1]
    dt1_2=dt_1[2]*dt_1[2]
    dt1_3=dt_1[1]*dt_1[2]/np.sqrt(2)*2
    dt1_4=dt_1[2]*dt_1[0]/np.sqrt(2)*2
    dt1_5=dt_1[0]*dt_1[1]/np.sqrt(2)*2
    dt_1i=[dt1_0,dt1_1,dt1_2,dt1_3,dt1_4,dt1_5]
    
    dt2_0=dt_2[0]*dt_2[0]
    dt2_1=dt_2[1]*dt_2[1]
    dt2_2=dt_2[2]*dt_2[2]
    dt2_3=dt_2[1]*dt_2[2]/np.sqrt(2)*2
    dt2_4=dt_2[2]*dt_2[0]/np.sqrt(2)*2
    dt2_5=dt_2[0]*dt_2[1]/np.sqrt(2)*2
    dt_2i=[dt2_0,dt2_1,dt2_2,dt2_3,dt2_4,dt2_5]
    
    return d_i,dt_1i,dt_2i

def plotSurface(i,j,X,Y,Z,mag): 
    # Hide grid lines
    #ax[i][j].grid(False)

    # # Hide axes ticks
    # ax[i][j].set_xticks([])
    ax.set_yticks([])
    # ax[i][j].set_zticks([])
    
    color_dimension = mag # change to desired fourth dimension
    minn, maxx = color_dimension.min(), color_dimension.max()
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    # plot
    #bound=np.amax(np.abs(mag))  
    bound=12
    # ax[i][j].plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
    # ax[i][j].set_xlim3d([-bound,bound])
    # ax[i][j].set_ylim3d([-bound,bound])
    # ax[i][j].set_zlim3d([-bound,bound])
    ax.plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
    ax.view_init(elev=0, azim=-90)
    # ax.set_xlim3d([0,bound])
    # ax.set_ylim3d([0,bound])
    # ax.set_zlim3d([0,bound])

    ax.set_xlim3d([-bound,bound])
    ax.set_ylim3d([-bound,bound])
    ax.set_zlim3d([-bound,bound])

###############################################################################

print(f" E: {round(  E_1,3)}\t,{round(  E_2,3)}\t,{round(  E_3,3)}\n")
print(f"nu: {round(nu_23,3)}\t,{round(nu_13,3)}\t,{round(nu_12,3)}\n")
print(f" G: {round( G_23,3)}\t,{round( G_13,3)}\t,{round( G_12,3)}\n")

#S=components(nu_23,nu_13,nu_12,G_23, G_13, G_12, E_1,  E_2,  E_3)
S=componentsC()
res=1 # degrees per line
nElev=round(180/res)
nAzim=round(360/res)
[d,dt_1,dt_2]=generateMesh(nElev,nAzim)
[d_i,dt_1i,dt_2i]=generateFlatTensor(d,dt_1,dt_2)

###############################################################################

E_inv=0                    
for i in range(6):
    for j in range(6):
        E_inv=E_inv+S[i][j]*d_i[i]*d_i[j]
        
E=1/E_inv
X = E * d[0]
Y = E * d[1]
Z = E * d[2]

plotSurface(0, 0,Y, Z,X, E)
ax[0][0].set_title("Elastic modulus")

# ###############################################################################

# K_inv=0    
# I=[1,1,1,0,0,0]                
# for i in range(6):
#     for j in range(6):
#         K_inv=K_inv+S[i][j]*I[i]*d_i[j]

# K=1/K_inv/3       
# X = K * d[0]
# Y = K * d[1]
# Z = K * d[2]

# plotSurface(0, 1, X, Y, Z, K)
# ax[0][1].set_title("Bulk modulus")

# ###############################################################################

# Et1_inv=0                    
# for i in range(6):
#     for j in range(6):
#         Et1_inv=Et1_inv+S[i][j]*d_i[i]*dt_1i[j]
        
# Et1=-Et1_inv
# X = Et1 * d[0]
# Y = Et1 * d[1]
# Z = Et1 * d[2]

# plotSurface(1, 0, X, Y, Z, Et1)
# ax[1][0].set_title("Poisson ratio t1")

# ###############################################################################

# Et2_inv=0                    
# for i in range(6):
#     for j in range(6):
#         Et2_inv=Et2_inv+S[i][j]*d_i[i]*dt_2i[j]
# Et2=-Et2_inv
# X = Et2 * d[0]
# Y = Et2 * d[1]
# Z = Et2 * d[2]

# plotSurface(1, 1, X, Y, Z, Et2)
# ax[1][1].set_title("Poisson ratio t2")

# ###############################################################################

#tri = mtri.Triangulation(angles[0].flatten(),angles[1].flatten())
# surf = ax.plot_trisurf(E_vec[:,0],E_vec[:,1],E_vec[:,2],triangles=tri.triangles,facecolor=plt.cm.coolwarm(E_mag))