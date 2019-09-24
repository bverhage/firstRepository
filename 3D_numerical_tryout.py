# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:58:36 2019

@author: billy
"""

# -*- coding: utf-8 -*-
"""
This is a numeric differential equation solver
"""
print('------ begin of code ------')
## preperation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

## world variables
graf=float(9.81)
## list of the active forces
Forces=[]

##functions
## functions for numeric diff.eq solving
def EF(Time,w,dt):
    ### Euler forward numerical integration method
    ### w_(n+1)=w_n+dt*f(tn,wn)
    ans=w[:,-1:]+dt*f(Time[-1],w[:,-1:]);
    return(ans)

def TZ(Time,w,dt):
    ### trapezodial numerical integration method
    ### w_(n+1)=w_n+dt*(f(t_n,w_n)+f(tn+dt,w*_(n+1)))/2
    ### with w*_(n+1)=w_n+dt*f(tn,wn)
    
    ans= w[:,-1:]+dt/2*(f(Time[-1],w[:,-1:])+f(Time[-1]+dt,w[:,-1:]+dt*f(Time[-1],w[:,-1:])));
    
    return(ans)
    
## diff.eq  
    
    
    ## the differential eq      
def f(t,w):
    #creating a place matrix
    x=np.delete(w,5,axis=0)
    x=np.delete(x,3,axis=0)
    x=np.delete(x,1,axis=0)
    #creating a speed matrix
    v=np.delete(w,4,axis=0)
    v=np.delete(v,2,axis=0)
    v=np.delete(v,0,axis=0)
    
    #the homogeneus part
    A=np.array([[0,1,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,0,0],
                [0,0,0,0,0,1],
                [0,0,0,0,0,0]])
    #the non homogenus part
    ans=A.dot(w)+FSUM(x,v)
    return(ans)      




## Forces
def Fg():
    #Grafity
    ans=np.zeros([6,1])
    ans[5]=-graf
    return(ans)
    
def Fwl(v):
    #air resitance
    Cw=float(0.75)
    Area=float(0.04)
    Rho=float(1.225)
    ans=-Cw*Area*Rho*v[:,-1:]*np.linalg.norm(v[:,-1:])
    ans = np.array([[0],
           ans[0],
           [0],
           ans[1],
           [0],
           ans[2]])
    return(ans)

def FSUM(x,v):
    ###non homogenius time independent forces
    ans=Fg()+Fwl(v)
    
    ##gives problesm with TZ rule
    Forces.append(ans)
    return(ans)

## variables
t0=0
Time=[t0]
tE=10
dt=0.5
w0=np.array([[0],
             [10],
             [0],
             [10],
             [100],
             [10]])
w=w0

## The excecution
while(Time[-1]<tE):
    #solution
    wn=EF(Time,w,dt)
    w=np.append(w,wn,axis=1)

    #time step
    Time.append(Time[-1]+dt)
    
    ##break if ground is reached 
    #if(w[4,-1]<0):break
    
    
## plots
    #creating a place matrix
    x=np.delete(w,5,axis=0)
    x=np.delete(x,3,axis=0)
    x=np.delete(x,1,axis=0)
    #creating a speed matrix
    v=np.delete(w,4,axis=0)
    v=np.delete(v,2,axis=0)
    v=np.delete(v,0,axis=0)

#for i in range(len(Time)):
#    Forces.append(FSUM(w[:,i:]))


fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax.scatter(x[0,:], x[1,:], x[2,:], c='r', marker='o')
fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
ax.scatter(Time, x[1,:], x[2,:], c='r', marker='o')
fig2 = plt.figure()
plt.plot(Time,x[2,:],'-')
#plt.subplot(312)
#plt.plot(Time,v[2,:],'-')
#plt.subplot(313)
#plt.plot(Time[0:-1],Forces,'-')
#plt.show
print('------ end of code ------')    


