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
    
def RK(Time,w,dt):
    ### Runge-Kutta integrtion method
    
    k1=dt*f(Time[-1],w[:,-1:])
    
    k2=dt*f(Time[-1]+dt/2,w[:,-1:]+k1/2)
    
    k3=dt*f(Time[-1]+dt/2,w[:,-1:]+k2/2)
    
    k4=dt*f(Time[-1]+dt,w[:,-1:]+k3)
    
    ans=w[:,-1:]+1/6*(k1+2*k2+2*k3+k4)
    
    return(ans)
    
## diff.eq  
    
    
    ## the differential eq      
def f(t,w):
    
    #creating a place matrix
    x=w[0:3]
    #creating a speed matrix
    x_dot=w[3:6]

    #creating a place matrix
    o=w[6:9]

    #creating a speed matrix
    o_dot=w[9:12]


    
    #interial matrix 
    I=np.array([[3,0,0],
                [0,1,0],
                [0,0,1]])
    
    
    A=np.block([
                [np.zeros((3,3)),np.identity(3)],
                [np.zeros((3,3)),np.zeros((3,3))]
                ])
    A=np.block([
                [A,np.zeros((6,6))],
                [np.zeros((6,6)),A]
                ])
    
    
    #the non homogenus part
    ans=A.dot(w)+Internal(o_dot,I)+External(I)+Fg()+Fwl(x_dot)

    
    return(ans)      




## Forces
def Internal(o_dot,I):

    
    #linear part of the rotional acceleration
    Internal=-1*np.linalg.inv(I).dot(np.transpose(np.cross(np.transpose(I.dot(o_dot)),np.transpose(o_dot))))

    ans=np.block([
                [np.zeros((9,1))],
                [Internal]
                ])
    return(ans)
    
def External(I):

    M=np.array([[0],
                [0],
                [0]])
    
    ans=np.block([
                [np.zeros((9,1))],
                [np.linalg.inv(I).dot(M)]
                ])
    return(ans)

def Fg():
    #Grafity
    ans=np.zeros([12,1])
    ans[5]=-graf
    return(ans)
    
def Fwl(v):
    #air resitance
    Cw=float(0.75)
    Area=float(0.04)
    Rho=float(1.225)
    ans=-Cw*Area*Rho*v[:,-1:]*np.linalg.norm(v[:,-1:])
    ans = np.block([
                [np.zeros((3,1))],
                [ans],
                [np.zeros((6,1))]
                ])
    return(ans)

## variables
t0=0
Time=[t0]
tE=20
dt=0.1
w0=np.array([[0],
             [0],
             [100],
             [0],
             [0],
             [10],
             [-1],
             [0],
             [1],
             [0.5],
             [1],
             [-1]])
w=w0

## The excecution
while(Time[-1]<tE):
    #solution
    wn=RK(Time,w,dt)
    w=np.append(w,wn,axis=1)

    #time step
    Time.append(Time[-1]+dt)
    

    
    
## plots
    #creating a place matrix
    x=w[0:3]
    #creating a speed matrix
    x_dot=w[3:6]

    #creating a place matrix
    o=w[6:9]

    #creating a speed matrix
    o_dot=w[9:12]

#for i in range(len(Time)):
#    Forces.append(FSUM(w[:,i:]))


#fig1 = plt.figure()
#ax = fig1.add_subplot(111, projection='3d')
#ax.scatter(o[0,:], o[1,:], o[2,:], c='r', marker='o')
#fig2 = plt.figure()
#ax = fig2.add_subplot(111, projection='3d')
#ax.scatter(Time, o[1,:], o[2,:], c='r', marker='o')
fig1 = plt.figure()
plt.subplot(211)
plt.plot(Time,o_dot[0,:],'-',color='black')
plt.plot(Time,o_dot[1,:],'-',color='blue')
plt.plot(Time,o_dot[2,:],'-',color='red')
plt.xlabel('Time(s)')
plt.ylabel('z')
plt.title('omega')
plt.legend(('p','q','r'))

plt.subplot(212)
o=np.mod(o,2*np.pi)
plt.plot(Time,o[0,:],'-',color='black')
plt.plot(Time,o[1,:],'-',color='blue')
plt.plot(Time,o[2,:],'-',color='red')

plt.title('orientation')
plt.legend(('theta','phi','xi'))
plt.show()

fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
ax.scatter(x[0,:], x[1,:], x[2,:], c='r', marker='o')

fig3 = plt.figure()
ax = fig3.add_subplot(121)
plt.plot(Time,x[2,:],'-')
plt.xlabel('Time(s)')
plt.ylabel('z')
plt.title('z')
ax = fig3.add_subplot(122)
plt.plot(Time,x_dot[2,:],'-')
plt.title('z_dot')
plt.ylabel('z_dot(m/s)')
plt.xlabel('Time(s)')
#plt.subplot(313)
#plt.plot(Time[0:-1],Forces,'-')
#plt.show
print('------ end of code ------')    


