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
print('Current assumprions:Flat earh')
## preperation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


## world variables
graf=float(9.81)


##functions
## functions for numeric diff.eq solving
def EF(Time,w,dt):
    ''' Euler forward numerical integration method
        w_(n+1)=w_n+dt*f(tn,wn)
        
        Computational easiest integration method and Analysticly worst
        
        input is a w matrix consisting of w=[w0,w1,w2,...,wn]
        with wi=[u1(i),u2(i),u3(i),...,um(i)]^T
        
        Returns: wn+1=[u1(n+1),u2(n+1),...,un(n+1)]^T'''
    
    ans=w[:,-1:]+dt*f(Time[-1],w[:,-1:]);
    return(ans)

def TZ(Time,w,dt):
    ''' Trapezodial numerical integration method
        w_(n+1)=w_n+dt*(f(t_n,w_n)+f(tn+dt,w*_(n+1)))/2
        with w*_(n+1)=w_n+dt*f(tn,wn)
        
        Averge on computuational time and averge on analytics
        
        input is a w matrix consisting of w=[w0,w1,w2,...,wn]
        with wi=[u1(i),u2(i),u3(i),...,um(i)]^T
        
        Returns: wn+1=[u1(n+1),u2(n+1),...,un(n+1)]^T'''
        
        
    
    ans= w[:,-1:]+dt/2*(f(Time[-1],w[:,-1:])+f(Time[-1]+dt,w[:,-1:]+dt*f(Time[-1],w[:,-1:])));
    
    return(ans)
    
def RK(Time,w,dt):
    ''' Runge-Kutta integrtion method
        w_(n+1)=w_n+1/6(k1+2k2+2k3+k4)
        
        with    k1=dt*f(tn,wn)
                k2=dt*f(tn+dt/2,w_n+k1/2) 
                k3=dt*f(tn+dt/2,w_n+k2/2)
                k4=dt*f(t_n+dt,w_n+k3)
                
        Computationaly the hardest method but analystically the best.
        
        input is a w matrix consisting of w=[w0,w1,w2,...,wn]
        with wi=[u1(i),u2(i),u3(i),...,um(i)]^T
        
        Returns: wn+1=[u1(n+1),u2(n+1),...,un(n+1)]^T'''
    
    k1=dt*f(Time[-1],w[:,-1:])
    
    k2=dt*f(Time[-1]+dt/2,w[:,-1:]+k1/2)
    
    k3=dt*f(Time[-1]+dt/2,w[:,-1:]+k2/2)
    
    k4=dt*f(Time[-1]+dt,w[:,-1:]+k3)
    
    ans=w[:,-1:]+1/6*(k1+2*k2+2*k3+k4)
    
    return(ans)
    
## diff.eq  
    
    
    ## the differential eq      
def f(t,w):
    ''' The differential equation that has to be solved.
        imputs: t is time array
                w matrix consisting of w=[w0,w1,w2,...,wn]
                with wi=[u1(i),u2(i),u3(i),...,um(i)]^T
                
        Returns:The next numerical approximation of the solution.
                wn+1=[u1(n+1),u2(n+1),...,un(n+1)]^T'''
    

    #Linear part of the differntial equation matrix 
    I=np.array([[5,0,0],
                [0,5,0],
                [0,0,1]])
    
    
    A=np.block([
                [np.zeros((3,3)),np.identity(3)],
                [np.zeros((3,3)),np.zeros((3,3))]
                ])
    A=np.block([
                [A,np.zeros((6,6))],
                [np.zeros((6,6)),A]
                ])
    
    # this is purly for comsetical reasons.
    #extracting a place matrix
    x=w[0:3]
    
    #exctracting a speed matrix
    x_dot=w[3:6]

    #exctrating a place matrix
    o=w[6:9]

    #exctracting a speed matrix
    o_dot=w[9:12]
    
    
    # this is experimental code
    cp=1
    #rcpcp is the vector form the centere of mass to the centre of pressure.
    #https://www.grc.nasa.gov/www/k-12/rocket/rktcp.html
    rcmcp=np.transpose(
            np.array(
                    [[np.cos(np.pi/2-o[2,-1])*np.sin(o[1,-1])*cp],
                    [np.cos(np.pi/2-o[2,-1])*np.cos(o[1,-1])*cp],
                    [np.cos(np.pi/2-o[2,-1])*cp]]))
    

    
    #the non linear part of the differential equation.
    ans=A.dot(w)+Rot_Internal(o_dot,I)+Rot_External(I,x_dot,rcmcp)+Lat_Fg()+Lat_Fwl(x_dot)

    return(ans)    




## Forces
def Rot_Internal(o_dot,I):

    
    #linear part of the rotional acceleration
    Internal=-1*np.linalg.inv(I).dot(np.transpose(np.cross(np.transpose(I.dot(o_dot)),np.transpose(o_dot))))

    ans=np.block([
                [np.zeros((9,1))],
                [Internal]
                ])
    return(ans)
    
def Rot_External(I,x_dot,rcmcp):


    M=np.transpose(np.cross(np.transpose(Lat_Fwl(x_dot)[3:6]),rcmcp))
    
    ans=np.block([
                [np.zeros((9,1))],
                [np.linalg.inv(I).dot(M)]
                ])
    return(ans)

def Lat_Fg():
    #Grafity
    ans=np.zeros([12,1])
    ans[5]=-graf
    return(ans)
    
def Lat_Fwl(v):
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

## Inital condition
# first time step
t0=0
Time=[t0]
#last time step
tE=20
#step time
dt=0.1

# inital conditions
# This has no coding purpose. 
# This is purly for clarivication

#initial place vector
x0=np.array([[0],
             [0],
             [80]])
    
#intial speed vector
x0_dot=np.array([[0],
                 [0],
                 [40]])

#inital orientation vector
o0=np.array([[0],
             [0],
             [1/2*np.pi]])  

#inital angular velocity vector
o0_dot=np.array([[0],
                 [0],
                 [0]]) 
    

#creating the initial w matrix
# with the earlier defined inital conditions
w0=np.block([
                [x0],
                [x0_dot],
                [o0],
                [o0_dot]
                ])
w=w0



##----------- The excecution --------------
while(Time[-1]<tE):
    
    #appling the numerical method over the differnential eqation f.
    
    wn=RK(Time,w,dt)
    
    #adding the next point to the numerical matrix w
    
    w=np.append(w,wn,axis=1)
    

    #going to the next time step.
    
    Time.append(Time[-1]+dt)
    
    #reitterating everything until tE
    
    #unless the z coordinate has become 0.
    if(w[2,-1]<0):
        break
    

    
    
##------------- The plots -----------------
    
#Creating these matrices is purley cosmetical.
    
#creating a place matrix
x=w[0:3]
#creating a speed matrix
x_dot=w[3:6]

#creating a place matrix
o=w[6:9]

#creating a speed matrix
o_dot=w[9:12]


#experimental code
cp=1

Rcmcp=np.array([[np.cos(np.pi/2-o[2])*np.sin(o[1])*cp],
                [np.cos(np.pi/2-o[2])*np.cos(o[1])*cp],
                [np.cos(np.pi/2-o[2])*cp]])


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
plt.legend(('p','q','r'))
plt.show()

fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
ax.scatter(x[0,:], x[1,:], x[2,:], c='r', marker='o')
ax.scatter(x[0,:]+Rcmcp[0,:], x[1,:]+Rcmcp[1,:], x[2,:]+Rcmcp[2,:], c='b', marker='*')

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


print('------ end of code ------')    


