from numpy import exp, cos, linspace
import matplotlib.pyplot as plt
import os, time, glob
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np

err_val=0.0001 #maximum relative error
maxiter=500 #maximum number of iteration 

#helper function to return the relative error
def relerror(old,new):
    i=np.where(abs(new) > 0)
    err=np.mean(((new[i]-old[i])/new[i])**2.)
    return err

#helper function to plot a 3d graph
def myplot3d(X,Y,data):
    fig = plt.figure(figsize=(5,5), dpi=100)
    ax = fig.gca(projection='3d')                      
    surf = ax.plot_wireframe(X,Y,data[:])
    
def solvePoisson_Jacobi(p,f,deltax):
    iter=0
    
    while True:
        pn=p.copy()
        iter=iter+1
        p[:]=0
        #boundary conditions
        pn[0,:] = 0
        pn[-1,:] = 0
        pn[:,0] = 0
        pn[:,-1] = 0
        #finite difference scheme
        p[1:-1,1:-1]=0.25*(deltax*deltax*f[1:-1,1:-1]+(pn[1:-1,2:]+pn[1:-1,0:-2]+pn[2:,1:-1]+pn[0:-2,1:-1]))
        #stop if relative error is below bound
        if relerror(pn,p)<err_val:
            break
        if iter>maxiter:
            break
    #print iter
    #print relerror(pn,p)
    return p

def solvePoisson_GaussSeidel(p,f,deltax):
    iter=0
    while True:
        pn=p.copy()
        p[:]=0.
        iter=iter+1
        
        #boundary conditions
        pn[0,:] = 0
        pn[-1,:] = 0
        pn[:,0] = 0
        pn[:,-1] = 0
        #finite difference scheme
        for i in range(1,p.shape[0]-1): 
            for j in range (1,p.shape[1]-1):
                p[i,j]=0.25*(deltax*deltax*f[i,j]+pn[i,j+1]+p[i,j-1]+pn[i+1,j]+p[i-1,j])
        #stop if relative error is below bound
        if relerror(pn,p)<err_val:
            break
        if iter>maxiter:
            break
    #print iter
    #print relerror(pn,p)
    return p

#@autojit
def solvePoisson_SOR(p,f,deltax):
    iter=0
    beta=1.12
    while True:
        pn=p.copy()
        iter=iter+1
        
        #boundary conditions
        pn[0,:] = 0
        pn[-1,:] = 0
        pn[:,0] = 0
        pn[:,-1] = 0
        #finite difference scheme
        for i in range(1,p.shape[0]-1): 
            for j in range (1,p.shape[1]-1):
                p[i,j]=(1-beta)*pn[i,j]+0.25*beta*(deltax*deltax*f[i,j]+pn[i,j+1]+p[i,j-1]+pn[i+1,j]+p[i-1,j])

                #stop if relative error is below bound
        if relerror(pn,p)<err_val:
            break
        if iter>maxiter:
            break
    #print iter
    #print relerror(pn,p)
    return p
def compute_CFD(maxstep=4000, n=71, L=1, Uwall=2.0, nu=0.05, deltat = 0.0003):
#    maxstep=4000 #number of steps
#    n=71 #grid
#    L=1.#length 
    deltax=L/(n-1)
    ps = np.zeros((n,n)) #Psi
    om = np.zeros((n,n)) #omega
    x = y = np.linspace(0, L, n)
#    Uwall=2.0 #velocity of the wall
#    nu=0.05 #kinematic viscosity
#    deltat=0.0003 #delta t

    CFL1=2.*nu*deltat/deltax/deltax #should be smaller than 1/2
    print ("CFL1={0}".format(CFL1))
    CFL2=Uwall*deltat/deltax #should be smaller than 1/2
    print ("CFL2={0}".format(CFL2))
    Re=L*Uwall/nu
    print ("Re={0}".format(Re))
    for tstep in range(maxstep):
        #Poisson Solver
        ps=solvePoisson_Jacobi(ps,om,deltax)
        #Boundary conditions for the vorticity
        omc=om.copy()
        omc[1:-1,0]=-2.0*ps[1:-1,1]/deltax/deltax #bottom wall
        omc[1:-1,-1]=-2.0*ps[1:-1,-2]-2.0*Uwall/deltax #top wall
        omc[0,1:-1]=-2.0*ps[1,1:-1]/deltax/deltax #left wall
        omc[-1,1:-1]=-2*ps[-2,1:-1]/deltax/deltax #right wall
        
        #Finite difference scheme for the vorticity transport equation: omega^n_ij=omega^n_ij+A+B
        A=-deltat/4./deltax/deltax*((ps[1:-1,2:]-ps[1:-1,0:-2])*(omc[2:,1:-1]-omc[0:-2,1:-1])+\
                                (ps[0:-2,1:-1]-ps[2:,1:-1])*(omc[1:-1,2:]-omc[1:-1,0:-2]))
        B=nu*deltat/deltax/deltax*(omc[2:,1:-1]+omc[1:-1,2:]-4*omc[1:-1,1:-1]+omc[0:-2,1:-1]+omc[1:-1,0:-2])
        om[1:-1,1:-1]=omc[1:-1,1:-1]+A+B

    #Now calculate the velocity vectors    
    u1=(ps[1:-1,2:]-ps[1:-1,0:-2])/2./deltax
    u2=-(ps[2:,1:-1]-ps[0:-2,1:-1])/2./deltax
    X, Y = np.meshgrid(x,y) #needed for plotting
    nn=3 #reduction of no of points for arrow plot (quiver) 
    fig = plt.figure(figsize=(12,12), dpi=200) #figure size
    ax1 = fig.add_subplot(1, 2, 1,aspect='equal') #two figures side-by-side
    ax2 = fig.add_subplot(1, 2, 2,aspect='equal')
    ax1.quiver(X.T[1:-1:nn,1:-1:nn],Y.T[1:-1:nn,1:-1:nn],u1[::nn,::nn],u2[::nn,::nn])#arrow plot
    ax1.imshow(om.T, interpolation='nearest',extent=[0,L,0,L],origin='lower')#color showing abs velocity
    uabs=np.sqrt(u1**2+u2**2)
    ax2.streamplot(X[1:-1,1:-1],Y[1:-1,1:-1],u1.T,u2.T,color=uabs.T, linewidth=2, cmap=plt.cm.autumn)#streamline plot

    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        # Remove old plot files
        for filename in glob.glob(os.path.join('static', '*.png')):
            os.remove(filename)
    # Use time since Jan 1, 1970 in filename in order make
    # a unique filename that the browser has not chached
    plotfile = os.path.join('static', str(time.time()) + '.png')
    plt.savefig(plotfile)
    return plotfile