#Libraries---------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
def main():
    #Functions---------------------------------------------------------------------
    #Animation
    def anim(k):
        fig=plt.figure(k)
        x_values=[i*10 for i in np.arange(-50,60,10)]
        ax=fig.add_subplot(111,projection="3d")
        plt.setp(ax,xticks=np.arange(0,110,10),xticklabels=x_values)
        plt.setp(ax,yticks=np.arange(0,110,10),yticklabels=x_values)
        ax.set_ylabel("Distance (km)")
        ax.set_xlabel("Distance (km)")

        x=np.arange(1,98,1)
        x,y =np.meshgrid(x,x)
        h_surface=[ax.plot_surface(y,y,h_dict[0],cmap=cm.coolwarm,rstride=1,cstride=1)]
        cmap=plt.colorbar(h_surface[0],label="Height (m)")

        def update_surface(frameNum,h_dict,h_surface):
            h_surface[0].remove()
            h_surface[0]=ax.plot_surface(x,y,h_dict[frameNum],cmap=cm.coolwarm)

        ani=animation.FuncAnimation(fig,update_surface,frames=range(0,timesteps,1),fargs=(h_dict,h_surface))
        #ani.save("water_tank.mp4",writer="ffmpeg",fps=15)
        plt.show()
    #1. Set model parameters-------------------------------------------------------
    dx=10000
    dy=10000
    dt=50
    timesteps=250
    idim=99
    h_0=0.1
    H=100
    sigma=50e3
    f=0
    g=9.81
    x_0=idim/2
    y_0=idim/2

    h_dict={}
    u_dict={}
    v_dict={}
    #2. Initialize model variables-------------------------------------------------
    h=np.zeros((idim,idim))
    u=np.zeros((idim,idim))
    v=np.zeros((idim,idim))

    htend=np.zeros((idim,idim))
    hold=np.zeros((idim,idim))
    hnew=np.zeros((idim,idim))

    utend=np.zeros((idim,idim))
    uold=np.zeros((idim,idim))
    unew=np.zeros((idim,idim))

    vtend=np.zeros((idim,idim))
    vold=np.zeros((idim,idim))
    vnew=np.zeros((idim,idim))
    #Initialize heightfield
    for x in range(0,idim):
        for y in range(0,idim):
            h[x][y]=H+h_0*np.exp(-((((x-x_0)*dx)**2)/(2*sigma**2)+(((y-y_0)*dy)**2)/(2*sigma**2)))

    h_dict[0]=h[1:98,1:98].copy()
    u_dict[0]=u[1:98,1:98].copy()
    v_dict[0]=v[1:98,1:98].copy()

    for k in range(1,timesteps+1):
        #3. Compute tendencies of state variables----------------------------------
        for i in range(1,idim-1):
            for j in range(1,idim-1):
                if i==idim-2:
                    if j==idim-2:
                        htend[i][j]=-u[i][j]*(h[1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][1]-u[i][j-1])/(2*dy)\
                                    -g*(h[1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][1]-h[i][j-1])/(2*dy)
                    elif j==1:
                        htend[i][j]=-u[i][j]*(h[1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][idim-2])/(2*dy)\
                                    -h[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][idim-2])/(2*dy)\
                                    -g*(h[1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][idim-2])/(2*dy)
                    else:
                        htend[i][j]=-u[i][j]*(h[1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][j-1])/(2*dy)\
                                    -g*(h[1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][j-1])/(2*dy)
                if i==1:
                    if j==idim-2:
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[idim-2][j])/(2*dx)\
                                    -v[i][j]*(h[i][1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -h[i][j]*(v[i][1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -v[i][j]*(u[i][1]-u[i][j-1])/(2*dy)\
                                    -g*(h[i+1][j]-h[idim-2][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[idim-2][j])/(2*dx)\
                                    -v[i][j]*(v[i][1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][1]-h[i][j-1])/(2*dy)
                    elif j==1:
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[idim-2][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][idim-2])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][idim-2])/(2*dy)\
                                    -g*(h[i+1][j]-h[idim-2][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[idim-2][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][idim-2])/(2*dy)
                    else:
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[idim-2][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[idim-2][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][j-1])/(2*dy)\
                                    -g*(h[i+1][j]-h[idim-2][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[idim-2][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][j-1])/(2*dy)
                else:
                    if (j==idim-2):
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][1]-u[i][j-1])/(2*dy)\
                                    -g*(h[i+1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][1]-h[i][j-1])/(2*dy)
                    if (j==1):
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][idim-2])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][idim-2])/(2*dy)\
                                    -g*(h[i+1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][idim-2])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][idim-2])/(2*dy)
                    else:
                        htend[i][j]=-u[i][j]*(h[i+1][j]-h[i-1][j])/(2*dx)\
                                    -v[i][j]*(h[i][j+1]-h[i][j-1])/(2*dy)\
                                    -h[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -h[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)
                        utend[i][j]=-u[i][j]*(u[i+1][j]-u[i-1][j])/(2*dx)\
                                    -v[i][j]*(u[i][j+1]-u[i][j-1])/(2*dy)\
                                    -g*(h[i+1][j]-h[i-1][j])/(2*dx)
                        vtend[i][j]=-u[i][j]*(v[i+1][j]-v[i-1][j])/(2*dx)\
                                    -v[i][j]*(v[i][j+1]-v[i][j-1])/(2*dy)\
                                    -g*(h[i][j+1]-h[i][j-1])/(2*dy)
        #4. Extrapolate in time----------------------------------------------------
        if k==1:
            for i in range(1,idim-1):
                for j in range(1,idim-1):
                    hnew[i][j]=h[i][j]+dt*htend[i][j]
                    unew[i][j]=u[i][j]+dt*utend[i][j]
                    vnew[i][j]=v[i][j]+dt*vtend[i][j]
        else:
            for i in range(1,idim-1):
                for j in range(1,idim-1):
                    hnew[i][j]=hold[i][j]+2*dt*htend[i][j]
                    unew[i][j]=uold[i][j]+2*dt*utend[i][j]
                    vnew[i][j]=vold[i][j]+2*dt*vtend[i][j]
        for i in range(1,idim-1):
            for j in range(1,idim-1):
                hold[i][j]=h[i][j]
                h[i][j]=hnew[i][j]
                uold[i][j]=u[i][j]
                u[i][j]=unew[i][j]
                vold[i][j]=v[i][j]
                v[i][j]=vnew[i][j]
        #5. Update boundary values-------------------------------------------------

        #6. Store computed values--------------------------------------------------
        h_dict[k]=h[1:98,1:98].copy()
        u_dict[k]=u[1:98,1:98].copy()
        v_dict[k]=v[1:98,1:98].copy()
        #7. Ready?-----------------------------------------------------------------
    anim(4)

if __name__ == "__main__":
    main()
