#
# module for doing tests of quad
# 
import numpy as np
from numpy.fft import rfft, irfft
import scipy.linalg as la

class QuadShit(object):
    def __init__ (self, N=1024, Nb=5, A2=1.0, mask=None):
        self.N=N    # Number of pixels in skewer
        self.Nfft=self.N/2+1
        self.Nb=Nb   # number of lowk bins
        self.A2=A2  # amplitude at z=2
        self.ramp=np.linspace(1.0,A2,self.N)
        if mask is None:
            self.mask=np.ones(N)
        else:
            self.mask=mask
        self.CI=1./(self.ramp**2) * self.mask
        
    def genskew (self):
        return np.random.normal(0.0, 1.0, self.N)*self.ramp

    def genskew_masked(self):
        return self.genskew()*self.mask
    
    def getPkIdiotic (self, Ng=1000):
        Pks=np.array([abs(rfft(self.genskew())**2)/self.N for i in range(Ng)])
        Pk=Pks.mean(axis=0)
        Pke=np.sqrt(Pks.var(axis=0)/Ng)
        ks=2*np.pi/1.0*np.arange(len(Pk))
        return ks,Pk,Pke


    def getSponez(self):
    ## return S prime, derivative in corr func for one redshift accross space
        Sp=[]
        xil=[]
        print "Creating Sp"
        N=self.N
        St=np.zeros((N,N))
        for i in range(self.Nb):
            # Andreu: I believe irfft expects as an input an array of N real numbers, not Nfft complex numbers
            # Although the code doesn't complain, so probably you can also do this...
            ar=np.zeros(self.Nfft, np.complex)
            if (i<self.Nb-1):
                kup=(i+1)**2
            else:
                kup=self.N/2
            kdown=i**2
            xip=(np.sin(kup*2*np.pi*np.arange(N)/N)-np.sin(kdown*2*np.pi*np.arange(N)/N))/(2*np.pi*np.arange(N)/N)
            xip[0]=(kup-kdown)
            xip/=(N/2.)
            xil.append(xip)
            M=np.zeros((N,N))
            for j in range(N):
                if j<N-1:
                    M[j,j:]=xip[:N-j]
                    M[j,:j]=xip[j:0:-1]
                else:
                    M[j,:]=xip[::-1]
            Sp.append(M)
            St+=M
        #Sp.append(np.diag(np.ones(N))-St)
        #xil.append(Sp[-1][0,:])
        return Sp, xil

    def getSp3(self):
    ### return Sprime, 3 redshift bins, 
        Sp1,xil=self.getSponez()
        ## We demand that Spl, Spm and Sph add to Sp1
        N=self.N
        ## we make window at twice resolution
        wl,wh=np.zeros(2*N),np.zeros(2*N)
        cosa=0.5*np.cos(2*np.pi*np.arange(2*N)/float(2*N))+0.5
        wl[:N]=cosa[:N]
        wh[N:]=cosa[N:]
        wm=(1-wl-wh)
        ## X below is the sum of indices, i.e. X[i,j]=i+j
        X=np.outer(np.ones(N,int),np.arange(N))+np.outer(np.arange(N),np.ones(N,int))
        ## the hocus pocus below make Wl[i,j]=wl[i+j]
        Xf=X.flatten()
        Wl=wl[Xf].reshape((N,N))
        Wm=wm[Xf].reshape((N,N))
        Wh=wh[Xf].reshape((N,N))
        Sp3=[]
        for W in [Wl,Wm,Wh]:
            for S in Sp1:
                Sp3.append(W*S)
        return Sp3,xil
    
        
    
    def getOQE(self, Sp, Ng=10):
        ## get OQE given some derivatives
        print "Getting FD"
        Nb=len(Sp)
        FD=np.zeros(Nb)
        for i in range(Ng):
            cid=self.CI*self.genskew_masked()
            ## FD = 1/2 * (Cinv d)^T S, (Cinv d)
            FD+=0.5*np.array([np.dot(np.dot(Sp[j],cid),cid) for j in range(Nb)])
        print "Getting SD"
        SD=np.zeros((Nb,Nb))
        ## this is an idiotic way, since we know CI is diagonal, but let's keep it simple
        ## and sure
        CI=np.diag(self.CI)
        for i in range(Nb):
            for j in range(i,Nb):
                # SD = 0.5 Tr(S, Cinv S, Cinv)
                #SD[i,j]=0.5*np.trace(np.dot(CI,np.dot(Sp[i],np.dot(CI,Sp[j]))))
                SD[i,j]=0.5*np.trace(np.dot(np.dot(CI,Sp[i]),np.dot(CI,Sp[j])))
                #SD[i,j]=0.5*np.trace(np.dot(np.dot(CI,Sp[i]),np.dot(CI,Sp[j]).T))
                SD[j,i]=SD[i,j]
        print SD[1,1]
        ## Ng is # of goes
        SD*=Ng
        return FD, SD

    def getPkfromD(self, FD, SD):
        SDI=la.inv(SD)
        Pk=np.dot(SDI,FD)
        err=np.sqrt(SDI.diagonal())
        ks=np.arange(len(Pk))
        chi2=np.dot(np.dot(SD,(Pk-1.0)),(Pk-1.0))
        chi2d=((Pk-1.0)**2/(err**2)).sum()
        print "chi2=",chi2, "chi2 diag=",chi2d
        return ks,Pk,err
                                   
        
        
        
    
