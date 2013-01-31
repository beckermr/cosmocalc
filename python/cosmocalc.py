import numpy as np
import copy
import scipy.interpolate
import scipy.integrate 
from numpy import linalg as la

class w0waCosmo:
    def __init__(self):
        self.om = 0.0
        self.ol = 0.0
        self.ob = 0.0
        self.onu = 0.0
        self.ok = 0.0
        self.h = 0.0
        self.s8 = 0.0
        self.ns = 0.0
        self.w0 = 0.0
        self.wa = 0.0
    
    def weff(self,a):
        if a != 1.0:
            return self.w0 + self.wa - self.wa*(a - 1.0)/np.log(a)
        else:
            return self.w0
    
    def hubble(self,a):
        return np.sqrt(self.om/a/a/a + self.ok/a/a + self.ol*np.exp(3.0*(1.0 + self.weff(a))))
    
class FLRWDistances:
    def __init__(self,cd):
        self.DH = 2997.92458
        self.cd = copy.deepcopy(cd)
        self.amin = 0.001
        self.amax = 1.0
        self.Na = 100
        if self.cd.ok != 0.0:
            self.DH_sqrtok = self.DH/np.sqrt(np.abs(self.cd.ok))
        else:
            self.DH_sqrtok = 0.0
        self.atab = np.arange(self.amin,self.amax,(self.amax-self.amin)/(self.Na-1.0),dtype='f8')
        self.cdtab = self.comvdist_exact(self.atab)
        self.cdinterp = scipy.interpolate.interp1d(self.atab,self.cdtab,kind='cubic',copy=False)
        self.ratab = self.atab[::-1]
        self.rcdtab = self.cdtab[::-1]
        self.acdinterp = scipy.interpolate.interp1d(self.rcdtab,self.ratab,kind='cubic',copy=False)
        
    def comvdist_int(self,a):
        return 1.0/a/a/self.cd.hubble(a)
    
    def comvdist_exact(self,avec):
        cd = []
        for a in avec:
            cd.append(scipy.integrate.quad(self.comvdist_int,a,1.0)[0]*self.DH)
        return cd
    
    def comvdist(self,a):
        return self.cdinterp(a)
    
    def acomvdist(self,cd):
        return self.acdinterp(cd)
    
    def angdist(self,a):
        if self.cd.ok > 0.0:
            return self.DH_sqrtok*np.sinh(self.comvdist(a)/self.DH_sqrtok)*a
        elif self.cd.ok < 0.0:
            return self.DH_sqrtok*np.sin(self.comvdist(a)/self.DH_sqrtok)*a
        else:
            return self.comvdist(a)*a
    
    def lumdist(self,a):
        if self.cd.ok > 0.0:
            return self.DH_sqrtok*np.sinh(self.comvdist(a)/self.DH_sqrtok)/a
        elif self.cd.ok < 0.0:
            return self.DH_sqrtok*np.sin(self.comvdist(a)/self.DH_sqrtok)/a
        else:
            return self.comvdist(a)/a

    def angdistdiff(self,a1,a2):
        if a1 < a2:
            if self.cd.ok > 0.0:
                return self.DH_sqrtok*np.sinh((self.comvdist(a1)-self.comvdist(a2))/self.DH_sqrtok)/a1;
            elif self.cd.ok < 0:
                return self.DH_sqrtok*np.sin((self.comvdist(a1)-self.comvdist(a2))/self.DH_sqrtok)/a1;
            else:
                return (self.comvdist(a1)-self.comvdist(a2))*a1;
        else:
            if self.cd.ok > 0.0:
                return self.DH_sqrtok*np.sinh((self.comvdist(a2)-self.comvdist(a1))/self.DH_sqrtok)/a2;
            elif self.cd.ok < 0:
                return self.DH_sqrtok*np.sin((self.comvdist(a2)-self.comvdist(a1))/self.DH_sqrtok)/a2;
            else:
                return (self.comvdist(a2)-self.comvdist(a1))*a2;
