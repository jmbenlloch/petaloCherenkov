"""
Scintillator
Abstract class for scintillators
"""
from abc import ABCMeta, abstractmethod
from Centella.physical_constants import *
import sys
from Util import *
from Xenon import *


nm = nanometer 
mol = mole
micron = micrometer

LXeRefractionIndex =[
[6.4, 1.58587, 0.0964027],
[6.6, 1.61513, 0.508607],
[6.8, 1.6505, 1.33957],
[7, 1.69447, 1.69005],
[7.2, 1.75124, 1.02138],
[7.4, 1.82865, 0.295683],
[7.6, 1.94333, 0.098714]
]
LysoScintSpectrum =[
[2.2543, 0.0643],
[2.4797, 0.1929],
[2.6436, 0.4],
[2.8700, 1.],
[3.0996, 0.4071]
]
############################################################
def sortRI(elem):
  """
  A helper function used to sort the hits. The hits are organises like this:
  (id, [x,y,z,A,t]): the function returns the time as key to sort the hit list
  """
  return elem[0]


class Scintillator(object):
  __metaclass__ = ABCMeta
  
  @abstractmethod
  def Name(self):
    """
    Name of scintillator 
    """
    pass

  @abstractmethod
  def EffectiveAtomicNumber(self):
    """
    EffectiveAtomicNumber
    """
    pass

  @abstractmethod
  def RadiationLength(self):
    """
    Radiation length/Density  
    """
    pass

  @abstractmethod     
  def RefractionIndex(self):
    """
    returns the refraction index to the scintillation light
    """
    pass

  @abstractmethod     
  def RefractionIndexL(self,lamda):
    """
    returns the refraction index to the scintillation light
    """
    pass

  @abstractmethod     
  def RefractionIndexBlue(self):
    """
    returns the refraction index to the blue scintillation light
    """
    pass

    

  @abstractmethod
  def DecayConstant(self):
    """
    returns the decay time of the scintillation light
    """
    pass

  @abstractmethod   
  def Density(self):
    """
    Density 
    """
    pass

  @abstractmethod
  def ScintillationWavelength(self):
    """
    average scintillation wavelength 
    """
    pass

  @abstractmethod
  def PhotoelectricFractionAt511KeV(self):
    """
    Photoelectric Xs
    """
    pass

  @abstractmethod
  def AttenuationAt511KeV(self,Z):
    """
    Attenuation of a beam of energy E 511 keV in a thickness Z
    """
    pass

  @abstractmethod      
  def EfficiencyAt511KeV(self,Z):
    """
    Fraction of gammas of energy E 511 keV interacting in thickness Z
    """
    pass

  #@abstractmethod      
  def PhotoelectricEfficiencyAt511KeV(self,Z):
    """
    Fraction of gammas of energy E 511 keV interacting in thickness Z
    """
    return self.EfficiencyAt511KeV(Z)*self.PhotoelectricFractionAt511KeV()

  @abstractmethod
  def ScintillationPhotonsAt511KeV(self):
    """
    Number of scintillation photons produced by a photon of energy 511 keV
    """
    pass

  def DisplayProperties(self):
    s= """
        Name = %s  Z = %d  
        Density = %7.4g g/cm3 X0= %7.4g cm
        Refraction Index = %7.2f 
        Decay constant = %7.2f ns
        Number of scintillation photons = %7.2f 
        Scintillation wavelength = %7.2f nm
        Photoelectric fraction at 511 keV = %7.2f 
        Attenuation per cm = %7.2f 
        Efficiency per cm = %7.2f 
        """%(self.Name(), self.EffectiveAtomicNumber(),
          self.Density()/(g/cm3),self.RadiationLength()/cm,
          self.RefractionIndex(), 
          self.DecayConstant()/ns,
          self.ScintillationPhotonsAt511KeV(),
          self.ScintillationWavelength()/nm, 
          self.PhotoelectricFractionAt511KeV(),
          self.AttenuationAt511KeV(1*cm),
          self.EfficiencyAt511KeV(1*cm))
    return s

  def __str__(self):
    return self.DisplayProperties()


class LXe(Scintillator):
  def __init__(self,wi=15.6*eV,ws=16.6*eV,lambdaScint=172*nm,rayleigh=36*mm,
      tau1=2.2*ns,tau2=27*ns,tau3=45*ns,
      rtau1=0.065,rtau2=0.935,rtau3=0.0,nUV=1.70,nBlue=1.4):
    
    self.Z = 54
    self.A = 131.29*g/mol
    self.T = 160*kelvin
    self.x0 = 8.48 * g/cm2
    self.rho = 3*g/cm3
    self.dedx=0.35*keV/micron  
    self.wi=wi
    self.ws=ws
    self.lambdaScint = lambdaScint
    self.rayleigh = rayleigh
    self.tau1=tau1
    self.tau2=tau2
    self.tau3=tau3
    self.rtau1=rtau1
    self.rtau2=rtau2
    self.rtau3=rtau3
    self.nUV = nUV
    self.nBlue=nBlue
    lxri =[]  #transform to nm
    for elem in LXeRefractionIndex:
      ene = elem[0]
      n = elem[1]
      f = elem[2]
      x =[(1240./ene)*nm,n,f]
      #print x[0]/nm
      lxri.append(x)

    self.LXRI = sorted(lxri, key=sortRI)
    #print self.LXRI

    l,n = self.AverageLamdaAndRI()
    self.AverageLamda = l
    self.AverageRefractionIndexUV = n


  def Name(self):
    """
    Interface: Name
    """
    return "LXE"

  def EffectiveAtomicNumber(self):
    """
    Interface: Xenon atomic number
    """
    return self.Z

  def RadiationLength(self):
    """
    Interface: X0/rho
    """
    return self.x0/self.rho

  def RefractionIndex(self):
    """
    Interface: Take the UV 
    """
    return self.AverageRefractionIndexUV

  def DecayConstant(self):
    """
    Interface: returns the decay time of the scintillation light
    In LXe one has 2 constants, return the weighted average
    """
    return (self.Lifetime(1)*self.LifetimeRatio(1)+
           self.Lifetime(2)*self.LifetimeRatio(2))/2.
           
  def Density(self):
    """
    Interface: Density of LXe 
    """
    return self.rho

  def ScintillationWavelength(self):
    """
    Interface: average scintillation wavelength 
    """
    return self.lambdaScint


  def PhotoelectricFractionAt511KeV(self):
    """
    Interface: PE fraction 
    """
    return self.PhotoelectricFraction(511*keV)
    
  def AttenuationAt511KeV(self,Z):
    """
    Interface: Attenuation of a beam of energy 511 keV in a thickness Z
    """
    return self.Attenuation(511*keV,Z)
  
  def EfficiencyAt511KeV(self,Z):
    """
    Interface. Fraction of gammas of energy 511 keV interacting in thickness E
    """
    return self.Efficiency(511*keV,Z)

  def ScintillationPhotonsAt511KeV(self):
    """
     Interface: Number of scintillation photons produced by a photon of energy 511 keV
    """
    return self.ScintillationPhotons(511*keV)

  def AverageLamdaAndRI(self):
    """
    Returns the average lamda and refraction index
    """
    l=0.
    n=0.
    w=0.
    for elem in self.LXRI:
      l+=elem[0]*elem[2]
      n+=elem[1]*elem[2]
      w+=elem[2]
    return (l/w,n/w)

  def RmsRI(self):
    """
    Returns the rms of the refraction index
    """ 
    ns=0.
    w=0.
    lw, nw = self.AverageLamdaAndRI()

    print "nw = %7.2f"%(nw)
    for elem in self.LXRI:

      ns+=elem[2]*(elem[1] - nw)**2
      w+=elem[2]
    
      print " ni = %7.4g, ni - nw = %7.4g, wi =%7.4g ns = %7.4g"%(
        elem[1],elem[1] - nw,elem[2],ns)

    N=len(self.LXRI)
    a = N*ns
    b = (N-1)*w
    sw = sqrt(a/b)

    print " N = %7.2f, ns = %7.2f, w = %7.2f, a = %7.2f b = %7.2f"%(
        N,ns,w,a,b)
    return sw

  def RefractionIndexL(self,lamda):
    """
    returns the refraction index
    """

    if lamda < self.LXRI[0][0]:
      return self.LXRI[0][1] 
    elif lamda > self.LXRI[6][0]:
      return self.LXRI[6][1]
    else:
      for i in xrange(len(self.LXRI)-1):
        elem = self.LXRI[i]
        x0 = elem[0]
        y0 = elem[1]
        elem = self.LXRI[i+1]
        x1 = elem[0]
        y1 = elem[1]
        if lamda >= x0 and lamda < x1:
          break
      return lin(lamda,x0,y0,x1,y1)

  def AtomicNumber(self):
    """
    Xenon atomic number
    """
    return self.Z

  def AtomicMass(self):
    """
    Xenon atomic mass
    """
    return self.A

  def X0(self):
    """
    X0 in gr/cm2
    """
    return self.x0

  def TemperatureAt1Bar(self):
    """
    LXe Temperature
    """
    return self.T
   
  def RefractionIndexUV(self):
    return self.AverageRefractionIndexUV

  def RefractionIndexBlue(self):
    return self.nBlue

  def Lifetime(self,i):
    """
    i ->(1,3) for the three lifetimes.
    """
    if i == 1: 
      return self.tau1
    elif i == 2: 
      return self.tau2
    elif i == 3: 
      return self.tau3
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def LifetimeRatio(self,i):
    """
    i ->(1,3) for the three lifetimes.
    """
    if i == 1: 
      return self.rtau1
    elif i == 2: 
      return self.rtau2
    elif i == 3: 
      return self.rtau3
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def Wi(self):
    """
    Energy needed to produce an ionization pair
    """
    return self.wi

  def Ws(self):
    """
    Energy needed to produce scintillation photons
    """
    return self.ws

  def Rayleigh(self):
    """
    Attenuation due to Rayleigh Scattering
    """
    return self.rayleigh

  def dEdX(self):
    return self.dedx
  
  def ComptonCrossSection(self,E):
    """
    Compton = Incoherent Scattering
    """
    return ScatterIncoherent(E)

  def PhotoelectricCrossSection(self,E):
    """
    Photoelectric Xs
    """
    return Photoelectric(E)

  def TotalCrossSection(self,E):
    """
    Total Xs
    """
    return TotalInteraction(E)

  def PhotoelectricFraction(self,E):
    """
    Interface: PE fraction 
    """
    return self.PhotoelectricCrossSection(E)/self.TotalCrossSection(E)

  def Attenuation(self,E,Z):
    """
    Attenuation of a beam of energy E in a thickness Z
    """
    return TransmittedBeam(E,Z,self.Density())
  
  def Efficiency(self,E,Z):
    """
    Fraction of gammas of energy E interacting in thickness E
    """
    return InteractionFraction(E,Z,self.Density())

  def GammaPathLength(self,E):
    """
    gamma path length in xenon 
    """
    xz = self.TotalCrossSection(E)*self.Density()
    return 1./xz

  def ScintillationPhotons(self,E):
    """
     Number of scintillation photons produced by a photon of energy E
    """
    return E/self.Ws()

  def SPhotonsAt511KeV(self,i):
    if i == 1: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(1)
    elif i == 2: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(2)
    elif i == 3: 
      return self.ScintillationPhotons(511*keV)*self.LifetimeRatio(3)
    else:
      print "index must be 1,2 or 3"
      sys.exit(0)

  def IonizationElectrons(self,E):
    """
     Number of ionization electrons produced by a photon of energy E
    """
    return E/self.Wi()

  def CostPerGram(self):
    """
    Cost per gram
    """
    return 1.0  #in euro
  
  def __str__(self):
    s= """
        Name = LXe  Z = %d  A = %7.4g g/mole 
        Temperature at atm pressure (1 bar) = %7.2f kelvin
        Density = %7.4g g/cm3 X0= %7.4g g/cm2 X1= %7.4g cm
        de/dx = %7.4g keV/cm Ws = %7.4g eV Wi = %7.4g eV
        Rayleigh Scattering = %7.2g cm
        Scintillation wavelength = %7.2f nm
        Refraction Index (UV) = %7.2f 
        Refraction Index (Blue) = %7.2f 
        Lifetimes:
        tau1 = %7.2f ns, ratio tau 1 = %7.2f
        tau2 = %7.2f ns, ratio tau 2 = %7.2f
        tau3 = %7.2f ns, ratio tau 3 = %7.2f 
        """%(self.AtomicNumber(),self.AtomicMass()/(g/mol),
          self.TemperatureAt1Bar()/kelvin,
          self.Density()/(g/cm3),
          self.X0()/(g/cm2),(self.X0()/self.Density())/cm,
          self.dEdX()/(keV/cm),self.Ws()/eV, self.Wi()/eV,
          self.Rayleigh()/cm, self.ScintillationWavelength()/nm, self.RefractionIndex(), 
          self.nBlue,
          self.tau1/ns, self.rtau1,self.tau2/ns,self.rtau2,
          self.tau3/ns, self.rtau3)

    return s 

class LYSO(Scintillator):
    def __init__(self,Z=54,rho=7.3*g/cm3, n=1.82, X0=1.16*cm, LambdaAtt = 0.87*(1./cm),
        LambdaPeak=420*nm, tau = 50*ns, PhotoFraction = 0.3, Nphot = 15000):
    
        """
        Represents lyso
        """  
        self.x0 = X0  
        self.Z = Z
        self.rho = rho
        self.tau=tau
        self.n = n
        self.mu=LambdaAtt
        self.lambdaScint = LambdaPeak
        self.tau = tau
        self.photoF = PhotoFraction
        self.Nphot = Nphot

        lysct =[]  #transform to nm
        for elem in LysoScintSpectrum:
            ene = elem[0]
            w =elem[1]
            x =[(1240./ene)*nm,w]
            lysct.append(x)

        self.LYSC = sorted(lysct, key=sortRI)
    
        print "scintillation spectrum"
        for elem in self.LYSC:
            print " lambda = %7.2f nm w= %7.2g"%(elem[0]/nm,elem[1])

        l = self.AverageLamda()
        
    def Name(self):
      """
      Interface: Name
      """
      return "LYSO"

    def EffectiveAtomicNumber(self):
      """
      Interface: Lyso atomic number
      """
      return self.Z

    def RadiationLength(self):
      """
      Interface: X0/rho
      """
      return self.x0/cm

    def RefractionIndex(self):
      """
      Interface: 
       """
      return self.n

    def RefractionIndexL(self,lamda):
      """
      returns the refraction index
      """
      return self.n

    def RefractionIndexBlue(self):
      """
      returns the refraction index
      """
      return self.n

    def DecayConstant(self):
      """
      Interface
      """
      return self.tau
           
    def Density(self):
      """
      Interface: Density of LXe 
      """
      return self.rho

    def ScintillationWavelength(self):
      """
      Interface: average scintillation wavelength 
      """
      return self.lambdaScint

    def PhotoelectricFractionAt511KeV(self):
        """
        Interface Photoelectric Xs
        """
        return self.photoF

    def AttenuationAt511KeV(self,Z):
        """
        Interface: Attenuation of a beam of energy E 511 keV in a thickness Z
        """
        return exp(-Z*self.mu)
  
    def EfficiencyAt511KeV(self,Z):
        """
        Interface: Fraction of gammas of energy E 511 keV interacting in thickness Z
        """
        return 1. - self.AttenuationAt511KeV(Z)

    def ScintillationPhotonsAt511KeV(self):
        """
        Number of scintillation photons produced by a photon of energy 511 keV
        """
        return self.Nphot
    
    def AverageLamda(self):
        """
        Returns the average lamda 
        """
        l=0.
        w=0.
        for elem in self.LYSC:
            l+=elem[0]*elem[1]
            w+=elem[1]
        return (l/w)

    def X0(self):
        """
        EffectiveAtomicNumber
        """
        return self.x0

    def Lifetime(self):
        """
        returns the lifetime
        """
        return self.tau 

    def __str__(self):
        s= """
            Name = LYSO  Z = %d  
            Density = %7.4g g/cm3 X0= %7.4g g/cm2 
            Scintillation wavelength = %7.2f nm
            Refraction Index (Blue) = %7.2f 
            Lifetime = %7.2f ns
            ScintillationPhotons =  %7.2f
            Attenuation in 1 cm =  %7.2f
            PhotoelectricFraction = %7.2f
            """%(self.EffectiveAtomicNumber(),
            self.Density()/(g/cm3),
            self.X0()/cm,
            self.ScintillationWavelength()/nm, self.RefractionIndex(),
            self.Lifetime(), self.ScintillationPhotonsAt511KeV(),
            self.AttenuationAt511KeV(1.*cm),self.PhotoelectricFractionAt511KeV())

        return s 


def testLxe():
  lxe = LXe()

  print lxe 
  print lxe.DisplayProperties() 

  for l in drange(150*nm,450*nm,10*nm):
    print """
    for lamda = %7.2f nm (%7.2f eV) n = %7.2f
    """%(l/nm, 1240./(l/nm), lxe.RefractionIndexL(l))

  l,n = lxe.AverageLamdaAndRI()
  print """
  Average lamda = %7.2f nm ; average n = %7.2f
  """%(l/nm, n)

  rmsri = lxe.RmsRI()
  print """
    Weighted rms of n = %7.2f 
    """%(rmsri)

  print """
    Dn/n = %7.2f 
    """%(rmsri/n)

  print "Efficiency for 511 keV photons" 

  for z in drange(1., 11., 1.):
    print """
       z = %7.2g cm LXe eff = %7.2g 
      """%(z,
           lxe.Efficiency(511*keV,z*cm))

  print """Photoelectric fraction 
              at 511 keV photons = %7.2g"""%(
           lxe.PhotoelectricCrossSection(511*keV)/
           lxe.TotalCrossSection(511*keV))
    
  print """
       Gamma path lenght in LXe for 511 keV photons = %7.2g cm 
      """%(
           lxe.GammaPathLength(511*keV)/cm)

    
  print """
        Number of scintillation photons Ns (511 keV, LXe)= %7.2g 
        with tau1 = %7.2f ns lifetime: = %7.2f
        with tau2 = %7.2f ns lifetime: = %7.2f
        with tau3 = %7.2f ns lifetime: = %7.2f
      """%(
           lxe.ScintillationPhotons(511*keV),lxe.Lifetime(1),
           lxe.SPhotonsAt511KeV(1),
           lxe.Lifetime(2),
           lxe.SPhotonsAt511KeV(2),
           lxe.Lifetime(3),
           lxe.SPhotonsAt511KeV(3)
           )
           
def testLyso():
  lyso = LYSO()

  print lyso 
  lyso.DisplayProperties() 

  print "Average Lamda = %7.2f"%(lyso.AverageLamda()/nm)

  for z in drange(1., 11., 1.):
    print """
       z = %7.2g cm LYSO eff = %7.2g 
      """%(z,
           lyso.EfficiencyAt511KeV(z*cm))

def plotLXe():
  Lambda=[]
  I=[]
  N=[]
  for elem in LXeRefractionIndex:
    ene = elem[0]
    n = elem[1]
    f = elem[2]
    Lambda.append(1240./ene)
    I.append(f)
    N.append(n)
  
  plt.plot(Lambda,I)
  plt.show()
  plt.plot(Lambda,N)
  plt.show()

def plotLYSO():
  Lambda=[]
  I=[]
  
  for elem in LysoScintSpectrum:
    ene = elem[0]
    f = elem[1]
    Lambda.append(1240./ene)
    I.append(f)
    
  
  plt.plot(Lambda,I)
  plt.show()


if __name__ == '__main__':
    
    #testLxe()
    #testLyso()
    plotLXe()
    plotLYSO()
    

        