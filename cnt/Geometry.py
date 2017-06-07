from abc import ABCMeta, abstractmethod
from Util import *

# class Cell:
#   __metaclass__ = ABCMeta
#   @abstractmethod
#   def Volume(self):
#     pass

#   @abstractmethod
#   def Sr(self):
#     pass

#   @abstractmethod
#   def L(self):
#     pass

#   @abstractmethod
#   def R(self):
#     pass

###########################################################
class XYZBox:

  def __init__(self,x,y,z):
    """
    Defines a box of dimensions x,y,z
    z is considered the longitudinal dimension 
    """
    self.x=x
    self.y=y
    self.z=z
        
  def V(self):
    return self.x*self.y*self.z

  def Sxy(self):
    return self.x*self.y

  def Sxz(self):
    return self.x*self.z

  def Syz(self):
    return self.x*self.z


  def __str__(self):
        
    s= """
        XYZBox:
        x = %7.2f y = %7.2f  z = %7.2f
        Sxy = %7.2f Sxz = %7.2f Syx = %7.2f 
        V = %7.2f
      """%(self.x, self.y, self.z, self.Sxy(),
        self.Sxz(), self.Syz(), self.V())
    return s

###########################################################
class Box(XYZBox):

  def __init__(self, boxCoord):
    """
        Defines a box in terms of the box coordinates in the reference system
        of G4. The box coordinates are defined by 8 vertices.
        Each vertex is an (x,y,z) point. 
        Define the box with the following convention:
        v1 = (xl,yb,zb), where xl = lefmost x coordinate,
                               yb = bottom y coordinate
                               zb = back z coordinate
        x-------x (xr,yt)
        |       |
        |       |
        x-------x            zb ----- zf
        (xl,yb)

        v2 = (xl,yt,zb), where xl = lefmost x coordinate,
                               yt = top y coordinate
                          

        v3 = (xr,yb,zb), where xr = rightmost x coordinate
        v4 = (xr,yt,zb)
        v5 = (xl,yb,zf), where zf = front z coordinate
        v6 = (xl,yt,zf)
        v7 = (xl,yt,zf)
        v8 = (xr,yb,zf)                           
 
    """
    self.boxCoord =boxCoord
    v1 = self.boxCoord[0] 
      
    self.xmin = v1[0]
    self.ymin = v1[1]
    self.zmin =v1[2]

    v8 = self.boxCoord[7] 
    
    self.xmax = v8[0]
    self.ymax = v8[1]
    self.zmax = v8[2]

    x = abs(self.xmax-self.xmin)
    y = abs(self.ymax-self.ymin)
    z = abs(self.zmax-self.zmin)

    XYZBox.__init__(self, x,y,z)

  def Active(self,coord):
    """
    returns true if coord = (x,y,z) is in the active volume of the box
    defined by xmin-xmax, ymin-ymax, zmin-zmax, false otherwise
    """

    x,y,z = coord
    box = False
    if x >= self.xmin and x <= self.xmax:
      if y >= self.ymin and y <= self.ymax:
        if z < 0:
          if z <= self.zmin and z >= self.zmax:
            box = True
        else:
          if z >= self.zmin and z <= self.zmax:
            box = True
    
    return box

  def __str__(self):
       
    s =  XYZBox.__str__(self)
    s+= """
        Box:
        xmin = %7.2f ymin = %7.2f zmin = %7.2f
        xmax = %7.2f ymax = %7.2f zmax = %7.2f
      """%(self.xmin, self.ymin, self.zmin,
           self.xmax, self.ymax, self.zmax)
    return s

###########################################################
def distance((x,y,z),(x0,y0,z0)):
  """
  distance between two points
  """
  return sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)

###########################################################
def extrapToZ(zc,(x0,y0,z0),(px,py,pz)):
  """
  Extrapolate to a plane located at z =zc
  V = (x0,y0,z0)
  P = (px,py,pz)
  The equation of the straight line in parametrics 
  x(t) = x0 + px*t
  y(t) = y0 + py*t
  z(t) = z0 + pz*t

  Extrap to z = zc -- > zc = z0 + pz* t --> t = (zc-z0)/pz
  and then:
  x(t) = x0 + [(zc-z0)/pz]*px = x0+ (px/pz)*(zc-z0)
  y(t) = y0 + [(zc-z0)/pz]*py = y0+ (py/pz)*(zc-z0)
  """
  x = x0+ (px/pz)*(zc-z0)
  y = y0+ (py/pz)*(zc-z0)

  return (x,y)

def extrapToX(xc,(x0,y0,z0),(px,py,pz)):
  """
  Extrapolate to a plane located at x =xc
  V = (x0,y0,z0)
  P = (px,py,pz)
  """
  
  z = z0+ (pz/px)*(xc-x0)
  y = y0+ (py/px)*(xc-x0)

  return (y,z)

def extrapToY(yc,(x0,y0,z0),(px,py,pz)):
  """
  Extrapolate to a plane located at x =xc
  V = (x0,y0,z0)
  P = (px,py,pz)
  """
  z = z0+ (pz/py)*(yc-y0)
  x = y0+ (px/py)*(yc-y0)

  return (x,z)

###########################################################
def pathInBox((x0,y0,z0), (px,py,pz), box):
    """
    Given an entry point (x0,y0,z0) (in face zmin) and a momentum (px,py,pz)
    compute the extrapolation to all possible exiting surfaces
    (zmax,xmin,xmax,ymin,ymax), and return the smallest one
    """
    
    z = box.zmax
    x,y = extrapToZ(z,(x0,y0,z0),(px,py,pz))
    d1 = distance((x,y,z),(x0,y0,z0))

    x = box.xmin
    y,z = extrapToX(x,(x0,y0,z0),(px,py,pz))
    d2 = distance((x,y,z),(x0,y0,z0))

    x = box.xmax
    y,z = extrapToX(x,(x0,y0,z0),(px,py,pz))
    d3 = distance((x,y,z),(x0,y0,z0))
    
    y = box.ymin
    x,z = extrapToX(y,(x0,y0,z0),(px,py,pz))
    d4 = distance((x,y,z),(x0,y0,z0))

    y = box.ymax
    x,z = extrapToX(y,(x0,y0,z0),(px,py,pz))
    d5 = distance((x,y,z),(x0,y0,z0))


    d = min(d1,d2,d3,d4,d5)

    return d

###########################################################
def propagateInBox((x0,y0,z0), (px,py,pz), d):
    """
    Given a photon with entry point (x0,y0,z0),normalized momentum (px,py,pz)
    propagate the photon to distance d inside the box 
    """
    x = x0 + px*d
    y = y0 + py*d
    z = z0 + pz*d

    return (x,y,z)
       

    

###########################################################
class Point3D:
  def __init__(self,x=0,y=0,z=0):
    """
    A 3D point  
    
    """
    
    self.x = x  # creation 
    self.y = y
    self.z = z
    
  def X(self):
    """
    x coordinate 
    """
    return self.x

  def Y(self):
    """
    Y coordinate 
    """
    return self.y

  def Z(self):
    """
    z coordinate 
    """
    return self.z

  def XYZ(self):
    """
    xyz coordinates
    """
    return (self.x,self.y,self.z)

  def __str__(self):

    s= """
        x = %7.2f  y = %7.2f  z = %7.2f  
      """%(self.X(), self.Y(), self.Z())

    return s

###########################################################
class Point4D(Point3D):
  def __init__(self,x=0,y=0,z=0, t=0):
    """
    A 4D point  
    
    """
    Point3D.__init__(self,x,y,z)
    self.t = t

  def T(self):
    """
    t coordinate 
    """
    return self.t

  def XYZT(self):
    """
    t coordinate 
    """
    return (self.x,self.y,self.z,self.t)

  def __str__(self):

    s= """
        x = %7.2f  y = %7.2f z = %7.2f  t = %7.2f  
      """%(self.X(), self.Y(), self.Z(), self.T())

    return s 


if __name__ == '__main__':
    
    boxCoord =[[-12.8, -12.8, -100],[-12.8, 12.8, -100],[12.8, -12.8, -100],  
    [12.8, 12.8, -100],[-12.8, -12.8, -130],[-12.8, 12.8, -130],
    [12.8, -12.8, -130],[12.8, 12.8, -130]]
    box = Box(boxCoord)
    fid = box.Active([10.,5., -110.])
    print fid
    
    print box