from math import *
from scipy import interpolate

def decompose(GCL1,GCL2,XX,YY):
  #XX,YY is the grid of 2D mesh
  # this function decompose the line (X1,Y1)->(X2,Y2) into different pieces
  if 1==2:# dont interpolate
  #if GCL1.hasE() and GCL2.hasE(): 
    Xmin=min(GCL1.value['X'],GCL2.value['X'])
    Xmax=max(GCL1.value['X'],GCL2.value['X'])
    Ymin=min(GCL1.value['Y'],GCL2.value['Y'])
    Ymax=max(GCL1.value['Y'],GCL2.value['Y'])
    GCL_out=[]
    GCL_out_tmp=[]
    GCL_tmp=[]
    distance2index={}
    distance=[]
    #select by X
    for x in XX:
      if (x>Xmin)and(x<Xmax):
        GCL_tmp=GCL2.clone()
        #inside
        f_tmp=(x-GCL1.value['X'])/(GCL2.value['X']-GCL1.value['X'])
        #X
        GCL_tmp.value['X']=x
        #Y
        GCL_tmp.value['Y']=(GCL2.value['Y']-GCL1.value['Y'])*f_tmp+GCL1.value['Y']
        #Z
        GCL_tmp.value['Z']=(GCL2.value['Z']-GCL1.value['Z'])*f_tmp+GCL1.value['Z']
        #E
        GCL_tmp.value['E']=(GCL2.value['E']-GCL1.value['E'])*f_tmp+GCL1.value['E']
        #compute distance
        d_tmp=(GCL_tmp.value['X']-GCL1.value['X'])**2+(GCL_tmp.value['Y']-GCL1.value['Y'])**2
        distance.append(d_tmp)
        distance2index[d_tmp]=len(distance)-1
        GCL_out_tmp.append(GCL_tmp)
    #select by y
    for y in YY:
      if (y>Ymin)and(y<Ymax):
        GCL_tmp=GCL2.clone()
        #inside
        f_tmp=(y-GCL1.value['Y'])/(GCL2.value['Y']-GCL1.value['Y'])
        #X
        GCL_tmp.value['X']=(GCL2.value['X']-GCL1.value['X'])*f_tmp+GCL1.value['X']
        #Y
        GCL_tmp.value['Y']=y
        #Z
        GCL_tmp.value['Z']=(GCL2.value['Z']-GCL1.value['Z'])*f_tmp+GCL1.value['Z']
        #E
        GCL_tmp.value['E']=(GCL2.value['E']-GCL1.value['E'])*f_tmp+GCL1.value['E']
        #compute distance
        d_tmp=(GCL_tmp.value['X']-GCL1.value['X'])**2+(GCL_tmp.value['Y']-GCL1.value['Y'])**2
        distance.append(d_tmp)
        distance2index[d_tmp]=len(distance)-1
        GCL_out_tmp.append(GCL_tmp)
    #sort by distance
    #GCL_out.append(GCL1.clone())
    distance.sort()
    for d in distance:
      GCL_out.append(GCL_out_tmp[distance2index[d]])
    GCL_out.append(GCL2.clone())
  else:
    GCL_out=[GCL2.clone()]
  return GCL_out
  
class GcodeLine:
  def __init__(self,line):
    line=line.replace('\r','').replace('\n','')
    temp=line.split()
    #print("length: ",len(temp))
    try:
      if temp[0] in ["G0", "G1"]:
        self.type="data"
        self.value={}
        for term in temp:
          if term[0] in ['G', 'F']:
            self.value[term[0]]=int(term[1:])
          elif term[0] in ['X', 'Y', 'Z', 'E']:
            self.value[term[0]]=float(term[1:])
          else:
            print ("Ignored : ",term)
          
      else:
        self.type="info"
        self.value=line
    except:
      self.type="info"
      self.value=line
      
  def toString(self,ff=1):
    if self.type=="info":
      res=self.value
    else:
      res=""
      for term in ['G', 'F', 'X', 'Y', 'Z', 'E']:
        try:
          if term in ['X', 'Y', 'Z', 'E']:
            res+=term+"%.3f "%(self.value[term])
          else:
            res+=term+str(self.value[term])+' '
        except:
          pass
      try:
        if res[-1]==" ":
          res=res[:-1]
      except:
        pass
    return res
    
  def hasZ(self):
    res=0
    try:
      if self.value['Z']!=None:
        res=1
    except:
      pass
    return res
    
  def setZ(self,Z):
    if self.type=="data":
      self.value['Z']=Z
  
  def hasXY(self):
    res=0
    try:
      if (self.value['X']!=None) and (self.value['Y']!=None):
        res=1
    except:
      pass
    return res
    
  def setXY(self,X,Y):
    if self.type=="data":
      self.value['X']=X
      self.value['Y']=Y
      
  def hasE(self):
    try: 
      self.value['E']
      res=1
    except:
      res=0
    return res
  
  def clone(self):
    return GcodeLine(self.toString())

class GcodeManager:
  def __init__(self):
    self.data=[]
    self.data_new=[]
    self.fn=None
    self.dZ=None
    self.ff=1.0
    self.Xmin=2000
    self.Xmax=-2000
    self.Ymin=2000
    self.Ymax=-2000
    self.Zmin=2000
    self.Zmax=-2000
    self.e=0.1
    self.notes=[]
    
  def read(self, gcodef):
    self.fn=gcodef.replace('.gcode','_gcm.gcode')
    f=open(gcodef,'r')
    line=f.readline()
    currentZ=None
    currentX=None
    currentY=None
    while line:
      #if line[0]==';':
      #  pass
      if 1==2:
        pass
      else:
        tmp=GcodeLine(line)
        #explicite Z
        if tmp.hasZ():
          currentZ=tmp.value['Z']
        else:
          #tmp.setZ(currentZ)
          pass
        #explicite X Y
        if tmp.hasXY():
          currentX=tmp.value['X']
          currentY=tmp.value['Y']
        else:
          #tmp.setXY(currentX,currentY)
          pass
        self.data.append(tmp)
      line=f.readline()
    f.close()
    self.getEnvelope()
    
  def getEnvelope(self):
    for line in self.data:
      if line.type=="data":
        try:
          #print(line.value)
          self.Xmin=min(self.Xmin,line.value['X'])
          self.Xmax=max(self.Xmax,line.value['X'])
          self.Ymin=min(self.Ymin,line.value['Y'])
          self.Ymax=max(self.Ymax,line.value['Y'])
          self.Zmin=min(self.Zmin,line.value['Z'])
          self.Zmax=max(self.Zmax,line.value['Z'])
        except:
          print(line.value)
    self.Xmin=floor(self.Xmin/10)*10
    self.Xmax=ceil(self.Xmax/10)*10
    self.Ymin=floor(self.Ymin/10)*10
    self.Ymax=ceil(self.Ymax/10)*10
    print("X range = %4d ->  %4d [ mm ]"%(self.Xmin,self.Xmax))
    print("Y range = %4d ->  %4d [ mm ]"%(self.Ymin,self.Ymax))
    print("Z range = %.2f -> %.2f [ mm ]"%(self.Zmin,self.Zmax))
    #leveling(self.Xmin,self.Xmax,self.Ymin,self.Ymax,1,5)
    
  def offset(self,dz):
    self.notes.append(";[GCM] vertical offset dZ=%.3f\n"%(dz))
    #offset to test
    if len(self.data_new)==0:
      self.data_new=self.data
    for line in self.data_new:
      if line.type=="data":
        try:
          line.value['Z']+=dz
        except:
          pass
    #end of test
  
  def offsetFeedrate(self,coef):
    self.notes.append(";[GCM] correct print speed by factor %.3f\n"%(coef))
    if len(self.data_new)==0:
      self.data_new=self.data
    for line in self.data_new:
      if line.type=="data":
        try:
          if (line.value['F']<7000):
            line.value['F']=int(line.value['F']*coef/10.)*10
        except:
          pass
  
  def offsetExtruder(self,coef):
    self.notes.append(";[GCM] correct extruder quantity by factor %.3f\n"%(coef))
    if len(self.data_new)==0:
      self.data_new=self.data
    for line in self.data_new:
      if line.type=="data":
        try:
          line.value['E']=line.value['E']*coef
        except:
          pass
    
  def correct(self,calibFile):
    #read calib file and construct dZ(x,y)
    f=open(calibFile,'r')
    lines=f.readlines()
    f.close()
    N=int(lines[0].replace('\r','').replace('\n','').replace('=',' ').split()[1])
    x=[]
    y=[]
    dz=[]
    for i in range(N):
      tmp=lines[5+i*4].replace('X','').replace('Y','').replace('\r','').replace('\n','').split()
      x.append(float(tmp[1]))
      y.append(float(tmp[2]))
      tmp=lines[6+i*4].replace('Z','').replace('\r','').replace('\n','').split()
      dz.append(float(tmp[1]))
    self.dZ=interpolate.interp2d(x,y,dz)
    
    #i=1,j=1 iterate hand in hand
    j=0
    prePos=None
    thisPos=None
    self.data_new=[]
    for line in self.data:
      if line.type=="info":
        self.data_new.append(line)
      else:
        # this is a data line
        if prePos==None:
          #no need to split
          tmp=line.clone()
          
          try:
            tmp.value['Z']+=self.dZ(tmp.value['X'],tmp.value['Y'])[0]
          except:
            print(tmp.value)
          self.data_new.append(tmp)
          prePos=line
        else:
          #split if needed
          for gcl in decompose(prePos,line,x,y):
            tmp=gcl.clone()
            try:
              tmp.value['Z']+=self.dZ(tmp.value['X'],tmp.value['Y'])[0]
            except:
              print(tmp.value)
            self.data_new.append(tmp)
        
  def export(self):
    f=open(self.fn,'w')
    for line in self.notes:
      f.write(line)
    for line in self.data_new:
      f.write(line.toString(self.ff)+"\n")
    f.close()

if __name__=='__main__':
  GM=GcodeManager()
  #GM.read("Sample.gcode")
  #GM.read("PortePortable.gcode")
  #GM.read("Badge_NG7F949.gcode")
  #GM.read("PortePortable20pc.gcode")
  #GM.read("Badge_Crapo.gcode")
  #GM.read("Badge_NG7F949_v4.gcode")
  GM.read("Crapolat_Prototype_3-0.gcode")
  #GM.offset(4.99+0.1)
  GM.offset(4.99+0.1)
  GM.offsetFeedrate(0.8)
  #GM.offsetExtruder(1.5)
  GM.offsetExtruder(2.0)
  GM.export()
