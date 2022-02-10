import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import time
from copy import deepcopy as copy


def defaultKeyVal(d,k,v):
  if k in d.keys():
    return d[k]
  else:
    return v
    

class dictClass(dict):
  """ 
  This class is essentially a subclass of dict
  with attribute accessors, one can see which attributes are available
  using the `keys()` method.
  """
  def __dir__(self):
      return self.keys()
    
  def __getattr__(self, name):
    try:
      return self[name]
    except KeyError:
      raise AttributeError(name)
  if dict==None:
    __setattr__ = {}.__setitem__
    __delattr__ = {}.__delitem__
  else:
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

  # def __repr__(self):
    # if self.keys():
      # L = list(self.keys())
      # L = [str(L[i]) for i in range(len(L))]
      # m = max(map(len, L)) + 1
      # f = ''
      # for k, v in self.items():
        # if isinstance(v,dict):
          # f = f + '\n'+str(k).rjust(m) + ': ' + repr(k) + ' class'
        # else:
          # f = f + '\n'+str(k).rjust(m) + ': ' + repr(v)
      # return f
    # else:
      # return self.__class__.__name__ + "()"

  def find_key(self,val):
    if val==None :
      return
    for k in self.keys():
      if self[k]==val:
        return k


def read_ascii(fid,format_id,npt=None):

    if format_id==1:
        tmp = pd.read_csv("ascii.#"+f"{fid:03d}",header=None,nrows=npt,skiprows=4,delim_whitespace=True)
        tmp.columns =["charge", "Amass", "x [cm]", "x' [rad]", "y [cm]", "y' [rad]", "phase [rad]", "beta", "zloss [cm]", "spin"]
    elif format_id==2:
        tmp = pd.read_csv("ascii.#"+f"{fid:03d}",header=None,nrows=npt,skiprows=4,delim_whitespace=True)
        tmp.columns =["charge", "Amass", "x [cm]", "x' [mrad]", "y [cm]", "y' [mrad]", "time [nsec]", "W-W0 [MeV/ u]", "zloss [cm]", "spin"]
    elif format_id==3:
        tmp = pd.read_csv("ascii.#"+f"{fid:03d}",header=None,nrows=npt,skiprows=2,delim_whitespace=True)
        tmp.columns =["Q", "A",  "Z",  "Zs [m]", "x [cm]", "x' [mrad]", "y [cm]", "y' [mrad]", "time [nsec]", "W [MeV/ u]", "spin"]

    return tmp
    
    
def read_dist(fname="read_dis.out",npt=None,delim_whitespace=True):   
    out = pd.read_csv(fname,header=None,nrows=npt,skiprows=4,delim_whitespace=delim_whitespace)
    out.columns =["x [cm]", "x' [rad]", "y [cm]", "y' [rad]", "phase [rad]", "beta", "spin"]
    with open(fname,"r") as f:
        for i in range(4):
            tmp = f.readline()
        beta_ref = tmp.split()[-2]
        if delim_whitespace:
            beta_ref = np.float(beta_ref)
        else:
            beta_ref = np.float(beta_ref[:-1]) 
    
    return out, beta_ref
    
    
def keV2beta(keV,Mu=931.49432e3):
    beta = (((keV/Mu)**2+2*(keV/Mu))/(1.+keV/Mu)**2)**0.5
    return beta
    
    
def beta2keV(beta,Mu=931.49432e3):
    gamma = (1./(1.-beta**2))**0.5
    keV = (gamma - 1)*Mu
    return keV  
    
    
def write_dist(dist,charge=17,beta_ref=None,keV_ref=None,fname="read_dis.dat"):
    if beta_ref==None:
        if keV_ref==None:
            raise ValueError("at least one of 'beta_ref' or 'beta_ref' must be specified.")
        beta_ref = keV2beta(keV_ref)
    if keV_ref==None:
        keV_ref = beta2keV(beta_ref)
    read_dist_header = [str(keV_ref)+" 1\n",
                        str(len(dist))+"\n",
                        str(charge)+"\n",
                        "0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 "+str(beta_ref)+" 0\n"]

    with open(fname,'w') as f:
        for i in range(4):
            f.write(read_dist_header[i])

        if isinstance(dist, pd.DataFrame):
            isspin = "spin" in dist.columns
            if isspin:
                try:
                    tmp = dist[["x [cm]", "x' [rad]", 
                                "y [cm]", "y' [rad]", 
                                "phase [rad]", 'beta', "spin"]].values
                except:
                    raise ValueError("invalid dataframe shape")
                for i in range(len(dist)):
                    tmp_str = ''
                    for val in tmp[i][:6]:
                        tmp_str += str(val)+' '
                    f.write(tmp_str+str(int(tmp[i][-1]))+'\n')
            else:
                try:
                    tmp = dist[["x [cm]", "x' [rad]", 
                                "y [cm]", "y' [rad]", 
                                "phase [rad]", 'beta']].values
                except:
                    raise ValueError("invalid dataframe shape")
                    
                for i in range(len(dist)):
                    tmp_str = ''
                    for val in tmp[i][:6]:
                        tmp_str += str(val)+' '
                    f.write(tmp_str+'0\n')
         
        else:
            tmp = dist
            for i in range(len(dist)):
                tmp_str = ''
                for val in tmp[i]:
                    tmp_str += str(val)+' '
                f.write(tmp_str+'0\n')
    
    
def gaussian_kde(x,y,npt=0):
    if npt!=0 and len(x)<npt:
        #print('cal ipt')
        #ipt = np.random.choice(np.arange(len(x)),size=npt,replace=False)
        #print('cal kde')
        #return stats.gaussian_kde([x[ipt],y[ipt]])
        return stats.gaussian_kde([x[:npt],y[:npt]])
    else:
        return stats.gaussian_kde([x,y])


def df_scatter_plot(df,column_x,column_y,s=4,density=True,npt=0,color="red",xlim=None,ylim=None):
    x=df[column_x].values
    y=df[column_y].values
    if density:
        kde = gaussian_kde(x,y,npt)
        color = kde.evaluate([x,y])
    plt.scatter(x,y,c=color,s=s)
    plt.xlabel(column_x)
    plt.ylabel(column_y)
    if xlim:
        plt.xlim(xlim)
    else:
        x99=np.quantile(x,0.99)
        x01=np.quantile(x,0.01)
        xr = x99-x01
        plt.xlim(x01-0.2*xr,x99+0.2*xr)
    if ylim:
        plt.ylim(ylim)
    else:
        x99=np.quantile(y,0.99)
        x01=np.quantile(y,0.01)
        xr = x99-x01
        plt.ylim(x01-0.2*xr,x99+0.2*xr)
    
     
def append_measurement(df,key,track):
    tmp = track.get_profile(key)
    tmp['name'] = key
    return df.append(tmp, ignore_index=True)    


def get_measurements(track,keys=None,sclinac_file=None):
    if sclinac_file==None and keys==None:
        raise ValueError("At least one of 'sclinac_file' or 'keys' needs to be specified")
        
    if keys==None:
        with open(sclinac_file,'r') as f:
            lattice = f.readlines()
        keys = []
        i=0
        for line in lattice:
            if 'stop' in line and line[0]!='!':
                break
            if '!== '  == line[:4]:
                tmp = line[4:].rstrip()
                keys.append(tmp)
                i+=1
                
    df = pd.DataFrame()
    for key in keys:
        df = append_measurement(df,key,track)
    df = df.set_index('name')
    return df
    
           
def adjust_trk_file(trk_file_in="track_ref.dat",trk_file_out="track.dat",**kwargs):
    
    def _replace_line(key,line,val):
        if key in line:
            ikey = line.find(key)
            i0 = line[ikey:].find("=") +ikey +1
            i1 = line[i0:].find(",") 
            if i1<0:
                i1 = line[i0:].find("\n")
            i1 += i0
            return line[:i0] + str(val) + line[i1:]
        else:
            return line
    

    keys = kwargs.keys()
    with open(trk_file_in,"r") as f:
        lines = f.readlines()
        i=0
        while i < len(lines):
            for key in keys:
                lines[i] = _replace_line(key,lines[i],kwargs[key])
            i+=1
                       
    with open(trk_file_out,"w") as f:
        f.writelines(lines)
        
        
def slice_lattice(pv_from=None,pv_to=None,
                  lattice_in="./lattices/sclinac_ref_allison2RFQ0.dat",
                  lattice_out="./lattices/sclinac_tmp.dat"):
    
    with open(lattice_in,"r") as f:
        lines = f.readlines()
    
    if pv_from==None:
        i_from = 0
    else:
        i_from = -1
    i_to = len(lines)
        
    for i,line in enumerate(lines):
        if i_from == -1:
            if pv_from in line:
                i_from = i
                if pv_to == None:
                    break
        else:
            if pv_to in line:
                i_to = i
                break
            
    with open(lattice_out,"w") as f:
        f.writelines(lines[i_from:i_to])


def get_lattice_pvs(lattice):
    
    with open(lattice,"r") as f:
        lines = f.readlines()

    pvs = []
    for line in lines:
        if line[:3] == "!==":
            pv = line[3:].rstrip().lstrip()
            if not pv=="":
                pvs.append(pv)
    return pvs
    
    
def getElem(elem_type,attributes=None):
    elem = dictClass({'type':elem_type})
    if elem_type == 'ascii':
        try:
            elem.format_id = attributes["format_id"]
        except:
            elem.format_id = 1
        try:
            elem.fname_id = attributes["fname_id"]
        except:
            elem.fname_id = 10000
    else:
        raise AttributeError(elem_type+'is not available type.')
        
    return elem


def _elem2str(elemDict): 
    if elemDict.type == 'ascii':
        elemStr = [str(elemDict.fname_id), "ascii", str(elemDict.format_id), "0"]
    else:
        raise AttributeError(elemDict.type+'is not available type.')
    
    elemStr.append("\n")
    
    return " ".join(elemStr)


def insert_elem_before(lattice_lines,pv,elem):
    for i,line in enumerate(lattice_lines):
        if pv in line:
            break
            
    lattice_lines.insert(i,_elem2str(elem))
    return lattice_lines