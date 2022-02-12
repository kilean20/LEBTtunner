from . import TRACKutil
from . import TRACKutil_kilean as kutil
# from .dictClass import dictClass

import numpy as np
import pandas as pd
import re

import requests
import os
import time
t0 = time.time()

class LEBT():
    def __init__(self, QperA=17/78, current=25., _dutyfactor=1, npt=2000, init_beam_arg=np.zeros(14),working_directory='./'):
        #for 78Kr17
        self._cS4 = -0.002891 
        self._scale = (17/78)/QperA
        self._current = current
        self._dutyfactor = _dutyfactor
        self._npt = npt
        if working_directory[-1]=="/":
            self._working_directory = working_directory
        else:
            self._working_directory = working_directory + "/"
        if not os.path.isdir(self._working_directory):
            os.makedirs(self._working_directory)
        self._download_TRACK_files()
            
        self._simulated = False
        self._measured_data = None
        self._name_family_simval = [['SOL_D0787','B',69.9*self._cS4*self._scale],
                                    ['SOL_D0802','B',74.0*self._cS4*self._scale],
                                    ['SOL_D0818','B',72.5*self._cS4*self._scale],
                                    #= L-LEBT
                                    ['SOL_D0951','B',78.8*self._cS4*self._scale],
                                    ['SOL_D0967','B',104.4*self._cS4*self._scale],
                                    ['SOL_D0982','B',146.0*self._cS4*self._scale],
                                    ['SOL_D0995','B',155.9*self._cS4*self._scale],
                                    #=== E-quads [V]
                                    #= U-LEBT
                                    ['QE_D0743','V',50.4*self._scale],
                                    ['QE_D0746','V',-1920.1*self._scale],
                                    ['QE_D0749','V',3164.5*self._scale],
                                    ['QE_D0767','V',-3290.*self._scale],
                                    ['QE_D0770','V',5248.*self._scale],
                                    ['QE_D0776','V',-4394.*self._scale],
                                    ['QE_D0780','V',2378*self._scale],
                                    #= V-LEBT
                                    ['QE_D0844','V',-1542*self._scale],
                                    ['QE_D0848','V',3981.*self._scale],
                                    ['QE_D0851','V',-3312.*self._scale],
                                    ['QE_D0871','V',2379.*self._scale],
                                    ['QE_D0874','V',-3547*self._scale],
                                    ['QE_D0878','V',798.*self._scale],
                                    ['QE_D0891','V',798*self._scale],
                                    ['QE_D0895','V',-3547.*self._scale],
                                    ['QE_D0898','V',2379.*self._scale],
                                    ['QE_D0918','V',-3312.*self._scale],
                                    ['QE_D0921','V',3981.*self._scale],
                                    ['QE_D0924','V',-1542*self._scale],
                                    #= U-LEBT
                                    ['DCHV_D0773','tm_xkick',0.0],
                                    ['DCHV_D0790','tm_xkick',0.0],
                                    ['DCHV_D0805','tm_xkick',0.0],
                                    ['DCHV_D0821','tm_xkick',0.0],
                                    ['DCHV_D0773','tm_ykick',0.0],
                                    ['DCHV_D0790','tm_ykick',0.0],
                                    ['DCHV_D0805','tm_ykick',0.0],
                                    ['DCHV_D0821','tm_ykick',0.0],
                                    #= V-LEBT
                                    ['DCHV_D0840','tm_xkick',0.0],
                                    ['DCHV_D0868','tm_xkick',0.0],
                                    ['DCHV_D0880','tm_xkick',0.0],
                                    ['DCHV_D0901','tm_xkick',0.0],
                                    ['DCHV_D0929','tm_xkick',0.0],
                                    ['DCHV_D0840','tm_ykick',0.0],
                                    ['DCHV_D0868','tm_ykick',0.0],
                                    ['DCHV_D0880','tm_ykick',0.0],
                                    ['DCHV_D0901','tm_ykick',0.0],
                                    ['DCHV_D0929','tm_ykick',0.0],
                                    #= L-LEBT
                                    ['DCHV_D0948','tm_xkick',0.0],
                                    ['DCHV_D0964','tm_xkick',0.0],
                                    ['DCHV_D0979','tm_xkick',0.0],
                                    ['DCHV_D0992','tm_xkick',0.0],
                                    ['DCHV_D0948','tm_ykick',0.0],
                                    ['DCHV_D0964','tm_ykick',0.0],
                                    ['DCHV_D0979','tm_ykick',0.0],
                                    ['DCHV_D0992','tm_ykick',0.0],
                                    #= Multi harmonic buncher
                                    ['MHB_D0987','phi0',-75.0],
                                    ['MHB_D0987','scl_fac0',11.515],
                                    ['MHB_D0987','phi1',-86.63],
                                    ['MHB_D0987','scl_fac1',8.907],
                                    ['MHB_D0987','phi2',3.58],
                                    ['MHB_D0987','scl_fac2',11.871], 
                                    #= MEBT Quads
                                    ['Q_D1057_D1060_D1062', 'scl_fac0', 13.05],
                                    ['Q_D1057_D1060_D1062', 'scl_fac1', -14.66],
                                    ['Q_D1057_D1060_D1062', 'scl_fac2', 14.26],
                                    ['Q_D1074_D1076_D1078', 'scl_fac0', 10.86],
                                    ['Q_D1074_D1076_D1078', 'scl_fac1', -13.68],
                                    ['Q_D1074_D1076_D1078', 'scl_fac2', 12.60],
                                    ['Q_D1095_D1098_D1100', 'scl_fac0', 10.33],
                                    ['Q_D1095_D1098_D1100', 'scl_fac1', -12.99],
                                    ['Q_D1095_D1098_D1100', 'scl_fac2', 10.03],
                                  ]

        self._init_beam(init_beam_arg)


    def _init_beam(self,arg=np.zeros(14)):
        """
        arg: represent beam parameter error. should be bounded [-1,1]. 
            will be normalized such that [-1,1] range corresponds to [-20%,20%] error in unit of cov_norm_ref.
        cov_norm_ref: scaling of the cov matrix. 
                    expected to be equal to (4D emittance)**0.25 
                                        or to np.linalg.det(cov_ref)**0.25 
        """
        beta_ref = 0.005075945
        cov_sqrt_ref = np.array([[ 2.86557295e-01,  6.54452658e-03, -1.34972639e-01,  2.70814331e-04],
                                [ 6.54452658e-03,  5.63021572e-03, -1.74878926e-03, -1.45392852e-03],
                                [-1.34972639e-01, -1.74878926e-03,  6.69979613e-01, -7.85738478e-04],
                                [ 2.70814331e-04, -1.45392852e-03, -7.85738478e-04,  2.11537307e-03]])
        cov_ref = cov_sqrt_ref.dot(cov_sqrt_ref.T)
        
        scale4 = np.array([cov_sqrt_ref[i,i] for i in range(4)])
        cov_sqrt = cov_sqrt_ref.copy()
        for i in range(4):
            for j in range(4):
                cov_sqrt[i,j] /= np.sqrt(scale4[i]*scale4[j])
                
        dcov_sqrt = np.zeros((4,4))
        dcov_sqrt[0,0] = arg[4]
        dcov_sqrt[1,1] = arg[5]
        dcov_sqrt[0,1] = dcov_sqrt[1,0] = arg[6]
        dcov_sqrt[2,2] = arg[7]
        dcov_sqrt[3,3] = arg[8]
        dcov_sqrt[2,3] = dcov_sqrt[3,2] = arg[9]
        dcov_sqrt[0,2] = dcov_sqrt[2,0] = arg[10]
        dcov_sqrt[0,3] = dcov_sqrt[3,0] = arg[11]
        dcov_sqrt[1,2] = dcov_sqrt[2,1] = arg[12]
        dcov_sqrt[1,3] = dcov_sqrt[3,1] = arg[13]
        dcov_sqrt *= 0.2*np.linalg.det(cov_sqrt)**0.25 # 20% error in unit of det[cov_sqrt]
        cov_sqrt += dcov_sqrt
        for i in range(4):
            for j in range(4):
                cov_sqrt[i,j] *= np.sqrt(scale4[i]*scale4[j])
        cov = cov_sqrt.dot(cov_sqrt.T)
        
        scale4 = np.array([cov_ref[i,i]**0.5 for i in range(4)])
        
        mean = 0.2*arg[:4]*scale4 # offset in unit of 20% of rms beam size
        pdata = np.zeros([self._npt,6])
        pdata[:,:4] = np.random.multivariate_normal(mean, cov, size=self._npt)
        pdata[:,4] = np.random.rand(self._npt)*2*np.pi-np.pi
        pdata[:,5] = beta_ref*(1+1e-4*np.random.randn(self._npt))
        
        kutil.write_dist(pdata,beta_ref=0.005075945,fname=self._working_directory+"read_dis.dat")
        self._simulated = False


    def _get_track(self, Win=12e3):
        trk_kwargs = {"npat":self._npt,
                      "Win":Win}
        kutil.adjust_trk_file("track_ref.dat",
                              "track_adj.dat",
                              **trk_kwargs)
        track = TRACKutil.exeClass(lat_file='sclinac_ref.dat', 
                                   trk_file='track_adj.dat', 
                                   fi_in_file='fi_in_ref.dat')
        for name,family,simval in self._name_family_simval:
            track.reconfigure(name,family,simval)
        return track

    def _runTRACK(self,batch=True):
        if self._simulated:
            return
        
        old_path = os.getcwd()
        os.chdir(self._working_directory)
        print("get TRACK")
        track = self._get_track()
        print("run TRACK")
        track.run(batch=batch)
        print("run done, old_path=",old_path)
        self._measured_data = kutil.get_measurements(track,
                                                     keys=["BCM_D0989",
                                                           "BCM_D1055",
                                                           "BPM_D1056",
                                                           "BPM_D1072",
                                                           "BPM_D1094",
                                                           "FC_D1102"])
        os.chdir(old_path)
        print("os.chdir(old_path).....?")
        self._simulated = True
        print("_runTRACK done")
        
        
    def _download_TRACK_files(self,path=None):
        from distutils.dir_util import copy_tree
  
        if path is None:
            from .package_path import package_path
            path = package_path+"/TRACK_data"
        else:
            assert type(path) is str
            if path[-1]=="/":
                path = path[:-1]
                
        from_directory = path
        to_directory = self._working_directory
        copy_tree(from_directory, to_directory)
            
            

        
    def _from_PV_to_name_family_Dnum(self,PVname):
        name = re.findall("\w+_D\d\d\d\d",PVname)[0]
        Dnum = int(re.findall("D\d\d\d\d",name)[0][1:])
        if 'DCH_' in name:
            name = name.replace("DCH_","DCHV_")
            family = "tm_xkick"
        elif 'DCV_' in name:
            name = name.replace("DCV_","DCHV_")
            family = "tm_ykick"
        
        else:
            family = re.findall("_D\d\d\d\d:\w+",PVname)[0][7:]

        return name,family,Dnum
            

    def _conversion_from_TRACK_to_PV(self,PVname):
        name,family,Dnum = self._from_PV_to_name_family_Dnum(PVname)
        
        for i in range(len(self._name_family_simval)):
            if name == self._name_family_simval[i][0] and family==self._name_family_simval[i][1]:
                simval = self._name_family_simval[i][2]
                flag_exist = True
                break
        if not flag_exist:
            raise ValueError("PVname not represent in TRACK simulation configure: ",PVname)

        if family=="tm_xkick":
            if Dnum in [709, 723, 773, 840, 868, 880, 901, 929]:
                return -simval/(3.484E-4*self._scale)
            elif Dnum in [790, 805, 821, 948, 964, 979, 992]:
                return -simval/(2.598E-4*self._scale)
        elif family =="tm_ykick":
            if Dnum in [709, 723, 773, 840, 868, 880, 901, 929]:
                return simval/(3.484E-4*self._scale)
            elif Dnum in [790, 805, 821, 948, 964, 979, 992]:
                return simval/(2.598E-4*self._scale)
        else:
            raise ValueError("not yet implemted name or family: ",name,family)


    def _conversion_from_PV_to_TRACK(self,PVname,PVval):
        name,family,Dnum = self._from_PV_to_name_family_Dnum(PVname)
        simval = PVval

        if family=="tm_xkick":
            if Dnum in [709, 723, 773, 840, 868, 880, 901, 929]:
                simval =  -3.484E-4*self._scale*PVval
            elif Dnum in [790, 805, 821, 948, 964, 979, 992]:
                simval =  -2.598E-4*self._scale*PVval
        elif family =="tm_ykick":
            if Dnum in [709, 723, 773, 840, 868, 880, 901, 929]:
                simval =  3.484E-4*self._scale*PVval
            elif Dnum in [790, 805, 821, 948, 964, 979, 992]:
                simval =  2.598E-4*self._scale*PVval
        else:
            raise ValueError("not yet implemted PV: ",PVname)
        
        return name,family,simval


    def caput(self,PVname,PVval):
        name,family,simval = self._conversion_from_PV_to_TRACK(PVname,PVval)
        
        for i in range(len(self._name_family_simval)):
            if name == self._name_family_simval[i][0] and family==self._name_family_simval[i][1]:
                self._name_family_simval[i][2] = simval
                flag_exist = True
                break
        if not flag_exist:
            self._name_family_simval.append([name,family,simval])
        
        self._simulated = False


    def caget(self,PVname):
        print("in caget")
        self._runTRACK()
        print("caget time",time.time()-t0)
        name,family,Dnum = self._from_PV_to_name_family_Dnum(PVname)
        if "BPM_" in name:
            if family == "XPOS_RD":
                return self._measured_data.loc[name].xcen*10
            elif family == "YPOS_RD":
                return self._measured_data.loc[name].ycen*10
            elif family == "PHASE_RD":
                if Dnum == 1056:
                    offset = 79.0
                elif Dnum == 1072:
                    offset = -26.7 -85.7
                elif Dnum == 1094:
                    offset = -19.4 -55.0
                else:
                    offset = 0
                return self._measured_data.loc[name].zcen + offset
            else:
                raise ValueError("not yet implemented PV for measurement: ",PVname)
                
        elif "BCM_" in name:
            if family == "AVGPK_RD":
                return self._measured_data.loc[name].np/self._npt*self._current
            elif family == "AVG_RD":
                return self._measured_data.loc[name].np/self._npt*self._current*self._duty
            else:
                raise ValueError("not yet implemented PV for measurement: ",PVname)
                
        elif "FC_" in name:
            if family == "PKAVG_RD":
                return self._measured_data.loc[name].np/self._npt*self._current
            elif family == "AVG_RD":
                return self._measured_data.loc[name].np/self._npt*self._current*self._duty
            else:
                raise ValueError("not yet implemented PV for measurement: ",PVname)
                
        elif "HVP_" in name:
            if family == "I_RD":
                return self._current
            else:
                raise ValueError("not yet implemented PV for measurement: ",PVname)

        else:
            try:
                return self._conversion_from_TRACK_to_PV(PVname)
            except:
                raise ValueError("not yet implemented PV for measurement: ",PVname)    
                
        print("caget done",time.time()-t0)
