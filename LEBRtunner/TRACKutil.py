import os
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import string
import numpy as np
from scipy.linalg import cholesky
import time
import shutil
from scipy.optimize import minimize
from IPython.display import display, clear_output
import subprocess
import errno


# plt.rcParams.update({
#     'font.size':20,
#     'font.family':'Times New Roman'
# })

class exeClass:
    def __init__(self, lat_file = 'sclinac.dat', trk_file = 'track.dat', fi_in_file=''):
        self.lat_file = lat_file
        self.trk_file = trk_file
        self.fi_in_file = fi_in_file
        
        #self.execute(lat_file)
        if os.path.exists(self.lat_file):
            f = open(lat_file,'r')
            self.lat = f.readlines()
            f.close()
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.lat_file)
        if os.path.exists(self.trk_file):
            f = open(trk_file,'r')
            self.trk = f.readlines()
            f.close()
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.trk_file)
    #################
    def run(self, eelem='', moment0 = [], batch = False, wait=False):
        self.execute(eelem, moment0, batch)
        return
    #################
    def execute(self, eelem='', moment0 = [], batch = False, wait=False):
        os.system('rm sclinac.dat')
        os.system('rm track.dat')
        os.system('rm beam.out')
        os.system('rm step.out')

        lat = self.lat
        if eelem !='':
            for (il,l) in enumerate(lat):
                if eelem in l:
                    lat.insert(il+1, '0 stop\n')
                    break

        f = open('sclinac.dat','w')
        f.writelines(lat)
        f.close()
        f = open('track.dat','w')
        f.writelines(self.trk)
        f.close()
        if len(self.fi_in_file)>0:
            shutil.copy(self.fi_in_file, 'fi_in.dat')

        if wait:
            os.system('track')
        else:
            p = subprocess.Popen(['track'])

        if batch == False:
            f1 = plt.figure(1, figsize=(12, 8))
            ht = self.draw_lattice(fig=f1)
            ht[0].legend(ncol=2)

        b = None
        while b != False:
            if batch == False:
                r1= self.get_data_step()
                ht[0].cla()
                ht[0].plot(r1['pos'],r1['xrms'],'b-', label = 'xrms')
                ht[0].plot(r1['pos'],r1['yrms'],'r-', label = 'yrms')
                ht[0].set_xlim(ht[-1].get_xlim())
                #ht[0].set_ticklabels( () )
                ht[0].xaxis.set_tick_params(labelsize=0)
                display(f1)
                clear_output(wait = True)
            #print(p.poll())
            #print(str(p.stdin)+'\t'+str(p.stdout)+'\r'+str(p.stderr))
            #plt.pause(0.5)
            time.sleep(0.5)
            b = p.poll()
        if batch == False:
            r1= self.get_data_step()
            ht[0].cla()
            ht[0].plot(r1['pos'],r1['xrms'],'b-', label = 'xrms')
            ht[0].plot(r1['pos'],r1['yrms'],'r-', label = 'yrms')
            ht[0].set_xlim(ht[-1].get_xlim())
            #ht[0].set_ticklabels( () )
            ht[0].xaxis.set_tick_params(labelsize=0)
            display(f1)
            clear_output(wait = True)
        
        time.sleep(1)

        '''
        fn = 'step.out'
        b = os.path.exists(fn)
        if  b == False:
            os.system('/usr/bin/track >/dev/null 2>&1')
            time.sleep(5)
        s = os.path.getsize(fn)
        while  s == 0:
            os.system('/usr/bin/track >/dev/null 2>&1')            
            time.sleep(5)

        fn = 'beam.out'
        b = os.path.exists(fn)
        while  b == False:
            os.system('/usr/bin/track >/dev/null 2>&1')            
            time.sleep(5)
        s = os.path.getsize(fn)
        while  s == 0:
            os.system('/usr/bin/track >/dev/null 2>&1')            
            time.sleep(5)
        '''


    #################
    def init_lattice(self, lat_file=''):
        if lat_file == '':
            lat_file  = self.lat_file
        f = open(lat_file,'r')
        self.lat = f.readlines()
        f.close()
    #################
    def get_lattice_info(self, lat_file = ''):
        if lat_file == '':
            lat_file = self.lat_file
        z = 0.0
        f = open(lat_file)
        a = f.readlines()
        f.close()
        
        zstart = 0.0
        for l in a:
            l = l.lstrip().rstrip()
            elem = l.split()
            if (len(l)==0) or (l[0] == '#') or (l[0] == '!') or (len(elem)==0):
                continue

            #elem = l.replace(' ', '').replace(';','').split(',')
            #print elem
            '''
            if ('L' in l) and ('aper' in l):
                L_loc_No = 0
                for i in range(len(elem)):
                    if 'L=' in elem[i]:
                        L_loc_No = i
                        break
                length = float(elem[L_loc_No][2:])
            else:
                length = 0.0
            '''
            
            zstart = z
            if ('updat' in l) or ('matrx' in l) or ('scrch' in l) or ('stop' in l):
                pass
            elif  'drift' in l:      #Drift
                length = float(elem[2])/100.
                z= zstart+ length
                #print str(z)+'\t'+l,
            elif  'cav' in l:  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                amp = elem[4]
                print('cav\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.)+'\tamp: '+str(amp))
                
            elif  'rfq' in l:  #Cavity
                length = float(elem[3])/100.
                z= zstart+ length
                print('rfq\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  'mhb4' in l:  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                print(elem[1]+'\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  'mhb' in l:  #Cavity
                length = float(elem[4])/100.
                z= zstart+ length
                print('mhb\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  ('extrc' in l):  #Cavity
                length = float(elem[3])/100.
                z= zstart+ length
                print('extrc\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  ('bunch' in l):  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                print('bunch\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  'sol' in l:  #Solenoid
                length = float(elem[3])/100.
                z= zstart+ length
                B = float(elem[2])*1e-4
                print(elem[1]+'\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.)+'\tB: '+str(B))
            elif  ('hbmag' in l): #M dipole
                length = float(elem[3])*math.pi*2*abs(float(elem[2]))/360./100.
                z= zstart+ length
                print('hbmag\t'+str(zstart)+'\t'+str(z))
            elif  ('bmag' in l) or ('dipo' in l):
                length = float(elem[2])/100.
                z= zstart+ length
                print('bmag\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif ('edp3d' in l): #dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                z= zstart+ length
                print('edp3d\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  ('mag3d' in l): #dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                z= zstart+ length
                print('mag3d\t'+str(zstart)+'\t'+str(z)+'\t'+str((zstart+z)/2.))
            elif  ('mult' in l): #M dipole
                length = float(elem[2])/100.
                length_eff = float(elem[3])/100.
                length_dif = length-length_eff
                sw = int(elem[8])
                z= zstart+ length
                Ra = float(elem[7])
                B2 = float(elem[4])*1e-4/(Ra*1e-2)
                B3 = float(elem[5])*1e-4/(Ra*Ra*1e-4)
                B4 = float(elem[6])*1e-4/(Ra*Ra*Ra*1e-6)
                #print(str(length)+'\t'+str(length_eff))
                if sw == 0:
                    print('mult\t'+str(zstart+length_dif/2.)+'\t'+str(zstart+length_dif/2.+length_eff)+'\tB2: '+str(B2)+'\tB3: '+str(B3)+'\tB4: '+str(B4))
                elif sw == 1: # 1st half
                    print('mult\t'+str(zstart+length_dif)+'\t'+str(zstart+length_eff)+'\tB2: '+str(B2)+'\tB3: '+str(B3)+'\tB4: '+str(B4))
                elif sw == 2: # 1st half
                    print('mult\t'+str(zstart)+'\t'+str(zstart+length_eff)+'\tB2: '+str(B2)+'\tB3: '+str(B3)+'\tB4: '+str(B4))
                    
                    #print str(z)+'\t'+l,
            elif ('quad' in l): # magnetic quad
                length = float(elem[3])/100.
                length_eff = float(elem[4])/100.
                length_dif = length-length_eff
                z= zstart+ length
                B=float(elem[2])*1e-4
                Ra = float(elem[5])*1e-2
                B2 = B/Ra
                print('quad\t'+str(zstart+length_dif/2.)+'\t'+str(z-length_dif/2.)+'\tB2: '+str(B2))
            elif ('mq3d' in l): # magnetic quad
                length = float(elem[2])/100.
                z= zstart+ length
                B=float(elem[3])*1e-4
                Ra = float(elem[4])*1e-2
                print('mq3d\t'+str(zstart)+'\t'+str(z)+'\tB2: '+str(B/Ra))
            elif ('eq3d' in l) or ('equad' in l):  #eQuad
                length = float(elem[3])/100.
                z= zstart+ length
                V = float(elem[2])
                print('eq3d\t'+str(zstart)+'\t'+str(z)+'\tV: '+str(V))
            elif ('bpm' in l): #FLAME
                z= zstart
            elif ('marker' in l): #FLAME
                z= zstart
                pass
            elif ('zero' in l) or ('corr' in l): #Corrector
                z= zstart
            elif ('slit' in l): # collimator slit
                z= zstart
                xgap = float(elem[2])*2
                ygap = float(elem[3])*2
                print('slit: '+str(z)+'\txgap: '+str(xgap)+'\tygap: '+str(ygap))
            elif ('mon' in l) or ('MON' in l):
                print('mon\t'+str(z))
            else:
                print(l),
        print('end\t'+str(z))

        return 
            
    #################
    def reconfigure(self, ename, pname, val):
        n = len(self.lat)
        ii = 0
        lat_tmp = []
        bDone = False
        while ii < n:
            l = self.lat[ii]
            lat_tmp.append(l)
            #print ename+'\t'+l,
            if ename in l:
                ii+=1
                l = self.lat[ii]
                while l[0] == '!':
                    lat_tmp.append(l)
                    ii += 1
                    l = self.lat[ii]
                
                el = l.split()
                
                if 'mhb4' in l:
                    if pname == 'scl_fac0':
                        el[6] = str(val)
                    elif pname == 'phi0':
                        el[4] = str(val)
                    elif pname == 'scl_fac1':
                        el[9] = str(val)
                    elif pname == 'phi1':
                        el[7] = str(val)
                    elif pname == 'scl_fac2':
                        el[12] = str(val)
                    elif pname == 'phi2':
                        el[10] = str(val)
                    l = '\t'.join(el)+'\n'
                    lat_tmp.append(l)
                    bDone = True
                elif 'mag3d' in l:
                    if pname == 'angle':
                        el[5] = str(val)
                    l = '\t'.join(el)+'\n'
                    lat_tmp.append(l)
                    bDone = True
                elif 'cav' in l:
                    if pname == 'scl_fac':
                        el[4] = str(val)
                        l = '\t '.join(el)+'\n'                    
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    lat_tmp.append(l)
                        
                    bDone = True
                elif 'sol3c' in l:
                    if pname == 'scl_fac0':
                        el[2] = str(val)
                    elif pname == 'scl_fac1':
                        el[7] = str(val)
                    elif pname == 'scl_fac2':
                        el[8] = str(val)
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    l = '\t '.join(el)+'\n'                    
                    lat_tmp.append(l)
                    bDone = True
                elif 'sol' in l:
                    if pname == 'B':
                        el[2] = str(10000.*val) # Tesla to gauss
                        l = '\t '.join(el)+'\n'                    
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    lat_tmp.append(l)
                    bDone = True
                elif 'quad' in l:
                    if pname == 'B2':
                        el[2] = str(val*float(el[5])*100.)
                        l = '\t '.join(el)+'\n'                    
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    lat_tmp.append(l)
                    #print ename+' '+str(val)+' '+str(el[5])+' '+str(el[2])
                    bDone = True
                elif 'eq3d' in l:
                    if (pname == 'V'): # equad
                        el[2] = str(val/1000.) # V to kV
                        l = '\t '.join(el)+'\n'                    
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    lat_tmp.append(l)
                    #print ename+' '+str(val)+' '+str(el[5])+' '+str(el[2])
                    bDone = True
                elif 'edp3d' in l:
                    if (pname == 'V'): # equad
                        el[2] = str(-val/1000.) # V to kV
                        l = '\t '.join(el)+'\n'                    
                    else:
                        print('unknown: '+ename+'\t'+pname)
                    lat_tmp.append(l)
                    #print ename+' '+str(val)+' '+str(el[5])+' '+str(el[2])
                    bDone = True
                elif 'corr' in l:# steerer in the unit of [rad]
                    if (pname == 'theta_x'):
                        if el[0] == '0':
                            el[7] = str(val*1e+3) # [rad] to [mrad]
                        else:
                            el = ['0','corr','0.0','0.0','0.0','0.0','0.0',str(val*1e+3),'0.0']
                    elif (pname == 'theta_y'):
                        if el[0] == '0':
                            el[8] = str(val*1e+3) # [rad] to [mrad]
                        else:
                            el = ['0','corr','0.0','0.0','0.0','0.0','0.0','0.0',str(val*1e+3)]
                    elif (pname == 'tm_xkick'):
                        if el[0] == '2':
                            el[7] = str(val)
                        else:
                            el = ['2','corr','0.0','0.0','0.0','0.0','0.0',str(val),'0.0']
                    elif (pname == 'tm_ykick'):
                        if el[0] == '2':
                            el[8] = str(val)
                        else:
                            el = ['2','corr','0.0','0.0','0.0','0.0','0.0','0.0',str(val)]
                    else:
                        print('Wrong fieldWrong field: '+str(ename)+'\t'+str(pname))

                    lat_tmp.append('\t '.join(el)+'\n')
                    #print ename+' '+str(val)+' '+str(el[5])+' '+str(el[2])
                    bDone = True

            ii+=1
        if bDone == False:
            print ('No element found: '+ename+'\t'+pname+'\t'+str(val))
                
        self.lat = lat_tmp
                    
                
    #################
    '''
    def reconfigure(self, ename, pname, val):
        n = len(self.lat)
        ii = 0
        lat_tmp = []
        bDone = False
        while ii < n:
            l = self.lat[ii]
            lat_tmp.append(l)
            #print ename+'\t'+l,
            if ename in l:
                if (pname == 'scl_fac'): # cavity
                    #print ename
                    ii+=1
                    l = self.lat[ii]
                    e1 = l.split()
                    e1[4] = str(val)
                    #l = string.join(e1,'\t')+'\n'
                    l = '\t '.join(e1)+'\n'                    
                    lat_tmp.append(l)
                    bDone = True
                elif (pname == 'B'): # Solenoid
                    ii+=1
                    l = self.lat[ii]
                    e1 = l.split()
                    e1[2] = str(10000.*val)
                    #l = string.join(e1,'\t')+'\n'
                    l = '\t '.join(e1)+'\n'                    
                    lat_tmp.append(l)
                    bDone = True
                elif (pname == 'B2'): # Quad
                    ii+=1
                    l = self.lat[ii]
                    e1 = l.split()
                    e1[2] = str(val*float(e1[5])*100.)
                    #l = string.join(e1,'\t')+'\n'
                    l = '\t '.join(e1)+'\n'                    
                    lat_tmp.append(l)
                    #print ename+' '+str(val)+' '+str(e1[5])+' '+str(e1[2])
                    #print l
                    bDone = True
                elif (pname == 'V'): # equad
                    ii+=1
                    l = self.lat[ii]
                    e1 = l.split()
                    e1[2] = str(val/1000.)
                    #l = string.join(e1,'\t')+'\n'
                    l = '\t '.join(e1)+'\n'                    
                    lat_tmp.append(l)
                    #print ename+' '+str(val)+' '+str(e1[5])+' '+str(e1[2])
                    #print l
                    bDone = True
                elif (pname == 'tm_xkick'): # steerer in the unit of T-m
                    ii+=1
                    l = self.lat[ii]
                    if 'sol3c' in l:
                        e1 = l.split()
                        e1[7] = str(val*1e+4)
                        #l = string.join(e1,'\t')+'\n'
                        l = '\t '.join(e1)+'\n'                    
                        lat_tmp.append(l)
                        #print ename+' '+str(val)+' '+str(e1[5])+' '+str(e1[2])
                        #print l
                        bDone = True
                elif (pname == 'tm_ykick'): # steerer in the unit of T-m
                    ii+=1
                    l = self.lat[ii]
                    if 'sol3c' in l:
                        e1 = l.split()
                        e1[8] = str(val*1e+4)
                        #l = string.join(e1,'\t')+'\n'
                        l = '\t '.join(e1)+'\n'                    
                        lat_tmp.append(l)
                        #print ename+' '+str(val)+' '+str(e1[5])+' '+str(e1[2])
                        #print l
                        bDone = True

            ii+=1
        if bDone == False:
            print ('No element found: '+ename)
                
        self.lat = lat_tmp
    '''
                    
                
    #################
    def getdata(self, fn = 'beam.out'):
        b = os.path.exists(fn)
        while  b == False:
            time.sleep(1)
            b = os.path.exists(fn)
        s = os.path.getsize(fn)
        while  s == 0:
            time.sleep(1)
            s = os.path.getsize(fn)
        
            
        res = self.get_data(fn)
        return res
    #################
    def get_data(self, fn = 'beam.out'):
        b = os.path.exists(fn)
        while  b == False:
            time.sleep(1)
            b = os.path.exists(fn)            
        s = os.path.getsize(fn)
        while  s == 0:
            time.sleep(1)
            s = os.path.getsize(fn)
        f = open(fn,'r')
        a = f.readlines()
        f.close()
        
        eidl = []
        typel = []
        sl = []
        Ek_refl = []
        dwwpl = []
        xrmsl,yrmsl = [],[]
        xmaxl,ymaxl = [],[]
        x99p5l,y99p5l = [],[]
        xp99p5l,yp99p5l = [],[]
        phirmsl,phimaxl,phiminl = [],[],[]
        dwwl = []
        xnemitl,xnemit90l,xnemitmaxl = [],[],[]
        ynemitl,ynemit90l,ynemitmaxl = [],[],[]
        znemitl,znemit90l,znemitmaxl = [],[],[]
        xcenl,xpcenl,ycenl,ypcenl,bgl,zcenl = [],[],[],[],[],[]
        xalphal,xbetal,yalphal,ybetal,zalphal,zbetal,nlossl,npl = [],[],[],[],[],[],[],[]

        for l in a[1:]:
            e = l.split()
            eidl.append(int(eval(e[0])))
            typel.append(e[1])
            sl.append(float(e[2])) # [m]
            Ek_refl.append(float(e[3])) # [MeV/u]
            dwwpl.append(float(e[4])) #
            xrmsl.append(float(e[5])*10.) # [mm]
            yrmsl.append(float(e[6])*10.) # [mm]
            xmaxl.append(float(e[7])*10.) # [mm]
            ymaxl.append(float(e[8])*10.) # [mm]
            phirmsl.append(float(e[9])) # [deg]
            phimaxl.append(float(e[10])) # [deg]
            dwwl.append(float(e[11])) # [rel.u.]
            xnemitl.append(float(e[12])*10./4.) # [mm-mrad]
            xnemit90l.append(float(e[13])*10.) # [mm-rad]
            xnemitmaxl.append(float(e[14])*10.) # [mm-rad]
            ynemitl.append(float(e[15])*10./4.) # [mm-rad]
            ynemit90l.append(float(e[16])*10.) # [mm-rad]
            ynemitmaxl.append(float(e[17])*10.) # [mm-rad]
            znemitl.append(float(e[18])/4.) # [kev/u*ns]
            znemit90l.append(float(e[19])) # [kev/u*ns]
            znemitmaxl.append(float(e[20])) # [kev/u*ns]
            xcenl.append(float(e[21])*10.) # [mm]
            xpcenl.append(float(e[22])) # [mrad]
            ycenl.append(float(e[23])*10.) # [mm]
            ypcenl.append(float(e[24])) #[mrad]
            bgl.append(float(e[25]))
            zcenl.append(float(e[26])) #[deg]
            xalphal.append(float(e[27]))
            xbetal.append(float(e[28])*10.) #[mm/mrad]
            yalphal.append(float(e[29]))
            ybetal.append(float(e[30])*10.) #[mm/mrad]
            zalphal.append(float(e[31]))
            zbetal.append(float(e[32])) #[deg/(%ofD_W/W)]
            nlossl.append(float(e[33]))
            npl.append(float(e[34]))


            if float(e[12])> 0.:
                x99p5l.append(float(e[5])*10.*math.sqrt((float(e[13])*10.)/(float(e[12])*10./4.)))
            if float(e[15])> 0.:
                y99p5l.append(float(e[6])*10.*math.sqrt((float(e[16])*10.)/(float(e[15])*10./4.)))

        n = len(sl)
        xpadvl,ypadvl = [0.],[0.]
        xadv_int = 0.
        yadv_int = 0.
        for i in range(1,n):
            ds = 1000.*(sl[i]-sl[i-1]) # [mm]
            #print str(sl[i])+'\t'+str(sl[i-1])+'\t'+str(ds)
            #print str(xbetal[i])+'\t'+str(xbetal[i-1])
            if xbetal[i] > 0. and ybetal[i] > 0.:
                xadv_int += 1./((xbetal[i]+xbetal[i-1])/2.)*ds*180./math.pi/1000.
                yadv_int += 1./((ybetal[i]+ybetal[i-1])/2.)*ds*180./math.pi/1000.
            xpadvl.append(xadv_int)
            ypadvl.append(yadv_int)

            
        res = {}
        res['eid'] = np.array(eidl)
        res['type'] = np.array(typel)
        res['s'] = np.array(sl)
        res['pos'] = np.array(sl)
        res['Ek_ref'] = np.array(Ek_refl)
        res['dwwp'] = np.array(dwwpl)
        res['xrms'] = np.array(xrmsl)
        res['yrms'] = np.array(yrmsl)
        res['xmax'] = np.array(xmaxl)
        res['ymax'] = np.array(ymaxl)
        res['phirms'] = np.array(phirmsl)
        res['phimax'] = np.array(phimaxl)
        res['dww'] = np.array(dwwl)
        res['xnemit'] = np.array(xnemitl)
        res['xnemit90'] = np.array(xnemit90l)
        res['xnemitmax'] = np.array(xnemitmaxl)
        res['ynemit'] = np.array(ynemitl)
        res['ynemit90'] = np.array(ynemit90l)
        res['ynemitmax'] = np.array(ynemitmaxl)
        res['znemit'] = np.array(znemitl)
        res['znemit90'] = np.array(znemit90l)
        res['znemitmax'] = np.array(znemitmaxl)
        res['xcen'] = np.array(xcenl)
        res['xpcen'] = np.array(xpcenl)
        res['ycen'] = np.array(ycenl)
        res['ypcen'] = np.array(ypcenl)
        res['bg'] = np.array(bgl)
        res['zcen'] = np.array(zcenl)
        res['xalpha'] = np.array(xalphal)
        res['xbeta'] = np.array(xbetal)
        res['yalpha'] = np.array(yalphal)
        res['ybeta'] = np.array(ybetal)
        res['zalpha'] = np.array(zalphal)
        res['zbeta'] = np.array(zbetal)
        res['nloss'] = np.array(nlossl)
        res['np'] = np.array(npl)
        res['xpadv'] = np.array(xpadvl)
        res['ypadv'] = np.array(ypadvl)
        res['x99p5'] = np.array(x99p5l)
        res['y99p5'] = np.array(y99p5l)
        
        return res
    #################
    def get_data_step(self, fn = 'step.out', bshow = False):
        m = 931.49432
        b = os.path.exists(fn)
        while  b == False:
            time.sleep(1)
            b = os.path.exists(fn)            
        s = os.path.getsize(fn)
        while  s == 0:
            time.sleep(1)
            s = os.path.getsize(fn)
        
        f = open(fn,'r')
        a = f.readlines()
        f.close()
        
        sl = []
        xcenl,ycenl = [],[]        
        xpcenl,ypcenl = [],[]
        zcenl,Ek_refl = [],[]        
        xrmsl,yrmsl = [],[]
        xprmsl,yprmsl = [],[]
        phirmsl,dwwrmsl = [],[]
        xmaxl,ymaxl = [],[]
        xpmaxl,ypmaxl = [],[]
        phimaxl,dwwmaxl = [],[]
        xnemitl,ynemitl,znemitl = [],[],[]

        for l in a[1:]:
            e = l.split()
            if len(e)<22:
                continue
            sl.append(float(e[0])/100.) # [m]
            xcenl.append(float(e[1])*10.) # [mm]
            xpcenl.append(float(e[2])) # [mrad]
            ycenl.append(float(e[3])*10.) # [mm]
            ypcenl.append(float(e[4])) #[mrad]
            zcenl.append(float(e[5])) #[deg]
            Ek_refl.append(float(e[6])) # [MeV/u]
            xrmsl.append(float(e[7])*10.) # [mrad]
            xprmsl.append(float(e[8])) # [mm]
            yrmsl.append(float(e[9])*10.) # [mm]
            yprmsl.append(float(e[10])) # [mrad]
            phirmsl.append(float(e[11])) # [deg]
            dwwrmsl.append(float(e[12])) # [rel.u.]
            xmaxl.append(float(e[13])*10.) # [mm]
            xpmaxl.append(float(e[14])*10.) # [mm]
            ymaxl.append(float(e[15])*10.) # [mm]
            ypmaxl.append(float(e[16])*10.) # [mm]
            phimaxl.append(float(e[17])) # [deg]
            dwwmaxl.append(float(e[18])) # [rel.u.]
            xnemitl.append(float(e[19])*10.) # [mm-mrad]
            ynemitl.append(float(e[20])*10.) # [mm-rad]
            znemitl.append(float(e[21])) # [kev/u*ns]

        n = len(sl)
        xpadvl,ypadvl = [0.],[0.]
        bgl = []
        xemitl,xalphal,xbetal,xgammal = [],[],[],[]
        yemitl,yalphal,ybetal,ygammal = [],[],[],[]
        xadv_int = 0.
        yadv_int = 0.

        ## Calculate initial parameters
        i = 0
        g = (m+Ek_refl[i])/m
        b = math.sqrt(1.-1./g/g)
        bg0 = b*g
        xemit0 = xnemitl[i]/bg0
        if bshow:
            print(xnemitl[i], end = '\t')
        if xemit0>0.:
            xbeta0 = xrmsl[i]*xrmsl[i]/xemit0
            if bshow:
                print(xbeta0, end = '\t')
            xgamma0 = xprmsl[i]*xprmsl[i]/xemit0
            if bshow:
                print(xgamma0, end = '\t')
            if xbeta0*xgamma0-1.0>0.:
                xalpha0 = math.sqrt(xbeta0*xgamma0-1.0)
            else:
                xalpha0 = 0.0
                if bshow:
                    print(xalpha0, end = '\n')
        else:
            xbeta0 = 1e-15
            xgamma0 = 1e-15
            xalpha0 = 1e-15
        yemit0 = ynemitl[i]/bg0
        if yemit0>0.:
            ybeta0 = yrmsl[i]*yrmsl[i]/yemit0
            ygamma0 = yprmsl[i]*yprmsl[i]/yemit0
            if ybeta0*ygamma0-1.0>0.:
                yalpha0 = math.sqrt(ybeta0*ygamma0-1.0)
            else:
                yalpha0 = 1e-15                
        else:
            ybeta0 = 1e-15
            ygamma0 = 1e-15
            yalpha0 = 1e-15
        bgl = [bg0]
        xemitl = [xemit0]
        xalphal = [xalpha0]
        xbetal = [xbeta0]
        xgammal = [xgamma0]
        yemitl = [yemit0]
        yalphal = [yalpha0]
        ybetal = [ybeta0]
        ygammal = [ygamma0]
        xpadvl = [xadv_int]
        ypadvl = [yadv_int]
        
        for i in range(1,n):
            g = (m+Ek_refl[i])/m
            b = math.sqrt(1.-1./g/g)
            bg1 = b*g
            xemit1 = xnemitl[i]/bg1
            if xemit1>0.:
                xbeta1 = xrmsl[i]*xrmsl[i]/xemit1
                xgamma1 = xprmsl[i]*xprmsl[i]/xemit1
                xalpha1 = math.sqrt(math.fabs(xbeta1*xgamma1-1.0))
            else:
                xbeta1 = 1e-15
                xgamma1 = 1e-15
                xalpha1 = 1e-15
            yemit1 = ynemitl[i]/bg1
            if yemit1>0.:
                ybeta1 = yrmsl[i]*yrmsl[i]/yemit1
                ygamma1 = yprmsl[i]*yprmsl[i]/yemit1
                #print('ybeta1 '+str(ybeta1)+', ygamma1 '+str(ygamma1)+', yemit1 '+str(yemit1)+', bg '+str(bg1))
                yalpha1 = math.sqrt(math.fabs(ybeta1*ygamma1-1.0))
            else:
                ybeta1 = 1e-15
                ygamma1 = 1e-15
                yalpha1 = 1e-15
            #print str(ybeta1)+'\t'+str(ygamma1)
            
            ds = 1000.*(sl[i]-sl[i-1]) # [mm]

            xadv_int += 1./((xbeta1+xbeta0)/2.)*ds*180./math.pi/1000.
            yadv_int += 1./((ybeta1+ybeta0)/2.)*ds*180./math.pi/1000.

            xalpha0 = xalpha1
            xbeta0 = xbeta1
            xgamma0 = xgamma1
            xemit0 = xemit1
            yalpha0 = yalpha1
            ybeta0 = ybeta1
            ygamma0 = ygamma1
            yemit0 = yemit1
            
            bgl.append(bg1)
            xemitl.append(xemit1)
            xalphal.append(xalpha1)
            xbetal.append(xbeta1)
            xgammal.append(xgamma1)
            yemitl.append(yemit1)
            yalphal.append(yalpha1)
            ybetal.append(ybeta1)
            ygammal.append(ygamma1)
            xpadvl.append(xadv_int)
            ypadvl.append(yadv_int)
        '''
        n = len(sl)
        xpadvl,ypadvl = [0.],[0.]
        xadv_int = 0.
        yadv_int = 0.
        for i in range(1,n):
            ds = 1000.*(sl[i]-sl[i-1]) # [mm]
            #print str(sl[i])+'\t'+str(sl[i-1])+'\t'+str(ds)
            #print str(xbetal[i])+'\t'+str(xbetal[i-1])
            xadv_int += 1./((xbetal[i]+xbetal[i-1])/2.)*ds*180./math.pi/1000.
            yadv_int += 1./((ybetal[i]+ybetal[i-1])/2.)*ds*180./math.pi/1000.
            xpadvl.append(xadv_int)
            ypadvl.append(yadv_int)
        '''

            
        res = {}
        res['s'] = np.array(sl)
        res['pos'] = np.array(sl)
        res['xcen'] = np.array(xcenl)
        res['xpcen'] = np.array(xpcenl)
        res['ycen'] = np.array(ycenl)
        res['ypcen'] = np.array(ypcenl)
        res['zcen'] = np.array(zcenl)
        res['Ek_ref'] = np.array(Ek_refl)
        res['xrms'] = np.array(xrmsl)
        res['yrms'] = np.array(yrmsl)
        res['xprms'] = np.array(xprmsl)
        res['yprms'] = np.array(yprmsl)
        res['phirms'] = np.array(phirmsl)
        res['dwwrms'] = np.array(dwwrmsl)
        res['xmax'] = np.array(xmaxl)
        res['xpmax'] = np.array(xpmaxl)
        res['ymax'] = np.array(ymaxl)
        res['ypmax'] = np.array(ypmaxl)
        res['phimax'] = np.array(phimaxl)
        res['dwwmax'] = np.array(dwwmaxl)
        res['xnemit'] = np.array(xnemitl)
        res['ynemit'] = np.array(ynemitl)
        res['znemit'] = np.array(znemitl)
        res['bg'] = np.array(bgl)
        res['xalpha'] = np.array(xalphal)
        res['xbeta'] = np.array(xbetal)
        res['xgamma'] = np.array(xgammal)
        res['yalpha'] = np.array(yalphal)
        res['ybeta'] = np.array(ybetal)
        res['ygamma'] = np.array(ygammal)
        res['xpadv'] = np.array(xpadvl)
        res['ypadv'] = np.array(ypadvl)


        return res
    #################
    def get_aperture(self, fn = '', soff = 0.0):
        if len(fn)==0:
            fn = self.lat_file
        print('lattice file: '+fn)
        f = open(fn,'r')
        a =  f.readlines()
        f.close()

        sl = []
        xl = []
        yl = []
        
        s = 0.0+soff
        for l in a:
            if l[0] == '!':
                continue
            #print(l)
            e = l.split()
            if 'drift' in l:
                ds = float(e[2])*0.01 # [m]
                xt = math.fabs(float(e[3]))*10 # [mm]
                yt = math.fabs(float(e[4]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

            elif 'cav' in l:
                ds = float(e[2])*0.01 # [m]
                xt = 20 # [mm]
                yt = 20 # [mm]
                
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)
                
            elif 'quad' in l:
                ds = float(e[3])*0.01 # [m]
                xt = math.fabs(float(e[5]))*10 # [mm]
                yt = math.fabs(float(e[5]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)
            elif 'bmag' in l:
                ds = float(e[2])*0.01 # [m]
                xt = float(e[6])*10
                yt = float(e[5])*10

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)
            elif 'eq3d' in l:
                ds = float(e[3])*0.01 # [m]
                xt = math.fabs(float(e[4]))*10 # [mm]
                yt = math.fabs(float(e[4]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

            elif 'edp3d' in l:
                if len(elem) == 8:
                    ds = float(e[2])/100. # [m]
                else:
                    ds = float(e[3])/100. # [m]
                
                xt = math.fabs(float(e[6]))*10/2. # [mm]
                yt = math.fabs(float(e[5]))*10/2. # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

            elif 'sol' in l:
                ds = float(e[3])*0.01 # [m]
                xt = math.fabs(float(e[4]))*10 # [mm]
                yt = math.fabs(float(e[4]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

            elif 'mult' in l:
                ds = float(e[2])*0.01 # [m]
                xt = math.fabs(float(e[7]))*10 # [mm]
                yt = math.fabs(float(e[7]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)

                s += ds
                sl.append(s)
                xl.append(xt)
                yl.append(yt)

            elif 'slit' in l:
                xt = math.fabs(float(e[2]))*10 # [mm]
                yt = math.fabs(float(e[3]))*10 # [mm]

                sl.append(s)
                xl.append(xt)
                yl.append(yt)
            elif 'shrt' in l:
                pass
            elif 'corr' in l:
                pass
            elif 'zero' in l:
                pass
            elif 'stop' in l:
                break
            else:
                e = l.split()
                if len(e)>1:
                    print('Undefined element: '+l)
            
        res = {
            'pos':np.array(sl),
            'x':np.array(xl),
            'y':np.array(yl)
        }
        return res
    #################
    def get_dist(self,fn = 'read_dis.out',fntrk = '',bshow=False, blost = False, allcs = False, bnumpy = True, ics = -1):
        npatl,qql,amassl = [],[],[]
        nqtot = 0

        b = os.path.exists(fn)
        while  b == False:
            time.sleep(1)
            b = os.path.exists(fn)            
        s = os.path.getsize(fn)
        while  s == 0:
            time.sleep(1)
            s = os.path.getsize(fn)
        if len(fntrk)>0:
            trk_file = fntrk
        else:
            trk_file = self.trk_file
        
        #= Read track.dat
        if len(trk_file) > 0:
            f = open(trk_file,'r')
            a = f.readlines()
            f.close()
            for l in a:
                l = l.replace(' ','').replace('\t','').replace('\n','')
                if '!' in l:
                    continue
                if 'nqtot=' in l:
                    e=l.replace('nqtot=','')
                    print(e)
                    nqtot = int(e)
                    print('nqtot: '+str(nqtot))
                elif 'npat=' in l:
                    e=l.replace('npat=','').split(',')
                    for i in range(nqtot):
                        npatl.append(int(e[i]))
                    #print('npat: '+str(npatl))
                elif 'amass=' in l:
                    e=l.replace('amass=','').split(',')
                    for i in range(nqtot):
                        amassl.append(float(e[i]))
                    #print('amass: '+str(amassl))
                elif ('qq=' in l) and ('flag' not in l):
                    e=l.replace('qq=','').split(',')
                    for i in range(nqtot):
                        qql.append(float(e[i]))
                    #print('qq: '+str(qql))
                elif ('freqb=' in l):
                    e=l.replace('freqb=','').replace('d','e').replace('D','E')
                    print(e)
                    freq = float(e)
                
        
        f = open(fn,'r')
        a = f.readlines()
        f.close()

        ioff = 0
        E = float(a[ioff].split()[0])
        ncs = int(a[ioff].split()[1])

        nparticle,csl = [],[]
        while len(nparticle)<ncs:
            ioff += 1
            for e in a[ioff].split():
                nparticle.append(int(e))
        while len(csl)<ncs:
            ioff += 1
            for e in a[ioff].split():
                csl.append(float(e))
        m0 = 931.494
        ioff += 1
        
        xl,xpl,yl,ypl,phil = [],[],[],[],[]
        xcenl,xpcenl,ycenl,ypcenl = [],[],[],[]
        xrmsl,xprmsl,yrmsl,yprmsl = [],[],[],[]
        phirmsl,Ermsl = [],[]
        phicenl,Ecenl = [],[]
        xemitl,xnemitl,xalphal,xbetal,xgammal = [],[],[],[],[]
        yemitl,ynemitl,yalphal,ybetal,ygammal = [],[],[],[],[]
        zemitl,zalphal,zbetal,zgammal = [],[],[],[]
        moment1l,emit4dl,nemit4dl = [],[],[]
        cxyl,cxpyl,cxypl,cxpypl = [],[],[],[]
        bg_refl = []
        flagl = []
        El,bgl = [],[]
        npart0l,npartl=[],[]
        #
        xall,xpall,yall,ypall,phiall,bgall,Eall,flagall = [],[],[],[],[],[],[],[]
        m = m0
        itot = 0
        for iq in range(ncs):
            istrt = ioff+sum(npatl[0:iq])+iq
            iend = ioff+sum(npatl[0:iq+1])+iq+1
            print(str(istrt)+'\t'+str(iend))
            xtl,xptl,ytl,yptl,phitl,Etl,flagtl = [],[],[],[],[],[],[]
            bgtl = []
            bg_ref = -1.0
            #print(istrt)
            for il,l in enumerate(a[istrt:iend]):
                e = l.replace(',',' ').replace('  ',' ').split()
                flag = int(e[6])
                if (blost == False) and (flag < 0):
                    continue
                #print(e)
                xtl.append(eval(e[0])*10) # [mm]
                xptl.append(eval(e[1])*1000) # [mrad]
                ytl.append(eval(e[2])*10) # [mm]
                yptl.append(eval(e[3])*1000) # [mrad]
                phitl.append(eval(e[4])*180./math.pi) # [deg]
                bg = eval(e[5])
                bgtl.append(bg)
                if il == 0:
                    bg_ref= bg
                p = bg*m
                E = math.sqrt(m*m+p*p)-m
                Etl.append(E)
                flagtl.append(flag)
                #
                xall.append(eval(e[0])*10) # [mm]
                xpall.append(eval(e[1])*1000) # [mrad]
                yall.append(eval(e[2])*10) # [mm]
                ypall.append(eval(e[3])*1000) # [mrad]
                phiall.append(eval(e[4])*180./math.pi) # [deg]
                Eall.append(E)
                flagall.append(flag)
                bgall.append(bg)
            #####
            npart0 = len(xtl)
            c11,c12,c22,c33,c34,c44 = 0.0,0.0,0.0,0.0,0.0,0.0
            c13,c14,c23,c24 = 0.0,0.0,0.0,0.0
            c55,c56,c66 = 0.0,0.0,0.0

            if npart0 > 0:
                xcen,xpcen = sum(xtl)/npart0,sum(xptl)/npart0
                ycen,ypcen = sum(ytl)/npart0,sum(yptl)/npart0
                phicen,Ecen = sum(phitl)/npart0,sum(Etl)/npart0
            else:
                xcen,xpcen = 0.0, 0.0
                ycen,ypcen = 0.0, 0.0
                phicen,Ecen = 0.0, 0.0
            npart = 0
            for ipar in range(npart0):
                if flagtl[ipar] != 0:
                    #print('Something wrong to calculate a beam matrix')
                    continue
                c11 += math.pow(xtl[ipar]-xcen,2)
                c22 += math.pow(xptl[ipar]-xpcen,2)
                c12 += (xtl[ipar]-xcen)*(xptl[ipar]-xpcen)
                c33 += math.pow(ytl[ipar]-ycen,2)
                c44 += math.pow(yptl[ipar]-ypcen,2)
                c34 += (ytl[ipar]-ycen)*(yptl[ipar]-ypcen)
                c13 += (xtl[ipar]-xcen)*(ytl[ipar]-ycen)
                c14 += (xtl[ipar]-xcen)*(yptl[ipar]-ypcen)
                c23 += (xptl[ipar]-xpcen)*(ytl[ipar]-ycen)
                c24 += (xptl[ipar]-xpcen)*(yptl[ipar]-ypcen)
                c55 += math.pow(phitl[ipar]-phicen,2)
                c66 += math.pow(Etl[ipar]-Ecen,2)
                c56 += (phitl[ipar]-phicen)*(Etl[ipar]-Ecen)
                npart += 1
                #print(str(npart)+'/'+str(npart0)+'\t'+str(npart/npart0*100)+' %')
                if itot%100000 == 1:
                    print('.',end = '')
                itot+=1
            if npart>0:
                c11 /= npart
                c22 /= npart
                c12 /= npart
                c33 /= npart
                c44 /= npart
                c34 /= npart
                c13 /= npart
                c14 /= npart
                c23 /= npart
                c24 /= npart
                c55 /= npart
                c56 /= npart
                c66 /= npart

            xrms = math.sqrt(c11)
            xprms = math.sqrt(c22)
            yrms = math.sqrt(c33)
            yprms = math.sqrt(c44)
            phirms = math.sqrt(c55)
            Erms = math.sqrt(c66)
            xemit = c11*c22-c12*c12
            if xemit > 0.0:
                xemit = math.sqrt(xemit)
            yemit = c33*c44-c34*c34                
            if yemit > 0.0:
                yemit = math.sqrt(yemit)
            zemit = c55*c66-c56*c56
            if zemit > 0.0:
                zemit = math.sqrt(zemit)
            xnemit = xemit*bg_ref
            ynemit = yemit*bg_ref

            if (xrms > 0.) and (yrms > 0.) and (npart>0):
                cxy = c13/(xrms*yrms)
                cxpy = c23/(xprms*yrms)
                cxyp = c14/(xrms*yprms)
                cxpyp = c24/(xprms*yprms)
            else:
                cxy = c13
                cxpy = c23
                cxyp = c14
                cxpyp = c24
            try:
                xalpha = -1*c12/xemit
                xbeta = c11/xemit
                xgamma = c22/xemit
                yalpha = -1*c34/yemit
                ybeta = c33/yemit
                ygamma = c44/yemit
            except:
                xalpha = 0.0
                xbeta = 0.0
                xgamma = 0.0
                yalpha = 0.0
                ybeta = 0.0
                ygamma = 0.0
            try:
                zalpha = -1*c56/zemit
                zbeta = c55/zemit
                zgamma = c66/zemit
            except:
                zalpha = 0.0
                zbeta = 0.0
                zgamma = 0.0
            #print(c14)
            #print(c24)
            #print(c34)
            #print(c44)

            moment1 = np.array([
                [c11,c12*1.0e-3,c13,c14*1.0e-3],
                [c12*1.0e-3,c22*1.0e-6,c23*1.0e-3,c24*1.0e-6],
                [c13,c23*1.0e-3,c33,c34*1.0e-3],
                [c14*1.0e-3,c24*1.0e-6,c34*1.0e-3,c44*1.0e-6]        
            ])

            if np.linalg.det(moment1) > 0.:
                emit4d = math.sqrt(np.linalg.det(moment1))
                nemit4d = emit4d*math.pow(bg_ref,2)
            else:
                emit4d = 0.0
                nemit4d = 0.0

            #####
            bg_refl.append(bg_ref)
            xl.append(np.array(xtl))
            xpl.append(np.array(xptl))
            yl.append(np.array(ytl))
            ypl.append(np.array(yptl))
            phil.append(np.array(phitl))
            bgl.append(np.array(bgtl))
            El.append(np.array(Etl))
            flagl.append(flagtl)
            #####
            xcenl.append(xcen)
            xpcenl.append(xpcen)
            ycenl.append(ycen)
            ypcenl.append(ypcen)
            phicenl.append(phicen)
            Ecenl.append(Ecen)
            xrmsl.append(xrms)
            xprmsl.append(xprms)
            yrmsl.append(yrms)
            yprmsl.append(yprms)
            phirmsl.append(phirms)
            Ermsl.append(Erms)
            #####
            xemitl.append(xemit)
            yemitl.append(yemit)
            zemitl.append(zemit)
            xnemitl.append(xnemit)
            ynemitl.append(ynemit)
            xalphal.append(xalpha)
            xbetal.append(xbeta)
            xgammal.append(xgamma)
            yalphal.append(yalpha)
            ybetal.append(ybeta)
            ygammal.append(ygamma)
            zalphal.append(zalpha)
            zbetal.append(zbeta)
            zgammal.append(zgamma)
            #####
            moment1l.append(moment1)
            emit4dl.append(emit4d)
            nemit4dl.append(nemit4d)
            #####
            cxyl.append(cxy)
            cxpyl.append(cxpy)
            cxypl.append(cxyp)
            cxpypl.append(cxpyp)
            #####
            #npart0l.append(npart0)
            npartl.append(npart)

        if allcs:
            npart0 = len(xall)
            c11,c12,c22,c33,c34,c44 = 0.0,0.0,0.0,0.0,0.0,0.0
            c13,c14,c23,c24 = 0.0,0.0,0.0,0.0
            c55,c56,c66 = 0.0,0.0,0.0
            
            xcenall,xpcenall = sum(xall)/len(xall),sum(xpall)/len(xpall)
            ycenall,ypcenall = sum(yall)/len(yall),sum(ypall)/len(ypall)
            phicenall,Ecenall = sum(phiall)/len(phiall),sum(Eall)/len(Eall)
            npart = 0
            for ipar in range(npart0):
                if flagall[ipar] != 0:
                    #print('Something wrong to calculate a beam matrix')
                    continue
                c11 += math.pow(xall[ipar]-xcen,2)
                c22 += math.pow(xpall[ipar]-xpcen,2)
                c12 += (xall[ipar]-xcen)*(xpall[ipar]-xpcen)
                c33 += math.pow(yall[ipar]-ycen,2)
                c44 += math.pow(ypall[ipar]-ypcen,2)
                c34 += (yall[ipar]-ycen)*(ypall[ipar]-ypcen)
                c13 += (xall[ipar]-xcen)*(yall[ipar]-ycen)
                c14 += (xall[ipar]-xcen)*(ypall[ipar]-ypcen)
                c23 += (xpall[ipar]-xpcen)*(yall[ipar]-ycen)
                c24 += (xpall[ipar]-xpcen)*(ypall[ipar]-ypcen)
                c55 += math.pow(phiall[ipar]-phicen,2)
                c66 += math.pow(Eall[ipar]-Ecen,2)
                c56 += (phiall[ipar]-phicen)*(Eall[ipar]-Ecen)
                npart += 1
                #print(str(npart)+'/'+str(npart0)+'\t'+str(npart/npart0*100)+' %')
                if itot%100000 == 1:
                    print('.',end = '')
                itot+=1
            c11 /= npart
            c22 /= npart
            c12 /= npart
            c33 /= npart
            c44 /= npart
            c34 /= npart
            c13 /= npart
            c14 /= npart
            c23 /= npart
            c24 /= npart
            c55 /= npart
            c56 /= npart
            c66 /= npart

            xrmsall = math.sqrt(c11)
            xprmsall = math.sqrt(c22)
            yrmsall = math.sqrt(c33)
            yprmsall = math.sqrt(c44)
            phirmsall = math.sqrt(c55)
            Ermsall = math.sqrt(c66)
            xemitall = math.sqrt(c11*c22-c12*c12)
            yemitall = math.sqrt(c33*c44-c34*c34)
            zemitall = math.sqrt(c55*c66-c56*c56)
            xnemitall = xemit*bg_ref
            ynemitall = yemit*bg_ref

            if (xrmsall > 0.) and (yrmsall > 0.):
                cxyall = c13/(xrmsall*yrmsall)
                cxpyall = c23/(xprmsall*yrmsall)
                cxypall = c14/(xrmsall*yprmsall)
                cxpypall = c24/(xprmsall*yprmsall)
            else:
                cxyall = c13
                cxpyall = c23
                cxypall = c14
                cxpypall = c24
            try:
                xalphaall = -1*c12/xemitall
                xbetaall = c11/xemitall
                xgammaall = c22/xemitall
                yalphaall = -1*c34/yemitall
                ybetaall = c33/yemitall
                ygammaall = c44/yemitall
            except:
                xalphaall = 0.0
                xbetaall = 0.0
                xgammaall = 0.0
                yalphaall = 0.0
                ybetaall = 0.0
                ygammaall = 0.0
            try:
                zalphaall = -1*c56/zemitall
                zbetaall = c55/zemitall
                zgammaall = c66/zemitall
            except:
                zalphaall = 0.0
                zbetaall = 0.0
                zgammaall = 0.0
            #print(c14)
            #print(c24)
            #print(c34)
            #print(c44)

            moment1all = np.array([
                [c11,c12*1.0e-3,c13,c14*1.0e-3],
                [c12*1.0e-3,c22*1.0e-6,c23*1.0e-3,c24*1.0e-6],
                [c13,c23*1.0e-3,c33,c34*1.0e-3],
                [c14*1.0e-3,c24*1.0e-6,c34*1.0e-3,c44*1.0e-6]        
            ])

            if np.linalg.det(moment1) > 0.:
                emit4dall = math.sqrt(np.linalg.det(moment1all))
                nemit4dall = emit4dall*math.pow(bg_ref,2)
            else:
                emit4dall = 0.0
                nemit4dall = 0.0
            
        if ncs == 1: #= single charge simulation
            res = {
                'bg_ref': bg_refl[0],
                'x': np.array(xl[0]),
                'xp': np.array(xpl[0]),
                'y': np.array(yl[0]),
                'yp': np.array(ypl[0]),
                'phi': np.array(phil[0]),
                'bg': np.array(bgl[0]),
                'E': np.array(El[0]),
                'flag':np.array(flagl[0]),
                'qq':qql[0],
                'amass':amassl[0],
                'nqtot':nqtot,
                'npat':npatl[0],
                'xcen':xcenl[0],
                'xpcen':xpcenl[0],
                'ycen':ycenl[0],
                'ypcen':ypcenl[0],
                'xrms':xrmsl[0],
                'xprms':xprmsl[0],
                'yrms':yrmsl[0],
                'yprms':yprmsl[0],
                'phirms':phirmsl[0],
                'Erms':Ermsl[0],
                'xemit':xemitl[0],
                'yemit':yemitl[0],
                'zemit':zemitl[0],
                'xnemit':xnemitl[0],
                'ynemit':ynemitl[0],
                'xalpha':xalphal[0],
                'xbeta':xbetal[0],
                'xgamma':xgammal[0],
                'yalpha':yalphal[0],
                'ybeta':ybetal[0],
                'ygamma':ygammal[0],
                'zalpha':zalphal[0],
                'zbeta':zbetal[0],
                'zgamma':zgammal[0],
                'moment1':moment1l[0],
                'emit4d':emit4dl[0],
                'nemit4d':nemit4dl[0],
                'cxy':cxyl[0],
                'cxpy':cxpyl[0],
                'cxyp':cxypl[0],
                'cxpyp':cxpypl[0],
                'npart0':npatl[0],
                'npart':npartl[0]
            }
        elif allcs: #= calculate all particles
            res = {
                'bg_ref': bg_refl,
                'x': np.array(xall),
                'xp': np.array(xpall),
                'y': np.array(yall),
                'yp': np.array(ypall),
                'phi': np.array(phiall),
                'bg': np.array(bgall),
                'E': np.array(Eall),
                'flag':np.array(flagall),
                'qq':qql,
                'amass':amassl,
                'nqtot':nqtot,
                'npat':npatl,
                'xcen':xcenall,
                'xpcen':xpcenall,
                'ycen':ycenall,
                'ypcen':ypcenall,
                'phicen':phicenall,
                'Ecen':Ecenall,
                'xrms':xrmsall,
                'xprms':xprmsall,
                'yrms':yrmsall,
                'yprms':yprmsall,
                'phirms':phirmsall,
                'Erms':Ermsall,
                'xemit':xemitall,
                'yemit':yemitall,
                'zemit':zemitall,
                'xnemit':xnemitall,
                'ynemit':ynemitall,
                'xalpha':xalphaall,
                'xbeta':xbetaall,
                'xgamma':xgammaall,
                'yalpha':yalphaall,
                'ybeta':ybetaall,
                'ygamma':ygammaall,
                'zalpha':zalphaall,
                'zbeta':zbetaall,
                'zgamma':zgammaall,
                'moment1':moment1all,
                'emit4d':emit4dall,
                'nemit4d':nemit4dall,
                'cxy':cxyall,
                'cxpy':cxpyall,
                'cxyp':cxypall,
                'cxpyp':cxpypall,
                'npart0':sum(npatl),
                'npart':sum(npartl)
            }
        elif ics>=0:
            res = {
                'bg_ref': bg_refl[ics],
                'x': np.array(xl[ics]),
                'xp': np.array(xpl[ics]),
                'y': np.array(yl[ics]),
                'yp': np.array(ypl[ics]),
                'phi': np.array(phil[ics]),
                'bg': np.array(bgl[ics]),
                'E': np.array(El[ics]),
                'flag':np.array(flagl[ics]),
                'qq':qql[ics],
                'amass':amassl[ics],
                'nqtot':nqtot,
                'npat':npatl[ics],
                'xcen':xcenl[ics],
                'xpcen':xpcenl[ics],
                'ycen':ycenl[ics],
                'ypcen':ypcenl[ics],
                'xrms':xrmsl[ics],
                'xprms':xprmsl[ics],
                'yrms':yrmsl[ics],
                'yprms':yprmsl[ics],
                'phirms':phirmsl[ics],
                'Erms':Ermsl[ics],
                'xemit':xemitl[ics],
                'yemit':yemitl[ics],
                'zemit':zemitl[ics],
                'xnemit':xnemitl[ics],
                'ynemit':ynemitl[ics],
                'xalpha':xalphal[ics],
                'xbeta':xbetal[ics],
                'xgamma':xgammal[ics],
                'yalpha':yalphal[ics],
                'ybeta':ybetal[ics],
                'ygamma':ygammal[ics],
                'zalpha':zalphal[ics],
                'zbeta':zbetal[ics],
                'zgamma':zgammal[ics],
                'moment1':moment1l[ics],
                'emit4d':emit4dl[ics],
                'nemit4d':nemit4dl[ics],
                'cxy':cxyl[ics],
                'cxpy':cxpyl[ics],
                'cxyp':cxypl[ics],
                'cxpyp':cxpypl[ics],
                'npart0':npatl[ics],
                'npart':npartl[ics]
            }
            
        elif allcs:
            res = {
                'bg_ref': bg_refl,
                'x': np.array(xall),
                'xp': np.array(xpall),
                'y': np.array(yall),
                'yp': np.array(ypall),
                'phi': np.array(phiall),
                'bg': np.array(bgall),
                'E': np.array(Eall),
                'flag':np.array(flagall),
                'qq':qql,
                'amass':amassl,
                'nqtot':nqtot,
                'npat':npatl,
                'xcen':xcenl,
                'xpcen':xpcenl,
                'ycen':ycenl,
                'ypcen':ypcenl,
                'phicen':phicenl,
                'Ecen':Ecenl,
                'xrms':xrmsl,
                'xprms':xprmsl,
                'yrms':yrmsl,
                'yprms':yprmsl,
                'phirms':phirmsl,
                'Erms':Ermsl,
                'xemit':xemitl,
                'yemit':yemitl,
                'zemit':zemitl,
                'xnemit':xnemitl,
                'ynemit':ynemitl,
                'xalpha':xalphal,
                'xbeta':xbetal,
                'xgamma':xgammal,
                'yalpha':yalphal,
                'ybeta':ybetal,
                'ygamma':ygammal,
                'zalpha':zalphal,
                'zbeta':zbetal,
                'zgamma':zgammal,
                'moment1':moment1l,
                'emit4d':emit4dl,
                'nemit4d':nemit4dl,
                'cxy':cxyl,
                'cxpy':cxpyl,
                'cxyp':cxypl,
                'cxpyp':cxpypl
            }
        else:
            res = {
                'bg_ref': bg_refl,
                'x': np.array(xl),
                'xp': np.array(xpl),
                'y': np.array(yl),
                'yp': np.array(ypl),
                'phi': np.array(phil),
                'bg': np.array(bgl),
                'E': np.array(El),
                'flag':np.array(flagl),
                'qq':qql,
                'amass':amassl,
                'nqtot':nqtot,
                'npat':npatl,
                'xcen':xcenl,
                'xpcen':xpcenl,
                'ycen':ycenl,
                'ypcen':ypcenl,
                'phicen':phicenl,
                'Ecen':Ecenl,
                'xrms':xrmsl,
                'xprms':xprmsl,
                'yrms':yrmsl,
                'yprms':yprmsl,
                'phirms':phirmsl,
                'Erms':Ermsl,
                'xemit':xemitl,
                'yemit':yemitl,
                'zemit':zemitl,
                'xnemit':xnemitl,
                'ynemit':ynemitl,
                'xalpha':xalphal,
                'xbeta':xbetal,
                'xgamma':xgammal,
                'yalpha':yalphal,
                'ybeta':ybetal,
                'ygamma':ygammal,
                'zalpha':zalphal,
                'zbeta':zbetal,
                'zgamma':zgammal,
                'moment1':moment1l,
                'emit4d':emit4dl,
                'nemit4d':nemit4dl,
                'cxy':cxyl,
                'cxpy':cxpyl,
                'cxyp':cxypl,
                'cxpyp':cxpypl,
                'npart0':sum(npatl),
                'npart':sum(npartl)
            }

        if bshow:
            if ncs == 1:
                print('bg_ref:\t'+str(bg_refl[0]))
                print('xcen:\t'+str(xcenl[0]))
                print('xpcen:\t'+str(xpcenl[0]))
                print('ycen:\t'+str(ycenl[0]))
                print('ypcen:\t'+str(ypcenl[0]))
                print('phicen:\t'+str(phicenl[0])+'\t[deg]')
                print('Ecen:\t'+str(Ecenl[0])+'\t[MeV/u]')
                print('xrms:\t'+str(xrmsl[0])+'\t[mm]')
                print('xprms:\t'+str(xprmsl[0])+'\t[mrad]')
                print('yrms:\t'+str(yrmsl[0])+'\t[mm]')
                print('yprms:\t'+str(yprmsl[0])+'\t[mrad]')
                print('phirms:\t'+str(phirmsl[0])+'\t[deg]')
                print('Erms:\t'+str(Ermsl[0])+'\t[MeV/u]')
                print('xemit:\t'+str(xemitl[0])+'\t[mm-mrad]')
                print('yemit:\t'+str(yemitl[0])+'\t[mm-mrad]')
                print('zemit:\t'+str(zemitl[0])+'\t[deg-MeV/u], '+str(zemitl[0]/180*math.pi)+'\t[rad-MeV/u], '+str(zemitl[0]*1000.*1e+9/360/freq)+'\t[ns-keV/u]')
                print('xnemit:\t'+str(xnemitl[0])+'\t[mm-mrad]')
                print('ynemit:\t'+str(ynemitl[0])+'\t[mm-mrad]')
                print('xalpha:\t'+str(xalphal[0]))
                print('xbeta:\t'+str(xbetal[0]))
                print('xgamma:\t'+str(xgammal[0]))
                print('yalpha:\t'+str(yalphal[0]))
                print('ybeta:\t'+str(ybetal[0]))
                print('ygamma:\t'+str(ygammal[0]))
                print('zalpha:\t'+str(zalphal[0]))
                print('zbeta:\t'+str(zbetal[0])+'\t[deg/MeVu], '+str(zbetal[0]*math.pi/180)+'\t[rad/MeVu]')
                print('zgamma:\t'+str(zgammal[0]))
                print('emit4d:\t'+str(emit4dl[0]))
                print('nemit4d:\t'+str(nemit4dl[0]))
                print('cxy:\t'+str(cxyl[0]))
                print('cxpy:\t'+str(cxpyl[0]))
                print('cxyp:\t'+str(cxypl[0]))
                print('cxpyp:\t'+str(cxpypl[0]))
                print('npart0:\t'+str(npatl[0]))
                print('npart:\t'+str(npartl[0]))
                for i1 in range(4):
                    l = ''
                    for i2 in range(4):
                        l+=str(moment1l[0][i1][i2])+'\t'
                    print(l)
            elif allcs:
                print('bg_ref:\t'+str(bg_refl))
                print('xcen:\t'+str(xcenall))
                print('xpcen:\t'+str(xpcenall))
                print('ycen:\t'+str(ycenall))
                print('ypcen:\t'+str(ypcenall))
                print('phicen:\t'+str(phicenall))
                print('Ecen:\t'+str(Ecenall))
                print('xrms:\t'+str(xrmsall))
                print('xprms:\t'+str(xprmsall))
                print('yrms:\t'+str(yrmsall))
                print('yprms:\t'+str(yprmsall))
                print('phirms:\t'+str(phirmsall))
                print('Erms:\t'+str(Ermsall))
                print('xemit:\t'+str(xemitall))
                print('yemit:\t'+str(yemitall))
                print('zemit:\t'+str(np.array(zemitall)*1000.*1e+9/360/freq)+'\t[ns-keV/u]')
                print('zemit:\t'+str(zemitall)+'\t[deg-MeVu]')
                print('xnemit:\t'+str(xnemitall))
                print('ynemit:\t'+str(ynemitall))
                print('xalpha:\t'+str(xalphaall))
                print('xbeta:\t'+str(xbetaall))
                print('xgamma:\t'+str(xgammaall))
                print('yalpha:\t'+str(yalphaall))
                print('ybeta:\t'+str(ybetaall))
                print('ygamma:\t'+str(ygammaall))
                print('zalpha:\t'+str(zalphaall))
                print('zbeta:\t'+str(zbetaall))
                print('zgamma:\t'+str(zgammaall))
                print('emit4d:\t'+str(emit4dall))
                print('nemit4d:\t'+str(nemit4dall))
                print('cxy:\t'+str(cxyall))
                print('cxpy:\t'+str(cxpyall))
                print('cxyp:\t'+str(cxypall))
                print('cxpyp:\t'+str(cxpypall))
                print('npart0:\t'+str(sum(npatl)))
                print('npart:\t'+str(sum(npartl)))
                print('moment1 [mm,mrad] (FLAME coordinate: [mm-rad]):')
                for i1 in range(4):
                    l = ''
                    for i2 in range(4):
                        l+=str(moment1all[i1][i2])+'\t'
                    print(l)
            else:
                print('bg_ref:\t'+str(bg_refl))
                print('xcen:\t'+str(xcenl))
                print('xpcen:\t'+str(xpcenl))
                print('ycen:\t'+str(ycenl))
                print('ypcen:\t'+str(ypcenl))
                print('phicen:\t'+str(phicenl))
                print('Ecen:\t'+str(Ecenl))
                print('xrms:\t'+str(xrmsl))
                print('xprms:\t'+str(xprmsl))
                print('yrms:\t'+str(yrmsl))
                print('yprms:\t'+str(yprmsl))
                print('phirms:\t'+str(phirmsl))
                print('Erms:\t'+str(Ermsl))
                print('xemit:\t'+str(xemitl))
                print('yemit:\t'+str(yemitl))
                print('zemit:\t'+str(np.array(zemitl)*1000.*1e+9/360/freq)+'\t[ns-keV/u]')
                print('zemit:\t'+str(zemitl)+'\t[deg-MeVu]')
                print('xnemit:\t'+str(xnemitl))
                print('ynemit:\t'+str(ynemitl))
                print('xalpha:\t'+str(xalphal))
                print('xbeta:\t'+str(xbetal))
                print('xgamma:\t'+str(xgammal))
                print('yalpha:\t'+str(yalphal))
                print('ybeta:\t'+str(ybetal))
                print('ygamma:\t'+str(ygammal))
                print('zalpha:\t'+str(zalphal))
                print('zbeta:\t'+str(zbetal))
                print('zgamma:\t'+str(zgammal))
                print('emit4d:\t'+str(emit4dl))
                print('nemit4d:\t'+str(nemit4dl))
                print('cxy:\t'+str(cxyl))
                print('cxpy:\t'+str(cxpyl))
                print('cxyp:\t'+str(cxypl))
                print('cxpyp:\t'+str(cxpypl))
                print('npart0:\t'+str(npatl))
                print('npart:\t'+str(npartl))

                print('moment1 [mm,mrad] (FLAME coordinate: [mm-rad]):')
                for ics in range(ncs):
                    for i1 in range(4):
                        l = ''
                        for i2 in range(4):
                            l+=str(moment1l[ics][i1][i2])+',\t'
                        print(l)
                    print('')

        #if bshow:
        #    return
        else:
            return res
    ##########
    def get_profile(self,ename):

        # get s at ename

        n = len(self.lat)
        s = 0
        ielem = 1
        etype = ''
        for i in range(n):
            line = self.lat[i].lstrip(' ')
            elem = line.split()
            if (line[0] == '#') or (line[0] == '!') or (len(elem)==0):
                if ename in line:
                    j = 1
                    while self.lat[i+j][0] == '!':
                        j+=1
                    elem = self.lat[i+j].split()
                    etype = elem[1]
                    break
            else:
                ielem += 1

        # print (ename+':\t'+etype+'\t'+str(ielem)+'\t'+str(s))


        res0 = self.get_data()

        sl = res0['pos']
        #print (str(s))
        #print (str(sl))
        n = len(sl)
        #ielem = 0
        for i in range(n):
            # print('**',i,res0['eid'][i],ielem,res0['type'][i],etype)
            #if math.fabs(sl[ielem]-s)<0.0001:
            if (res0['eid'][i] == ielem) and (res0['type'][i] == etype):
                iline = i
                break

            #ielem += 1

        res = {}
        res['eid'] = res0['eid'][iline]
        res['type'] = res0['type'][iline]
        res['s'] = res0['s'][iline]
        res['pos'] = res0['s'][iline]
        res['Ek_ref'] = res0['Ek_ref'][iline]
        res['dwwp'] = res0['dwwp'][iline]
        res['xrms'] = res0['xrms'][iline]
        res['yrms'] = res0['yrms'][iline]
        res['xmax'] = res0['xmax'][iline]
        res['ymax'] = res0['ymax'][iline]
        res['phirms'] = res0['phirms'][iline]
        res['phimax'] = res0['phimax'][iline]
        res['dww'] = res0['dww'][iline]
        res['xnemit'] = res0['xnemit'][iline]
        res['xnemit90'] = res0['xnemit90'][iline]
        res['xnemitmax'] = res0['xnemitmax'][iline]
        res['ynemit'] = res0['ynemit'][iline]
        res['ynemit90'] = res0['ynemit90'][iline]
        res['ynemitmax'] = res0['ynemitmax'][iline]
        res['znemit'] = res0['znemit'][iline]
        res['znemit90'] = res0['znemit90'][iline]
        res['znemitmax'] = res0['znemitmax'][iline]
        res['xcen'] = res0['xcen'][iline]
        res['xpcen'] = res0['xpcen'][iline]
        res['ycen'] = res0['ycen'][iline]
        res['ypcen'] = res0['ypcen'][iline]
        res['bg'] = res0['bg'][iline]
        res['zcen'] = res0['zcen'][iline]
        res['xalpha'] = res0['xalpha'][iline]
        res['xbeta'] = res0['xbeta'][iline]
        res['yalpha'] = res0['yalpha'][iline]
        res['ybeta'] = res0['ybeta'][iline]
        res['zalpha'] = res0['zalpha'][iline]
        res['zbeta'] = res0['zbeta'][iline]
        res['nloss'] = res0['nloss'][iline]
        res['np'] = res0['np'][iline]
        res['xpadv'] = res0['xpadv'][iline]
        res['ypadv'] = res0['ypadv'][iline]
        
        return res

    ##########
    '''
    def get_profile(self,ename):

        # get s at ename

        n = len(self.lat)
        s = 0
        ielem = 1
        etype = ''
        for i in range(n):
            line = self.lat[i].lstrip(' ')
            elem = line.split()
            if (line[0] == '#') or (line[0] == '!') or (len(elem)==0):
                if ename in line:
                    j = 1
                    while self.lat[i+j][0] == '!':
                        j+=1
                    elem = self.lat[i+j].split()
                    
                    if  'drift' in elem:      #Drift
                        length = float(elem[2])/100.
                        s += length
                    elif  'cav' in elem:  #Cavity
                        length = float(elem[2])/100.
                        s += length
                    elif  'rfq' in elem:  #Cavity
                        length = float(elem[3])/100.
                        s += length
                    elif 'mhb4' in elem:
                        length = float(elem[2])/100.
                        s += length
                    elif  'mhb' in elem:  #Cavity
                        length = float(elem[4])/100.
                        s += length
                    elif  'bunch' in elem:  #Cavity
                        length = float(elem[2])/100.
                        s += length
                    elif  'sol' in elem:  #Solenoid
                        length = float(elem[3])/100.
                        s += length
                    elif  ('bmag' in elem) or ('dipo' in elem):
                        length = float(elem[2])/100.
                        s += length
                    elif ('edp3d' in elem): # dipoles
                        if len(elem) == 8:
                            length = float(elem[2])/100.
                        else:
                            length = float(elem[3])/100.
                        s += length
                    elif  ('mag3d' in elem): #M dipole
                        if len(elem) == 8:
                            length = float(elem[2])/100.
                        else:
                            length = float(elem[3])/100.
                        s += length
                    elif  ('mult' in elem): #M dipole
                        length = float(elem[2])/100.
                        s += length
                    elif ('quad' in elem): # magnetic quad
                        length = float(elem[3])/100.
                        s += length
                    elif ('mq3d' in elem): # magnetic quad
                        length = float(elem[2])/100.
                        s += length
                    elif ('eq3d' in elem) or ('equad' in elem):  #eQuad
                        length = float(elem[3])/100.
                        s += length
                    elif ('bpm' in elem): #FLAME
                        pass
                    elif ('marker' in elem): #FLAME
                        pass
                    elif ('zero' in elem) or ('corr' in elem): #Corrector
                        pass
                    elif ('shrt' in elem): #Corrector
                        pass
                    #print str(s)+'\t'+str(elem)
                    break

            elif ('updat' in line) or ('matrx' in line) or ('scrch' in line) or ('stop' in line):
                ielem += 1
            elif  'drift' in line:      #Drift
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif  'cav' in line:  #Cavity
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif  'rfq' in line:  #Cavity
                length = float(elem[3])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif 'mhb4' in elem:
                length = float(elem[2])/100.
                s += length
                ielem += 1
            elif  'mhb' in line:  #Cavity
                length = float(elem[4])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif  'bunch' in line:  #Cavity
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif  'sol' in line:  #Solenoid
                length = float(elem[3])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif  ('bmag' in line) or ('dipo' in line):
                length = float(elem[2])/100.
                s += length
                ielem += 1
            elif ('edp3d' in line): #M dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                s += length
                ielem += 1
            elif  ('mag3d' in line): #M dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                s += length
                ielem += 1
            elif  ('mult' in line): #M dipole
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
                ielem += 1
            elif ('quad' in line): # magnetic quad
                length = float(elem[3])/100.
                s += length
                ielem += 1
            elif ('mq3d' in line): # magnetic quad
                length = float(elem[2])/100.
                s += length
                ielem += 1
            elif ('eq3d' in line) or ('equad' in line):  #eQuad
                length = float(elem[3])/100.
                s += length
                ielem += 1
            elif ('bpm' in line) or ('mon' in line): #FLAME
                ielem += 1
            elif ('marker' in line): #FLAME
                ielem += 1
            elif ('shrt' in line): #FLAME
                ielem += 1
            elif ('zero' in line) or ('corr' in line): #Corrector
                ielem += 1
            else:
                ielem += 1
                print (line),

        print (ename+':\t'+str(ielem)+'\t'+str(s))


        res0 = self.get_data()

        sl = res0['pos']
        #print (str(s))
        #print (str(sl))
        n = len(sl)
        #ielem = 0
        for i in range(n):
            #if math.fabs(sl[ielem]-s)<0.0001:
            if res0['eid'][i] == ielem:
                iline = i
                break

            #ielem += 1

        res = {}
        res['eid'] = res0['eid'][iline]
        res['type'] = res0['type'][iline]
        res['s'] = res0['s'][iline]
        res['pos'] = res0['s'][iline]
        res['Ek_ref'] = res0['Ek_ref'][iline]
        res['dwwp'] = res0['dwwp'][iline]
        res['xrms'] = res0['xrms'][iline]
        res['yrms'] = res0['yrms'][iline]
        res['xmax'] = res0['xmax'][iline]
        res['ymax'] = res0['ymax'][iline]
        res['phirms'] = res0['phirms'][iline]
        res['phimax'] = res0['phimax'][iline]
        res['dww'] = res0['dww'][iline]
        res['xnemit'] = res0['xnemit'][iline]
        res['xnemit90'] = res0['xnemit90'][iline]
        res['xnemitmax'] = res0['xnemitmax'][iline]
        res['ynemit'] = res0['ynemit'][iline]
        res['ynemit90'] = res0['ynemit90'][iline]
        res['ynemitmax'] = res0['ynemitmax'][iline]
        res['znemit'] = res0['znemit'][iline]
        res['znemit90'] = res0['znemit90'][iline]
        res['znemitmax'] = res0['znemitmax'][iline]
        res['xcen'] = res0['xcen'][iline]
        res['xpcen'] = res0['xpcen'][iline]
        res['ycen'] = res0['ycen'][iline]
        res['ypcen'] = res0['ypcen'][iline]
        res['bg'] = res0['bg'][iline]
        res['zcen'] = res0['zcen'][iline]
        res['xalpha'] = res0['xalpha'][iline]
        res['xbeta'] = res0['xbeta'][iline]
        res['yalpha'] = res0['yalpha'][iline]
        res['ybeta'] = res0['ybeta'][iline]
        res['zalpha'] = res0['zalpha'][iline]
        res['zbeta'] = res0['zbeta'][iline]
        res['nloss'] = res0['nloss'][iline]
        res['np'] = res0['np'][iline]
        res['xpadv'] = res0['xpadv'][iline]
        res['ypadv'] = res0['ypadv'][iline]
        
        return res
    '''
    #################
    ##########
    def get_profile_step(self,ename):

        # get s at ename
        '''
        if fnlat == '':
            fnlat = self.lat_file
        f= open(fnlat,'r')
        a = f.readlines()
        f.close()
        n = len(a)
        '''

        n = len(self.lat)
        s = 0
        for i in range(n):
            line = self.lat[i]
            elem = line.split()
            if (line[0] == '#') or (line[0] == '!') or (len(elem)==0):
                if ename in line:
                    j = 1
                    while self.lat[i+j][0] == '!':
                        j+=1
                    elem = self.lat[i+j].split()
                    
                    if  'drift' in elem:      #Drift
                        length = float(elem[2])/100.
                        s += length
                    elif  'cav' in elem:  #Cavity
                        length = float(elem[2])/100.
                        s += length
                    elif  'rfq' in elem:  #Cavity
                        length = float(elem[3])/100.
                        s += length
                    elif  'mhb' in elem:  #Cavity
                        length = float(elem[4])/100.
                        s += length
                    elif  'bunch' in elem:  #Cavity
                        length = float(elem[2])/100.
                        s += length
                    elif  'sol' in elem:  #Solenoid
                        length = float(elem[3])/100.
                        s += length
                    elif  ('bmag' in elem) or ('dipo' in elem):
                        length = float(elem[2])/100.
                        s += length
                    elif ('edp3d' in elem): # dipoles
                        if len(elem) == 8:
                            length = float(elem[2])/100.
                        else:
                            length = float(elem[3])/100.
                        s += length
                    elif  ('mag3d' in elem): #M dipole
                        if len(elem) == 8:
                            length = float(elem[2])/100.
                        else:
                            length = float(elem[3])/100.
                        s += length
                    elif  ('mult' in elem): #M dipole
                        length = float(elem[2])/100.
                        s += length
                    elif ('quad' in elem): # magnetic quad
                        length = float(elem[3])/100.
                        s += length
                    elif ('mq3d' in elem): # magnetic quad
                        length = float(elem[2])/100.
                        s += length
                    elif ('eq3d' in elem) or ('equad' in elem):  #eQuad
                        length = float(elem[3])/100.
                        s += length
                    elif ('bpm' in elem): #FLAME
                        pass
                    elif ('marker' in elem): #FLAME
                        pass
                    elif ('zero' in elem) or ('corr' in elem): #Corrector
                        pass
                    elif ('shrt' in elem): #Corrector
                        pass
                    #print str(s)+'\t'+str(elem)
                    break

            elif ('updat' in line) or ('matrx' in line) or ('scrch' in line) or ('stop' in line):
                pass
            elif  'drift' in line:      #Drift
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  'cav' in line:  #Cavity
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  'rfq' in line:  #Cavity
                length = float(elem[3])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  'mhb' in line:  #Cavity
                length = float(elem[4])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  'bunch' in line:  #Cavity
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  'sol' in line:  #Solenoid
                length = float(elem[3])/100.
                s += length
                #print str(s)+'\t'+line,
            elif  ('bmag' in line) or ('dipo' in line): #M dipole
                length = float(elem[2])/100.
                s += length
            elif  'edp3d' in line: #E dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                s += length
            elif  ('mag3d' in line): #M dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                s += length
            elif  ('mult' in line): #M dipole
                length = float(elem[2])/100.
                s += length
                #print str(s)+'\t'+line,
            elif ('quad' in line): # magnetic quad
                length = float(elem[3])/100.
                s += length
            elif ('mq3d' in line): # magnetic quad
                length = float(elem[2])/100.
                s += length
            elif ('eq3d' in line) or ('equad' in line):  #eQuad
                length = float(elem[3])/100.
                s += length
            elif ('bpm' in line): #FLAME
                pass
            elif ('marker' in line): #FLAME
                pass
            elif ('shrt' in line): #FLAME
                pass
            elif ('zero' in line) or ('corr' in line): #Corrector
                pass
            else:
                print (line),

        #print ename+':\t'+str(s)


        res0 = self.get_data_step()

        sl = res0['pos']
        #print (str(s))
        #print (str(sl))
        n = len(sl)
        ielem = 0
        for i in range(n):
            if math.fabs(sl[ielem]-s)<0.0001:
                break

            ielem += 1

        res = {}
        res['s'] = res0['s'][ielem]
        res['pos'] = res0['s'][ielem]
        res['Ek_ref'] = res0['Ek_ref'][ielem]
        res['dwwrms'] = res0['dwwrms'][ielem]
        res['xrms'] = res0['xrms'][ielem]
        res['yrms'] = res0['yrms'][ielem]
        res['xmax'] = res0['xmax'][ielem]
        res['ymax'] = res0['ymax'][ielem]
        res['phirms'] = res0['phirms'][ielem]
        res['phimax'] = res0['phimax'][ielem]
        res['dwwrms'] = res0['dwwrms'][ielem]
        res['xnemit'] = res0['xnemit'][ielem]
        res['ynemit'] = res0['ynemit'][ielem]
        res['znemit'] = res0['znemit'][ielem]
        res['xcen'] = res0['xcen'][ielem]
        res['xpcen'] = res0['xpcen'][ielem]
        res['ycen'] = res0['ycen'][ielem]
        res['ypcen'] = res0['ypcen'][ielem]
        res['bg'] = res0['bg'][ielem]
        res['zcen'] = res0['zcen'][ielem]
        res['xalpha'] = res0['xalpha'][ielem]
        res['xbeta'] = res0['xbeta'][ielem]
        res['yalpha'] = res0['yalpha'][ielem]
        res['ybeta'] = res0['ybeta'][ielem]
        res['xpadv'] = res0['xpadv'][ielem]
        res['ypadv'] = res0['ypadv'][ielem]
        
        return res

    #################
    def set_initial_moment(self, moment0 = {}, moment1 = {}):
        trk = []
        for l in self.trk:
            for par in ['x00', 'xp00', 'y00', 'yp00', 'ph00', 'dww00']:
                if (par in l) and (par in  moment0.keys()):
                    if (par == 'x00') or (par == 'y00'):
                        val = moment0[par]*100. # [m] to [cm]
                    elif par == 'ph00':
                        val = moment0[par]*180./np.pi # [rad] to [deg]
                    else:
                        val = moment0[par]
                    e = l.split(',')
                    n = len(e)
                    for j in range(n):
                        if par in e[j]:
                            e[j] = par+' = '+str(val)
                    l = ', '.join(e)+'\n'
            '''
            if ('x00' in l) and ('x' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'x00' in e[j]:
                        e[j] = 'x00 = '+str(moment0['x']/10.)
                l = ', '.join(e)+'\n'
            if ('xp00' in l) and ('xp' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'xp00' in e[j]:
                        e[j] = 'xp00 = '+str(moment0['xp'))
                l = ', '.join(e)+'\n'
            if ('y00' in l) and ('y' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'y00' in e[j]:
                        e[j] = 'y00 = '+str(moment0['y']/10.)
                l = ', '.join(e)+'\n'
            if ('yp00' in l) and ('yp' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'yp00' in e[j]:
                        e[j] = 'yp00 = '+str(moment0['yp']))
                l = ', '.join(e)+'\n'
            if ('ph00' in l) and ('phi' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'ph00' in e[j]:
                        e[j] = 'ph00 = '+str(moment0['phi'])
                l = ', '.join(e)+'\n'
            if ('dww00' in l) and ('E' in moment0.keys()):
                e = l.split(',')
                n = len(e)
                for j in range(n):
                    if 'dww00' in e[j]:
                        e[j] = 'dww00 = '+str(moment0['E'])
                l = ', '.join(e)+'\n'
            '''
            for par in ['alfax', 'betax', 'epsnx', 'alfay', 'betay', 'epsny']:
                if (par in l) and (par in moment1.keys()):
                    if 'alfa' in par:
                        val = moment1[par]
                    elif 'beta' in par:
                        val = moment1[par]/10. # [m/mrad] to [cm/rad]
                    elif 'epsn' in par:
                        val = moment1[par]*100*6 # [m-rad] to 6x[cm-mrad]
                    e = l.split(',')
                    n = len(e)
                    for j in range(n):
                        if par in e[j]:
                            e[j] = par+' = '+str(val)
                    l = ', '.join(e)+'\n'
            trk.append(l)
        self.trk = trk
        
    #################
    def draw_lattice(self, fig, fnlat='', smin = -1, smax = -1, nplot = 1, bshow = False):

        if len(fnlat)==0:
            fnlat = self.lat_file
            
        figw,figh = fig.get_figwidth(),fig.get_figheight()
        hlatheight = 0.04*figw/figh
        hplotheight= 0.8-(hlatheight+0.01)
        
        
        if fnlat == '':
            fnlat = self.lat_file
        hlat = fig.add_axes((0.1, 0.1, 0.8, hlatheight))
        
        hlat.set_yticklabels(())
        z = 0.0
        file= open(fnlat)
        zstart = 0
        for line in iter(file):
            line = line.lstrip().rstrip()
            elem = line.split()
            if (len(line)==0) or (line[0] == '#') or (line[0] == '!') or (len(elem)==0):
                continue

            #elem = line.replace(' ', '').replace(';','').split(',')
            #print elem
            '''
            if ('L' in line) and ('aper' in line):
                L_loc_No = 0
                for i in range(len(elem)):
                    if 'L=' in elem[i]:
                        L_loc_No = i
                        break
                length = float(elem[L_loc_No][2:])
            else:
                length = 0.0
            '''
            
            zstart = z
            if ('updat' in line) or ('matrx' in line) or ('scrch' in line) or ('stop' in line):
                pass
            elif  'drift' in line:      #Drift
                length = float(elem[2])/100.
                z= zstart+ length
                #print str(z)+'\t'+line,
            elif  'cav' in line:  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                if bshow:
                    print('cav\t'+str(zstart)+'\t'+str(z))
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  'rfq' in line:  #Cavity
                length = float(elem[3])/100.
                z= zstart+ length
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  ('mhb4' in line):  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                if bshow:
                    print('mhb4\t'+str(zstart)+'\t'+str(z))
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  'mhb' in line:  #Cavity
                length = float(elem[4])/100.
                z= zstart+ length
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  ('extrc' in line):  #Cavity
                length = float(elem[3])/100.
                z= zstart+ length
                if bshow:
                    print('extrc\t'+str(zstart)+'\t'+str(z))
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
            elif  ('bunch' in line):  #Cavity
                length = float(elem[2])/100.
                z= zstart+ length
                if bshow:
                    print('bunch\t'+str(zstart)+'\t'+str(z))
                
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.4),length, 0.8,
                    facecolor = 'Blue',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  'sol3c' in line:  # 3 magnetic field component with a field map
                length = float(elem[3])/100.
                z= zstart+ length
                if bshow:
                    print('sol3c\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, 0.-0.15),length, 0.3,
                    facecolor = 'Red',
                    edgecolor=None,
                    alpha=0.8,
                ))
            elif  'sol' in line:  #Solenoid
                length = float(elem[3])/100.
                z= zstart+ length
                if bshow:
                    print('sol\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, 0.-0.15),length, 0.3,
                    facecolor = 'Red',
                    edgecolor=None,
                    alpha=0.8,
                ))
                #print str(z)+'\t'+line,
            elif  ('hbmag' in line): #M dipole
                length = float(elem[3])*math.pi*2*abs(float(elem[2]))/360./100.
                z= zstart+ length
                #print('hbmag\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.5),length, 1,
                    facecolor = 'Green',
                    edgecolor=None,
                    alpha=0.5,
                ))
                #print str(z)+'\t'+line,
            elif  ('bmag' in line) or ('dipo' in line): #dipole
                length = float(elem[2])/100.
                z= zstart+ length
                #print('bmag\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.5),length, 1,
                    facecolor = 'Green',
                    edgecolor=None,
                    alpha=0.5,
                ))
                #print str(z)+'\t'+line,
            elif  ('edp3d' in line):
                if len(elem) == 9:
                    length = float(elem[3])/100.                    
                else:
                    length = float(elem[2])/100.
                z= zstart+ length
                #print('bmag\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.5),length, 1,
                    facecolor = 'Green',
                    edgecolor=None,
                    alpha=0.5,
                ))
                #print str(z)+'\t'+line,
            elif  ('mag3d' in line): #dipole
                if len(elem) == 8:
                    length = float(elem[2])/100.
                else:
                    length = float(elem[3])/100.
                z= zstart+ length
                #print('mag3d\t'+str(zstart)+'\t'+str(z))
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -0.5),length, 1,
                    facecolor = 'Green',
                    edgecolor=None,
                    alpha=0.5,
                ))
            elif  ('mult' in line): #M dipole
                length = float(elem[2])/100.
                length_eff = float(elem[3])/100.
                length_dif = length-length_eff
                sw = int(elem[8])
                z= zstart+ length
                #print('mult\t'+str(zstart+length_dif/2.)+'\t'+str(z-length_dif/2.))
                #print(str(length)+'\t'+str(length_eff))
                if sw == 0:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart+length_dif/2., -0.5),length_eff, 1,
                        facecolor = 'slategray',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                elif sw == 1: # 1st half
                    hlat.add_patch(ptc.Rectangle(
                        (zstart+length_dif, -0.5),length_eff, 1,
                        facecolor = 'slategray',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                elif sw == 2: # 1st half
                    hlat.add_patch(ptc.Rectangle(
                        (zstart, -0.5),length_eff, 1,
                        facecolor = 'slategray',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                    
                    #print str(z)+'\t'+line,
            elif ('quad' in line): # magnetic quad
                length = float(elem[3])/100.
                length_eff = float(elem[4])/100.
                length_dif = length-length_eff
                z= zstart+ length
                #print('quad\t'+str(zstart+length_dif/2.)+'\t'+str(z-length_dif/2.))
                B=float(elem[2])
                if B >=0.:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart+length_dif/2., -0.2),length_eff, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                else:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart+length_dif/2., -0.6),length_eff, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                #print str(z)+'\t'+line,
            elif ('mq3d' in line): # magnetic quad
                length = float(elem[2])/100.
                z= zstart+ length
                #print('mq3d\t'+str(zstart)+'\t'+str(z))
                B=float(elem[3])
                '''
                for e in elem:
                    if 'B2=' in e:
                        e2 = e.split('=')
                        B=float(e2[1])
                '''
                if B >=0.:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart, -0.2),length, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                else:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart, -0.6),length, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
            elif ('eq3d' in line) or ('equad' in line):  #eQuad
                length = float(elem[3])/100.
                z= zstart+ length
                if bshow:
                    print('eq3d\t'+str(zstart)+'\t'+str(z))
                V=0.0
                '''
                for e in elem:
                    if 'V=' in e:
                        e2 = e.split('=')
                        V=float(e2[1])
                '''
                V = float(elem[2])
                if V >=0.:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart, -0.2),length, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
                else:
                    hlat.add_patch(ptc.Rectangle(
                        (zstart, -0.6),length, 0.8,
                        facecolor = 'Purple',
                        edgecolor=None,
                        alpha=0.5,
                    ))
            elif ('bpm' in line): #FLAME
                z= zstart
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -1),0.05, 1,
                    facecolor = 'Green',
                    edgecolor=None,
                    alpha=0.5,
                ))
            elif ('marker' in line): #FLAME
                z= zstart
                pass
            elif ('zero' in line) or ('corr' in line): #Corrector
                z= zstart
                hlat.add_patch(ptc.Rectangle(
                    (zstart, -1),0.05, 2,
                    facecolor = 'Orange',
                    edgecolor=None,
                    alpha=0.5,
                ))
            elif ('mon' in line) or ('MON' in line):
                pass
                #print('mon\t'+str(z))
            elif ('ascii' in line):
                pass
            elif ('shrt' in line):
                pass
            elif ('filtr' in line):
                pass
            else:
                print (line),
        #print('end\t'+str(z))
        file.close
        hlat.set_autoscale_on(False)

        if smin>0:
            z1=smin
        else:
            z1=0.0
        if smax>=0:
            z2=smax
        else:
            z2=z
        if z2>=z1:
            hlat.axis((z1,z2,-1,1))
        else:
            hlat.axis((z2,z1,-1,1))
        hlat.set_xlim(z1, z2)
        #print('smin: '+str(smin)+'\tsmax: '+str(smax))
        #print('z1: '+str(z1)+'\tz2: '+str(z2))
        hlat.plot([z1,z2],[0,0],'--k',linewidth=0.5)
        #print str(z1)+'\t'+str(z2)

        hist = []
        if nplot == 1:
            #ht = fig.add_axes((0.1, 0.21, 0.8, 0.69))
            ht = fig.add_axes((0.1, 0.1+hlatheight+0.01, 0.8, hplotheight))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)
        elif nplot == 2:
            hlataxes = (hplotheight-0.01)/2.
            ht = fig.add_axes((0.1, 0.1+hlatheight+hlataxes+0.02, 0.8, hlataxes))
            #ht = fig.add_axes((0.1, 0.55, 0.8, 0.34))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)
            
            ht = fig.add_axes((0.1, 0.1+hlatheight+0.01, 0.8, hlataxes))
            #ht = fig.add_axes((0.1, 0.21, 0.8, 0.34))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)
        elif nplot == 3:
            hlataxes = (hplotheight-0.02)/3.
            #ht = fig.add_axes((0.1, 0.68, 0.8, 0.22))
            ht = fig.add_axes((0.1, hlatheight+0.1+2*hlataxes+0.03, 0.8, hlataxes))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)
            
            #ht = fig.add_axes((0.1, 0.45, 0.8, 0.22))
            ht = fig.add_axes((0.1, 0.1+hlatheight+hlataxes+0.02, 0.8, hlataxes))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)

            #ht = fig.add_axes((0.1, 0.22, 0.8, 0.22))
            ht = fig.add_axes((0.1, 0.1+hlatheight+0.01, 0.8, hlataxes))
            ht.set_xticklabels( () )
            ht.set_xlim(z1,z2)
            hist.append(ht)
        hist.append(hlat)
        return hist

    #################
    def gen_dist_long_accep(self, q = 18, A = 40, n = 10000, phimax = 180, phimin = -180., Emax = 22, Emin = 18, fname = 'read_dis_tmp.dat'):
        deg = math.pi/180.
        m  = 931.49432
        
        El = np.random.uniform(Emin, Emax, n)
        phil = np.random.uniform(phimin*deg, phimax*deg, n)

        gl = (El+m)/m
        bl = np.sqrt(1.-1./gl/gl)
        bgl = bl*gl
        #Ecen = np.average(El)
        Ecen = (Emin+Emax)/2.
        #bgcen = np.average(bgl)
        gcen = (Ecen+m)/m
        bcen = math.sqrt(1.-1./gcen/gcen)
        bgcen = bcen*gcen
        #phicen = np.average(phil)
        phicen = 0.0
        
        f = open(fname,'w')
        f.write(str(Ecen*1e+3)+'\t1\n')
        f.write(str(n)+'\n')
        f.write(str(q)+'\n')
        f.write('0.0\t0.0\t0.0\t0.0\t'+str(phicen)+'\t'+str(bgcen)+'\t0\n')
        for i in range(n):
            f.write('0.0\t0.0\t0.0\t0.0\t'+str(phil[i])+'\t'+str(bgl[i])+'\t0\n')
        f.close()
        
        
    #################        
    def gen_dist_trans_accep(self, ran, q = 18, A = 40, n = 10000, fname = 'read_dis_tmp.dat', seed = 1):
        deg = math.pi/180.
        m  = 931.49432

        np.random.seed(seed)
        xl = np.random.uniform(ran['xmin'], ran['xmax'], n)/10. # [mm] -> [cm]
        xpl = np.random.uniform(ran['xpmin'], ran['xpmax'], n)/1000. # [mrad] -> [rad]
        yl = np.random.uniform(ran['ymin'], ran['ymax'], n)/10. # [mm] -> [cm]
        ypl = np.random.uniform(ran['ypmin'], ran['ypmax'], n)/1000. # [mrad] -> [rad]
        El = np.random.uniform(ran['Emin'], ran['Emax'], n) # [MeV/u]
        phil = np.random.uniform(ran['phimin']*deg, ran['phimax']*deg, n)# [deg] -> [rad]

        gl = (El+m)/m
        bl = np.sqrt(1.-1./gl/gl)
        bgl = bl*gl
        #Ecen = np.average(El)
        Ecen = (ran['Emin']+ran['Emax'])/2.
        #bgcen = np.average(bgl)
        gcen = (Ecen+m)/m
        bcen = math.sqrt(1.-1./gcen/gcen)
        bgcen = bcen*gcen
        #phicen = np.average(phil)
        phicen = 0.0
        
        f = open(fname,'w')
        f.write(str(Ecen*1e+3)+'\t1\n')
        f.write(str(n)+'\n')
        f.write(str(q)+'\n')
        f.write('0.0\t0.0\t0.0\t0.0\t'+str(phicen)+'\t'+str(bgcen)+'\t0\n')
        for i in range(n):
            f.write(str(xl[i])+'\t'+str(xpl[i])+'\t'+str(yl[i])+'\t'+str(ypl[i])+'\t'+str(phil[i])+'\t'+str(bgl[i])+'\t0\n')
        f.close()
        
        
    #################
    def fit_gauss(self, p, bshow = False):
        wl = [1]*3
        nbins = len(p[0])

        xl = np.array([(p[1][i] + p[1][i+1])/2. for i in range(nbins)])
        yl = np.array(p[0])

        xcen = sum(xl*yl)/sum(yl)
        ysum = np.sum(yl)
        rms = np.sqrt(1/ysum*np.sum(np.power(xl-xcen, 2)))
        f0l = [np.max(yl), xcen, rms]
        ###
        def get_chi2(f1l):
            fl = [f0l[i]+f1l[i]*wl[i] for i in range(3)]
            y1l = fl[0]*np.exp(-1*np.power((xl-fl[1])/fl[2],2))
            chi2 = np.sqrt(np.sum(np.power(yl-y1l, 2)))

            return chi2
        ###
        inp = 'y'
        chi2last = get_chi2([0]*3)
        while inp != 'n':
            res = minimize(get_chi2, [0]*3, method = 'Nelder-Mead', tol = 1e-7)
            for i in range(3):
                f0l[i] += wl[i]*res.x[i]
            chi2 = get_chi2([0]*3)
            if chi2>=chi2last:
                inp='n'
            #print(chi2)
            chi2last = chi2

        y1l = f0l[0]*np.exp(-1*np.power((xl-f0l[1])/f0l[2],2))
            
        d = {
            'height':f0l[0],
            'cen':f0l[1],
            'sig':f0l[2],
            'res':f0l,
            'chi2':chi2,
            'x':xl,
            'y':y1l
            }
        if bshow:
            print('Height:\t'+str(f0l[0])+'\tcentroid:\t'+str(f0l[1])+'\tsigma:\t'+str(f0l[2])+'\tchi2:\t'+str(d['chi2']))
        return d
        
    #################        
    #################        
    def TRACK_dist_gen6D(self, PNum=[10000],
                         Ekmevu=[0.012],  # accelerated kinetic energy in eV
                         kEsprd=[0.0005],
                         A=[40],
                         Q=[9],
                         M=[39.91],
                         emitx=[0.066 * 10 ** (-5)],
                         alphax=[0],
                         betax=[0.065],
                         emity=[0.066 * 10 ** (-5)],
                         alphay=[0],
                         betay=[0.065],
                         emitz=[0.066 * 10 ** (-5)],
                         alphaz=[0],
                         betaz=[0.065],
                         rho13=[0],
                         rho14=[0],
                         rho15=[0],
                         rho16=[0],
                         rho23=[0],
                         rho24=[0],
                         rho25=[0],
                         rho26=[0],
                         rho35=[0],
                         rho36=[0],
                         rho45=[0],
                         rho46=[0],
                         sigmal = [[1.866103774559999984e+00,4.972447301369999984e-04,3.959648994790000015e-03,1.718973260559999869e-05,2.654352962009999886e-04,1.675530925640000032e-05,0.000000000000000000e+00,
                                    4.972447301369999984e-04,5.744434998460000419e-06,-5.751544512039999873e-06,2.142679787250000069e-08,-1.692103655169999977e-07,2.685342374179999817e-09,0.000000000000000000e+00,
                                    3.959648994790000015e-03,-5.751544512039999873e-06,1.562692086209999998e+00,2.281577507259999920e-04,4.501222748789999984e-04,-1.480766713469999916e-05,0.000000000000000000e+00,
                                    1.718973260559999869e-05,2.142679787250000069e-08,2.281577507259999920e-04,7.126095562939999942e-06,-3.558792535480000069e-07,-5.507375936640000302e-08,0.000000000000000000e+00,
                                    2.654352962009999886e-04,-1.692103655169999977e-07,4.501222748789999984e-04,-3.558792535480000069e-07,1.063777264700000108e-03,1.840365217029999917e-06,0.000000000000000000e+00,
                                    1.675530925640000032e-05,2.685342374179999817e-09,-1.480766713469999916e-05,-5.507375936640000302e-08,1.840365217029999917e-06,1.122641260759999906e-06,0.000000000000000000e+00,
                                    0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00]
                               ],
                         moment0 = [0,0,0,0,0,0],
                         Flag_DCbeam=0,
                         Flag_initG=1,
                         gene_type=0):



        ####################################################################################################################
        sigama_z = 0.01  # rad
        mu_z = 0.0000
        sigama_zp = 2500.0 * 2  # eV
        mu_zp = 0.0000
        
        fo = open("partcl.data", "w")
        fo2 = open("read_dis.dat", "w")

        fo.write('#x[cm], xp[rad], y[cm], yp[rad], z_phase[rad], zp(betagamma)\n')

        #fo2.write('%f, %d \n' % float(Ekmevu[0] * 1000.0), len(Ekmevu))
        #fo2.write('%d\n' % int(PNum[0] - 1))
        #fo2.write('%d\n' % Q)
        print (Ekmevu[0], Ekmevu)
        
        Ekline = str(float(Ekmevu[0] * 1000.0)) + " " + str(len(Ekmevu))
        PNumline = ""
        Qline =""

        for i in range(len(Q)):
            PNumline    += str(int(PNum[i] - 1)) + " "
            Qline	 += str(Q[i]) + " "

        fo2.write(Ekline + "\n")
        fo2.write(PNumline + "\n")
        fo2.write(Qline + "\n")

        # multi charge states cal.
        for beamNO in range(0, len(Ekmevu), 1):
            m0 = 9.3149432e+8*M[beamNO]/A[beamNO]  # [eV]
            print('beamNO '+str(beamNO)+'\tM '+str(M[beamNO])+'\tA '+str(A[beamNO]))
            #mc2 = np.array(A) * m0  # eV
            if not(sigmal == 'None'):
                #sigmaM = np.array(sigmal[beamNO]).reshape(7,7)[0:6,0:6]
                if Flag_DCbeam == 1 or Flag_DCbeam == 3:
                    print('DC beam generation')
                    sigmaM = np.array(sigmal[beamNO]).reshape(7,7)[0:4,0:4]
                    print (sigmaM)
                    self.det = np.linalg.det(sigmaM)
                    print ("det (4D):", self.det)
                else:
                    sigmaM = np.array(sigmal[beamNO]).reshape(7,7)[0:6,0:6]
                    print (sigmaM)
                    self.det = np.linalg.det(sigmaM)
                    print ("det (6D):", self.det)


                
            # determinant

            if (self.det < 0):
                print ("det<0, then zero coupling is set")
                # determinant should be positive determinant of
                sigma11 = sigmaM[0][0]
                sigma12 = sigmaM[0][1]
                sigma21 = sigmaM[1][0]
                sigma22 = sigmaM[1][1]

                sigma33 = sigmaM[2][2]
                sigma34 = sigmaM[2][3]
                sigma43 = sigmaM[3][2]
                sigma44 = sigmaM[3][3]

                if Flag_DCbeam == 0:
                    sigma55 = sigmaM[4][4]
                    sigma56 = sigmaM[4][5]
                    sigma65 = sigmaM[5][4]
                    sigma66 = sigmaM[5][5]
                else:
                    sigma55 = 1.0
                    sigma56 = 0.0
                    sigma65 = 0.0
                    sigma66 = 1e-12

                sigma13 = 0.0
                sigma14 = 0.0
                sigma15 = 0.0
                sigma16 = 0.0
                
                sigma23 = 0.0
                sigma24 = 0.0
                sigma25 = 0.0
                sigma26 = 0.0
                
                sigma35 = 0.0
                sigma36 = 0.0
                sigma45 = 0.0
                sigma46 = 0.0

                sigma31 = sigma13
                sigma41 = sigma14
                sigma51 = sigma15
                sigma61 = sigma16
				
                sigma32 = sigma23
                sigma42 = sigma24
                sigma52 = sigma25
                sigma62 = sigma26
				
                sigma53 = sigma35
                sigma63 = sigma36
                sigma54 = sigma45
                sigma64 = sigma46		
                
                sigmaM = np.array([
                    [sigma11, sigma12, sigma13, sigma14, sigma15, sigma16],
                    [sigma21, sigma22, sigma23, sigma24, sigma25, sigma26],
                    [sigma31, sigma32, sigma33, sigma34, sigma35, sigma36],
                    [sigma41, sigma42, sigma43, sigma44, sigma45, sigma46],
                    [sigma51, sigma52, sigma53, sigma54, sigma55, sigma56],
                    [sigma61, sigma62, sigma63, sigma64, sigma65, sigma66],
                ])

            # inverse matrix
            sigmaM = np.array(sigmal[beamNO]).reshape(7,7)[0:6,0:6]
            print('sigmaM: '+str(sigmaM))
            '''
            inv_sigmaM = np.linalg.inv(sigmaM)

            #print ("det (6D):", self.det)
            print ("inv_sigmaM:", inv_sigmaM)
            '''

            # https://lucidfrontier45.wordpress.com/2011/07/23/multivariate_gauss/

            i = 0

            xrms, xprms = 0, 0
            yrms, yprms = 0, 0
            xyrms, xyprms = 0, 0
            xpyrms, xpyprms = 0, 0
                
            mean = moment0
            cov = sigmaM

            if gene_type == 0:  # python module only for Gaussian
                if Flag_DCbeam == 1 or Flag_DCbeam == 3:
                    x, xp, y, yp, = np.random.multivariate_normal(mean[0:4], cov[0:4,0:4], PNum[beamNO] + 1, check_valid='warn', tol = 1e-20).T
                else:
                    x, xp, y, yp, phase, dEkmevu = np.random.multivariate_normal(mean, cov, PNum[beamNO] + 1, check_valid='warn', tol = 1e-20).T

            elif gene_type == 1:  # Cholesky module
                v = cholesky(cov)
                np.random.seed(beamNO)
                cholesky_rand_init = np.random.randn(PNum[beamNO] + 1, 6)
                    
                #cholesky_rand_init = np.random.randn(1, 4)
                Y = np.dot(cholesky_rand_init, v)
                x, xp, y, yp, phase, dEkmevu = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4], Y[:, 5]
                    
            elif gene_type == 2:  # KV distribution with cholesky module (IMAPCT module copy in Distribution.f90)
                v = cholesky(cov)

                # initialization
                cholesky_rand_init = [[0 for i in range(0,6,1)] for i in range(0,PNum[beamNO] + 1,1)]
                
                twopi = 4*np.arcsin(1.0)

                for i in range(0,PNum[beamNO]+1):
                        
                    r1, r2, r3 = np.random.uniform(0, 1), np.random.uniform(0, 1), np.random.uniform(0, 1)
                    r4 = np.sqrt(r1)
                    r5 = np.sqrt(1.0 - r1)
                    r2 = r2 * twopi
                    r3 = r3 * twopi

                    x1 = 2 * r4 * np.cos(r2)
                    x2 = 2 * r4 * np.sin(r2)
                    x3 = 2 * r5 * np.cos(r3)
                    x4 = 2 * r5 * np.sin(r3)

                    r5 = np.random.uniform(0, 1)
                    r6 = np.random.uniform(0, 1)
                    
                    x5 = r5*np.sqrt(3)
                    x6 = r6*np.sqrt(3)

                    cholesky_rand_init[i] = [x1, x2, x3, x4, x5, x6]

                    Y = np.dot(cholesky_rand_init, v)
                    x, xp, y, yp, phase, dEkmevu = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4], Y[:, 5]

            elif gene_type == 3:  # 6D waterbag cholesky module
                v = cholesky(cov)

                # initialization
                cholesky_rand_init = [[0 for i in range(0,6,1)] for i in range(0,PNum[beamNO] + 1,1)]

                i = 0
                while i < PNum[beamNO]+1:

                    rnd0, rnd1, rnd2, rnd3, rnd4, rnd5 = np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)
                    if(rnd0**2 + rnd1**2 + rnd2**2 + rnd3**2 + rnd4**2 + rnd5**2 <= 1):
                        cholesky_rand_init[i] = [rnd0, rnd1, rnd2, rnd3, rnd4, rnd5]
                        i = i + 1

                        if i%1000000 == 0:
                            print (str(i)+' generated')
                Y = np.dot(cholesky_rand_init, v)*np.sqrt(8)

                x, xp, y, yp, phase, dEkmevu = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4], Y[:, 5]


            spin = 0


            """
            x, y : mm => cm
            xp, yp: rad => rad
            phase: rad => rad
            deKmevu: MeV/u => MeV/u
            """

            if Flag_DCbeam == 1:
                phase = [np.random.uniform(-np.pi, np.pi) for i in range(PNum[beamNO] + 1)]
                dEkmevu = [0.0 for i in range(PNum[beamNO] + 1)]
            elif Flag_DCbeam == 3:
                phase = [0.0 for i in range(PNum[beamNO] + 1)]
                dEkmevu = [0.0 for i in range(PNum[beamNO] + 1)]

                
            x, xp, y, yp, phase, dEkmevu = np.array(x)*0.1, np.array(xp), np.array(y)*0.1, np.array(yp), np.array(phase), np.array(dEkmevu)

            Ekmevul = Ekmevu[beamNO] + dEkmevu

            #Ekdgammal = Ekmevul/931.49432
            print('m0: '+str(m0))
            Ekdgammal = Ekmevul/(m0*1e-6)
            Ekgammal = 1 + Ekdgammal
            Ekgammabetal = np.sqrt(Ekdgammal * Ekdgammal + 2.0 * Ekdgammal)
                
            i = 0
            while i < PNum[beamNO]:
                i += 1
            i = 0
            while i < PNum[beamNO]:
                i += 1

                if i%1000000 == 0:
                    print (str(i)+' written')
                # x[cm], xp[rad], y[cm], yp[rad],phase [rad], kE, spin

                x[1], xp[1], y[1], yp[1], phase[1] = 0, 0, 0, 0, 0

                fo.write(
                    '%.3e, %.3e, %.3e, %.3e, %.3e, %.6e, %d\n' % (x[i], xp[i], y[i], yp[i], phase[i], Ekmevul[i], spin))
                    
                fo2.write(
                    '%.3e, %.3e, %.3e, %.3e, %.3e, %.6e, %d\n' % (x[i], xp[i], y[i], yp[i], phase[i], Ekgammabetal[i], spin))

                if i > -1:
                    xrms += x[i] * x[i]
                    xprms += xp[i] * xp[i]
                    yrms += y[i] * y[i]
                    yprms += yp[i] * yp[i]

                    xyrms += x[i] * y[i]
                    xyprms += x[i] * yp[i]
                    xpyrms += xp[i] * y[i]
                    xpyprms += xp[i] * yp[i]


            PNum0 = PNum[beamNO]

            xrms = np.sqrt(xrms / PNum0)
            xprms = np.sqrt(xprms / PNum0)
            yrms = np.sqrt(yrms / PNum0)
            yprms = np.sqrt(yprms / PNum0)

            xyrms = xyrms / PNum0
            xyprms = xyprms / PNum0
            xpyrms = xpyrms / PNum0
            xpyprms = xpyprms / PNum0
                    
            rho13_ans = xyrms / (xrms * yrms)
            rho14_ans = xyprms / (xrms * yprms)
            rho23_ans = xpyrms / (xprms * yrms)
            rho24_ans = xpyprms / (xprms * yprms)


            print ("xrms:", xrms," cm")
            print ("yrms:", yrms," cm")
            print ("xprms:", xprms)
            print ("yprms:", yprms)
            
            print ("rho13:", rho13, ", rho13_ans:", rho13_ans)
            print ("rho14:", rho14, ", rho14_ans:", rho14_ans)
            print ("rho23:", rho23, ", rho23_ans:", rho23_ans)
            print ("rho24:", rho24, ", rho24_ans:", rho24_ans)

        fo.truncate()
        fo.close()

        self.generate_partcldataG(fname='partcl.data',Flag_initG = Flag_initG)
    ##############
    def generate_partcldataG(self, fname='partcl.data',Flag_initG = 1):
        #print ('aaaa')
        # plot#############################################################################
        if Flag_initG:
            #print ('aaaa')
            fileDatatmp = np.array(np.loadtxt(fname, comments='#', delimiter=','))

            fileDatatmp12 = fileDatatmp[(fileDatatmp[:, 5] < 0.005) & (fileDatatmp[:, 5] > 0.0042)]
            fileDatatmp11 = fileDatatmp[(fileDatatmp[:, 5] < 0.0042) & (fileDatatmp[:, 5] > 0.004)]
            fileDatatmp10 = fileDatatmp[(fileDatatmp[:, 5] < 0.004) & (fileDatatmp[:, 5] > 0.0035)]
            fileDatatmp9 = fileDatatmp[(fileDatatmp[:, 5] < 0.0035) & (fileDatatmp[:, 5] > 0.0032)]
            fileDatatmp8 = fileDatatmp[(fileDatatmp[:, 5] < 0.0032) & (fileDatatmp[:, 5] > 0.0029)]
            fileDatatmp7 = fileDatatmp[(fileDatatmp[:, 5] < 0.0029) & (fileDatatmp[:, 5] > 0.0025)]
            fileDatatmp6 = fileDatatmp[(fileDatatmp[:, 5] < 0.0025) & (fileDatatmp[:, 5] > 0.002)]
            fileDatatmp5 = fileDatatmp[(fileDatatmp[:, 5] < 0.002)  & (fileDatatmp[:, 5] > 0.00175)]
            fileDatatmp4 = fileDatatmp[(fileDatatmp[:, 5] < 0.00175)& (fileDatatmp[:, 5] > 0.0012)]
            fileDatatmp3 = fileDatatmp[(fileDatatmp[:, 5] < 0.0012) & (fileDatatmp[:, 5] > 0.001)]

            fileDatatmp3 = fileDatatmp


            # 6D distribution
            axisfont = {'size': 25, 'color': '0', 'weight': 'normal', }  # 'serif', 'name':'arial'

            plt.figure(figsize=(15, 10))
            plt.figure(1)

            plt.title('logit')
            plt.subplot(221)
            plt.grid(b=True, which='major', color='0.7', linestyle='-')
            plt.xlabel('X [mm]', **axisfont)
            plt.ylabel('Xp [mrad]', **axisfont)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.axis([-10, 10, -100, 100])
            #plt.scatter(10 * fileDatatmp[:, 0], 1000 * fileDatatmp[:, 1], linewidth=2.0, color='indigo', alpha=1,
            #s=0.2)
            plt.scatter(10 * fileDatatmp3[:, 0], 1000 * fileDatatmp3[:, 1], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp4[:, 0], 1000 * fileDatatmp4[:, 1], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp5[:, 0], 1000 * fileDatatmp5[:, 1], linewidth=2.0, color='r', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp6[:, 0], 1000 * fileDatatmp6[:, 1], linewidth=2.0, color='c', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp7[:, 0], 1000 * fileDatatmp7[:, 1], linewidth=2.0, color='m', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp8[:, 0], 1000 * fileDatatmp8[:, 1], linewidth=2.0, color='y', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp9[:, 0], 1000 * fileDatatmp9[:, 1], linewidth=2.0, color='k', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp10[:, 0], 1000 * fileDatatmp10[:, 1], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp11[:, 0], 1000 * fileDatatmp11[:, 1], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp12[:, 0], 1000 * fileDatatmp12[:, 1], linewidth=2.0, color='r', alpha=1, s=0.2)


            plt.subplot(222)
            plt.grid(b=True, which='major', color='0.7', linestyle='-')
            plt.xlabel('Y [mm]', **axisfont)
            plt.ylabel('Yp [mrad]', **axisfont)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.axis([-10, 10, -100, 100])
            #plt.scatter(10 * fileDatatmp[:, 2], 1000 * fileDatatmp[:, 3], linewidth=2.0, color='indigo', alpha=1,
            #s=0.2)
            
            plt.scatter(10 * fileDatatmp3[:, 2], 1000 * fileDatatmp3[:, 3], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp4[:, 2], 1000 * fileDatatmp4[:, 3], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp5[:, 2], 1000 * fileDatatmp5[:, 3], linewidth=2.0, color='r', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp6[:, 2], 1000 * fileDatatmp6[:, 3], linewidth=2.0, color='c', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp7[:, 2], 1000 * fileDatatmp7[:, 3], linewidth=2.0, color='m', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp8[:, 2], 1000 * fileDatatmp8[:, 3], linewidth=2.0, color='y', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp9[:, 2], 1000 * fileDatatmp9[:, 3], linewidth=2.0, color='k', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp10[:, 2], 1000 * fileDatatmp10[:, 3], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp11[:, 2], 1000 * fileDatatmp11[:, 3], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp12[:, 2], 1000 * fileDatatmp12[:, 3], linewidth=2.0, color='r', alpha=1, s=0.2)

            
            plt.subplot(223)
            plt.grid(b=True, which='major', color='0.7', linestyle='-')
            plt.xlabel('X [mm]', **axisfont)
            plt.ylabel('Y [mm]', **axisfont)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.axis('equal')
            plt.axis([-10, 10, -10, 10])
            #plt.scatter(10 * fileDatatmp[:, 0], 10 * fileDatatmp[:, 2], linewidth=2.0, color='indigo', alpha=1,
            #s=0.2)

            plt.scatter(10 * fileDatatmp3[:, 0], 10 * fileDatatmp3[:, 2], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp4[:, 0], 10 * fileDatatmp4[:, 2], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp5[:, 0], 10 * fileDatatmp5[:, 2], linewidth=2.0, color='r', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp6[:, 0], 10 * fileDatatmp6[:, 2], linewidth=2.0, color='c', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp7[:, 0], 10 * fileDatatmp7[:, 2], linewidth=2.0, color='m', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp8[:, 0], 10 * fileDatatmp8[:, 2], linewidth=2.0, color='y', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp9[:, 0], 10 * fileDatatmp9[:, 2], linewidth=2.0, color='k', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp10[:, 0], 10 * fileDatatmp10[:, 2], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp11[:, 0], 10 * fileDatatmp11[:, 2], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(10 * fileDatatmp12[:, 0], 10 * fileDatatmp12[:, 2], linewidth=2.0, color='r', alpha=1, s=0.2)

            # equivalent:
            # plt.scatter(x,y,s=80, c=z, marker=None, verts=verts)

            plt.subplot(224)
            plt.grid(b=True, which='major', color='0.7', linestyle='-')
            plt.xlabel('Phase [rad]', **axisfont)
            plt.ylabel('dE [MeV/u]', **axisfont)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            #plt.scatter(fileDatatmp[:, 4], fileDatatmp[:, 5], linewidth=2.0, color='indigo', alpha=1, s=0.2)
            plt.scatter(fileDatatmp3[:, 4], fileDatatmp3[:, 5], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(fileDatatmp4[:, 4], fileDatatmp4[:, 5], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(fileDatatmp5[:, 4], fileDatatmp5[:, 5], linewidth=2.0, color='r', alpha=1, s=0.2)
            plt.scatter(fileDatatmp6[:, 4], fileDatatmp6[:, 5], linewidth=2.0, color='c', alpha=1, s=0.2)
            plt.scatter(fileDatatmp7[:, 4], fileDatatmp7[:, 5], linewidth=2.0, color='m', alpha=1, s=0.2)
            plt.scatter(fileDatatmp8[:, 4], fileDatatmp8[:, 5], linewidth=2.0, color='y', alpha=1, s=0.2)
            plt.scatter(fileDatatmp9[:, 4], fileDatatmp9[:, 5], linewidth=2.0, color='k', alpha=1, s=0.2)
            plt.scatter(fileDatatmp10[:, 4], fileDatatmp10[:, 5], linewidth=2.0, color='b', alpha=1, s=0.2)
            plt.scatter(fileDatatmp11[:, 4], fileDatatmp11[:, 5], linewidth=2.0, color='g', alpha=1, s=0.2)
            plt.scatter(fileDatatmp12[:, 4], fileDatatmp12[:, 5], linewidth=2.0, color='r', alpha=1, s=0.2)
            

            plt.tight_layout()
            plt.savefig('A_phase plot.png', dpi=80)

            # plt.show()
            plt.clf()  # Clear the figure for 


#######################################
if __name__ == '__main__':
    track = exeClass()
    track.TRACK_dist_gen6D(PNum=[10001],
                           Ekmevu=[138.8],  # accelerated kinetic energy in MeV
                           kEsprd=[0.0005,0.0005,0.0005],
                           A=[238],
                           Q=[90],
                           M=[238],
                           emitx=[0.066 * 10 ** (-5)],
                           alphax=[0],
                           betax=[0.065],
                           emity=[0.066 * 10 ** (-5)],
                           alphay=[0],
                           betay=[0.065],
                           emitz=[0.066 * 10 ** (-5)],
                           alphaz=[0],
                           betaz=[0.065],
                           rho13=[0],
                           rho14=[0],
                           rho15=[0],
                           rho16=[0],
                           rho23=[0],
                           rho24=[0],
                           rho25=[0],
                           rho26=[0],
                           rho35=[0],
                           rho36=[0],
                           rho45=[0],
                           rho46=[0],
                           sigmal = p['moment1'].reshape(1,7,7),
                           Flag_DCbeam=0, # 1: DC beam, 2: bunched, 3: point beam in phi-E space
                           Flag_initG=1,
                           gene_type=0)

    sys.exit()
