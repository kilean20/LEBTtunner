__all__ = ['_decision_PVs', '_target_ref', 'LEBT_tunner']


import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from datetime import datetime
import os
from .dictClass import dictClass
from .package_path import package_path
from .utils import is_notebook
if is_notebook():
    from IPython.display import display, clear_output


# PVs for beam current measure
_current_PVs = ['FE_LEBT:BCM_D0989:AVGPK_RD',
                'FE_MEBT:BCM_D1055:AVGPK_RD',
                'FE_MEBT:FC_D1102:PKAVG_RD',
                'FE_ISRC1:HVP_D0679:I_RD' ]


# target
#=== 2022/01/17 36Ar10+ ===
target  = { 'FE_MEBT:BPM_D1056:XPOS_RD':0.0,
            'FE_MEBT:BPM_D1056:YPOS_RD':0.0,
            'FE_MEBT:BPM_D1056:PHASE_RD':78.98,
            
            'FE_MEBT:BPM_D1072:XPOS_RD':0.0,
            'FE_MEBT:BPM_D1072:YPOS_RD':0.0,
            'FE_MEBT:BPM_D1072:PHASE_RD':-26.71,
            
            'FE_MEBT:BPM_D1094:XPOS_RD':0.0,
            'FE_MEBT:BPM_D1094:YPOS_RD':0.0,
            'FE_MEBT:BPM_D1094:PHASE_RD':-19.41 }

target = pd.DataFrame(target,index=['target'])


weights = { 'FE_MEBT:BPM_D1056:XPOS_RD':1.0,
            'FE_MEBT:BPM_D1056:YPOS_RD':1.0,
            'FE_MEBT:BPM_D1056:PHASE_RD':1.0,
            
            'FE_MEBT:BPM_D1072:XPOS_RD':0.4,
            'FE_MEBT:BPM_D1072:YPOS_RD':0.4,
            'FE_MEBT:BPM_D1072:PHASE_RD':1.0,
            
            'FE_MEBT:BPM_D1094:XPOS_RD':0.2,
            'FE_MEBT:BPM_D1094:YPOS_RD':0.2,
            'FE_MEBT:BPM_D1094:PHASE_RD':1.0 }

weights = pd.DataFrame(weights,index=['weights'])
_target_ref = pd.concat([target,weights])
_target_ref.loc["weights"] = _target_ref.loc["weights"]/_target_ref.loc["weights"].sum()

# Decision PVs
_decision_PVs = ["FE_LEBT:DCH_D0979:I_CSET",
                 "FE_LEBT:DCV_D0979:I_CSET",
                 "FE_LEBT:DCH_D0992:I_CSET",
                 "FE_LEBT:DCV_D0992:I_CSET"]


# LEBT_tunner
class LEBT_tunner():
    def __init__(self, 
                 target_ref = None, 
                 decision_PVs = None, 
                 current_PVs = None, 
                 display_plots = False,
                 QperA = 17/78, 
                 virtual = False, virtual_kwarg={}):
        '''
        decision_PVs: python list of PV for tunning (i.e. for decision making)
        target_ref: DataFrame specifying optimization targets of BPM PVs and weghts for each PV
        '''  
        def briefer(PV):
            name = re.findall("\w+_D\d\d\d\d",PV)[0]
            family = re.findall("_D\d\d\d\d:\w+_",PV)[0][6:-1]
            return name+family
        
        if decision_PVs == None:
            decision_PVs = _decision_PVs
        for PV in decision_PVs:
            assert PV[-4:] == "CSET", "invalid decision_PVs"
        self.decision_PVs = decision_PVs
        self.decision_PVs_brief = [briefer(PV) for PV in decision_PVs]
        
        if target_ref == None:
            target_ref = _target_ref
        target_ref.loc["weights"] = target_ref.loc["weights"]/target_ref.loc["weights"].sum()
        self.target_ref = target_ref
        self.target_PVs = list(target_ref.columns)
        self.target_PVs_brief = [briefer(PV) for PV in target_ref.columns]
        
        if current_PVs == None:
            current_PVs = _current_PVs
        tmp = ['FE_MEBT:FC_D1102:PKAVG_RD',
               'FE_ISRC1:HVP_D0679:I_RD' ]
        self.current_PVs = [PV for PV in current_PVs if PV not in tmp]
        # flowing order of beam current PVs list are important. They are used in loss function for good transmission
        self.current_PVs = self.current_PVs + ['FE_MEBT:FC_D1102:PKAVG_RD',
                                               'FE_ISRC1:HVP_D0679:I_RD' ]
        self.current_PVs_brief = [briefer(PV) for PV in self.current_PVs]
        
        self.QperA = QperA
        
        self.virtual = virtual
        if virtual:
            from . import TRACK_virtual_accelerator as _virtual_AC
            self.virtual_AC = _virtual_AC.LEBT(QperA = self.QperA, **virtual_kwarg)
            self.virtual_AC._simulated = True
            self.caget = self.virtual_AC.caget
            self.caput = self.virtual_AC.caput
        else:
            from phantasy import caget,caput
            self.caget = caget
            self.caput = caput
           
        # check if normalization info (made from archived data) present for decision_PVs
        LEBT_statistical_info = pd.read_json(package_path+'/LEBT_statistical_info.json')
        for PV in self.decision_PVs:
            if not PV+"*Q/A" in LEBT_statistical_info.columns:
                raise ValueError("the decision PV (",PV ,") is not present in 'LEBT_statistical_info.json' file. Use other devision PVs ")
        
        self.decision_stat = LEBT_statistical_info[[PV+"*Q/A" for PV in self.decision_PVs]]/self.QperA
        def mapper(tmp):
            return tmp.replace("*Q/A","")
        self.decision_stat.rename(mapper,axis="columns",inplace=True)
        self.target_stat = LEBT_statistical_info[self.target_PVs]
        
        self.x0 = self.read_normalized_decision_values()
        
        self.history = dictClass()
        self.history.x_head = self.decision_PVs_brief
        self.history.x = []
        self.history.y_head= self.target_PVs_brief + self.current_PVs_brief
        self.history.y = []
        self.history.y_std = []
        self.history.loss = []
        self.history.loss_runtime = []
        
        self.date = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if display_plots:
            self.init_plot()
            self.callbacks = [self.re_plot]
        else:
            self.callbacks = []

        
        
    def is_signal_stationary(self,data):
        '''
        check if a time series data is stationary 
        (assume time interval is short so that signal is stationary unless a linear fit is succesful with non-negligble slope )
        '''
        #t0=time.time()
        T,N = data.shape
        residue = np.zeros(N)
        a1 = np.zeros(N)
        x = np.linspace(0, 1, T)
        #plt.figure(figsize=(16,4))
        for i in range(N):
            #if i < int(N/4)*4:
                #plt.subplot(int(N/4),4,i+1)
                #plt.plot(data[:,i])            
            a1[i],a0 = np.polyfit(x,data[:,i],deg=1)
            residue[i] = np.std(data[:,i]-a1[i]*x-a0)
        #plt.show()
        noise = np.std(data,axis=0)
        if np.all(noise < 1e-12):
            is_stationary = True
        else:
            isFitTrusted = residue <= 0.9*noise
            isSlopeLarge = np.abs(a1) > 1.1*residue
            #t1=time.time()
            #print("isFitTrusted: ",isFitTrusted)
            #print("isSlopeLarge: ",isSlopeLarge)
            #print("isStaionary: ",np.logical_not( isFitTrusted *isSlopeLarge ))
            is_stationary = np.sum(np.logical_not( isFitTrusted *isSlopeLarge )) > 0.8*len(isFitTrusted)  # 80% of the PVs are stationary
        return is_stationary 
    
    
    def read_normalized_decision_values(self):
        x0 = np.array([self.caget(pv.replace("CSET","RSET")) for pv in self.decision_PVs])
        if self.virtual:
            x0 = x0*self.decision_stat.loc["standard_deviation"].values
        else:
            x0 = (x0-self.decision_stat.loc["average"].values)*self.decision_stat.loc["standard_deviation"].values
        return x0
    
    
    def wait_RD_reach_CSET(self, PVs, Refs, Tolerances):
        time.sleep(0.1)
        V = np.array([self.caget(PV.replace("CSET","RD")) for PV in PVs])
        iitr = 0
        while(np.max(np.abs(V-Refs)/Tolerances) > 1.) :
            time.sleep(0.1)
            V = np.array([self.caget(PV.replace("CSET","RD")) for PV in PVs])
            iitr +=1
            if iitr>20:
                print("current / voltage ramping not yet stabilized after 2 secconds. Ignoreing CSET vs RD stablization... ")
                break
        return
    
    
    def wait_RD_stablize_then_read(self, PVs, min_t=1.2, max_t=4):
        '''
        mint_t : minimum time window (in second) for averageging
        '''
        t0 = time.time()
        dt = 0.1  # measure period. need to verify with the maximum measurable frequency
        n = int(min_t/dt)
        V = np.zeros([n,len(PVs)])
        for i in range(n):
            V[i,:] = [self.caget(PV) for PV in PVs]
            time.sleep(dt) 
    
        while(time.time()-t0 < max_t):
            if self.is_signal_stationary(V):
                break
            # empty 20% of the head of time series data and re-fill new time data in the tail
            m = int(0.2*n)+1
            V[:n-m,:] = V[m:,:]
            for j in range(m):
                V[n-m+j,:] = [self.caget(PV) for PV in PVs]
                time.sleep(dt)            
       
        #t1 = time.time()
        #m = int(max_t - t1)
        #if m > 1:
        #    concatV = np.zeros([m,len(PVs)])
        #    for i in range(m):
        #        concatV[i,:] = [self.caget(PV) for PV in PVs]
        #        time.sleep(dt) 
        #    V = np.vstack((V,concatV))
        
        return np.mean(V,axis=0) #, np.std(V,axis=0)
        
    
    def set_n_measure(self,decision_values):
        self.history.x.append(decision_values)
        for i,PV in enumerate(self.decision_PVs):
            self.caput(PV, decision_values[i])
        # wailt CSET ramping complete
        self.wait_RD_reach_CSET(self.decision_PVs, decision_values, 0.1*self.decision_stat.loc["standard_deviation"].values)
        # wait target RD stablize and measure 
#         RD_mean, RD_std = self.wait_RD_stablize_then_read(self.target_PVs + self.current_PVs , min_t=2, max_t=5)
        RD_mean = self.wait_RD_stablize_then_read(self.target_PVs + self.current_PVs , min_t=1.2, max_t=4)
        self.history.y.append(list(RD_mean))
#         self.history.y_std.append(RD_std)
        return RD_mean
    
    
    def loss_func(self,normalized_decision_values,weight_on_transmission=0.1):
        t0 = time.time()
        if self.virtual:
            x = normalized_decision_values/self.decision_stat.loc["standard_deviation"].values 
        else:
            x = normalized_decision_values/self.decision_stat.loc["standard_deviation"].values + self.decision_stat.loc["average"].values 
        y_measure = self.set_n_measure(x)
        if self.virtual and np.all(y_measure==0):
            loss = np.inf
        else:
            y_true = self.target_ref.loc["target"].values
            y_std = self.target_stat.loc["standard_deviation"].values
            y_weight = self.target_ref.loc["weights"].values
            
            n_target = len(self.target_PVs)
            loss = np.sum(np.sqrt(((y_true - y_measure[:n_target])/y_std)**2)*y_weight)
            loss += weight_on_transmission*y_measure[-2]/y_measure[-1]
        self.history.loss.append(loss)

        # call callback functions if any
        for func in self.callbacks:
            func()
        t1 = time.time()
        self.history.loss_runtime.append(t1-t0)
        return loss 
        
        
    def init_plot(self):
        self.plot = dictClass()
        self.plot.fig, self.plot.subs = plt.subplots(5, 1, figsize=(12, 14))
        
        # subplot[0] -- decision parameters
        # self.plot.subs[0].set_title('decision parameters')
        self.plot.subs[0].yaxis.grid()
        self.plot.sub0 = [
                           self.plot.subs[0].plot([],[],label=label)[0]
                           for label in self.decision_PVs_brief
                          ]
        self.plot.subs[0].legend(ncol=int(len(self.decision_PVs_brief)/2))
        
        # subplot[1] -- X, Y
#         self.plot.subs[1].set_title('BPM position')
        self.plot.subs[1].yaxis.grid()
        self.plot.subs[1].set_ylabel('$x,y (mm)$')
#         self.plot.subs[1].set_ylim(-3,3)
        self.plot.sub1 = [
                            self.plot.subs[1].plot([],[],label=label)[0] 
                            for label in self.target_PVs_brief if "POS" in label
                         ]
        self.plot.subs[1].legend(ncol=3)
        
        # subplot[2] -- phase differences
#         self.plot.subs[2].set_title('BPM phase differences')
        self.plot.subs[2].yaxis.grid()
        self.plot.subs[2].set_ylabel('$\Delta \phi (deg)$')
#         self.plot.subs[2].set_ylim(-3,3)
        self.plot.sub2 = [
                            self.plot.subs[2].plot([],[],label=label)[0]
                            for label in self.target_PVs_brief if "PHASE" in label
                         ]
        self.plot.subs[2].legend(ncol=3)
        
        # subplot[3] -- beam current 
#         self.plot.subs[3].set_title('Beam Current (Peak)')
        self.plot.subs[3].yaxis.grid()
        self.plot.subs[3].set_ylabel('(uA)')
        self.plot.sub3 = [
                            self.plot.subs[3].plot([],[],label=label)[0]
                            for label in self.current_PVs_brief 
                         ]
        self.plot.subs[3].legend(ncol=int(len(self.current_PVs_brief)/2))
        
        # subplot[4] -- loss
#         self.plot.subs[4].set_title('loss')
        self.plot.subs[4].yaxis.grid()
        self.plot.subs[4].set_ylabel('loss')
        self.plot.sub4 = self.plot.subs[4].plot([],[],label="loss")[0]
        self.plot.subs[4].legend()
        
        self.plot.fig.tight_layout()
        if is_notebook():
            display(self.plot.fig)
            # clear_output(wait = True)
        else:
            self.plot.fig.show()
            self.plot.fig.canvas.draw()
        
        
    def re_plot(self):
        if is_notebook():
            clear_output(wait = True)
        niter = len(self.history.loss)
        x = np.arange(niter)
        y = np.array(self.history.x)
        for i,handle in enumerate(self.plot.sub0):
            handle.set_data(x,y[:,i])
            
        k = 0
        y = np.array(self.history.y)
        for i,pv in enumerate(self.target_PVs):
            if "POS" in pv:
                self.plot.sub1[k].set_data(x,y[:,i])
                k += 1
                
        k = 0
        for i,pv in enumerate(self.target_PVs):
            if "PHASE" in pv:
                self.plot.sub2[k].set_data(x,y[:,i]-self.target_ref[pv].target)
                k += 1
                
        n = len(self.target_PVs)
        for i,pv in enumerate(self.current_PVs):
            self.plot.sub3[i].set_data(x,y[:,i+n])
            
        self.plot.sub4.set_data(x,self.history.loss) 
        
        for i in range(5):
            self.plot.subs[i].relim()
            self.plot.subs[i].autoscale()
        
        if is_notebook():
            display(self.plot.fig)
            # clear_output(wait = True)
        else:
            self.plot.fig.canvas.draw()
            self.plot.fig.canvas.flush_events()
        if not os.path.isdir("./data"):
            os.mkdir("./data")
        self.plot.fig.savefig('./data/LEBT_trajectory_scan_'+self.date+'.png', dpi=(144), bbox_inches='tight')
        
        
    def plot_results(self,fname):
        niter = len(self.history.loss)
        x = np.arange(niter)
        
        plot = dictClass()
        plot.fig, plot.subs = plt.subplots(5, 1, figsize=(15, 18))       
        
        # subplot[0] -- decision parameters
        plot.subs[0].set_title('decision parameters')
        plot.subs[0].yaxis.grid()
        y = np.array(self.history.x)
        for i,label in enumerate(self.decision_PVs_brief):
            plot.subs[0].plot(x,y[:,i],label=label)
        plot.subs[0].legend(ncol=int(len(self.decision_PVs_brief)/2))
        
        # subplot[1] -- X, Y
#         plot.subs[1].set_title('BPM position')
        plot.subs[1].yaxis.grid()
        plot.subs[1].set_ylabel('$x,y (mm)$')
#         plot.subs[1].set_ylim(-3,3)
        y = np.array(self.history.y)
        for i,label in enumerate(self.target_PVs_brief):
            if "POS" in label:
                plot.subs[1].plot(x,y[:,i],label=label)
        plot.subs[1].legend(ncol=3)
        
        # subplot[2] -- phase differences
#         plot.subs[2].set_title('BPM phase differences')
        plot.subs[2].yaxis.grid()
        plot.subs[2].set_ylabel('$\Delta \phi (deg)$')
#         plot.subs[2].set_ylim(-3,3)
        for i,label in enumerate(self.target_PVs_brief):
            if "PHASE" in label:
                plot.subs[2].plot(x,y[:,i]-self.target_ref[self.target_PVs[i]].target,label=label)
        plot.subs[2].legend(ncol=3)
        
        # subplot[3] -- beam current 
#         plot.subs[3].set_title('Beam Current (Peak)')
        plot.subs[3].yaxis.grid()
        plot.subs[3].set_ylabel('(uA)')
        n = len(self.target_PVs)
        for i,label in enumerate(self.current_PVs_brief):
            plot.subs[3].plot(x,y[:,i+n],label=label)
        plot.subs[3].legend(ncol=int(len(current_PVs)/2))
        
        # subplot[4] -- loss
#         plot.subs[4].set_title('loss')
        plot.subs[4].yaxis.grid()
        plot.subs[4].set_ylabel('loss')
        plot.subs[4].plot(x,self.history.loss,label="loss")[0]
        plot.subs[4].legend()
        
        plot.fig.tight_layout()
        plot.fig.savefig(fname, dpi=(144), bbox_inches='tight')
