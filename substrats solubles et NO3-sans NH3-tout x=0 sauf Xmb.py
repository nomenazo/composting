#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from pandas import Series,DataFrame
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from math import exp,expm1
from math import sqrt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from math import log


# In[2]:


pip install xlrd


# In[6]:


class Compostage:
    def __init__(self, data_path, data_T, technology, aer, duration, MB, TB, MA, TA, MF, TF, Xa, xtemp):
        
        self.technology = technology
        self.aer =aer
        self.duration = duration
        self.MB = MB
        self.TB = TB
        self.MA = MA
        self.TA = TA
        self.MF = MF
        self.TF = TF
        self.Xa = Xa
        self.xtemp = xtemp
        
        with pd.ExcelFile(data_path) as f:
            
            #variables d'état 
            composition = pd.read_excel(data_path, 'Variables', index_col=0).fillna(0.0)
            mass = [composition[c] for c in composition.columns]
            self.mass = mass
            self.va=mass[0]
            
            #Paramètres cinétiques
            parameters = pd.read_excel(data_path, 'kinetics', index_col=0).fillna(0.0)
            Valeurs = [parameters[c] for c in parameters.columns]
            self.Valeurs = Valeurs
            self.k = Valeurs[0] 
            
            #paramètres stochiométriques
            stoech = pd.read_excel(data_path, 'stoe', index_col=0).fillna(0.0)
            val = [stoech[c] for c in stoech.columns]
            self.val=val
            self.st = val[0]
            
            #biomass yield
            biomass = pd.read_excel(data_path,'Xyield', index_col=0).fillna(0.0)
            self.Y_S = 1/biomass['Y(X/S)']
            self.Y_O2 = 1/biomass['Y(X/O2)']
            self.Y_NH3 = 1/biomass['Y(X/NH3)']
            self.Y_CO2 = 1/biomass['Y(X/CO2)']
            self.Y_W = 1/biomass['Y(X/H2O)']
            
            
            
            #paramètres d'opérations
            oper = pd.read_excel(data_path, 'tech', index_col=0).fillna(0.0)
            self.cond = oper['u']
            self.area = oper['A']
            self.V = oper['Vreacteur']
            
            
            #aeration
            aeration = pd.read_excel(data_path, 'aer', index_col=0).fillna(0.0)
            self.air = aeration['Qair']
            
            #paramètres régionaux
            reg = pd.read_excel(data_path, 'reg', index_col=0).fillna(0.0)
            self.Temp = reg['Ta'] #température ambiant
        
        with pd.ExcelFile(data_T) as e:
            temp = pd.read_excel(data_T, 'Temp', index_col=0)
            self.Ti = [temp[c] for c in temp.columns]
            self.Tint =self.xtemp + np.array(self.Ti[0])
    
    def U(self):
        if self.technology == 'home composting,wood':
            return self.cond['hc wood']
        elif self.technology == 'home composting,plastic':
            return self.cond['hc plastic']
        elif self.technology == 'community composting':
            return self.cond['CC']
        elif self.technology =='industrial composting':
            return self.cond['IC']
        elif self.technology =='test_fitting':
            return self.cond['test']
        
    
    def A(self):
        if self.technology == 'home composting,wood':
            return self.area['hc wood']
        elif self.technology == 'home composting,plastic':
            return self.area['hc plastic']
        elif self.technology == 'community composting':
            return self.area['CC']
        elif self.technology =='industrial composting':
            return self.area['IC']
        elif self.technology =='test_fitting':
            return self.area['test']
        
    def Vreact(self):
        if self.technology == 'home composting,wood':
            return self.V['hc wood']
        elif self.technology == 'home composting,plastic':
            return self.V['hc plastic']
        elif self.technology == 'community composting':
            return self.V['CC']
        elif self.technology =='industrial composting':
            return self.area['IC']
        elif self.technology =='test_fitting':
            return self.V['test']
        
    
    def Qair(self):
        if self.aer == 'Passive aeration':
            return self.air['Passive']
        elif self.aer == 'Forced aeration':
            return self.air['Forced']
    
    def resolution(self):
        def system(t, Y):
            self.C=Y[0]
            self.P=Y[1]
            self.L=Y[2]
            self.H=Y[3]
            self.CE=Y[4]
            self.LG=Y[5]
            self.Xi=Y[6]
            self.Sc=Y[7]
            self.Sp=Y[8]
            self.Sl=Y[9]
            self.Sh=Y[10]
            self.Slg=Y[11]
            self.Xmb=Y[12]
            self.Xtb=Y[13]
            self.Xma=Y[14]
            self.Xta=Y[15]
            self.Xmf=Y[16]
            self.Xtf=Y[17]
            self.Xa = Y[18]
            self.Xdb=Y[19]
            self.CO2=Y[20]
            #self.NH3=Y[21]
            self.NO3=Y[21]
            self.N2=Y[22]
            self.N2O = Y[23]
            self.CH4gen = Y[24]
            self.CH4oxi = Y[25]
            self.CH4emis = Y[26]
            #self.T= Y[27]
            self.W = Y[27]
            self.Sn = Y[28]
            #self.NH3 = Y[29]
            
            #self.T =  (self.Tint[int(t)])
            
            #I. Growth limiting functions 
            #I.1. Substrate
            #a.Carbohydrates
            self.fSc=(self.Sc/self.W)/(self.k['ks']+(self.Sc/self.W)) #np.exp(logsumexp #variable (W) has to be changed to self.H2O
            #b.Proteins
            self.fSp=(self.Sp/self.W)/(self.k['ks']+(self.Sp/self.W))
            #c.Lipids
            self.fSl=(self.Sl/self.W)/(self.k['ks']+(self.Sl/self.W))                                               
            #d.Hemicelluloses
            self.fSh=(self.Sh/self.W)/(self.k['ks']+(self.Sh/self.W))
            #e.Lignin
            self.fSlg=(self.Slg/self.W)/(self.k['ks']+(self.Slg/self.W))
            
            #f.nitrate
            self.fNO3= (self.NO3/self.W)/(self.k['kNO3']+(self.NO3/self.W)) #np.exp(logsumexp
            
            #g.ammonium for heterotrophic activities
           # self.fSn= (self.Sn/self.W)/(self.k['kNH4']+(self.Sn/self.W))
            
            #h. ammonium for nitrification
            #self.fSn_nit = (self.Sn/self.W)/((10**((0.051*self.T)-7.158)) +(self.Sn/self.W))
            
            #I.2. Substrate availability
            ##for bacteries:
            #self.faB_C = self.Sc/(self.Sc+self.Sl) ###availability of Sc
            #print (self.Sc/(self.Sc+self.Sl))
            #self.faB_L =self.Sl/(self.Sc+self.Sl)  ###availability of Sc
            #for actinomycetes
            #self.faA_C =self.Sc/(self.Sc+self.Sl+self.Sh) #np.exp(logsumexp
            #self.faA_L= self.Sl/(self.Sc+self.Sl+self.Sh)
            #self.faA_H=self.Sh/(self.Sc+self.Sl+self.Sh)
            ##for fungi
            #self.faF_C =self.Sc/(self.Sc+self.Sl+self.Sh+self.Slg)
            #self.faF_L=self.Sl/(self.Sc+self.Sl+self.Sh+self.Slg)
            #self.faF_H= self.Sh/(self.Sc+self.Sl+self.Sh+self.Slg)
            #self.faF_LG= self.Slg/(self.Sc+self.Sl+self.Sh+self.Slg)
            
            
             #I.3. Temperature
                        
            #I.3.1. For heterotrophic activities
            #self.flT1 = (((self.T - 44) *(self.T - 5.1)*(self.T - 5.1))/((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 35.4))-((35.4- 44)*(35.4+5.1-(2*(self.T)))))))
            #print ((((self.T - 44) *(self.T - 5.1)*(self.T - 5.1))/((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 35.4))-((35.4- 44)*(35.4+5.1-(2*(self.T))))))))
                                    
           # flT1 = np.exp(logsumexp(((self.T -317.15) *((self.T -278.15)*(self.T -278.15))))/\
                                  #  ((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 308.55))-((35.4- 44)*(308.55+278.15-(2*self.T))))))
            #self.flT2 = ((self.T - 65.5) *(self.T - 30.8)*(self.T - 30.8))/((57.2 - 30.8)*(((57.2 - 30.8)*(self.T - 57.2))-((57.2- 65.5)*(57.2+30.8-(2*(self.T))))))                              
                                    
          #  flT2 = np.exp(logsumexp(((self.T - 338.15) *((self.T - 303.95)*(self.T - 303.95))))/\
                                   # ((57.2 - 30.8)*(((57.2 - 30.8)*(self.T - 330.35))-((57.2- 65.5)*(330.35+303.95-(2*self.T))))))
            
            #I.3.2. For denitrification
            #self.fT_denit = exp((self.T - 20)*log(2.1)/10)
            
            #II. Aerobic degradation
            #II.1-Hydrolysis phase
            v1 = self.k['kh1C']*self.Xmb*(self.C/((self.k['khS']*self.Xmb)+self.C))
            v2 = self.k['kh2P']* self.Xmb*(self.P/((self.k['khS']*self.Xmb)+self.P)) #(np.exp(logsumexp
            v3 = self.k['kh3L']*self.Xmb*(self.L/((self.k['khS']*self.Xmb)+self.L))
            v4 = self.k['kh4C']*self.Xtb*(self.C/((self.k['khS']*self.Xtb)+self.C))
            v5= self.k['kh5P']*self.Xtb*(self.P/((self.k['khS']*self.Xtb)+self.P))
            v6 = self.k['kh6L']*self.Xtb*(self.L/((self.k['khS']*self.Xtb)+self.L))
            v7 = self.k['kh7H']*self.Xta*(self.H/((self.k['khS']*self.Xta)+self.H))
            v8 = self.k['kh8CE']*self.Xtf*(self.CE/((self.k['khS']*self.Xtf)+self.CE))
            v9 = self.k['kh9LG']*self.Xtf*(self.LG/((self.k['khS']*self.Xtf)+self.LG))
            v10 = self.k['kh10H']*self.Xma*(self.H/((self.k['khS']*self.Xma)+self.H))
            v11 = self.k['kh11CE']*self.Xmf*(self.CE/((self.k['khS']*self.Xmf)+self.CE))
            v12 = self.k['kh12LG']*self.Xmf*(self.LG/((self.k['khS']*self.Xmf)+self.LG))
            
            #II.2-Microorganisms growth
            #2-1-Mesophilic Bacterias growth
            v13 = self.k['µmb']*self.Xmb  *self.fSc #*self.fSn  *self.flT1*self.faB_C
            v14 = self.k['µmb']*self.Xmb  * self.fSp  # *self.flT1
            v15 = self.k['µmb']*self.Xmb  *self.fSl  # *self.fSn  *self.flT1 *self.faB_L
            #2-2-Thermophilic bacterias growth
            v16 = self.k['µtb']*self.Xtb  * self.fSc  # *self.fSn   *self.flT2 *self.faB_C 
            v17 = self.k['µtb']*self.Xtb *self.fSp   #*self.flT2
            v18 = self.k['µtb']*self.Xtb *self.fSl  #*self.fSn    *self.flT2  *self.faB_L
            #2-3-Mesophilic Actinomycetes growth
            v19 = self.k['µma']*self.Xma   * self.fSc   #*self.fSn     *self.flT1  *self.faA_C
            v20 = self.k['µma']*self.Xma    *self.fSp     #*self.flT1
            v21 = self.k['µma']*self.Xma  *self.fSl  #*self.fSn     *self.flT1   *self.faA_L
            v22 = self.k['µma']*self.Xma   *self.fSh  #*self.fSn     *self.flT1   *self.faA_H
            #2-4-Thermophilic actinomycetes growth
            v23 = self.k['µta']*self.Xta  *self.fSc  #*self.fSn    *self.flT2  *self.faA_C
            v24 = self.k['µta']*self.Xta  *self.fSp  #*self.flT2
            v25 = self.k['µta']*self.Xta  *self.fSl  #*self.fSn      *self.flT2   *self.faA_L
            v26 = self.k['µta']*self.Xta  *self.fSh  #*self.fSn     *self.flT2    *self.faA_H
            #2-5-Mesophilic fungis growth
            v27 = self.k['µmf']*self.Xmf  * self.fSc  # *self.fSn    *self.flT1  *self.faF_C
            v28 = self.k['µmf']*self.Xmf  * self.fSp   # *self.flT1
            v29 = self.k['µmf']*self.Xmf  * self.fSl # *self.fSn  *self.flT1   *self.faF_L 
            v30 = self.k['µmf']*self.Xmf  *  self.fSh  #*self.fSn  *self.flT1   *self.faF_H 
            v31 = self.k['µmf']*self.Xmf  * self.fSlg  #*self.fSn  *self.flT1  *self.faF_LG
            #2-6- Thermophilic fungis growth
            v32 = self.k['µtf']*self.Xtf  * self.fSc  #*self.fSn      *self.flT2  *self.faF_C
            v33 = self.k['µtf']*self.Xtf   *self.fSp   # *self.flT2
            v34 = self.k['µtf']*self.Xtf   *self.fSl   #*self.fSn      *self.flT2  *self.faF_L 
            v35 = self.k['µtf']*self.Xtf    *self.fSh  #*self.fSn       *self.flT2  *self.faF_H
            v36 = self.k['µtf']*self.Xtf   *self.fSlg   #*self.fSn       *self.flT2 *self.faF_LG
            
             #II.3-Heterotroph Microorganisms death
            v37 = self.k['bmb']*self.Xmb
            v38 = self.k['btb']*self.Xtb
            v39 = self.k['bma']*self.Xma
            v40 = self.k['bta']*self.Xta
            v41 = self.k['bmf']*self.Xmf
            v42 = self.k['btf']*self.Xtf
            
            #II.4-Lysis of microorganisms
            v43 = self.k['kdec']*self.Xdb
            
            #III. Nitrification-denitrification
            #III.1-Nitrification
            v44 = self.k['µa']*self.Xa #* self.fSn_nit
            
            #III.2-Denitrification
            v45 = self.k['pmaxdenit']*self.NO3 *self.fNO3  #*self.fT_denit 
            
            #III.3-Autotroph microorganisms death
            v46 = self.k['ba']*self.Xa 
            
            #III.4. Emission of ammoniac
        
            #He_NH3 = exp(160.559 - (8621.06/(self.T+273))-(25.6767*log(self.T+273))+(0.035388*(self.T+273)))
            
            
            #NH3gas= (0.9*self.Sn * He_NH3)/ (0.08206 * (self.T+273))
            
            #Vgas = self.Vreact() - (self.va['TM']/self.k['rho_waste']) #0.020 
            #v47 = (self.Qair()*NH3gas)/(self.k['rho_air']*Vgas)
            
            #Global equations
            self.dC_dt= -v1 -v4
            self.dP_dt =   -v5 -v2 + (1-self.k['fi'])*(v43+v46)
            self.dL_dt = -v3-v6
            self.dH_dt= -v7-v10
            self.dCE_dt= -v8-v11
            self.dLG_dt= -v9-v12
            self.dXi_dt= self.k['fi']*(v43+v46)
            self.dSc_dt =  v1 +v4 +v11 +v8 - self.Y_S['r13']*(v13+v16+v19+v23+v27+v32)    
            self.dSp_dt =   v2  +v5 -self.Y_S['r13']*(v14+v17+v20+v24+v28+v33)
            self.dSl_dt = v3+v6-self.Y_S['r13']*(v15+v18+v21+v25+v29+v34)
            self.dSh_dt = v7+v10-self.Y_S['r13']*(v22+v26+v30+v35)
            self.dSlg_dt = v9+v12-self.Y_S['r13']*(v31+v36)
            self.dXmb_dt = v13+ v14 + v15 -v37
            self.dXtb_dt = v16 + v17 + v18 - v38
            self.dXma_dt = v19+v20 + v21 +v22- v39
            self.dXta_dt = v23 + v24 + v25 +v26- v40
            self.dXmf_dt = v27 + v28 + v29 +v30+ v31 - v41
            self.dXtf_dt = v32 + v33 + v34 +v35+v36 - v42
            
            self.dXdb_dt = v37+v38+v39+v40+v41+v42-v43
            
            self.dCO2_dt = (self.Y_CO2['r13'])*v13+(self.Y_CO2['r14'])*v14+(self.Y_CO2['r15'])*v15+            (self.Y_CO2['r16'])*v16+(self.Y_CO2['r17'])*v17+(self.Y_CO2['r18'])*v18+(self.Y_CO2['r19'])*v19+            (self.Y_CO2['r20'])*v20+(self.Y_CO2['r21'])*v21+(self.Y_CO2['r22'])*v22+(self.Y_CO2['r23'])*v23+            (self.Y_CO2['r24'])*v24+(self.Y_CO2['r25'])*v25+(self.Y_CO2['r26'])*v26+(self.Y_CO2['r27'])*v27+            (self.Y_CO2['r28'])*v28+(self.Y_CO2['r29'])*v29+(self.Y_CO2['r30'])*v30+(self.Y_CO2['r31'])*v31+            (self.Y_CO2['r32'])*v32+(self.Y_CO2['r33'])*v33+(self.Y_CO2['r34'])*v34+(self.Y_CO2['r35'])*v35+            (self.Y_CO2['r36'])*v36
            
            self.dXa_dt = (self.st['l']*v44)- v46
            
            self.dNO3_dt = (self.st['n']*v44)- v45
            
            self.dN2O_dt = self.k['pmaxdenit']*v45
            self.dN2_dt = (1-self.k['pmaxdenit'])*v45
            self.dCH4gen_dt = (((self.k['YCH4,Sc'])*(v1+v4+v8+v11))+((self.k['YCH4,Sp'])*(v2+v5))+((self.k['YCH4,Sl'])*(v3+v6))+                               ((self.k['YCH4,Sh'])*(v7+v10))+((self.k['YCH4,Slg'])*(v9+v12)))*(1/(1+(self.k['η']*500)))

            self.dCH4oxi_dt = self.k['vmax']*((self.CH4gen/self.W)/(self.k['Km']+(self.CH4gen/self.W)))*(500/(self.k['Ko2,ch4']+500))
                                                         
            self.dCH4emis_dt = self.dCH4gen_dt - self.dCH4oxi_dt
            
            self.dSn_dt = (-self.st['a']*v13)+((4-self.st['b'])*v14)-(self.st['c']*v15)-(self.st['a']*v16)+((4-self.st['b'])*v17)            -(self.st['c']*v18)-(self.st['a']*v19)+((4-self.st['b'])*v20)-(self.st['c']*v21)-(self.st['d']*v22)-(self.st['a']*v23)            +((4-self.st['b'])*v24)-(self.st['c']*v25)-(self.st['d']*v26)-(self.st['e']*v27)+((4-self.st['f'])*v28)-(self.st['g']*v29)            -(self.st['h']*v30)-(self.st['i']*v31)-(self.st['e']*v32)+((4-self.st['f'])*v33)-(self.st['g']*v34)-(self.st['h']*v35)            -(self.st['i']*v36)-v44 #-v47
            
            #self.dNH3_dt = v47
            
            self.dW_dt = (self.Y_W['r13'])*v13+(self.Y_W['r14'])*v14+(self.Y_W['r15'])*v15+            (self.Y_W['r16'])*v16+(self.Y_W['r17'])*v17+(self.Y_W['r18'])*v18+(self.Y_W['r19'])*v19+            (self.Y_W['r20'])*v20+(self.Y_W['r21'])*v21+(self.Y_W['r22'])*v22+(self.Y_W['r23'])*v23+            (self.Y_W['r24'])*v24+(self.Y_W['r25'])*v25+(self.Y_W['r26'])*v26+(self.Y_W['r27'])*v27+            (self.Y_W['r28'])*v28+(self.Y_W['r29'])*v29+(self.Y_W['r30'])*v30+(self.Y_W['r31'])*v31+            (self.Y_W['r32'])*v32+(self.Y_W['r33'])*v33+(self.Y_W['r34'])*v34+(self.Y_W['r35'])*v35+            (self.Y_W['r36'])*v36+(self.st['m']*v44)
            
            return [self.dC_dt,self.dP_dt,self.dL_dt,self.dH_dt,self.dCE_dt,self.dLG_dt,self.dXi_dt,self.dSc_dt,                    self.dSp_dt,self.dSl_dt,self.dSh_dt,self.dSlg_dt,self.dXmb_dt, self.dXtb_dt,self.dXma_dt,self.dXta_dt,                   self.dXmf_dt,self.dXtf_dt,self.dXa_dt,self.dXdb_dt,self.dCO2_dt,self.dNO3_dt,self.dN2O_dt,self.dN2_dt,                   self.dCH4gen_dt,self.dCH4oxi_dt,self.dCH4emis_dt,self.dW_dt, self.dSn_dt] #,self.dNH3_dt] #, self.flT1, self.flT2] #self.dT_dt,
        
        
       
        self.solution = solve_ivp(system, [0, self.duration], [self.va['G'],self.va['P'],self.va['L'],                                                self.va['HE'],self.va['CE'],self.va['LG'], self.va['Xi'],                                                    0,0,0,0,0,                                                     self.MB,self.TB,self.MA,self.TA,self.MF,                                                     self.TF,self.Xa, self.va['Xdb'], 0,                                                     0,0,0,0,0,0, self.va['W'], 0.000], method='RK45') #, t_eval = np.arange(0,self.duration,1)


   
        #for i in range(len(self.solution.y)):
            #self.solution.y[i][self.solution.y[i] < 0] = 0 
        
        return (self.solution.y[0])
    
    def emissions(self):
        variablesE = ['CO2','NO3', 'N2', 'N2O','CH4'] #,'NH3']
        functionsE= [self.solution.y[20], self.solution.y[21], self.solution.y[22], self.solution.y[23], self.solution.y[26]] #,self.solution.y[29]]
        num_variablesE = len(variablesE)
        
        for i in range(num_variablesE): 
            print ('total emission of', variablesE[i], ':' , trapz(functionsE[i], self.solution.t))
    
    def courbes(self):
        variablesC = ['C','Sc', 'Xmb', 'Xtb', 'Xdb','CO2']
        functions = [self.solution.y[0], self.solution.y[7], self.solution.y[12], self.solution.y[13], self.solution.y[19],                     self.solution.y[20], self.solution.y[26]]
        num_variablesC = len(variablesC)
        
        #for i in range(num_variables): 
            #print ('total emission of', variables[i], ':' , trapz(functions[i], self.solution.t))
        
        for i in range(num_variablesC):
            plt.subplot(3, 2, i + 1)
            plt.plot(self.solution.t, functions[i], label=variablesC[i])
            plt.xlabel('Temps (h)')
            plt.ylabel(f'{variablesC[i]} (kg/kg)')
            plt.legend()
            plt.grid(True)
            #plt.title(f'{variables[i]}')
            #plt.figure(figsize=(12, 4))
        
        
        plt.tight_layout()
        
        plt.show()
                                 


# In[7]:


Data1 = 'Data.xlsx'
Data2 = 'Data HC.xlsx'
Data3 = 'Data for fitting en h.xlsx'
Data4 = 'T Margaritis h.xlsx'
Data5 = 'Temp Margaritis+5.xlsx'


# In[8]:


biomass = pd.read_excel('Data for fitting en h.xlsx','Xyield', index_col=0).fillna(0.0)


# In[9]:


essai1 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 1000, 1e-1, 0, 0,0,0,0,0,0)


# In[49]:


essai1.resolution()


# In[158]:


plt.plot(essai1.solution.t[0:10],essai1.solution.y[9][0:10])


# In[11]:


essai7 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 1000, 1e-7,0, 0,0,0,0,0,0)


# In[12]:


essai7.resolution()


# In[13]:


essai7.solution


# In[26]:


plt.plot(essai7.solution.t, essai7.solution.y[12])


# In[15]:


essai7.courbes()


# In[16]:


essai7.solution.y[1]


# In[20]:


i = 5

plt.plot(essai7.solution.y[i])
plt.xlabel('Temps (h)')
plt.ylabel('sc(t)')
plt.legend()
plt.grid(True)
plt.title('Courbe pour Sc')


# In[21]:


essai7.emissions()


# In[29]:


essai3 = Compostage(Data3, Data5, 'test_fitting', 'Passive aeration', 2000, 0.001, 0, 0,0,0,0,0,0)


# In[30]:


essai3.resolution()


# In[26]:


essai3.solution


# In[48]:


plt.plot(essai3.solution.t,essai3.solution.y[12])


# In[33]:


essai3.courbes()


# In[42]:


essai4 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 1000, 0.0001, 0.0001, 0.0001,0.0001,0.0001,0.0001,0.0001,0)


# In[43]:


essai4.resolution()


# In[44]:


essai4.solution


# In[45]:


essai4.emissions()


# In[46]:


essai4.courbes()


# In[57]:


plt.plot(essai4.solution.t, essai4.solution.y[16])


# In[1171]:


essai5 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 480, 1e-5, 1e-5, 1e-5,1e-5,1e-5,1e-5,1e-5,0)


# In[1172]:


essai5.resolution()


# In[1173]:


essai5.emissions()


# In[1371]:


plt.plot(essai5.solution.y[7])


# In[1182]:


essai6 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 480, 1e-6, 1e-6, 1e-6,1e-6,1e-6,1e-6,1e-6,0)


# In[1183]:


essai6.resolution()


# In[1184]:


essai6.emissions()


# In[1376]:


plt.plot(essai6.solution.y[8])


# In[1205]:


essai8 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 480, 1e-8, 1e-8, 1e-8,1e-8,1e-8,1e-8,1e-8,0)


# In[1207]:


essai8.resolution()


# In[1208]:


essai8.emissions()


# In[1231]:


plt.plot(essai8.solution.y[20])


# In[1329]:


essai9 = Compostage(Data3, Data4, 'test_fitting', 'Passive aeration', 480, 1e-9, 1e-9, 1e-9,1e-9,1e-9,1e-9,1e-9,0)


# In[1330]:


essai9.resolution()


# In[1331]:


essai9.emissions()


# In[1332]:


plt.plot(essai9.solution.y[0])


# In[1340]:


plt.plot(essai9.solution.y[17])


# In[1334]:


plt.plot(essai9.solution.y[1])


# In[1342]:


plt.plot(essai3.solution.y[7][0:350])


# In[ ]:




