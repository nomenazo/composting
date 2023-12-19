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
from scipy.interpolate import interp1d
from scipy.optimize import minimize


# In[2]:


pip install xlrd


# In[56]:


class Compostage:
    def __init__(self, data_path, data_CO2, technology, aer, duration, MB, TB, MA, TA, MF, TF):
        
        self.technology = technology
        self.aer =aer
        self.duration = duration #durée de la simulation en h-1
        self.MB = MB #valeur initiale des bactéries mésophiles
        self.TB = TB #valeur initiale des bactéries thermophiles
        self.MA = MA #valeur initiale des actinomycètes mésophiles
        self.TA = TA #valeur initiale des actinomycètes thermophiles
        self.MF = MF #valeur initiale des fungi mésophiles
        self.TF = TF #valeur initiale des fungi thermophiles
       
        
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
            self.Y_S = 1/biomass['Y(X/S)'] #rendement de la biomasse sur le substrat kgX/kgS
            self.Y_O2 = 1/biomass['Y(X/O2)'] #rendement de la biomasse sur l'oxygène' kgX/kgO2
            self.Y_NH3 = 1/biomass['Y(X/NH3)'] #rendement de la biomasse sur le NH3 kgX/kgNH3
            self.Y_CO2 = 1/biomass['Y(X/CO2)'] #rendement de la biomasse sur le CO2 kgX/kgCO2
            self.Y_W = 1/biomass['Y(X/H2O)'] #rendement de la biomasse sur l'eau kgX/kgW
            self.Y_NO3=1/biomass['Y(X/NO3-)'] #rendement de la biomasse sur NO3 kgX/kgNO3
            
            
            
            #paramètres d'opérations
            oper = pd.read_excel(data_path, 'tech', index_col=0).fillna(0.0)
            self.cond = oper['u']
            self.area = oper['A']
           # self.V = oper['Vreacteur']
            
            
            #aeration
            aeration = pd.read_excel(data_path, 'aer', index_col=0).fillna(0.0)
            self.air = aeration['Qair']
            
            #paramètres régionaux
            reg = pd.read_excel(data_path, 'reg', index_col=0).fillna(0.0)
            self.Temp = reg['Ta'] #température ambiant
        
        with pd.ExcelFile(data_CO2) as e:
            co = pd.read_excel(data_CO2, 'Data', index_col=0)
            self.gas = [co[c] for c in co.columns]
            self.CO2 = np.array(self.gas[0]) 
    
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
        
    #def Vreact(self):
        #if self.technology == 'home composting,wood':
            #return self.V['hc wood']
        #elif self.technology == 'home composting,plastic':
            #return self.V['hc plastic']
      #  elif self.technology == 'community composting':
       #     return self.V['CC']
       # elif self.technology =='industrial composting':
       #     return self.area['IC']
       # elif self.technology =='test_fitting':
        #    return self.V['test']
        
    
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
            
            self.Xdb=Y[18]
            self.CO2=Y[19]
            self.CH4gen = Y[20]
            self.CH4oxi = Y[21]
            self.CH4emis = Y[22]
            self.W = Y[23]
            self.T= Y[24]
            #self.O2 = Y[30]
            #self.NH3= Y[29]
        
            
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
            
            
             #I.3. Temperature
                        
            #I.3.1. For heterotrophic activities
            #self.flT1 = (((self.T - 44) *(self.T - 5.1)*(self.T - 5.1))/((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 35.4))-((35.4- 44)*(35.4+5.1-(2*(self.T)))))))
            #print ((((self.T - 44) *(self.T - 5.1)*(self.T - 5.1))/((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 35.4))-((35.4- 44)*(35.4+5.1-(2*(self.T))))))))
            
            
            
            self.flT1 =1 # ((self.T -317.15) *((self.T -278.15)*(self.T -278.15)))/\
                                  #((35.4 - 5.1)*(((35.4 - 5.1)*(self.T - 308.55))-((35.4- 44)*(308.55+278.15-(2*self.T)))))
            #self.flT2 = ((self.T - 65.5) *(self.T - 30.8)*(self.T - 30.8))/((57.2 - 30.8)*(((57.2 - 30.8)*(self.T - 57.2))-((57.2- 65.5)*(57.2+30.8-(2*(self.T))))))                              
                                    
            self.flT2 =1 #((self.T - 338.15) *((self.T - 303.95)*(self.T - 303.95)))/\
                                   #((57.2 - 30.8)*(((57.2 - 30.8)*(self.T - 330.35))-((57.2- 65.5)*(330.35+303.95-(2*self.T)))))
            
            #modèle Higgins
            #self.flT = ((self.T - (71.6+273)) *(self.T - (5+273))*(self.T - (5+273)))/((58.6 - 5)*(((58.6 - 5)*(self.T - (58.6+273)))-((58.6- 71.6)*(58.6+273+5+273-(2*(self.T))))))
            
            
            #II. Aerobic degradation
            #II.1-Hydrolysis phase
            v1 = 0.0716*self.Xmb*(self.C/((self.k['khS']*self.Xmb)+self.C)) 
            v2 = 0.0676* self.Xmb*(self.P/((self.k['khS']*self.Xmb)+self.P)) 
            v3 = 0.0901*self.Xmb*(self.L/((self.k['khS']*self.Xmb)+self.L)) 
            v4 = 0.0001*self.Xtb*(self.C/((self.k['khS']*self.Xtb)+self.C)) 
            v5= 0.0451*self.Xtb*(self.P/((self.k['khS']*self.Xtb)+self.P))  
            v6 = 0.0568*self.Xtb*(self.L/((self.k['khS']*self.Xtb)+self.L)) 
            v7 = 0.009*self.Xta*(self.H/((self.k['khS']*self.Xta)+self.H))
            v8 = 0.0097*self.Xtf*(self.CE/((self.k['khS']*self.Xtf)+self.CE)) 
            v9 = 0.0140*self.Xtf*(self.LG/((self.k['khS']*self.Xtf)+self.LG)) 
            v10 = 0.009*self.Xma*(self.H/((self.k['khS']*self.Xma)+self.H)) 
            v11 = 0.0094*self.Xmf*(self.CE/((self.k['khS']*self.Xmf)+self.CE)) 
            v12 = 0.0140*self.Xmf*(self.LG/((self.k['khS']*self.Xmf)+self.LG))
            
            #II.2-Microorganisms growth
            #2-1-Mesophilic Bacterias growth
            v13 = self.k['µmb']*self.Xmb  *self.fSc *self.flT1 #*self.fSn  *self.faB_C
            v14 = self.k['µmb']*self.Xmb  * self.fSp  *self.flT1
            v15 = self.k['µmb']*self.Xmb  *self.fSl *self.flT1 # *self.fSn  *self.flT1 *self.faB_L
            #2-2-Thermophilic bacterias growth
            v16 = self.k['µtb']*self.Xtb  * self.fSc *self.flT2 # *self.fSn   *self.flT2 *self.faB_C 
            v17 = self.k['µtb']*self.Xtb *self.fSp *self.flT2  #*self.flT2
            v18 = self.k['µtb']*self.Xtb *self.fSl *self.flT2 #*self.fSn    *self.flT2  *self.faB_L
            #2-3-Mesophilic Actinomycetes growth
            v19 = self.k['µma']*self.Xma   * self.fSc*self.flT1   #*self.fSn     *self.flT1  *self.faA_C
            v20 = self.k['µma']*self.Xma    *self.fSp  *self.flT1  #*self.flT1
            v21 = self.k['µma']*self.Xma  *self.fSl *self.flT1 #*self.fSn     *self.flT1   *self.faA_L
            v22 = self.k['µma']*self.Xma   *self.fSh *self.flT1 #*self.fSn     *self.flT1   *self.faA_H
            #2-4-Thermophilic actinomycetes growth
            v23 = self.k['µta']*self.Xta  *self.fSc*self.flT2  #*self.fSn    *self.flT2  *self.faA_C
            v24 = self.k['µta']*self.Xta  *self.fSp *self.flT2 #*self.flT2
            v25 = self.k['µta']*self.Xta  *self.fSl *self.flT2 #*self.fSn      *self.flT2   *self.faA_L
            v26 = self.k['µta']*self.Xta  *self.fSh *self.flT2 #*self.fSn     *self.flT2    *self.faA_H
            #2-5-Mesophilic fungis growth
            v27 = self.k['µmf']*self.Xmf  * self.fSc *self.flT1  # *self.fSn    *self.flT1  *self.faF_C
            v28 = self.k['µmf']*self.Xmf  * self.fSp  *self.flT1 # *self.flT1
            v29 = self.k['µmf']*self.Xmf  * self.fSl*self.flT1 # *self.fSn  *self.flT1   *self.faF_L 
            v30 = self.k['µmf']*self.Xmf  *  self.fSh *self.flT1 #*self.fSn  *self.flT1   *self.faF_H 
            v31 = self.k['µmf']*self.Xmf  * self.fSlg *self.flT1 #*self.fSn  *self.flT1  *self.faF_LG
            #2-6- Thermophilic fungis growth
            v32 = self.k['µtf']*self.Xtf  * self.fSc *self.flT2 #*self.fSn      *self.flT2  *self.faF_C
            v33 = self.k['µtf']*self.Xtf   *self.fSp *self.flT2  # *self.flT2
            v34 = self.k['µtf']*self.Xtf   *self.fSl *self.flT2  #*self.fSn      *self.flT2  *self.faF_L 
            v35 = self.k['µtf']*self.Xtf    *self.fSh *self.flT2 #*self.fSn       *self.flT2  *self.faF_H
            v36 = self.k['µtf']*self.Xtf   *self.fSlg  *self.flT2 #*self.fSn       *self.flT2 *self.faF_LG
            
             #II.3-Heterotroph Microorganisms death
            v37 = self.Xmb*self.k['bmb']
            v38 = self.Xtb*self.k['btb']
            v39 = self.Xma*self.k['bma']
            v40 = self.Xta*self.k['bta']
            v41 = self.Xmf*self.k['bmf']
            v42 = self.Xtf*self.k['btf']
            
            #II.4-Lysis of microorganisms
            v43 = self.k['kdec']*self.Xdb
            
            #Global equations
            self.dC_dt= -v1 -v4
            self.dP_dt =   -v5 -v2 + (1-self.k['fi'])*(v43) #+v46)
            self.dL_dt = -v3-v6
            self.dH_dt= -v7-v10
            self.dCE_dt= -v8-v11
            self.dLG_dt= -v9-v12
            self.dXi_dt= self.k['fi']*(v43)
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
            
            self.dCH4gen_dt = (((self.k['YCH4,Sc'])*(v1+v4+v8+v11))+((self.k['YCH4,Sp'])*(v2+v5))+((self.k['YCH4,Sl'])*(v3+v6))+                               ((self.k['YCH4,Sh'])*(v7+v10))+((self.k['YCH4,Slg'])*(v9+v12)))*(1/(1+(self.k['η']*500)))

            self.dCH4oxi_dt = self.k['vmax']*((self.CH4gen/self.W)/(self.k['Km']+(self.CH4gen/self.W)))*(500/(self.k['Ko2,ch4']+500))
                                                         
            self.dCH4emis_dt = self.dCH4gen_dt - self.dCH4oxi_dt
            
            #self.dSn_dt = (-self.Y_NH3['r13']*v13)+(self.Y_NH3['r14']*v14)-(self.Y_NH3['r15']*v15)-(self.Y_NH3['r16']*v16)+(self.Y_NH3['r17']*v17)\
            #-(self.Y_NH3['r18']*v18)-(self.Y_NH3['r19']*v19)+(self.Y_NH3['r20']*v20)-(self.Y_NH3['r21']*v21)-(self.Y_NH3['r22']*v22)-(self.Y_NH3['r23']*v23)\
            #+(self.Y_NH3['r24']*v24)-(self.Y_NH3['r25']*v25)-(self.Y_NH3['r26']*v26)-(self.Y_NH3['r27']*v27)+(self.Y_NH3['r28']*v28)-(self.Y_NH3['r29']*v29)\
            #-(self.Y_NH3['r30']*v30)-(self.Y_NH3['r31']*v31)-(self.Y_NH3['r32']*v32)+(self.Y_NH3['r33']*v33)-(self.Y_NH3['r34']*v34)-(self.Y_NH3['r35']*v35)\
            #-(self.Y_NH3['r36']*v36)-self.Y_S['r44']*v44  #-v47
            
            #self.dNH3_dt = v47
            
            self.dW_dt = (self.Y_W['r13'])*v13+(self.Y_W['r14'])*v14+(self.Y_W['r15'])*v15+            (self.Y_W['r16'])*v16+(self.Y_W['r17'])*v17+(self.Y_W['r18'])*v18+(self.Y_W['r19'])*v19+            (self.Y_W['r20'])*v20+(self.Y_W['r21'])*v21+(self.Y_W['r22'])*v22+(self.Y_W['r23'])*v23+            (self.Y_W['r24'])*v24+(self.Y_W['r25'])*v25+(self.Y_W['r26'])*v26+(self.Y_W['r27'])*v27+            (self.Y_W['r28'])*v28+(self.Y_W['r29'])*v29+(self.Y_W['r30'])*v30+(self.Y_W['r31'])*v31+            (self.Y_W['r32'])*v32+(self.Y_W['r33'])*v33+(self.Y_W['r34'])*v34+(self.Y_W['r35'])*v35+            (self.Y_W['r36'])*v36 #+(self.Y_W['r44'])*v44
            
            #self.dO2_dt = -(-(self.Y_O2['r13'])*v13-(self.Y_O2['r14'])*v14-(self.Y_O2['r15'])*v15-\
            #(self.Y_O2['r16'])*v16-(self.Y_O2['r17'])*v17-(self.Y_O2['r18'])*v18-(self.Y_O2['r19'])*v19-\
            #(self.Y_O2['r20'])*v20-(self.Y_O2['r21'])*v21-(self.Y_O2['r22'])*v22-(self.Y_O2['r23'])*v23-\
            #(self.Y_O2['r24'])*v24-(self.Y_O2['r25'])*v25-(self.Y_O2['r26'])*v26-(self.Y_O2['r27'])*v27-\
            #(self.Y_O2['r28'])*v28-(self.Y_O2['r29'])*v29-(self.Y_O2['r30'])*v30-(self.Y_O2['r31'])*v31-\
            #(self.Y_O2['r32'])*v32-(self.Y_O2['r33'])*v33-(self.Y_O2['r34'])*v34-(self.Y_O2['r35'])*v35-\
            #(self.Y_O2['r36'])*v36 -(self.Y_O2['r44'])*v44) #+O2in
            
            #température
            #biological temperature
            Qbio = self.k['hbio']*(((v13+v14+v15+v16+v17+v18+v19+v20+v21+v22+v23+v24+v25+v26)/0.113)+                             ((v27+v28+v29+v30+v31+v32+v33+v34+v35+v36)/0.247))*self.va['TM'] / self.k['Yo2']
            

            
            #Qbio = (self.k['hbio']*(O2cons) * self.va['TM']) / 0.032
            
            #Heat loss by conduction
            Qcond =0 # (self.T - 325)*self.A()*2 #*self.U()
            #Convctive heat loss
            Qconv = 1*( self.T- 293)*self.Qair() #supposons que pour calculer l'enthalpie de l'air sortant et entrant, on prend Ca=1kJ/K.kg
            
            
            self.dT_dt = (Qbio-Qcond-Qconv)/(self.va['TM']*2) #self.va['Cp'])
            

            
            return [self.dC_dt,self.dP_dt,self.dL_dt,self.dH_dt,self.dCE_dt,self.dLG_dt,self.dXi_dt,self.dSc_dt,                    self.dSp_dt,self.dSl_dt,self.dSh_dt,self.dSlg_dt,self.dXmb_dt, self.dXtb_dt,self.dXma_dt,self.dXta_dt,                   self.dXmf_dt,self.dXtf_dt,self.dXdb_dt,self.dCO2_dt,                   self.dCH4gen_dt,self.dCH4oxi_dt,self.dCH4emis_dt,self.dW_dt, self.dT_dt]#,self.dNH3_dt] #,self.dO2_dt
        
       
        dt = 1
            
        self.solution = solve_ivp(system, [0, self.duration], [self.va['G'],self.va['P'],self.va['L'],                                                self.va['HE'],self.va['CE'],self.va['LG'], self.va['Xi'],                                                    0,0,0,0,0,                                                     self.MB,self.TB,self.MA,self.TA,self.MF,                                                     self.TF, self.va['Xdb'], 0,                                                     0,0,0,self.va['W'],294], t_eval = np.arange(0,self.duration,dt),                                      method='RK45') #
            
        for i in range(len(self.solution.y)):
            self.solution.y[i][self.solution.y[i] < 0] = 0
            
        return (self.solution.y[0])    
        
       
        
    
    def emissions(self):
        variablesE = ['CO2','NO3', 'N2', 'N2O','CH4' ,'NH3']
        functionsE= [self.solution.y[20], self.solution.y[21], self.solution.y[22], self.solution.y[23], self.solution.y[26],self.solution.y[29]]
        num_variablesE = len(variablesE)
        
        for i in range(num_variablesE): 
            print ('total emission of', variablesE[i], ':' , functionsE[i][-1])
    
    def biomass(self):
        variables = ['Xmb','Xtb', 'Xma', 'Xta','Xmf','Xtf']
        functions = [self.solution.y[12], self.solution.y[13], self.solution.y[14], self.solution.y[15], self.solution.y[16],                     self.solution.y[17]]
        num_variablesC = len(variables)
        
        #for i in range(num_variables): 
            #print ('total emission of', variables[i], ':' , trapz(functions[i], self.solution.t))
        
        for i in range(num_variablesC):
            plt.subplot(3,2, i + 1)
            plt.plot(self.solution.t, functions[i], label=variables[i])
            plt.xlabel('Temps (h)')
            plt.ylabel(f'{variables[i]} (kg/kg)')
            plt.legend()
            plt.grid(True)
            
            #plt.title(f'{variables[i]}')
            #plt.figure(figsize=(12, 4))
        
        
        plt.tight_layout()
        
        plt.show()
    
    def insolubles(self):
        substrates = ['C','P', 'L', 'H', 'CE','LG']
        functions_substrates = [self.solution.y[0], self.solution.y[1], self.solution.y[2], self.solution.y[3], self.solution.y[4],self.solution.y[5]]
                     
        num_substrates = len(substrates)
        
        
        for i in range(num_substrates):
            plt.subplot(3, 2, i + 1)
            plt.plot(self.solution.t, functions_substrates[i], label=substrates[i])
            plt.xlabel('Temps (h)')
            plt.ylabel(f'{substrates[i]} (kg/kg)')
            plt.legend()
            plt.grid(True)
    
    def solubles(self):
        sol = ['Sc','Sp', 'Sl', 'Sh', 'Slg']
        functions_sol = [self.solution.y[7], self.solution.y[8], self.solution.y[9], self.solution.y[10],self.solution.y[11]]
        
                     
        num_sol = len(sol)
        
        #for i in range(num_variables): 
            #print ('total emission of', variables[i], ':' , trapz(functions[i], self.solution.t))
        
        for i in range(num_sol):
            plt.subplot(3, 2, i + 1)
            plt.plot(self.solution.t, functions_sol[i], label=sol[i])
            plt.xlabel('Temps (h)')
            plt.ylabel(f'{sol[i]} (kg/kg)')
            plt.legend()
            plt.grid(True)
    
    def CO2emission(self):
        plt.plot(self.solution.y[19], label = 'CO2 model')
        plt.plot(np.array(self.gas[0]), label = 'CO2 experimental')
        plt.xlabel('Temps (h)')
        plt.ylabel('CO2 (kg/kg)')
        plt.legend()
        plt.grid(True)


# In[84]:


Data1 = 'Data Komilis FWGW.xlsx'
Data2 = 'cO2 Komilis.xlsx'
Data3 = 'Data Komilis GW.xlsx'
Data4 = 'Data for fitting en h.xlsx'
Data5 = 'Data random.xlsx'
Data6 = 'Data random 2.xlsx'


# In[58]:


essai0 = Compostage(Data1, Data2, 'test_fitting', 'Passive aeration',1000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[62]:


essai0.resolution()


# In[63]:


essai0.solution


# In[20]:


essai0.insolubles()


# In[28]:


essai0.solution.y[4]


# In[21]:


essai0.solubles()


# In[22]:


essai0.biomass()


# In[23]:


essai0.CO2emission()


# In[65]:


essai1 = Compostage(Data1, Data2, 'test_fitting', 'Passive aeration',2000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[66]:


essai1.resolution()


# In[67]:


essai1.CO2emission()


# In[27]:


essai1.biomass()


# In[50]:


plt.plot(essai1.solution.y[24])


# In[29]:


essai1.insolubles()


# In[30]:


essai1.solubles()


# In[34]:


essaiGW = Compostage(Data3, Data2, 'test_fitting', 'Passive aeration',2000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[35]:


essaiGW.resolution()


# In[59]:


Margaritis = Compostage(Data4, Data2, 'test_fitting', 'Passive aeration',2000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[60]:


Margaritis.resolution()


# In[61]:


Margaritis.CO2emission()


# In[76]:


plt.plot(Margaritis.solution.y[19])
plt.plot(essai1.solution.y[19])
plt.plot(random.solution.y[19])


# In[45]:


Margaritis.insolubles()


# In[46]:


Margaritis.solubles()


# In[47]:


Margaritis.biomass()


# In[49]:


plt.plot(Margaritis.solution.y[24])


# In[55]:


plt.plot(Margaritis.solution.y[19])


# In[70]:


random = Compostage(Data5, Data2, 'test_fitting', 'Passive aeration',2000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[71]:


random.resolution()


# In[72]:


random.solution


# In[75]:


random.CO2emission()


# In[77]:


random.biomass()


# In[78]:


random.insolubles()


# In[80]:


essai1.insolubles()


# In[81]:


Margaritis.solution.y[19][-1]


# In[82]:


essai1.solution.y[19][-1]


# In[83]:


random.solution.y[19][-1]


# In[86]:


random2 = Compostage(Data6, Data2, 'test_fitting', 'Passive aeration',2000, 0.001, 0.001, 0.001,0.001,0.001,0.001)


# In[87]:


random2.resolution()


# In[89]:


random2.CO2emission()


# In[90]:


random2.solution.y[19][-1]


# In[95]:


plt.plot(random2.solution.y[24])


# In[91]:


random2_biomass = Compostage(Data6, Data2, 'test_fitting', 'Passive aeration',2000, 0.0001, 0.0001, 0.0001,0.0001,0.0001,0.0001)


# In[92]:


random2_biomass.resolution()


# In[93]:


random2_biomass.CO2emission()


# In[96]:


plt.plot(random2_biomass.solution.y[24])


# In[97]:


random2_biomass2 = Compostage(Data6, Data2, 'test_fitting', 'Passive aeration',2000, 0.01, 0.01, 0.01,0.01,0.01,0.01)


# In[98]:


random2_biomass2.resolution()


# In[100]:


random2_biomass2.CO2emission()


# In[101]:


random2_biomass2.biomass()


# In[102]:


random2_biomass2.insolubles()


# In[103]:


random2_biomass2.solubles()


# In[ ]:




