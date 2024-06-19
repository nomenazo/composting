function dydt=compostfitall(t,y); %fitting with more variables than compostfit (CH4, O2)
global thetag;
%1er essai de fitting par la constante de décès
%bd=thetag
kinparam = thetag;


%insolubles substrates
C=y(1);	% carbohydrate kg/kgTM
P=y(2);
L=y(3);
H=y(4);
CE=y(5);
LG=y(6);
Xi=y(7);
%solubles substrates
Sc=y(8);	% carbohydrate soluble kg/kgTM
Sp=y(9);
Sl=y(10);
Sh=y(11);
Slg=y(12);
%biomass
Xmb=y(13); % biomass kg/kgTM
Xtb=y(14);
Xma=y(15);
Xta=y(16);
Xmf=y(17);
Xtf=y(18);
Xdb=y(19);
CO2=y(20); %CO2 kg/kgTM
W=y(21);   %H2O kg/kgTM
T=y(22);
%O2cons=y(23); %oxygen consommé kg/kgTM

%O2trans=y(25); %oxygène transféré dans la phase liquide kg/kgTM
O2diss=y(23); %dissolved oxygen kg/kgTM
%O2trans=y(24); 
%O2in=y(25);%oxygen entrant dans la phase gazeuse par le débit d'air kg/kgTM
CH4gen=y(24);
CH4oxi=y(25);
CH4=y(26);

TM= 0.6; %kg total matter
Ta = 285; %External temperature

%kh = [0.0001 0.0378 0.2991 2.8873e-05 0.0271 0.0153 0.009 0.0025 0.0078 0.009 0.0025 0.0313] %hydrolysis constant fitted when initial values of microorganisms are 0.001 (h-1)
%kh = [0.0336 0.0161 0.0076 0.0173 0.0333 0.0082 0.0090 0.0089 0.0086 0.0090 0.0089 0.0086] %hydrolysis constant fitted when initial values of microorganisms are 0.01 (h-1)
%kinparam = [0.04 0.02 0.01 0.02 0.04 0.01 0.009 0.007 0.007 0.009 0.007 0.007 0.2 0.18 0.1 0.12 0.1 0.1 0.03 0.02 0.01 0.015 0.01 0.01]; %kh,mu,bd not fitted
%mu = [0.2 0.18 0.1 0.12 0.1 0.1]; %specific growth rate (h-1)
%bd = [0.03 0.02 0.01 0.015 0.01 0.01]; %death rate (h-1)
K=[6.2e-5 1e-4 0.2 0.0025 7e-9 ]; %kinetic parameters
KT = [440 0.072 1 2 0.09 1]; %parameters for temperature module
Yx_s= 0.35; %biomass yield on substrate kgX/kgS
Yx_co2 = [0.445717506 0.293234476	0.165284084	0.445717506	0.293234476	0.165284084	0.445717506	0.293234476	0.165284084	0.139609659	0.445717506	0.293234476	0.165284084	0.139609659	0.41509283	0.254264706	0.160882526	0.136456284	0.19653097	0.41509283	0.254264706	0.160882526	0.136456284	0.19653097
]; %Yield coeff of biomass on CO2 kgX/kgCO2
Yx_h = [0.716485507	1.791988467	0.380836347	0.716485507	1.791988467	0.380836347	0.716485507	1.791988467	0.380836347	0.329440099	0.716485507	1.791988467	0.380836347	0.329440099	0.830451489	2.728498673	0.410802056	0.351627849	0.625936213	0.830451489	2.728498673	0.410802056	0.351627849	0.625936213
]; %Yield coeff of biomass on H2O kgX/kgH2O
Yx_O2mol = [0.173553719	0.098678414	0.044403196	0.173553719	0.098678414	0.044403196	0.173553719	0.098678414	0.044403196	0.054361283	0.173553719	0.098678414	0.044403196	0.054361283	0.07678245	0.044286279	0.019267472	0.024607076	0.027217967	0.07678245	0.044286279	0.019267472	0.024607076	0.027217967];
%yield coeff biomass on O2 molX/molO2
Yx_o2 = [0.61286157	0.34845815	0.156798785	0.61286157	0.34845815	0.156798785	0.61286157	0.34845815	0.156798785	0.191963281	0.61286157	0.34845815	0.156798785	0.191963281	0.592664534	0.341834717	0.148720799	0.189935871	0.210088681	0.592664534	0.341834717	0.148720799	0.189935871	0.210088681	0.040522541
]; %Yield coeff biomass on O2 kgX/kgO2



Ks = K(1); %substrate saturation for monod kinetics (kgSc/dm3W)
Khs = K(2); %coefficient de saturation pour les cinétiques de contois
fi = K(3); %Proportion of dead biomass recycled to inert materials
kdec = K(4); %microorganisms decomposition constant
kO2=K(5); %oxygen saturation for heterotrophic activities (kgO2/l)

%Methane module parameters
Ych4_Sc = 0.267; %kgCH4/kgSC
Ych4_Sp = 0.375 ; %kgCH4/kgSp
Ych4_Sl = 0.707; %kgCH4/kgSl
Ych4_Sh = 0.284 ; 
Ych4_Slg = 0.535 ;
eta = 4e5 ; %L/mol %sensitivity of methangogenesis to inhibition by oxygen
Vmax = 5.35e-4; %kgCH4/kgTM.h vitesse maximale d'oxydation de methane
km = 0.72 ; %kg/l Michaelis constant for methane oxidation
Kch4_O2 = 0.033; %mol/l Michelis constant for oxygen in methane oxidation



hbio = KT(1); %chaleur dégagée par mol d'oxygène consommée (kJ/mol d'O2)
Qair = KT(2); %débit d'air par aération passive (kg/h) 0.072 value from Rasapoo, equivalent to 0.6 l/m3.kg 
Ca = KT(3); %Capacité calorifique de l'air sec (kJ/K.kg)
Cw = KT(4); %Capacité calorique des biodéchets (kJ/K.kg)
U = KT(5); %heat transfer coefficient of wall (kJ/m2.K.h) %%valeur dans de Guardia 2012 : 7W/m2.C = 0.09 kJ/m2.h.K
A = KT (6); %surface area of heat conduction (m2)

%hydrolysis constant
kh1C = kinparam(1);	
kh2P = kinparam(2);
kh3L = kinparam(3);	
kh4C = kinparam(4);	
kh5P = kinparam(5);	
kh6L = kinparam(6);	
kh7H = kinparam(7);	
kh8CE = kinparam(8);	
kh9LG = kinparam(9);	
kh10H = kinparam(10);	
kh11CE = kinparam(11);	
kh12LG = kinparam(12);

%growth rate
mmb = kinparam(13);
mtb = kinparam(14);
mma = kinparam(15);
mta = kinparam(16);
mmf = kinparam(17);
mtf = kinparam(18);

%death rate
bmb = kinparam(19);
btb = kinparam(20);
bma = kinparam(21)*1;
bta = kinparam(22)*1;
bmf = kinparam(23)*1;
btf = kinparam(24)*1;

%Inverse of yield of biomass on CO2
Ymb_c_c = 1/Yx_co2(1) ; %MB on Sc
Ymb_p_c = 1/Yx_co2(2) ; %MB on Sp
Ymb_l_c = 1/Yx_co2(3) ; %MB on Sl
Ytb_c_c = 1/Yx_co2(4) ; %TB on Sc
Ytb_p_c = 1/Yx_co2(5) ; %TB on Sp
Ytb_l_c = 1/Yx_co2(6) ; %TB on Sl
Yma_c_c = 1/Yx_co2(7) ; %MA on Sc
Yma_p_c = 1/Yx_co2(8) ; %MA on Sp
Yma_l_c = 1/Yx_co2(9) ; %MA on Sl
Yma_h_c = 1/Yx_co2(10) ; %MA on Sh
Yta_c_c = 1/Yx_co2(11) ; %TA on Sc
Yta_p_c = 1/Yx_co2(12) ; %TA on Sp
Yta_l_c = 1/Yx_co2(13) ; %TA on Sl
Yta_h_c = 1/Yx_co2(14) ; %TA on Sh
Ymf_c_c = 1/Yx_co2(15) ; %MF on Sc
Ymf_p_c = 1/Yx_co2(16) ; %MF on Sp
Ymf_l_c = 1/Yx_co2(17) ; %MF on Sl
Ymf_h_c = 1/Yx_co2(18) ; %MF on Sh
Ymf_lg_c = 1/Yx_co2(19) ; %MF on Slg
Ytf_c_c = 1/Yx_co2(20) ; %TF on Sc
Ytf_p_c = 1/Yx_co2(21) ; %TF onSp
Ytf_l_c = 1/Yx_co2(22) ; %TF on Sl
Ytf_h_c = 1/Yx_co2(23) ; %TF on Sh
Ytf_lg_c = 1/Yx_co2(24) ; %TF on Slg

%Inverse of yield of biomass on H2O
Ymb_c_h = 1/Yx_h(1) ; %MB on Sc
Ymb_p_h = 1/Yx_h(2) ; %MB on Sp
Ymb_l_h = 1/Yx_h(3) ; %MB on Sl
Ytb_c_h = 1/Yx_h(4) ; %TB on Sc
Ytb_p_h = 1/Yx_h(5) ; %TB on Sp
Ytb_l_h = 1/Yx_h(6) ; %TB on Sl
Yma_c_h = 1/Yx_h(7) ; %MA on Sc
Yma_p_h = 1/Yx_h(8) ; %MA on Sp
Yma_l_h = 1/Yx_h(9) ; %MA on Sl
Yma_h_h = 1/Yx_h(10) ; %MA on Sh
Yta_c_h = 1/Yx_h(11) ; %TA on Sc
Yta_p_h = 1/Yx_h(12) ; %TA on Sp
Yta_l_h = 1/Yx_h(13) ; %TA on Sl
Yta_h_h = 1/Yx_h(14) ; %TA on Sh
Ymf_c_h = 1/Yx_h(15) ; %MF on Sc
Ymf_p_h = 1/Yx_h(16) ; %MF on Sp
Ymf_l_h = 1/Yx_h(17) ; %MF on Sl
Ymf_h_h = 1/Yx_h(18) ; %MF on Sh
Ymf_lg_h = 1/Yx_h(19) ; %MF on Slg
Ytf_c_h = 1/Yx_h(20) ; %TF on Sc
Ytf_p_h = 1/Yx_h(21) ; %TF onSp
Ytf_l_h = 1/Yx_h(22) ; %TF on Sl
Ytf_h_h = 1/Yx_h(23) ; %TF on Sh
Ytf_lg_h = 1/Yx_h(24) ; %TF on Slg

%Inverse of yield of biomass on O2
Ymb_c_o2 = 1/Yx_o2(1) ; %MB on Sc
Ymb_p_o2 = 1/Yx_o2(2) ; %MB on Sp
Ymb_l_o2 = 1/Yx_o2(3) ; %MB on Sl
Ytb_c_o2 = 1/Yx_o2(4) ; %TB on Sc
Ytb_p_o2 = 1/Yx_o2(5) ; %TB on Sp
Ytb_l_o2 = 1/Yx_o2(6) ; %TB on Sl
Yma_c_o2 = 1/Yx_o2(7) ; %MA on Sc
Yma_p_o2 = 1/Yx_o2(8) ; %MA on Sp
Yma_l_o2 = 1/Yx_o2(9) ; %MA on Sl
Yma_h_o2 = 1/Yx_o2(10) ; %MA on Sh
Yta_c_o2 = 1/Yx_o2(11) ; %TA on Sc
Yta_p_o2 = 1/Yx_o2(12) ; %TA on Sp
Yta_l_o2 = 1/Yx_o2(13) ; %TA on Sl
Yta_h_o2 = 1/Yx_o2(14) ; %TA on Sh
Ymf_c_o2 = 1/Yx_o2(15) ; %MF on Sc
Ymf_p_o2 = 1/Yx_o2(16) ; %MF on Sp
Ymf_l_o2 = 1/Yx_o2(17) ; %MF on Sl
Ymf_h_o2 = 1/Yx_o2(18) ; %MF on Sh
Ymf_lg_o2 = 1/Yx_o2(19) ; %MF on Slg
Ytf_c_o2 = 1/Yx_o2(20) ; %TF on Sc
Ytf_p_o2 = 1/Yx_o2(21) ; %TF onSp
Ytf_l_o2 = 1/Yx_o2(22) ; %TF on Sl
Ytf_h_o2 = 1/Yx_o2(23) ; %TF on Sh
Ytf_lg_o2 = 1/Yx_o2(24) ; %TF on Slg



%limitation by the substrate
fSc = (Sc/W)/(Ks+(Sc/W));
fSp=(Sp/W)/(Ks+(Sp/W));
fSl=(Sl/W)/(Ks+(Sl/W));
fSh=(Sh/W)/(Ks+(Sh/W));
fSlg=(Slg/W)/(Ks+(Slg/W));

%limitation function by temperature

Ti = T ;
Tmax1 = 355; %319
Tmin1 =  0; %280.1;
Topt1 =  322; %310.4;
%Tmax2 = 340.5; %355;
%Tmin2 = 305.8; %0;
%Topt2 = 332.2; %322;
fT1 = ((Ti-Tmax1)*(Ti-Tmin1)^2)/((Topt1-Tmin1)*((Topt1-Tmin1) *(Ti-Topt1)-(Topt1-Tmax1)*(Topt1+Tmin1-2*Ti))) ;
%fT2 = ((Ti-Tmax2)*(Ti-Tmin2)^2)/((Topt2-Tmin2)*((Topt2-Tmin2) *(Ti-Topt2)-(Topt2-Tmax2)*(Topt2+Tmin2-2*Ti))) ; %Oudart parameters

%fT = 0.01;

%Limitation function by the oxygen
fO2 = (O2diss/W)/(kO2+(O2diss/W));

%hydrolisis rate
v1 = kh1C*Xmb*(C/((Khs*Xmb)+C));
v2 =kh2P* Xmb*(P/((Khs*Xmb)+P)) ;
v3 = kh3L *Xmb*(L/((Khs*Xmb)+L));
v4 = kh4C*Xtb*(C/((Khs*Xtb)+C));
v5= kh5P*Xtb*(P/((Khs*Xtb)+P));
v6 = kh6L*Xtb*(L/((Khs*Xtb)+L));
v7 = kh7H*Xta*(H/((Khs*Xta)+H));
v8 = kh8CE*Xtf*(CE/((Khs*Xtf)+CE));
v9 = kh9LG*Xtf*(LG/((Khs*Xtf)+LG));
v10 = kh10H*Xma*(H/((Khs*Xma)+H));
v11 = kh11CE*Xmf*(CE/((Khs*Xmf)+CE));
v12 = kh12LG*Xmf*(LG/((Khs*Xmf)+LG));

%Biomass growth reaction
% Mesophilic Bacterias growth
v13 = mmb*Xmb  *fSc *fT1 *fO2;
v14 = mmb*Xmb  * fSp *fT1*fO2 ;
v15 = mmb*Xmb  *fSl *fT1*fO2;
%Thermophilic bacterias growth
v16 = mtb*Xtb  * fSc *fT1*fO2; %*fT2;
v17 = mtb*Xtb *fSp *fT1*fO2;%*fT2;
v18 = mtb*Xtb *fSl *fT1*fO2;%*fT2;
%Mesophilic Actinomycetes growth
v19 = mma*Xma   * fSc *fT1*fO2;
v20 = mma*Xma    *fSp *fT1*fO2;
v21 = mma*Xma  *fSl *fT1*fO2;
v22 = mma*Xma   *fSh *fT1*fO2;
 %Thermophilic actinomycetes growth
v23 = mta*Xta  *fSc *fT1*fO2;%*fT2;
v24 = mta*Xta  *fSp*fT1*fO2;% *fT2;
v25 = mta*Xta  *fSl *fT1*fO2;%*fT2;
v26 = mta*Xta  *fSh*fT1*fO2; %*fT2;
%Mesophilic fungis growth
v27 = mmf*Xmf  *fSc *fT1*fO2;
v28 = mmf*Xmf  * fSp *fT1*fO2;
v29 = mmf*Xmf  * fSl *fT1*fO2;
v30 = mmf*Xmf  *  fSh *fT1*fO2;
v31 = mmf*Xmf  * fSlg *fT1*fO2; 
%Thermophilic fungis growth
v32 = mtf*Xtf  *fSc *fT1*fO2;%*fT2;
v33 = mtf*Xtf *fSp *fT1*fO2;%*fT2;
v34 = mtf*Xtf   *fSl *fT1*fO2;%*fT2; 
v35 = mtf*Xtf    *fSh*fT1*fO2;%*fT2;
v36 = mtf*Xtf   *fSlg *fT1*fO2;%*fT2;

%Heterotroph Microorganisms death
v37 = bmb*Xmb ;
v38 = btb*Xtb ;
v39 = bma*Xma ;
v40 = bta*Xta;
v41 = bmf*Xmf;
v42 = btf*Xtf;

%lysis of microorganisms
v43 = kdec*Xdb;

%Global equations
dCdt= -v1 -v4;
dPdt =   -v5 -v2 + (1-fi)*(v43);
dLdt = -v3-v6;
dHdt= -v7-v10;
dCEdt = -v8-v11;
dLGdt= -v9-v12;
dXidt= fi*(v43);
dScdt =  v1 +v4 +v11 +v8 - (1/Yx_s)*(v13+v16+v19+v23+v27+v32);    
dSpdt =  v2  +v5 -(1/Yx_s)*(v14+v17+v20+v24+v28+v33);
dSldt = v3+v6-(1/Yx_s)*(v15+v18+v21+v25+v29+v34);
dShdt = v7+v10-(1/Yx_s)*(v22+v26+v30+v35);
dSlgdt = v9+v12-(1/Yx_s)*(v31+v36);
dXmbdt = v13+ v14 + v15 -v37;
dXtbdt = v16 + v17 + v18 - v38;
dXmadt = v19+v20 + v21 +v22- v39;
dXtadt = v23 + v24 + v25 +v26- v40;
dXmfdt = v27 + v28 + v29 +v30+ v31 - v41;
dXtfdt = v32 + v33 + v34 +v35+v36 - v42;           
dXdbdt = v37+v38+v39+v40+v41+v42-v43;
            
dCO2dt = (Ymb_c_c)*v13+(Ymb_p_c)*v14+(Ymb_l_c)*v15+(Ytb_c_c)*v16+(Ytb_p_c)*v17+(Ytb_l_c)*v18+(Yma_c_c)*v19+(Yma_p_c)*v20+(Yma_l_c)*v21+...
(Yma_h_c)*v22+(Yta_c_c)*v23+(Yta_p_c)*v24+(Yta_l_c)*v25+(Yta_h_c)*v26+(Ymf_c_c)*v27+(Ymf_p_c)*v28+(Ymf_l_c)*v29+...
(Ymf_h_c)*v30+(Ymf_lg_c)*v31+(Ytf_c_c)*v32+(Ytf_p_c)*v33+(Ytf_l_c)*v34+(Ytf_h_c)*v35+(Ytf_lg_c)*v36;
            
dWdt = (Ymb_c_h)*v13+(Ymb_p_h)*v14+(Ymb_l_h)*v15+(Ytb_c_h)*v16+(Ytb_p_h)*v17+(Ytb_l_h)*v18+(Yma_c_h)*v19+(Yma_p_h)*v20+(Yma_l_h)*v21+...
    (Yma_h_h)*v22+(Yta_c_h)*v23+(Yta_p_h)*v24+(Yta_l_h)*v25+(Yta_h_h)*v26+(Ymf_c_h)*v27+(Ymf_p_h)*v28+(Ymf_l_h)*v29+...
    (Ymf_h_h)*v30+(Ymf_lg_h)*v31+(Ytf_c_h)*v32+(Ytf_p_h)*v33+(Ytf_l_h)*v34+(Ytf_h_h)*v35+(Ytf_lg_h)*v36; 

%O2cons = ((Ymb_c_o2)*v13+(Ymb_p_o2)*v14+(Ymb_l_o2)*v15+(Ytb_c_o2)*v16+(Ytb_p_o2)*v17+(Ytb_l_o2)*v18+(Yma_c_o2)*v19+(Yma_p_o2)*v20+(Yma_l_o2)*v21+...
    %(Yma_h_o2)*v22+(Yta_c_o2)*v23+(Yta_p_o2)*v24+(Yta_l_o2)*v25+(Yta_h_o2)*v26+(Ymf_c_o2)*v27+(Ymf_p_o2)*v28+(Ymf_l_o2)*v29+...
   % (Ymf_h_o2)*v30+(Ymf_lg_o2)*v31+(Ytf_c_o2)*v32+(Ytf_p_o2)*v33+(Ytf_l_o2)*v34+(Ytf_h_o2)*v35+(Ytf_lg_o2)*v36); %kg/kgTM.h

%dO2indt = ((Qair * 16e-3*210)/(1.293*22.4)); %-O2trans; %(Qair*MO2*pO2/rhoair*V_GP) %V_GP=volume molaire des GP %kg/kgTM.h

xO2 = ((O2diss/16e-3 )/ W)/55.56; %fraction molaire
Vreactor = 0.002 ; %m3
Vwaste = 0.6/350; %m3%mass/rhobiowaste
Vgas = Vreactor - Vwaste; % m3
rhoair = 1.2; %kg/m3
khO2 = 1.4e-3; %mol/m3.Pa Henry constant for oxygen, mais devra être en fonction de la T

R = 8.134; %Pa.m3/mol.K
He =  4.39751e9; %101325 *exp(66.7354-(8747.55/T)-24.4526*log(T/100)); %Pa %
%klO2 = 1.584e-6*W; %1e-7 ; %kg.h-1.Pa-1.kgTM-1 ; Oxygen liquid–gas mass transfer constant

%O2trans = 1e-7 * ( He*xO2 - ((O2in/16e-3)*R*T*TM/Vgas)); %O2 transferred to liquid phase kg/kgTM %kg/kgTM.h
%O2trans = 1e-7 * ( He*xO2 - ((((Qair *16e-3*210)/(1.293*22.4))/16e-3)*R*T*TM/Vgas)); 
%1e-7 * ( He*xO2 - ((((Qair * 16e-3*210)/(1.293*22.4))/16e-3)*R*T*TM/Vgas))-
O2gas = 3.33e-5;%kgO2/kgTM en considérant que la concentration de l'oxygène dans la phase gazeuse est maintenue à 15% tout au long du compostage
%dO2dissdt = klO2 * ( He*xO2 - ((O2gas/16e-3)*R*T*TM/Vgas));
 
%Théorie des deux films
dO2liqdt = khO2*(0.79/22.4)*(Qair/rhoair)*R*T/Vgas; %mol/m3
%dO2dissdt=(khO2*(0.79/22.4)*(Qair/rhoair)*R*T/Vgas)*W*1e-3*32e-3
dO2dissdt = (khO2*(0.79/22.4)*(Qair/rhoair)*R*T/Vgas)*W*1e-3*32e-3-((Ymb_c_o2)*v13+(Ymb_p_o2)*v14+(Ymb_l_o2)*v15+(Ytb_c_o2)*v16+(Ytb_p_o2)*v17+(Ytb_l_o2)*v18+(Yma_c_o2)*v19+(Yma_p_o2)*v20+(Yma_l_o2)*v21+...
(Yma_h_o2)*v22+(Yta_c_o2)*v23+(Yta_p_o2)*v24+(Yta_l_o2)*v25+(Yta_h_o2)*v26+(Ymf_c_o2)*v27+(Ymf_p_o2)*v28+(Ymf_l_o2)*v29+...
  (Ymf_h_o2)*v30+(Ymf_lg_o2)*v31+(Ytf_c_o2)*v32+(Ytf_p_o2)*v33+(Ytf_l_o2)*v34+(Ytf_h_o2)*v35+(Ytf_lg_o2)*v36);



%dO2dt= O2in - O2cons;

%temperature module
Qbio = hbio *(((v13/Yx_O2mol(1))+(v14/Yx_O2mol(2))+(v15/Yx_O2mol(3))+(v16/Yx_O2mol(4))+(v17/Yx_O2mol(5))+(v18/Yx_O2mol(6))+(v19/Yx_O2mol(7))+(v20/Yx_O2mol(8))+...
    (v21/Yx_O2mol(9))+(v22/Yx_O2mol(10))+(v23/Yx_O2mol(11))+(v24/Yx_O2mol(12))+(v25/Yx_O2mol(13))+(v26/Yx_O2mol(14)))/0.113+((v27/Yx_O2mol(15))+(v28/Yx_O2mol(16))+(v29/Yx_O2mol(17))+...
    (v30/Yx_O2mol(18))+(v31/Yx_O2mol(19))+(v32/Yx_O2mol(20))+(v33/Yx_O2mol(21))+(v34/Yx_O2mol(22))+(v35/Yx_O2mol(23))+(v36/Yx_O2mol(24)))/0.247)* TM ;

Qconv = Ca*( T- Ta)*Qair ;

Qcond = (T - Ta)*U*A;

dTdt= (Qbio-Qconv-Qcond)/(TM*Cw) ;

%Methane module
dCH4gendt= (Ych4_Sc*(v1+v4+v8+v11 )+Ych4_Sp*(v2+v5 )+Ych4_Sl*(v3+v6 )+Ych4_Sh* (v7+v10 )+Ych4_Slg*(v9+v12 ))*(1/(1+eta*(O2diss/(W*32e-3)))); %kg/kgTM.h


dCH4oxidt = Vmax*((CH4gen/W)/(km+(CH4gen/W)))*((O2diss/(W*32e-3))/(Kch4_O2+(O2diss/(W*32e-3))));


dCH4dt= dCH4gendt - dCH4oxidt ;

%C=max(0,C);
%Sc=max(0,Sc);
%Xmb=max(0,Xmb);
%LG = max(0,LG);
%CE = max(0,CE);
O2diss = max(0,O2diss);
%CH4gen = max(0,CH4gen);
%CH4oxi=max(0,CH4oxi);
%CH4=max(0,CH4);

%y(y<0)=0;

dydt = [dCdt,dPdt,dLdt,dHdt,dCEdt,dLGdt,dXidt, dScdt,dSpdt,dSldt,dShdt,dSlgdt, dXmbdt,dXtbdt,dXmadt,dXtadt,dXmfdt,dXtfdt, dXdbdt, dCO2dt, dWdt,dTdt,dO2dissdt,dCH4gendt, dCH4oxidt,dCH4dt]'; %,dO2transdt,dO2indt]' %]';