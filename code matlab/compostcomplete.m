function dydt=compostcomplete(t,y);

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
%CH4gen=y(22);
%CH4oxi=y(23);
%CH4emis=y(24);
Xa=y(22);
NO3=y(23);
N2O=y(24);
N2=y(25);
%Sn=y(29);
%T=y(30);

%kh = [0.04 0.02 0.01 0.02 0.04 0.01 0.009 0.007 0.007 0.009 0.007 0.007];%hydrolysis constant%literature values
kh = [0.0716 0.0676 0.0901 0.0001 0.0451 0.0568 0.0090 0.0097 0.0140 0.0090 0.0094 0.0140] %2nd fitting
mu = [0.2 0.18 0.1 0.12 0.1 0.1]; %specific growth rate (h-1)
bd = [0.03 0.02 0.01 0.015 0.01 0.01]; %death rate (h-1) (values from
%literature)
%bd = [0.0170 0.0131 1.0776 1.8710 0.0075 0.0073] %values from 1st fitting
K=[6.2e-5 1e-4 0.2 0.0025]; %kinetic parameters
Yx_s= 0.35; %biomass yield on substrate
Yx_co2 = [0.445717506 0.293234476	0.165284084	0.445717506	0.293234476	0.165284084	0.445717506	0.293234476	0.165284084	0.139609659	0.445717506	0.293234476	0.165284084	0.139609659	0.41509283	0.254264706	0.160882526	0.136456284	0.19653097	0.41509283	0.254264706	0.160882526	0.136456284	0.19653097
]; %Yield coeff of biomass on CO2
Yx_h = [0.716485507	1.791988467	0.380836347	0.716485507	1.791988467	0.380836347	0.716485507	1.791988467	0.380836347	0.329440099	0.716485507	1.791988467	0.380836347	0.329440099	0.830451489	2.728498673	0.410802056	0.351627849	0.625936213	0.830451489	2.728498673	0.410802056	0.351627849	0.625936213
]; %Yield coeff of biomass on H2O

kXa = [0.03 0.0083 0.042 0.2] %literature values of kinetic parameters for nitrfication/denitrification
Y_NO3 = 1/0.0390553; %inverse of autotrop biomass yield on NO3


Ks = K(1); %substrate saturation for monod kinetics (kgSc/dm3W)
Khs = K(2); %coefficient de saturation pour les cinétiques de contois
fi = K(3); %Proportion of dead biomass recycled to inert materials
kdec = K(4); %microorganisms decomposition constant


%hydrolysis constant
kh1C = kh(1)*1.4;	
kh2P = kh(2)*1.4;
kh3L = kh(3)*1.4;	
kh4C = kh(4)*1.4;	
kh5P = kh(5)*1.4;	
kh6L = kh(6)*1.4;	
kh7H = kh(7);	
kh8CE = kh(8);	
kh9LG = kh(9);	
kh10H = kh(10);	
kh11CE = kh(11);	
kh12LG = kh(12);

%growth rate
mmb = mu(1);
mtb = mu(2);
mma = mu(3);
mta = mu(4);
mmf = mu(5);
mtf = mu(6);

%death rate
bmb = bd(1);
btb = bd(2);
bma = bd(3)*1;
bta = bd(4)*1;
bmf = bd(5)*1;
btf = bd(6)*1;

%literature values of kinetic parameters for nitrfication/denitrification
ma = kXa(1);
ba = kXa(2);
pmaxdenit = kXa(3);
pN2Odenit =kXa(4);

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

%limitation by the substrate
fSc = (Sc/W)/(Ks+(Sc/W));
fSp=(Sp/W)/(Ks+(Sp/W));
fSl=(Sl/W)/(Ks+(Sl/W));
fSh=(Sh/W)/(Ks+(Sh/W));
fSlg=(Slg/W)/(Ks+(Slg/W));

%hydrolisis rate
v1 = kh1C*Xmb*(C/((Khs*Xmb)+C));
v2 =kh2P* Xmb*(P/((Khs*Xmb)+P)); 
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
v13 = mmb*Xmb  *fSc ;
v14 = mmb*Xmb  * fSp ;
v15 = mmb*Xmb  *fSl;
%Thermophilic bacterias growth
v16 = mtb*Xtb  * fSc;
v17 = mtb*Xtb *fSp ;
v18 = mtb*Xtb *fSl;
%Mesophilic Actinomycetes growth
v19 = mma*Xma   * fSc;
v20 = mma*Xma    *fSp;
v21 = mma*Xma  *fSl;
v22 = mma*Xma   *fSh ;
 %Thermophilic actinomycetes growth
v23 = mta*Xta  *fSc;
v24 = mta*Xta  *fSp;
v25 = mta*Xta  *fSl;
v26 = mta*Xta  *fSh;
%Mesophilic fungis growth
v27 = mmf*Xmf  *fSc;
v28 = mmf*Xmf  * fSp;
v29 = mmf*Xmf  * fSl;
v30 = mmf*Xmf  *  fSh;
v31 = mmf*Xmf  * fSlg; 
%Thermophilic fungis growth
v32 = mtf*Xtf  *fSc ;
v33 = mtf*Xtf *fSp;
v34 = mtf*Xtf   *fSl; 
v35 = mtf*Xtf    *fSh;
v36 = mtf*Xtf   *fSlg ;

%Heterotroph Microorganisms death
v37 = bmb*Xmb ;
v38 = btb*Xtb ;
v39 = bma*Xma ;
v40 = bta*Xta;
v41 = bmf*Xmf;
v42 = btf*Xtf;

%lysis of microorganisms
v43 = kdec*Xdb ;

%Nitrification
v44 = ma*Xa %*fSn_nit ;
            
%Denitrification
v45 = pmaxdenit*NO3 %*fNO3 ; %*self.fT_denit 

% Autotroph microorganisms death
v46 = ba*Xa; 

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

dXadt = (v44) - v46;

dNO3dt = Y_NO3*v44 - v45;
            
dN2Odt = pN2Odenit*v45;

dN2dt = (1-pN2Odenit)*v45;

%C=max(0,C);
%Sc=max(0,Sc);
%Xmb=max(0,Xmb);
LG = max(0,LG);
CE = max(0,CE);

dydt = [dCdt,dPdt,dLdt,dHdt,dCEdt,dLGdt,dXidt, dScdt,dSpdt,dSldt,dShdt,dSlgdt, dXmbdt,dXtbdt,dXmadt,dXtadt,dXmfdt,dXtfdt, dXdbdt, dCO2dt, dWdt, dXadt, dNO3dt, dN2Odt, dN2dt]';