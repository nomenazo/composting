clear;

%initialisation des variables

%y0=[0.14420 0.03928 0.03018 0.03199 0.07666 0.07061 0.01423 0 0 0 0 0 5e-4 ...
%5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.55223 293 0 0 0 1e-4 0 0 0 0 0.00026 0 0.003 5e-2];
%%Data Italy+WC 3:1

%y0=[0.1709 0.04655 0.03577 0.01422 0.036024 0.03177 0.01288 0 0 0 0 0 5e-4 ...
%5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.6360 293 0 0 0 1e-4 0 0 0 0 0.00031 0 0.003 5e-2]; %Italy+WC 8:1

%y0=[0.13852 0.03307 0.03812 0.03199 0.07553 0.07036 0.01211 0 0 0 0 0 5e-4 ...
%5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.5554 293 0 0 0 1e-4 0 0 0 0 0.00031 0 0.003 5e-2]; %aFinland+WC 3:1 ; O2gas: 0.21*1.42*0.4*(((3*1.5)/2)*6)/100000

y0=[0.16417 0.03919 0.04518 0.01422 0.03468 0.03147 0.01274 0 0 0 0 0 5e-4 ...
5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.6397 293 0 0 0 1e-4 0 0 0 0 0.00021 0 0.003 5e-2 0 0 0 0]; %aFinland+WC 8:1


tspan(1)=0;
for i=2:180; %91
    tspan(i)=tspan(i-1)+24;
end;

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14);

[t,y]=ode15s('windrow_vf6',[tspan],[y0],options);



C=y(:,1);
P=y(:,2);
L=y(:,3);
H=y(:,4);
CE=y(:,5);
LG=y(:,6);
Xi=y(:,7);
Sc=y(:,8);
Sp=y(:,9);
Sl=y(:,10);
Sh=y(:,11);
Slg=y(:,12);
Xmb=y(:,13);
Xtb=y(:,14);
Xma=y(:,15);
Xta=y(:,16);
Xmf=y(:,17);
Xtf=y(:,18);
Xdb=y(:,19);
CO2=y(:,20);
W=y(:,21);
T=y(:,22);
CH4gen=y(:,23);
CH4oxi=y(:,24);
CH4=y(:,25);
Xa=y(:,26);
NO3=y(:,27);
N2O=y(:,28);
N2=y(:,29);
NH3=y(:,30);
NH4=y(:,31);
Wvap = y(:,32);
O2diss = y(:,33);
O2gaz = y(:,34);
O2cons = y(:,35);
O2out = y(:,36);
air = y(:,37);
NH3gaz = y(:,38);

%y = max(0,y);
CH4 = max(0,CH4);
CH4gen = max(0,CH4gen);
CH4oxi = max(0,CH4oxi);

NH4 = max(0,NH4);

CCO2= (0.012/0.044) * CO2;
OCO2 = (0.032/0.044) * CO2;

CCH4 = (0.012/0.016)*CH4;

CCH4gen = (0.012/0.016)*CH4gen;

HCH4gen = (0.004/0.016)*CH4gen;

HH20comp = (0.002/0.018)*W;
OH20comp = (0.016/0.018)*W;

HH20vap = (0.002/0.018)*Wvap;
OH20vap = (0.016/0.018)*Wvap;

NNH3 = (0.014/0.017)*NH3;
HNH3 = (0.003/0.017)*NH3;

HNH4 = (0.004/0.018)*NH4;

NNH3gaz=(0.014/0.017)*NH3gaz;
HNH3gaz = (0.003/0.017)*NH3gaz;


NN2 = N2;
NNH4 = (0.014/0.018)*NH4;
NNO3 = (0.014/0.062)*NO3;

NN2O = (0.028/0.044)*N2O;
ON2O = (0.016/0.044)*N2O;

OO2diss = O2diss;
OO2gas = O2gaz;

%Cinit =0.1627;
Cinit = (C(1)*(6*0.012/0.180)) + (P(1)*(16*0.012/0.352)) + (L(1)*(25*0.012/0.393))+...
    (H(1)*(10*0.012/0.282))+(CE(1)*(6*0.012/0.180))+(LG(1)*(20*0.012/0.366)) +(Xmb(1)*(5*0.012/0.113))+(Xtb(1)*(5*0.012/0.113))+(Xma(1)*(5*0.012/0.113))+...
    (Xta(1)*(5*0.012/0.113))+(Xmf(1)*(10*0.012/0.247))+(Xtf(1)*(10*0.012/0.247))+(Xa(1)*(5*0.012/113));

Ninit = (P(1)*(4*0.014/0.352)) + ((Xmb(1)*(1*0.014/0.113))+(Xtb(1)*(1*0.014/0.113))+(Xma(1)*(1*0.014/0.113))+...
    (Xta(1)*(1*0.014/0.113))+(Xmf(1)*(1*0.014/0.247))+(Xtf(1)*(1*0.014/0.247))+(Xa(1)*(1*0.014/0.113)));
%Finland 8_1
% 0.1925; %F_3_1
% 0.1927; %It_3_1 
%0.1623; %It_8_1 
% 0.9797

%

%Ninit =0.0062; %F_8_1 
% 0.0053; %F_3_1 
% 0.0062; %It_3_1; 
% 0.0074; %It_8_1 
  


Hinit=0.096 ; %F_8_1
Oinit = 0.705+5e-2+0.003+21.3425



Ccompost= Cinit - CCO2 - cumsum(CCH4);
Cemit = (CCO2+cumsum(CCH4))/Cinit; %cumsum(CCH4)
Cemitcum = Cemit(end)

Hemit = (HCH4gen + HH20vap+ HH20comp+HNH3+HNH4+HNH3gaz)/Hinit;
Hemitcum = Hemit(end)

Ncompost = Ninit - NNH3- NN2 - NN2O-NNO3-NNH3gaz;
Nemit = (NNH3+NN2+NN2O+NNH3gaz)/Ninit;
Nmin = (NNO3+NNH4)/Ninit;
Nemitcum = Nemit(end)
Nmincum = Nmin(end)

Oemit = (OCO2+ON2O+OH20comp+OH20vap+O2out+OO2gas)/Oinit;
Oemitcum = Oemit(end);
Oemitcum

O2in = (OO2diss)/Oinit;
O2incum = O2in(end);
O2incum

CNratio = Ccompost/Ncompost;
%y(y<0)=0;








%%final composition
%emissions
CH4tot= cumsum(CH4);

CH4final = CH4tot(end)
CO2final= CO2(end)
%CH4final=cumsum(CH4)(end)
NH3final = NH3(end)
N2Ofinal=N2O(end)
N2final= N2(end)
Wevapfinal = Wvap(end)

Ccomfinal = Ccompost(end); %par rapport àTM
Ncomfinal = Ncompost(end);

N2O0 = 2.2556e-04;
CH40 = 0.655e-3;
CO20 = 0.2675;
NH30 = 6.8949e-04;

output_CH4 = ((CH4final - CH40) / CH40) * 100;
output_CO2 = ((CO2final - CO20) / CO20) * 100;
output_N2O = ((N2Ofinal - N2O0) / N2O0) * 100;
output_NH3 = ((NH3final - NH30) / NH30) * 100;

%masse finale de compost
Mcompost = (C(end) +P(end)+ L(end)+ H(end)+LG(end)+CE(end)+Sc(end)+Sp(end)+...
    Sl(end)+Sh(end)+Slg(end)+Xmb(end)+Xtb(end)+Xma(end)+Xta(end)+Xmf(end)+Xtf(end)+...
    Xdb(end)+Xi(end)+NH4(end)+W(end)+NO3(end)+O2diss(end)+Xa(end)+NNH3(end)) * 0.0809;


McompostTM = Mcompost/0.0809;
McompostTM

O2gaz = O2gaz(end)
O2gaz

MO = (C(end) +P(end)+ L(end)+ H(end)+LG(end)+CE(end)+Sc(end)+Sp(end)+...
    Sl(end)+Sh(end)+Slg(end)+Xmb(end)+Xtb(end)+Xma(end)+Xta(end)+Xmf(end)+Xtf(end)+...
    Xdb(end))/McompostTM
MO


MS = (Mcompost - (W(end)*0.0809))/Mcompost
MSfinal = MS(end);
MSfinal

Cprod = (Ccompost*0.0809)/Mcompost; %par rapport à la masse de compost
Cprodfinal = Cprod(end);
Cprodfinal
Nprod = (Ncompost*0.0809)/Mcompost;
Nprodfinal = Nprod(end);
Nprodfinal

Wvapfinal = cumsum(Wvap);

%%calcul des éléments
%calcul des éléments dans matières inertes formées C1.48H1.5O0.85N0.12
Xi_init = 0.01274; %kg/kgTM
Xi_formed = Xi(end) - Xi_init; %kg/kgTM
M_Xi = 0.03454 ; %kg/mol

M_element = [0.012, 0.014, 0.016, 0.001]; %masse molaire
n_element = [1.48, 0.12, 0.85, 1.5];

C_Xi =(Xi_formed*M_element(1)*n_element(1))/(M_Xi *Cinit);
N_Xi =(Xi_formed*M_element(2)*n_element(2))/(M_Xi *Ninit);
O_Xi =(Xi_formed*M_element(3)*n_element(3))/(M_Xi *Oinit);
H_Xi =(Xi_formed*M_element(4)*n_element(4))/(M_Xi *Hinit);

%calcul des éléments dans la matière organique du compost
compounds = {'C', 'P', 'L', 'H', 'CE', 'LG', 'Sc', 'Sp', 'Sl', 'Sh', 'Slg'};
Mmol = [0.180, 0.352, 0.393, 0.282, 0.180, 0.366, 0.180, 0.352, 0.393, 0.282, 0.366]; %masse molaire
atomC = [6, 16, 25, 10, 6, 20, 6, 16, 25, 10, 20];
atomH = [12, 24, 45, 18, 12, 30, 12, 24, 45, 18, 30];
atomO = [6, 5, 3, 9, 6, 6, 6, 5, 3, 9, 6];
atomN = [0, 4, 0,0,0,0,0,4,0,0,0];

% Initialisation du carbone organique
Corg = 0;
Horg = 0;
Oorg = 0;
Norg = 0;

% Boucle pour calculer la masse des élémnets pour chaque composé
for i = 1:length(compounds)
    end_value = eval([compounds{i}, '(end)']); % Récupérer la valeur finale du composé
    C_value = ((end_value / Mmol(i)) * atomC(i) * 12e-3) / Cinit;
    H_value = ((end_value / Mmol(i)) * atomH(i) * 1e-3) / Hinit;
    O_value = ((end_value / Mmol(i)) * atomO(i) * 16e-3) / Oinit;
    N_value = ((end_value / Mmol(i)) * atomN(i) * 14e-3) / Ninit;
    Corg = Corg + C_value;
    Horg = Horg + H_value;
    Oorg = Oorg + O_value;
    Norg = Norg + N_value;
end;

Corg
Horg
Oorg
Norg

%Calcul des éléments dans la biomasse vivante et morte du compost
biomass = {'Xmb', 'Xtb', 'Xma', 'Xta', 'Xmf', 'Xtf', 'Xa', 'Xdb'};
Mmol_X = [0.113, 0.113, 0.113, 0.113, 0.247, 0.247, 0.113, 0.113];
atomC_X = [5, 5, 5, 5, 10, 10, 5, 5];
atomH_X = [7, 7, 7, 7, 17, 17, 7, 7];
atomO_X = [2, 2,2, 2, 6, 6, 2, 2];
atomN_X = [1,1,1,1,1,1,1,1];

% Initialisation de la somme de carbone dans la biomasse
Cbio = 0;
Hbio = 0;
Obio = 0;
Nbio = 0;

% Boucle pour calculer la masse de carbone pour chaque composé X
for i = 1:length(biomass)
    end_value_X = eval([biomass{i}, '(end)']); % Récupérer la valeur finale du composé
    C_value_X = ((end_value_X / Mmol_X(i)) * atomC_X(i) * 12e-3) / Cinit;
    H_value_X = ((end_value_X / Mmol_X(i)) * atomH_X(i) * 1e-3) / Hinit;
    O_value_X = ((end_value_X / Mmol_X(i)) * atomO_X(i) * 16e-3) / Oinit;
    N_value_X = ((end_value_X / Mmol_X(i)) * atomN_X(i) * 14e-3) / Ninit;
    Cbio = Cbio + C_value_X;
    Hbio = Hbio + H_value_X;
    Obio = Obio + O_value_X;
    Nbio = Nbio + N_value_X;
end
Cbio
Hbio
Obio
Nbio


Ccheck = Cemitcum + Corg +Cbio + C_Xi %0.05669C_Xi/Cinit
Hcheck = Hemitcum + Horg + Hbio + H_Xi
Ocheck = Oemitcum + Oorg + Obio + O2incum + O_Xi
Ncheck = Nemitcum + Norg + Nbio + Nmincum + N_Xi

Ncompostcheck = (Norg +Nbio +Nmincum + 0.1974)*Ninit;

Nmineral = (Nmincum)/(Norg +Nbio +Nmincum+ 0.1974);

%valeurs initiales des variables
%compound
%%Carbon
C_compound_init = ((C(1)*(6*0.012/0.180)) + (P(1)*(16*0.012/0.352)) + (L(1)*(25*0.012/0.393))+...
    (H(1)*(10*0.012/0.282))+(CE(1)*(6*0.012/0.180))+(LG(1)*(20*0.012/0.366)))/Cinit;

C_biomass_init = ((Xmb(1)*(5*0.012/0.113))+(Xtb(1)*(5*0.012/0.113))+(Xma(1)*(5*0.012/0.113))+...
    (Xta(1)*(5*0.012/0.113))+(Xmf(1)*(10*0.012/0.247))+(Xtf(1)*(10*0.012/0.247)))/Cinit;

%%Nitrogen
N_compound_init = (P(1)*(4*0.014/0.352))/Ninit;

N_biomass_init = ((Xmb(1)*(1*0.014/0.113))+(Xtb(1)*(1*0.014/0.113))+(Xma(1)*(1*0.014/0.113))+...
    (Xta(1)*(1*0.014/0.113))+(Xmf(1)*(1*0.014/0.247))+(Xtf(1)*(1*0.014/0.247)))/Ninit;



%%sankey diagram
%for C
C_link = {'Ccompound', 'Cinit', C_compound_init; 'Cbioinit', 'Cinit', C_biomass_init;'Cinit','Cemit', Cemitcum ; 'Cinit', 'Corg', Corg ; 'Cinit', 'Cbiomass', Cbio;
    'Cinit', 'Cinert', C_Xi};

N_link = {'Ncompound', 'Ninit', N_compound_init; 'Nbioinit', 'Ninit', N_biomass_init;'Ninit','Nemit', Nemitcum ; 'Ninit', 'Norg', Norg ;
    'Ninit', 'Ninert', N_Xi; 'Ninit', 'Nmin', Nmincum; 'Ninit', 'Nbiomass', Nbio; 'Norg','Ncompost', Norg; 'Nbiomass', 'Ncompost', Nbio;
    'Nmin','Ncompost', Nmincum; 'Ninert','Ncompost', N_Xi };


SK.RenderingMethod='interp';  
% 修改对齐方式(Set alignment)
% 'up'/'down'/'center'(default)
SK.Align='center';
% 修改文本位置(Set text location)
% 'left'(default)/'right'/'top'/'center'/'bottom'
SK.LabelLocation='top';
SK.Sep=.1;
% 开始绘图(Start drawing)

SK_C=SSankey(C_link(:,1),C_link(:,2), C_link(:,3));

SK_N=SSankey(N_link(:,1),N_link(:,2), N_link(:,3));

%SK_C.draw()
SK_N.draw()




