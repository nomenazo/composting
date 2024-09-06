clear;

%initialisation des variables
%y0=[0.1012 0.0202 0.0075 0 0.0044 0.0008 0.010 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.86 293 0 0 0 1e-4 0 0 0 0 3e-5];%O2init %kg/kgTM %2.6954e-4

y0=[0.2247 0.0329 0.0162 4e-3 0.0034 0.0008 0.010 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.58 293 0 0 0 1e-4 0 0 0 0 3e-5 0 ] ; %biochemical composition from article


tspan(1)=0;
for i=2:50
    tspan(i)=tspan(i-1)+24;
end;

%tspan = [0 432];

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14);

%options = odeset('Events', @events);
%[t,y]=ode45_with_corrections('compostfitted',[tspan],[y0]);
[t,y]=ode15s('compostfitted_Peng',[tspan],[y0],options);




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





%y = max(0,y);
CH4 = max(0,CH4);

NH4 = max(0,NH4);

CCO2= (0.012/0.044) * CO2;

CCH4 = (0.012/0.016)*CH4;

NNH3 = (0.014/0.017)*NH3;


NN2 = N2;

NN2O = (0.028/0.044)*N2O;

Cinit= 0.1222; 
Ninit = 5e-3;

Ccompost= Cinit - CCO2 - cumsum(CCH4);

Ncompost = Ninit - NNH3- NN2 - NN2O;

CNratio = (Ccompost/Ncompost);
%y(y<0)=0;

texp(1)=0;
for i=2:19
    texp(i)=texp(i-1)+24;
end;

dataCO2= [0	0.009822223	0.022605962	0.035887098	0.043670426	0.059216425	0.075787411	0.087725741	0.100288143	0.110657032...
    0.1115704	0.113534804	0.121567152	0.130538493	0.13973814	0.148585943	0.155123658	0.160385789	0.165531836] ; %CO2

dataNH3 = [0	3.00805E-11	3.28907E-06	1.31573E-05	3.61749E-05	9.62537E-05	0.000172563	0.000271929	0.000501063	0.000560614	...
    0.000616367	0.000665734	0.000741557	0.000768856	0.00082108	0.000842463	0.000865504	0.000882548	0.000892408]; 

dataN2O = [0	4.53993E-05	5.4045E-05	5.96148E-05	6.39638E-05	7.14527E-05	7.70891E-05	8.17623E-05	8.42516E-05	8.8259E-05	...
    8.87089E-05	8.90272E-05	8.93473E-05	8.96674E-05	9.01802E-05	9.05643E-05	9.12308E-05	9.16794E-05	9.22556E-05]; 





subplot(3,1,1);
plot(t,NH3,'r-',texp,dataNH3,'k.');
xlabel('time (h)');
ylabel('NH3 (kg/kgTM)');
subplot(3,1,2);
plot(t,CO2,'r-',texp,dataCO2,'k.')
xlabel('time (h)');
ylabel('CO2 (kg/kgTM)');
subplot(3,1,3);
plot(t,N2O,'r-', texp,dataN2O,'k.')
xlabel('time (h)');
ylabel('N2O(kg/kgTM)');
%subplot(4,1,4);
%plot(t,CH4)
%xlabel('time (h)');
%ylabel('CH4(kg/kgTM)');
%subplot(5,1,5);
%plot(t,T)
%xlabel('time (h)');
%ylabel('T');


