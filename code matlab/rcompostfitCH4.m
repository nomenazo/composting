function z=rcompostfitCH4(kCH4,t);
global thetag;

thetag = kCH4;

%initialisation des variables
y0=[0.1208 0.023 0.0102 0.63e-3 0.00349 0.0008 0.009 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.60 293 0 0 0 1e-4 0 0 0 0 3e-5];%O2init %kg/kgTM %2.6954e-4
%y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-2 5e-2 1e-4 1e-4 5e-6 5e-6 0 0 0.710953 293 0 0 0 1e-5 0 0 0 0 0.00036];%O2init %kg/kgTM %2.6954e-4



tspan(1)=0;
for i=2:36
    tspan(i)=tspan(i-1)+24;
end;

%tspan = [0 2200];

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14);

%options = odeset('Events', @events);
%[t,y]=ode45_with_corrections('compostfitted',[tspan],[y0]);
[t,y]=ode15s('compostfitch4',[tspan],[y0],options);




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






%y = max(0,y);
CH4 = max(0,CH4);

NH4 = max(0,NH4);

CCO2= (0.012/0.044) * CO2;

CCH4 = (0.012/0.016)*CH4;

NNH3 = (0.014/0.017)*NH3;


NN2 = N2;

NN2O = (0.028/0.044)*N2O;

Cinit= 0.1288; 
Ninit = 0.008/0.6;

Ccompost= Cinit - CCO2 - cumsum(CCH4);
Cemit = (CCO2+cumsum(CCH4))/Cinit;

Ncompost = Ninit - NNH3- NN2 - NN2O;
Nemit = (NNH3+NN2+NN2O)/Ninit;

CNratio = Ccompost/Ncompost;

y(y<0)=0;

simulCH4 = [CH4(1), CH4(2), CH4(3), CH4(4), CH4(5), CH4(6), CH4(7), CH4(8), CH4(9), CH4(10), CH4(11), CH4(12), ... 
            CH4(13), CH4(14), CH4(15), CH4(16), CH4(17), CH4(18), CH4(19), CH4(20), CH4(21), CH4(22), CH4(23), CH4(24), ...
            CH4(25), CH4(26), CH4(27), CH4(28), CH4(29), CH4(30), CH4(31), CH4(32), CH4(33), CH4(34), CH4(35), CH4(36)];




z= [simulCH4];
