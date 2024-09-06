function z=rcompostfitvf(kh,t);
global thetag;
%thetag=bd;
thetag = kh;

%initialisation des variables
y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.710953 293 0 0 0 1e-4 0 0 0 0 0.00036];%O2init %kg/kgTM %2.6954e-4
%y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-2 5e-2 1e-4 1e-4 5e-6 5e-6 0 0 0.710953 293 0 0 0 1e-5 0 0 0 0 0.00036];%O2init %kg/kgTM %2.6954e-4



tspan(1)=0;
for i=2:96
   tspan(i)=tspan(i-1)+24;
end;

%tspan = [0 2200];

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14);

%options = odeset('Events', @events);
%[t,y]=ode45_with_corrections('compostfitted',[tspan],[y0]);
[t,y]=ode15s('compostfitvf',[tspan],[y0],options)




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

simulCO2 = [CO2(1), CO2(2), CO2(3), CO2(4), CO2(5), CO2(6), CO2(7), CO2(8), CO2(9), CO2(10), CO2(11), CO2(12), ...
    CO2(13), CO2(14), CO2(15), CO2(16), CO2(17), CO2(18), CO2(19), CO2(20), CO2(21), CO2(22), CO2(23), CO2(24), ...
    CO2(25), CO2(26), CO2(27), CO2(28), CO2(29), CO2(30), CO2(31), CO2(32), CO2(33), CO2(34), CO2(35), CO2(36), ...
    CO2(37), CO2(38), CO2(39), CO2(40), CO2(41), CO2(42), CO2(43), CO2(44), CO2(45), CO2(46), CO2(47), CO2(48), ...
    CO2(49), CO2(50), CO2(51), CO2(52), CO2(53), CO2(54), CO2(55), CO2(56), CO2(57), CO2(58), CO2(59), CO2(60), ...
    CO2(61), CO2(62), CO2(63), CO2(64), CO2(65), CO2(66), CO2(67), CO2(68), CO2(69), CO2(70), CO2(71), CO2(72), ...
    CO2(73), CO2(74), CO2(75), CO2(76), CO2(77), CO2(78), CO2(79), CO2(80), CO2(81), CO2(82), CO2(83), CO2(84), ...
    CO2(85), CO2(86), CO2(87), CO2(88), CO2(89), CO2(90), CO2(91), CO2(92), CO2(93), CO2(94), CO2(95), CO2(96)];

%q = cumsum(simulCO2);
simulCCO2 = (0.012/0.044)*simulCO2;

z= [simulCCO2];
