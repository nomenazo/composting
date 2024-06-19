function z=rcompostfitN(kinXa,t);
global thetag;
thetag = kinXa;

%initialisation des variables
y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 5e-6 5e-6 0 0 0.710953 293 0 0 0 1e-5 0 0 0 0 0.00036]; %O2init %kg/kgTM %2.6954e-4


tspan(1)=0;
for i=2:91
   tspan(i)=tspan(i-1)+24;
end;

%tspan = [0 5000];

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14); %'Maxstep', 1e-15,

%options = odeset('Events', @events);
%[t,y]=ode45_with_corrections('compostfitted',[tspan],[y0]);
[t,y]=ode15s('compostfit2N',[tspan],[y0],options);




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


CCO2= (0.012/0.044) * CO2;
NNH3 = (0.014/0.017)*NH3;

%y = max(0,y);
%CH4gen = max(0,CH4gen);
y(y<0)=0;

simulNH3 = [NH3(1), NH3(2), NH3(3), NH3(4), NH3(5), NH3(6), NH3(7), NH3(8), NH3(9), NH3(10), NH3(11), NH3(12), ...
    NH3(13), NH3(14), NH3(15), NH3(16), NH3(17), NH3(18), NH3(19), NH3(20), NH3(21), NH3(22), NH3(23), NH3(24), ...
    NH3(25), NH3(26), NH3(27), NH3(28), NH3(29), NH3(30), NH3(31), NH3(32), NH3(33), NH3(34), NH3(35), NH3(36), ...
    NH3(37), NH3(38), NH3(39), NH3(40), NH3(41), NH3(42), NH3(43), NH3(44), NH3(45), NH3(46), NH3(47), NH3(48), ...
    NH3(49), NH3(50), NH3(51), NH3(52), NH3(53), NH3(54), NH3(55), NH3(56), NH3(57), NH3(58), NH3(59), NH3(60), ...
    NH3(61), NH3(62), NH3(63), NH3(64), NH3(65), NH3(66), NH3(67), NH3(68), NH3(69), NH3(70), NH3(71), NH3(72), ...
    NH3(73), NH3(74), NH3(75), NH3(76), NH3(77), NH3(78), NH3(79), NH3(80), NH3(81), NH3(82), NH3(83), NH3(84), ...
    NH3(85), NH3(86), NH3(87), NH3(88), NH3(89), NH3(90), NH3(91)];


simulNNH3 = (0.014/0.017)*simulNH3;

z= [simulNH3];