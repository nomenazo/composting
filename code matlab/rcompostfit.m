function z=rcompostfit(kh,t);
global thetag;
%thetag=bd;
thetag = kh;

%initialisation des variables
y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 0 0 0.710953];

tspan(1)=0;
for i=2:42
   tspan(i)=tspan(i-1)+24;
end;

[t,y]=ode15s('compostfit',[tspan],[y0]);


%[t,y]=ode15s('compost',[0 400],y0);
%[t,y]=ode15s('compostfit',[0 1000],y0);


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

simulCO2=[CO2(1)	CO2(2)	CO2(3)	CO2(4)	CO2(5)	CO2(6)	CO2(7)	CO2(8)	CO2(9)	CO2(10) ...
    CO2(11)	CO2(12)	CO2(13)	CO2(14)	CO2(15)	CO2(16)	CO2(17)	CO2(18)	CO2(19)	CO2(20)	CO2(21)	CO2(22) ...
    CO2(23)	CO2(24)	CO2(25)	CO2(26)	CO2(27)	CO2(28)	CO2(29)	CO2(30)	CO2(31)	CO2(32)	CO2(33)	CO2(34)	CO2(35) ...
    CO2(36)	CO2(37)	CO2(38)	CO2(39)	CO2(40)	CO2(41)	CO2(42)];

z=[simulCO2];