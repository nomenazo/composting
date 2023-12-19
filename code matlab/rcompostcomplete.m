clear;


y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 0 0 0.710953 0.001 0 0 0];
%y0=[0.143648615 0 1 0 0 0.721716301]; %Xmb en mg/kg;

%[t,y]=ode15s('compost',[0 400],y0);
[t,y]=ode45('compostcomplete',[0 2000],y0);


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
Xa=y(:,22);
NO3=y(:,23);
N2O=y(:,24);
N2=y(:,25);

subplot(4,1,1);
plot(t,Xa);
xlabel('Time (h)');
ylabel('Xa');
subplot(4,1,2);
plot(t,NO3);
xlabel('Time (h)');
ylabel('NO3');
subplot(4,1,3);
plot(t,N2O);
xlabel('Time (h)');
ylabel('N2O');
subplot(4,1,4);
plot(t,N2);
xlabel('Time (h)');
ylabel('N2');

