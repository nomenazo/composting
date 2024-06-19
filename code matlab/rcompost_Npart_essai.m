clear;

%initialisation des variables
y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 5e-6 5e-6 0 0 0.710953 293 0 0 0 1e-5 0 0 0 0 0.00036];%O2init %kg/kgTM %2.6954e-4
%y0=[0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-8 5e-8 1e-4 1e-4 5e-6 5e-6 0 0 0.710953 293 0 0 0 1e-5 0 0 0 0 0.00036];%O2init %kg/kgTM %2.6954e-4



%tspan(1)=0;
%for i=2:42
 %  tspan(i)=tspan(i-1)+24;
%end;

tspan = [0 2000];

options = odeset( 'RelTol',1e-12,'AbsTol',1e-14); %'Maxstep', 1e-15,

%options = odeset('Events', @events);
%[t,y]=ode45_with_corrections('compostfitted',[tspan],[y0]);
[t,y]=ode15s('compostfitted_Npart',[tspan],[y0],options)




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

dataCO2= [0	0.010658307	0.012064092	0.020488396	0.023870271	0.030051839	0.032267393  ...
    0.034299089	0.038006534	0.039503369	0.040739904	0.042901283	0.044411015	0.045248496 ...
    0.047254084	0.048646905	0.049735599	0.049812045	0.051551822	0.051782538	0.053578481 ...
    0.054737495	0.05584828	0.056941472	0.057660209	0.058723218	0.059132392	0.059721854 ...
    0.060876461	0.062067494	0.062106234	0.062848158	0.063068093	0.063059189	0.063047242 ...
    0.063277004	0.063419717	0.063739185	0.06404356	0.064904292	0.065162834	0.064998938];

%subplot(1,1,1);
plot(t,CCO2,'r-') %,t,dataCO2,'k.');
%title('Production of CO2');
%xlabel('time (h)');
%ylabel('CO2-C (kg/kgTM)');




%subplot(7,1,2);
%plot(t,CCO2);
%xlabel('Time (h)');
%ylabel('CCO2 (kg/kgTM)');
%subplot(7,1,3);
%plot(t,O2diss);
%xlabel('Time (h)');
%ylabel('O2diss (kg/kgTM)');
%subplot(7,1,4);
%plot(t,P);
%xlabel('Time (h)');
%ylabel('P (kg/kgTM)');
%subplot(7,1,5);
%plot(t,L);
%xlabel('Time (h)');
%ylabel('L (kg/kgTM)');
%subplot(5,1,5);
%plot(t,O2in);
%xlabel('Time (h)');
%ylabel('O2in (kg/kgTM)');