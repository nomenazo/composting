clear;
%% ce code se différencie de fitcompostcumul par les unités de C-CO2 qui sont en kg/kgTM au lieu de kg/kgsec


t=[0	24	48	72	96	120	144	168	192	216	240	264	288	312	336	360	384	408 ...
    432	456	480	504	528	552	576	600	624	648	672	696	720	744	768	792	816	840	864 ...
    888	912	936	960	984];	

dataCO2= [0	0.010658307	0.012064092	0.020488396	0.023870271	0.030051839	0.032267393  ...
    0.034299089	0.038006534	0.039503369	0.040739904	0.042901283	0.044411015	0.045248496 ...
    0.047254084	0.048646905	0.049735599	0.049812045	0.051551822	0.051782538	0.053578481 ...
    0.054737495	0.05584828	0.056941472	0.057660209	0.058723218	0.059132392	0.059721854 ...
    0.060876461	0.062067494	0.062106234	0.062848158	0.063068093	0.063059189	0.063047242 ...
    0.063277004	0.063419717	0.063739185	0.06404356	0.064904292	0.065162834	0.064998938];  ...


kh0 = [0.0001 0.0378*0.5 0.2991 2.8873e-05 0.0271*0.5 0.0153 0.009 0.0025 0.0078 0.009 0.0025 0.0313];

LB=[0 0 0 0 0 0 0 0 0 0 0 0]';
UB=[inf inf inf inf inf inf inf inf inf inf inf inf]';

%options = optimoptions('lsqcurvefit', ...
                     %  'StepTolerance', 1e-20, ...
                     %  , ...
                     %  'FunctionTolerance', 1e-8, ...
                     %  'OptimalityTolerance', 1e-8, ...
                     %  'Display', 'iter');
options = optimoptions('lsqcurvefit','StepTolerance', 1e-15,'MaxFunctionEvaluations', 50000, 'Maxiterations',35000,'FunctionTolerance', 1e-10, 'Display', 'iter');

[kh]=lsqcurvefit('rcompostfitcumul',kh0,t,dataCO2,LB,UB, options)
yfit=rcompostfitcumul(kh,t);
CO2fit=yfit(1,:);
%nP=Pfit-dataP;
%nSt=Stfit-dataSt;
%dP=dataP-(1./10).*sum(dataP);
%dSt=dataSt-(1./10)*sum(dataSt);
%R2P=1-(sum(nP.^2)./sum(dP.^2));
%R2St=1-(sum(nSt.^2)./sum(dSt.^2));;

subplot(1,1,1);
plot(t,CO2fit,'r-',t,dataCO2,'k.');
title('Production of CO2');
xlabel('time (h)');
ylabel('CO2-C (kg/kgTM)');
