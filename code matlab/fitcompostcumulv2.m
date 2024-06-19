clear;
%% ce code se différencie de fitcompostcumul par les unités de C-CO2 qui sont en kg/kgTM au lieu de kg/kgsec


t=[0 24	48	72	96	120	144	168	192	216	240	264	288	312	336	360	384	408 ...
    432	456	480	504	528	552	576	600	624	648	672	696	720	744	768	792	816	840	864 ...
    888	912	936	960	984 1008	1032	1056	1080	1104	1128	1152 ...
    1176	1200	1224	1248	1272	1296	1320	1344	1368	1392 ...
    1416	1440	1464	1488	1512	1536	1560	1584	1608	1632 ...
    1656	1680	1704	1728	1752	1776	1800	1824	1848	1872 ...
    1896	1920	1944	1968	1992	2016	2040	2064	2088	2112 ...
    2136	2160	2184	2208	2232	2256	2280];	
%dataCO2=[0	0.030825453	0.038069825	0.054948653	0.07457083	0.08570688	0.100401451	0.118543452	0.135975422	0.146154358	0.152217452	0.16051648	0.166380789	0.172417868	0.176910715	0.184742316 ...
    %0.187609565	0.188741664	0.191862974	0.195344996	0.200357006	0.205940935	0.210422809	0.213240323	0.216394076	0.219748911	0.221979687	0.224193883	0.228225611	0.232231571	0.233427551	0.234536665 ...
    %0.236583391	0.236792542	0.236792542	0.237443205	0.238141787	0.238394479	0.239642021	0.242166207	0.243471514	0.244584];

dataCO2= [0	0.010658307	0.012064092	0.020488396	0.023870271	0.030051839	0.032267393  ...
    0.034299089	0.038006534	0.039503369	0.040739904	0.042901283	0.044411015	0.045248496 ...
    0.047254084	0.048646905	0.049735599	0.049812045	0.051551822	0.051782538	0.053578481 ...
    0.054737495	0.05584828	0.056941472	0.057660209	0.058723218	0.059132392	0.059721854 ...
    0.060876461	0.062067494	0.062106234	0.062848158	0.063068093	0.063059189	0.063047242 ...
    0.063277004	0.063419717	0.063739185	0.06404356	0.064904292	0.065162834	0.064998938  ...
    0.066034943	0.06645967	0.066818461	0.066716552	0.067741576	0.068010284	0.067975709	...
    0.068245469	0.06884277	0.069147466	0.069376341	0.069798145	0.06985545	0.069630769	...
    0.070680557	0.070905549	0.070939253	0.070902227	0.071318185	0.071257984	0.071469734	...
    0.071710832	0.071979559	0.072027407	0.072160354	0.072529026	0.072921378	0.072858746	...
    0.073311255	0.073726422	0.073782105	0.073784744	0.073763294	0.073971598	0.073904728	...
    0.074438034	0.074721419	0.07484848	0.075038274	0.075373127	0.075577342	0.075553739	...
    0.075661505	0.075922374	0.075913674	0.075927796	0.07593392	0.076546655	0.076825088	...
    0.076796612	0.07708573	0.077181414	0.077353359	0.076898372];

%dataCO2 = [0	0.030825453	0.007244372	0.016878828	0.019622177	0.011136051	0.01469457	0.018142001	0.01743197	0.010178936	0.006063094	0.008299028	0.005864309	0.006037079	0.004492846	0.007831602	...
  %  0.002867248	0.001132099	0.00312131	0.003482022	0.005012011	0.005583929	0.004481874	0.002817514	0.003153752	0.003354835	0.002230776	0.002214196	0.004031728	0.00400596	0.00119598	0.001109114	...
  %  0.002046726	0.000209151	0	0.000650663];	%0.000698582	0.000252692	0.001247541	0.002524187	0.001305307	0.001112485];


%kh0 = [0.04*1.2 0.02*0.5 0.01*5e-2 0.02 0.04*0.5 0.01*5e-2 0.009 0.007 0.007 0.009 0.007 0.007];
%kh0 = [0.0001*1.2  0.0378*0.5  0.2991*5e-2  2.8873e-05 0.0271*0.5 0.0153*5e-2 0.009 0.0025 0.0078 0.009 0.0025 0.0313];
kh0 = [0.04 0.02 0.01 0.02 0.04 0.01 0.009 0.007 0.007 0.009 0.007 0.007];

LB=[0 0 0 0 0 0 0 0 0 0 0 0]';
UB=[inf inf inf inf inf inf inf inf inf inf inf inf]';

%options = optimoptions('lsqcurvefit', ...
                     %  'StepTolerance', 1e-20, ...
                     %  , ...
                     %  'FunctionTolerance', 1e-8, ...
                     %  'OptimalityTolerance', 1e-8, ...
                     %  'Display', 'iter');
options = optimoptions('lsqcurvefit','StepTolerance', 1e-15,'MaxFunctionEvaluations', 50000, 'Maxiterations',35000,'FunctionTolerance', 1e-15, 'Display', 'iter');

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
