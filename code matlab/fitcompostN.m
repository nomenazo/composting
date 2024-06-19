clear;

t=[0 24	48	72	96	120	144	168	192	216	240	264	288	312	336	360	384	408 ...
    432	456	480	504	528	552	576	600	624	648	672	696	720	744	768	792	816	840	864 ...
    888	912	936	960	984 1008	1032	1056	1080	1104	1128	1152 ...
    1176	1200	1224	1248	1272	1296	1320	1344	1368	1392 ...
    1416	1440	1464	1488	1512	1536	1560	1584	1608	1632 ...
    1656	1680	1704	1728	1752	1776	1800	1824	1848	1872 ...
    1896	1920	1944	1968	1992	2016	2040	2064	2088	2112 ...
    2136	2160	2184];	

dataNH3 = [9.50523E-05	8.3184E-05	8.5979E-05	9.59316E-05	0.000113998	0.000148685	0.000275354	0.000453294	0.000663354 ...
    0.000877706	0.001094094	0.001477452	0.001856072	0.002462657	0.003061782	0.003358049	0.003850929	0.004161586	0.00448449 ...
    0.00491856	0.005432504	0.005863488	0.006036313	0.006277233	0.006532356	0.006828961	0.007146364	0.007286909	0.007411799 ...
    0.00759065	0.007764798	0.007871677	0.007871931	0.007945459	0.008094275	0.008162866	0.008215585	0.008415265	0.008606776 ...
    0.008647248	0.008756188	0.008901962	0.009115249	0.009276016	0.009311398	0.009375908	0.009428734	0.009494994	0.009580957	...
    0.009629212	0.009666273	0.009696949	0.00973921	0.009769989	0.009781303	0.009807082	0.009830799	0.009861079	0.00987839 ...
    0.009883663	0.009893537	0.009899967	0.009941319	0.009968676	0.009983343	0.010006289	0.009984702	0.009981679	0.010014646	...
    0.01003884	0.010053253	0.010052459	0.010067748	0.010092195	0.010103562	0.010108216	0.010127495	0.010150719	0.010153972	...
    0.010160904	0.010176637	0.010192993	0.010204598	0.010211709	0.010217376	0.010233124	0.01024688	0.010250307	0.010252842	0.01027212	0.010298132];

kinXa0 = 0.003 ;

LB = 0';
UB = inf';

options = optimoptions('lsqcurvefit','StepTolerance', 1e-50,'MaxFunctionEvaluations', 50000, 'Maxiterations',35000,'FunctionTolerance', 1e-15,'OptimalityTolerance', 1e-20, 'Display', 'iter');

[kinXa]=lsqcurvefit('rcompostfitN',kinXa0,t,dataNH3,LB,UB, options)
yfit=rcompostfitN(kinXa,t);
NH3fit=yfit(1,:);

subplot(1,1,1);
plot(t,NH3fit,'r-') ; %,t,dataNH3,'k.');
title('Production of NH3');
xlabel('time (h)');
ylabel('NH3-N (kg/kgTM)');