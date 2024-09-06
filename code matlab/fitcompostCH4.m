clear;
%% ce code se différencie de fitcompostcumul par les unités de C-CO2 qui sont en kg/kgTM au lieu de kg/kgsec


%t=[0 24	48	72	96	120	144	168	192	216	240	264	288	312	336	360	384	408 ...
   % 432	456	480	504	528	552	576	600	624	648	672	696	720	744	768	792	816	840	];

t(1)=0;
for i=2:36
    t(i)=t(i-1)+24;
end;

dataCH4d = [0 1.95155E-06 6.48274E-06 2.60566E-05 6.56842E-05 5.55233E-05 8.78856E-05 9.20516E-05 0.000114093 0.000103701 ...
8.80915E-05 0.000102455 0.000100069 8.53698E-05 7.47548E-05 6.26267E-05 5.13218E-05 2.99409E-05 2.60347E-05 2.16979E-05 ...
1.9065E-05 1.99399E-05 9.11593E-06 6.94192E-06 4.40578E-06 1.7358E-06 1.42009E-06 1.40639E-06 1.63831E-06 1.43084E-06 ...
1.18318E-06 1.39761E-06 1.27638E-06 2.11566E-06 2.59064E-06 2.59241E-06];


kCH40  = 100*[0.267 0.375 0.707 0.284 0.535]; % 0.72 0.033];



LB=[0 0 0 0 0]';
UB=[inf inf inf inf inf]';

%options = optimoptions('lsqcurvefit', ...
                     %  'StepTolerance', 1e-20, ...
                     %  , ...
                     %  'FunctionTolerance', 1e-8, ...
                     %  'OptimalityTolerance', 1e-8, ...
                     %  'Display', 'iter');
options = optimoptions('lsqcurvefit','StepTolerance', 1e-25,'MaxFunctionEvaluations',...
    50000, 'Maxiterations',35000,'FunctionTolerance', 1e-30, 'OptimalityTolerance', 1e-16,'Display', 'iter');




[kCH4]=lsqcurvefit('rcompostfitCH4',kCH40,t,dataCH4d,LB,UB, options) 

%[kh]=lsqcurvefit('rcompostfitNH3',kh0,t,dataNH3,LB,UB, options)
yfit=rcompostfitCH4(kCH4,t);
CH4fit=yfit(1,:);


subplot(1,1,1);
plot(t,CH4fit,'r-',t,dataCH4d,'k.');
title('Production of Ch4');
xlabel('time (h)');
ylabel('Ch4(kg/kgTM)');
