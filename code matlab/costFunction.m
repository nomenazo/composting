function cost = costFunction(eta)
    % Données observées (exemple, à remplacer par vos données)
    dataCH4d = [0 1.95155E-06 6.48274E-06 2.60566E-05 6.56842E-05 5.55233E-05 8.78856E-05 9.20516E-05 0.000114093 0.000103701 ...
8.80915E-05 0.000102455 0.000100069 8.53698E-05 7.47548E-05 6.26267E-05 5.13218E-05 2.99409E-05 2.60347E-05 2.16979E-05 ...
1.9065E-05 1.99399E-05 9.11593E-06 6.94192E-06 4.40578E-06 1.7358E-06 1.42009E-06 1.40639E-06 1.63831E-06 1.43084E-06 ...
1.18318E-06 1.39761E-06 1.27638E-06 2.11566E-06 2.59064E-06 2.59241E-06];
    
    tspan(1)=0;
    for i=2:36
        tspan(i)=tspan(i-1)+24;
    end;

    % Paramètres initiaux (remplacez par les valeurs appropriées)
    y0=[0.1208 0.023 0.0102 0.63e-3 0.00349 0.0008 0.009 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.60 293 0 0 0 1e-4 0 0 0 0 3e-5 0]

    % Options pour l'ODE solver
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);

    % Résoudre l'ODE avec le paramètre eta
    [T,Y] = ode45(@(t,y) fitch4(t, y, eta), tspan, y0, options);

    % Extraire les données simulées correspondantes (par exemple, la concentration de CH4)
    simulatedData = Y(:,25);  % Exemple pour CH4

    % Calculer la fonction de coût (erreur quadratique moyenne)
    cost = mean((dataCH4d - simulatedData).^2);
end
