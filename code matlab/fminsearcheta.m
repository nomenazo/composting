% Initialiser la valeur de eta
eta_initial = 4e5 * 0.02;

% Utiliser fminsearch pour trouver la valeur optimale de eta
optimal_eta = fminsearch(@costFunction, eta_initial);

% Afficher le résultat
disp(['Valeur optimale de eta : ', num2str(optimal_eta)]);
