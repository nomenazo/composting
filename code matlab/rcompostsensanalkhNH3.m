clear;

% Initialisation des variables
y0 = [0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.710953 293 0 0 0 1e-4 0 0 0 0 0.00036];

% Définition de xO2_values avec linspace
%Yxa_nh4_values = linspace(0.01, 0.9, 100); % Exemple de valeurs de xO2
%ma_values = linspace(0.0001,0.01,100);
%ma_values = linspace(1e-6,0.01,100);
%kO2nit_values = linspace(1e-6,0.01,100);
Qair_values = linspace(0.001,0.03,100);

% Définition de tspan avec linspace
tspan = linspace(0, 2200, 91); % De t=0 à t=2200 en 91 points

% Options pour ode15s
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Pré-allocation pour stocker les résultats de CH4
%output_CH4 = zeros(length(Yxa_nh4_values), 1);

% Boucle sur les valeurs de xO2 pour l'analyse de sensibilité de CH4
for i = 1:length(Qair_values)
    params = Qair_values(i);
    [~, y] = ode15s(@(t, y) compostsensanalysiskhNH3(t, y, params), tspan, y0, options);

    y(y < 0) = 0;

    NH3final = y(end, 30);

    
    output_results(i) = NH3final; %CH4_cumulative(end);

end

% Affichage des résultats de la sensibilité de CH4 cumulée par rapport à xO2
figure;
plot(Qair_values, output_results, 'o-');
xlabel('Valeur de xO2');
%ylabel('CH4 Cumulé final');
ylabel('NH3cumul (kg/kgTM)');
title('Sensibilité de NH3 par rapport à Qair');
grid on;