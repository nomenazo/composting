clear;

% Initialisation des variables
y0 = [0.203655 0.049914 0.025361 0 0.001704 0.000159 0.010254 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.710953 293 0 0 0 1e-4 0 0 0 0 0.00036];

% Définition de xO2_values avec linspace
xO2_values = linspace(0.08, 0.24, 10); % Exemple de valeurs de xO2

% Définition de tspan avec linspace
%tspan = linspace(0, 2200, 91); % De t=0 à t=2200 en 91 points

tspan(1)=0;
for i=2:91
    tspan(i)=tspan(i-1)+24;
end;

% Options pour ode15s
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Pré-allocation pour stocker les résultats de CH4
output_CH4 = zeros(length(xO2_values), 1);

% Boucle sur les valeurs de xO2 pour l'analyse de sensibilité de CH4
for i = 1:length(xO2_values)
    % Résolution des ODE avec ode15s pour chaque valeur de xO2
    params = xO2_values(i);
    [~, y] = ode15s(@(t, y) compostsensanalysis(t, y, params), tspan, y0, options);

    y(y < 0) = 0;

    % Récupération de CH4 à la fin de la simulation
    CH4_cumulative = cumsum(y(:,25)); %cumtrapz(tspan, y(:, 25)); 

    %CO2_final = y(end, 20);


    % Stockage des résultats (CH4 cumulée à la fin de la simulation)
    output_results(i) =  CH4_cumulative(end); %CO2_final;

end

% Affichage des résultats de la sensibilité de CH4 cumulée par rapport à xO2
figure;
plot(xO2_values, output_results, 'o-');
xlabel('Valeur de xO2');
ylabel('CH4 Cumulé final');
%ylabel('CO2 final (kg/kgTM');
%title('Sensibilité de CO2 Cumulé par rapport à xO2');
grid on;