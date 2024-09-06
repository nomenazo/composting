clear;

% Initialisation des variables
y0=[0.1208 0.023 0.0102 0.63e-3 0.00349 0.0008 0.009 0 0 0 0 0 5e-4 5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.60 293 0 0 0 1e-4 0 0 0 0 3e-5] ;
%Yang waste composition

tspan(1)=0;
for i=2:36
    tspan(i)=tspan(i-1)+24;
end;

% Définition de xO2_values avec linspace
FAS_values = linspace(0.1, 0.5, 4); % Exemple de valeurs de xO2


% Options pour ode15s
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Pré-allocation pour stocker les résultats de CH4
output_CH4 = zeros(length(FAS_values), 1);

% Boucle sur les valeurs de xO2 pour l'analyse de sensibilité de CH4
for i = 1:length(FAS_values)
    % Résolution des ODE avec ode15s pour chaque valeur de xO2
    params = FAS_values(i);
    [~, y] = ode15s(@(t, y) sensanalch4(t, y, params), tspan, y0, options);

    y(y < 0) = 0;

    % Récupération de CH4 à la fin de la simulation
    %CH4_cumulative = cumsum(y(:,25)); %cumtrapz(tspan, y(:, 25)); 

    %CO2_final = y(end, 20);

    %NH3_final = y(end, 30);

    N2O_final = y(end, 28);


    % Stockage des résultats (CH4 cumulée à la fin de la simulation)
    output_results(i) = N2O_final; %NH3_final; % CO2_final; %CH4_cumulative(end); %

end

% Affichage des résultats de la sensibilité de CH4 cumulée par rapport à xO2
figure;
plot(FAS_values, output_results, 'o-');
xlabel('Valeur de FAS');
%ylabel('CH4 Cumulé final');
ylabel('N2O final (kg/kgTM');
%title('Sensibilité de CO2 Cumulé par rapport à xO2');
grid on;