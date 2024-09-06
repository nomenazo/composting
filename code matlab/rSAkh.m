clear;

% Initialisation des variables
y0=[0.13852 0.03307 0.03812 0.03199 0.07553 0.07036 0.01211 0 0 0 0 0 5e-4 ...
5e-4 1e-4 1e-4 1e-6 1e-6 0 0 0.5554 293 0 0 0 1e-4 0 0 0 0 0.00031 0]; %aFinland+WC 3:1

%Values of variables
N2O0 = 2.2553e-4;
CH40 = 0.9214e-4;
CO20 = 0.2655;
NH30 = 0.0018;


tspan(1)=0;
for i=2:180; %91
    tspan(i)=tspan(i-1)+24;
end;

% Initial values of parameters
%kh0 = [0.0293    0.1508    1.4563e-07    0.0053    0.1731    0.0182    0.0090 0.0070    0.0068    0.0090    0.0070    0.0068];
%mu0 = [0.2 0.18 0.1 0.12 0.1 0.1 0.03];

%bd0 = [0.03 0.02 0.01 0.015 0.01 0.01 0.0083];
%K=[6.2e-5 1e-4 0.2 0.0025 0.0064 0.0608e-5];
KCH40 = [0.267 0.375 0.707 0.284 0.535 4e5 5.35e-4 0.72 0.033];
%Variation percentage

var = [-0.75, -0.25, 0.25, 0.75];


% Options pour ode15s
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-18);

% Pré-allocation pour stocker les résultats
output_CH4 = zeros(length(KCH40), length(var),1);
output_CO2 = zeros(length(KCH40), length(var),1);
output_N2O = zeros(length(KCH40), length(var),1);
output_NH3 = zeros(length(KCH40), length(var),1);



for i = 1:length(KCH40)
    for j = 1:length(var)
        KCH4 = KCH40;
        % positive values of var
        KCH4(i)=KCH40(i)*(1+var(j));
        [~, y_var_plus] = ode15s(@(t, y) SAKCH4(t, y, KCH4), tspan, y0, options);
        y_var_plus(y_var_plus < 0) = 0;

        % Stockage des résultats
        CH4endpos = cumsum(y_var_plus(:,25));
        output_CH4(i,j,1)= ((CH4endpos(end)-CH40)/CH40)*100;
        output_CO2(i,j,1) = ((y_var_plus(end, 20)-CO20)/CO20)*100;
        output_N2O(i,j,1)= ((y_var_plus(end, 28)-N2O0)/N2O0)*100;
        output_NH3(i,j,1) = ((y_var_plus(end, 30)-NH30)/NH30)*100;

    end
end

% Affichage des résultats
% Indice pour remplir la matrice
idx = 1;

% Boucle principale sur les paramètres
% Nombre de paramètres et de variations
n_params = length(KCH40);
n_vars = length(var);
n_gases = 4; % N2O, CH4, CO2, NH3

% Initialisation de la matrice pour stocker les résultats
results_matrix = zeros(n_params, n_vars * n_gases);

% Remplissage de la matrice
for i = 1:n_params
    % Chaque gaz remplit une partie de la matrice
    results_matrix(i, 1:n_vars) = output_N2O(i, :, 1);
    results_matrix(i, n_vars+1:n_vars*2) = output_CH4(i, :, 1);
    results_matrix(i, n_vars*2+1:n_vars*3) = output_CO2(i, :, 1);
    results_matrix(i, n_vars*3+1:n_vars*4) = output_NH3(i, :, 1);
end

% Ajout des labels pour les lignes
%row_labels = {'kh1C', 'kh2P', 'kh3L', 'kh4C', 'kh5P', 'kh6L', 'kh7H', 'kh8CE', 'kh9LG', 'kh10H', 'kh11CE', 'kh12LG'};



% Afficher la table
disp(results_matrix);
writematrix(results_matrix, 'var_KCH4.xlsx');

%%for i = 1:length(kh0)
    %fprintf('\n');
    
    % Boucle sur les gaz
    %fprintf('N2O:\n');
   %% for j = 1:length(var)
      %%  fprintf('%f\n', output_N2O(i, j, 1));
   %% end
    
    %fprintf('CH4:\n');
    %%for j = 1:length(var)
     %%   fprintf('%f\n',  output_CH4(i, j, 1));
   %% end
    
    %fprintf('CO2:\n');
   %% for j = 1:length(var)
     %%   fprintf('%f\n', output_CO2(i, j, 1));
  %%  end
    
    %fprintf('NH3:\n');
   %% for j = 1:length(var)
     %%   fprintf('%f\n', output_NH3(i, j, 1));
  %%  end
    
  %%  fprintf('\n'); % Saut de ligne pour séparer les paramètres
%%end
