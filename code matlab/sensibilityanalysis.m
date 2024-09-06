function sensitivity_analysis(params)
    % Définir les valeurs de paramètres pour l'analyse de sensibilité
    num_params = length(params);
    delta_params = params * 0.5; % Variation de 50% des paramètres

    % Initialiser les matrices pour stocker les résultats de sensibilité
    sensitivity_matrix = zeros(length(delta_params), length(params));

    % Boucle sur les paramètres pour calculer la sensibilité
    for i = 1:num_params
        perturbed_params = params;
        perturbed_params(i) = params(i) + delta_params(i);

        % Résoudre les ODE avec les paramètres perturbés
        [~, y_perturbed] = ode45(@(t,y) compostsensanalysis(t, y, perturbed_params), [0 10], [1; 0]);

        % Calculer la variation relative des sorties par rapport au paramètre i
        sensitivity_matrix(:, i) = (y_perturbed(:,1) - y(:,1)) ./ (params(i) * 0.01);
    end

    % Afficher les résultats de l'analyse de sensibilité
    fprintf('Matrice de sensibilité : \n');
    disp(sensitivity_matrix);
end
