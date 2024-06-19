function [value, isterminal, direction] = nonnegativeEvents(t, y)
    % Définir les conditions des événements
    value = y;  % Nous voulons vérifier si y devient négatif
    isterminal = zeros(size(y));  % Continuer l'intégration même si une condition est atteinte
    direction = zeros(size(y));  % Détecter les zéros dans toutes les directions
end