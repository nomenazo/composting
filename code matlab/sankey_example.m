% Données d'entrée
source = [1 1 2 3 3];  % Noeuds sources
target = [2 3 3 4 5];  % Noeuds cibles
weight = [10 5 15 5 10];  % Poids (épaisseur des liens)

% Création du diagramme de Sankey
figure;
ssankey(source, target, weight);

