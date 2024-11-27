%data = [randi(3,[20 1]), randi(5,[20 1]), randi(2,[20 1])]
data = [0.5 1; 1.5; 0.2 1.3]



% Options pour personnaliser le diagramme de Sankey
options.color_map = 'parula';         % Colormap
options.flow_transparency = 0.2;      % Opacité des chemins de flux
options.bar_width = 40;               % Largeur des blocs de catégories
options.show_perc = true;             % Afficher les pourcentages sur les blocs
options.text_color = [1 1 1];         % Couleur du texte pour les pourcentages
options.show_layer_labels = true;     % Afficher les noms des couches sous le graphique
options.show_cat_labels = true;       % Afficher les étiquettes des catégories sur les blocs
options.show_legend = true;           % Afficher la légende avec les noms des catégories

% Tracer le diagramme de Sankey
figure;
%plotSankeyFlowChart(data,options)
plotSankeyFlowChart(data, options);
