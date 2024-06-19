
function [t, y] = ode45_with_corrections(~, tspan, y0)
    t = tspan(1);
    y = y0(:).'; % Assurez-vous que y0 est un vecteur ligne
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    
    while t(end) < tspan(end)
        % Intégrer un petit pas de temps
        [t_step, y_step] = ode45('compostfitted', [t(end) min(t(end) + 0.1, tspan(end))], y(end, :), options);
        
        % Appliquer la correction pour forcer les valeurs négatives à zéro
        y_step(y_step < 0) = 0;
        
        % Accumuler les résultats
        t = [t; t_step(2:end)];
        y = [y; y_step(2:end, :)];
    end
end