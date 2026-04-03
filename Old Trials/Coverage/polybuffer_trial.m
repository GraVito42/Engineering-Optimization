% 1. Pulizia e conversione esplicita
R = double(R); 

% 2. Crea il polyshape in modo esplicito e semplice
% Rimuoviamo eventuali nan o duplicati prima
g_clean = [g(1, ~isnan(g(1,:))); g(2, ~isnan(g(2,:)))];
ps = polyshape(g_clean(1,:), g_clean(2,:), 'Simplify', true);

% 3. Esegui il buffer con controllo del tipo
% Se continua a dare errore, prova a definire il buffer senza 'lines' 
% (a volte polybuffer preferisce lavorare su poligoni chiusi)
buffer_ps = polybuffer(ps, 'lines', R);

% 4. BUFFER: Usiamo il parametro 'lines' in modo esplicito
% Se l'errore persiste, è probabile che ci siano punti sovrapposti
% Puliamo il polyshape prima del buffer
ps = rmmissing(ps); 
buffer_ps = polybuffer(ps, 'lines', R);

% 5. Plot
figure; hold on; axis equal; grid on;
plot(buffer_ps, 'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(g(1,:), g(2,:), 'r', 'LineWidth', 2);
title('Inviluppo tramite Polybuffer');

% 6. Area
area_totale = area(buffer_ps);
fprintf('L''area reale dell''unione dei cerchi è: %.4f\n', area_totale);