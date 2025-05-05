% Simulatore di orbite con legge di gravità F = G*m1*m2/r^n
% Studio della possibilità di vita in universi con diversa legge di gravità

function main()
    close all; 
    clc; 

    % definisco le costanti fisiche (in unità delsistema internazionale)
    costG = 6.67430e-11; 
    massaSole = 1.989e30; 
    massaTerra = 5.972e24; 
    unitaAU = 1.496e11; 
    
    % definisco i limiti minimi e massimi per n trovati numericamente a
    % parte
    nMinimo = 1.904;
    nMassimo = 2.76393328;
    
    % prendo l'nput dell'utente e ne controllo la validità
    n = input(['Inserisci il valore dell''esponente n per la legge di gravità F = G*m1*m2/r^n\n' ...
              'Il valore deve essere compreso tra ' num2str(nMinimo) ' e ' num2str(nMassimo) ':\n']);
    
    while n < nMinimo || n > nMassimo
        fprintf('Valore di n non valido. Deve essere compreso tra %.6f e %.6f\n', nMinimo, nMassimo);
        n = input('Inserisci un nuovo valore di n: ');
    end
    
    % stabilisco dei parametri simulazione
    massa1 = massaSole;      
    massa2 = massaTerra;    
    
    % calcolo parametri orbitali
    [momentoAngolare, velocitaOrbitale] = calcolaMomentoAngolareCompatibile(costG, massa1, massa2, n, unitaAU);
    [raggioMinimo, energiaMinima] = disegnaDiagrammaEnergia(costG, massa1, massa2, n, momentoAngolare);
    raggioOrbita = raggioMinimo;  
    energia = energiaMinima;  
    [perielio, afelio] = trovaEstremiOrbitali(costG, massa1, massa2, momentoAngolare, energia, n);
    periodo = calcolaPeriodoOrbitale(costG, massa1, massa2, perielio, afelio, n);
    
    % visualizzazione risultati
    mostraRisultati(costG, massa1, massa2, momentoAngolare, raggioOrbita, energia, perielio, afelio, periodo, n, velocitaOrbitale);
    
    % simulazione ed animazione dell'orbita
    simulaOrbita(costG, massa1, massa2, momentoAngolare, raggioOrbita, energia, n);
end

function [momentoAngolare, velocitaOrbitale] = calcolaMomentoAngolareCompatibile(costG, massa1, massa2, n, raggioDesiderato)
    % funzione che calcola un momento angolare compatibile con un'orbita circolare
    % ricordimo che per un'orbita circolare, la forza centrifuga deve bilanciare la forza gravitazionale
    
    % formula derivata dall'equilibrio tra forza gravitazionale e centrifuga:
    % F_g = F_c => G*M*m/r^n = m*v^2/r => v^2 = G*M/r^(n-1)
    momentoAngolare = massa2 * sqrt(costG * massa1 * raggioDesiderato^(3-n));
    
    % velocità orbitale necessaria per orbita circolare con dato n
    velocitaOrbitale = sqrt(costG * massa1 / raggioDesiderato^(n-1));
   
end

function [raggioMinimo, energiaMinima] = disegnaDiagrammaEnergia(costG, massa1, massa2, n, momentoAngolare)
    unitaAU = 1.496e11; 
    intervalloDiRaggio = linspace(0.1*unitaAU, 2*unitaAU, 1000);
    
    % calcolo il potenziale efficace:
    % U_eff = potenziale centrifugo + potenziale gravitazionale
    if n ~= 1
        potenzEffettivo = (momentoAngolare^2)./(2*massa2*(intervalloDiRaggio.^2)) + (costG*massa2*massa1)./(1-n)*(intervalloDiRaggio.^(1-n));
    else
        % caso speciale per n=1, l'andamento passa a logaritmico
        potenzEffettivo = (momentoAngolare^2)./(2*massa2*(intervalloDiRaggio.^2)) - costG*massa2*massa1*log(intervalloDiRaggio);
    end
    
    % plot energia
    figure('Name', ['Diagramma Energia-Distanza (n = ' num2str(n) ')'], 'Position', [100, 100, 800, 600]);
    plot(intervalloDiRaggio/unitaAU, potenzEffettivo, 'LineWidth', 2);
    [energiaMinima, indice] = min(potenzEffettivo);
    raggioMinimo = intervalloDiRaggio(indice);
    
    hold on;
    plot(raggioMinimo/unitaAU, energiaMinima, 'ro', 'MarkerSize', 10);
    plot(intervalloDiRaggio/unitaAU, ones(size(intervalloDiRaggio))*energiaMinima, 'g--');
    
    xlabel('Distanza [AU]');
    ylabel('Energia Potenziale Efficace [J]');
    title(['Energia potenziale efficace con n = ' num2str(n)]);
    grid on;
end

function simulaOrbita(costG, massa1, massa2, momentoAngolare, raggioOrbita, energia, n)
    % verifico la stabilità secondo il teorema di Bertrand
    orbitaStabile = isOrbitaBertrandStabile(n);
    
    % definisco alcuni parametri
    unitaAU = 1.496e11; 
    secondiAnno = 1/(365.25*86400);
    lunghezzaScia = 500;
    
    % calcoo i limiti assi iniziali
    [perielio, afelio] = trovaEstremiOrbitali(costG, massa1, massa2, momentoAngolare, energia, n);
    raggioMassimoIniziale = max([perielio, afelio])/unitaAU;
    limiteAsseIniziale = 1.2 * raggioMassimoIniziale;
    
    % inizializzazo figura
    fig = figure('Name', ['Simulazione Orbita con n = ' num2str(n)], 'Position', [100, 100, 800, 800]);
    hold on;
    grid on;
    axis equal;
    axis([-limiteAsseIniziale limiteAsseIniziale -limiteAsseIniziale limiteAsseIniziale]);
    
    % elementi grafici
    hStella = plot(0, 0, 'yo', 'MarkerSize', 20, 'MarkerFaceColor', 'yellow');
    hPianeta = plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'cyan');
    hScia = plot([], [], 'b-', 'LineWidth', 1.5);
    
    % testo 
    hTesto = text(0.05, 0.95, '', 'Units', 'normalized', 'FontSize', 12);
    hStabilita = text(0.5, 0.05, '', 'Units', 'normalized',...
                      'HorizontalAlignment', 'center', 'FontSize', 14,...
                      'FontWeight', 'bold', 'Color', 'red');

    % genero traiettoria in base alla stabilità
    if orbitaStabile
        [x, y, t, raggi] = generaOrbitaStabile(costG, massa1, raggioOrbita, n, 5000);
        set(hStabilita, 'String', 'ORBITA STABILE (Teorema di Bertrand)', 'Color', 'green');
    else
        [x, y, t, raggi] = generaOrbitaInstabile(costG, massa1, massa2, momentoAngolare, energia, n, raggioOrbita, 5);
        set(hStabilita, 'String', 'ORBITA INSTABILE (Teorema di Bertrand)', 'Color', 'red');
    end

    % creo animazione
    raggioMassimoDinamico = raggioMassimoIniziale;
    datiScia = [];
    for i = 1:length(t)
        if ~isvalid(fig), break; end % Interrompi se la figura è chiusa
        
        % aggiorno posizione pianeta
        posizioneAttuale = [x(i)/unitaAU, y(i)/unitaAU];
        set(hPianeta, 'XData', posizioneAttuale(1), 'YData', posizioneAttuale(2));
        
        % aggiorna scia <--- sta parte ancora non riesco a farla andare
        datiScia = [datiScia; posizioneAttuale];
        if size(datiScia,1) > lunghezzaScia
            datiScia = datiScia(end-lunghezzaScia:end,:);
        end
        set(hScia, 'XData', datiScia(:,1), 'YData', datiScia(:,2));
        
        % aggiorno lo zoom in modo dinamico se il pianeta si allontana
        raggioAttuale = norm(posizioneAttuale);
        if raggioAttuale > raggioMassimoDinamico
            raggioMassimoDinamico = raggioAttuale;
            nuovoLimite = 1.2 * raggioMassimoDinamico;
            axis([-nuovoLimite nuovoLimite -nuovoLimite nuovoLimite]);
        end
        
        % aggiorna informazioni
        tempoTrascorso = t(i) * secondiAnno;
        testoStato = sprintf(['Tempo: %.2f anni\n'...
                              'Distanza: %.3f AU\n'...
                              'Velocità: %.2f km/s'],...
                              tempoTrascorso, raggioAttuale,...
                              sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2)/(t(i+1)-t(i))/1000);
        set(hTesto, 'String', testoStato);
        
        drawnow;
        pause(0.001);
    end
end

function [x, y, t, raggi] = generaOrbitaStabile(costG, massa1, raggioOrbita, n, numPunti)
    % calcola il periodo orbitale in modo accurato in base alla legge gravitazionale
    if n == 2
        periodo = 2*pi*sqrt(raggioOrbita^3/(costG*massa1)); % legge di Newton (n=2)
    elseif n == 3
        periodo = 2*pi*sqrt(raggioOrbita^3/(costG*massa1)); % oscillatore armonico (n=3) (mai possibile ma questa parte l'avevo già iniziata a fare da prima che trovassimo dei limiti di n)
    else 
        periodo = 2*pi*sqrt(raggioOrbita^n/(costG*massa1)); % approssimazione generica
    end
    
    theta = linspace(0, 4*pi, numPunti); % 2 periodi completi
    if n == 2
        e = 0.0167; % impongo un eccentricità simile a Terra-Sole
        r = raggioOrbita*(1 - e^2)./(1 + e*cos(theta));
    else
        r = raggioOrbita*ones(size(theta)); % orbita circolare
    end
    
    x = r.*cos(theta);
    y = r.*sin(theta);
    t = linspace(0, 2*periodo, numPunti);
    raggi = sqrt(x.^2 + y.^2);
end

function [x, y, t, raggi] = generaOrbitaInstabile(costG, massa1, massa2, momentoAngolare, energia, n, raggioOrbita, numPeriodi)
    % per le orbite instabili, risolviamo numericamente le equazioni del moto
    numPeriodi = 10; % aumentato per vedere meglio l'instabilità
    periodoApprossimato = 2*pi*sqrt(raggioOrbita^n/(costG*massa1));
    
    % configura opzioni per l'integratore numerico
    opzioni = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,...
                    'Events', @(t,y) eventi(t,y,raggioOrbita));
    
    % risolvo numericamente le equazioni del moto passando alle coordinate polari
    [t,y] = ode45(@(t,y) odeFunzione(t,y,costG,massa1,massa2,n), [0 numPeriodi*periodoApprossimato],...
                 [raggioOrbita; 0; 0; momentoAngolare/(massa2*raggioOrbita)], opzioni);
    
    % converto le coordinate da polari a cartesiane
    x = y(:,1).*cos(y(:,2));
    y = y(:,1).*sin(y(:,2));
    t = t(1:end-1);
    raggi = sqrt(x.^2 + y.^2);
    
    % campiono solo per rendere l'animazione più fluida
    if length(t) > 5000
        idx = round(linspace(1, length(t), 5000));
        x = x(idx);
        y = y(idx);
        t = t(idx);
    end
end

function dydt = odeFunzione(t, y, costG, massa1, massa2, n)
    % equazioni del moto in coordinate polari
    % y(1) = r, y(2) = theta, y(3) = velocità radiale, y(4) = velocità angolare
    r = y(1);
    theta = y(2);
    vr = y(3);
    omega = y(4);
    
    % prevengo possibili divisioni per zero (parte scritta dopo 5 ore di studio filato potrei aver scritto delle assurdità)
    if r <= 0
        r = eps;
    end
    
    dydt = zeros(4,1);
    dydt(1) = vr;                                  % dr/dt = v_r
    dydt(2) = omega;                               % dθ/dt = omega
    dydt(3) = r*omega^2 - costG*massa1/r^(n-1);    % dvr/dt = r·ω² - G·M/r^(n-1)
    dydt(4) = -2*vr*omega/r;                       % dω/dt = -2·v_r·ω/r (conservazione momento angolare)
end

function [valore, terminale, direzione] = eventi(t, y, raggioIniziale)
    % funzione di controllo per interrompere l'integrazione nei seguenti casi:
    % - il pianeta si allontana troppo (10 volte il raggio iniziale)
    % - il pianeta si avvicina troppo al centro (1% del raggio iniziale)
    valore = [y(1) - 10*raggioIniziale;    % criterio raggio massimo
             y(1) - 0.01*raggioIniziale];  % criterio raggio minimo
    terminale = [1; 1];                    % interrompo integrazione quando raggiunti
    direzione = [1; -1];                   % direzione dell'attraversamento
end

function stabile = isOrbitaBertrandStabile(n)
    % in reatlà il teorema di Bertrand stabilisce che solo per n=2 (legge inversa del quadrato)
    % e n=3 (oscillatore armonico) le orbite chiuse sono stabili per piccole perturbazioni
    stabile = (abs(n - 2) < 1e-6) || (abs(n - 3) < 1e-6);
end

function [perielio, afelio] = trovaEstremiOrbitali(costG, massa1, massa2, momentoAngolare, energia, n)
    unitaAU = 1.496e11; 
    
    % funzione che calcola E - U_eff. I punti di intersezione con zero
    % rappresentano i raggi in cui la velocità radiale è nulla (punti di inversione)
    if n ~= 1
        f = @(r) energia - (momentoAngolare^2./(2*massa2*(r.^2)) + (costG*massa2*massa1)./(1-n)*(r.^(1-n)));
    else
        f = @(r) energia - (momentoAngolare^2./(2*massa2*(r.^2)) - costG*massa2*massa1*log(r));
    end
    
    intervalloRaggio = linspace(0.01*unitaAU, 10*unitaAU, 10000);
    valoriF = arrayfun(f, intervalloRaggio);
    cambiSegno = find(diff(sign(valoriF)) ~= 0);
    
    perielio = NaN;
    afelio = NaN;
    
    try
        if length(cambiSegno) < 2
            % se non troviamo due intersezioni evidenti, proviamo diversi punti iniziali
            raggioMinPotenziale = (momentoAngolare^2 / (costG*massa1*massa2))^(1/n);
            puntiProva = [0.1*unitaAU, 0.5*unitaAU, raggioMinPotenziale, 1.5*unitaAU, 5*unitaAU];
            
            for i = 1:length(puntiProva)
                try
                    rProva = fzero(f, puntiProva(i), optimset('Display', 'off'));
                    if rProva > 0
                        if isnan(perielio) || rProva < perielio
                            perielio = rProva;
                        end
                        if isnan(afelio) || rProva > afelio
                            afelio = rProva;
                        end
                    end
                catch
                    % ignora errori di convergenza
                end
            end
        else
            %  in questo caso abbiamo trovato almeno due cambi di segno (due intersezioni)
            perielio = fzero(f, intervalloRaggio(cambiSegno(1)));
            afelio = fzero(f, intervalloRaggio(cambiSegno(end)));
        end
        
        % controllo che perielio < afelio
        if perielio > afelio
            [perielio, afelio] = deal(afelio, perielio);
        end
    catch
        % non c'è andata bene, usiamo valori predefiniti
        perielio = unitaAU;
        afelio = unitaAU;
    end
end

function periodo = calcolaPeriodoOrbitale(costG, massa1, massa2, perielio, afelio, n)
    unitaAU = 1.496e11; 
    
    % Se non abbiamo trovato perielio e afelio validi
    if isnan(perielio) || isnan(afelio)
        periodo = 2 * pi * sqrt(unitaAU^n / (costG*massa1));
        return;
    end
    
    % aprrossiamiamo il periodo usando il raggio medio
    raggioMedio = (perielio + afelio)/2;
    
    %  generalizzazione della terza legge di Keplero: T² ∝ R^n
    periodo = 2 * pi * sqrt(raggioMedio^n / (costG*massa1));
end

function mostraRisultati(costG, massa1, massa2, momentoAngolare, raggioOrbita, energia, perielio, afelio, periodo, n, velocitaOrbitale)
    unitaAU = 1.496e11; 
    anno = 31557600; 
    
    fprintf('\n=== RISULTATI PER n = %.6f ===\n', n);
    fprintf('Momento angolare: %.3e kg·m²/s\n', momentoAngolare);
    fprintf('Velocità orbitale: %.2f km/s\n', velocitaOrbitale/1000);
    fprintf('Raggio orbita: %.3f AU\n', raggioOrbita/unitaAU);
    fprintf('Perielio: %.3f AU\n', perielio/unitaAU);
    fprintf('Afelio: %.3f AU\n', afelio/unitaAU);
    fprintf('Periodo orbitale: %.2f anni\n', periodo/anno);
    fprintf('Rapporto T²/r^n: %.3e\n', periodo^2/raggioOrbita^n);
    fprintf('4π²/(GM): %.3e\n', 4*pi^2/(costG*massa1));
    
    % verifichiamo la stabilità secondo Bertrand
    if abs(n - 2) < 1e-6 || abs(n - 3) < 1e-6
        fprintf('L''orbita è STABILE secondo il teorema di Bertrand.\n');
    else
        fprintf('L''orbita NON è stabile secondo il teorema di Bertrand.\n');
    end
end