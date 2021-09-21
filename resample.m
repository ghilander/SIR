function idx = resample(wk)

    % Approccio alternativo proposto da Lennart Svensson:
    Ns = length(wk);  % Ricavo il numero delle particelle
    
    % somma comulativa dei pesi stando attenti al round-off
    edges = min([0 cumsum(wk)'],1); 
    edges(end) = 1;  
    
    %uniforme u1~U[0,1]
    u1 = rand/Ns;
    
    % Cerco l'intervallo in cui si trova il campione e recupero gli indici
    [~, idx] = histc(u1:1/Ns:1, edges);
    
    return;
end