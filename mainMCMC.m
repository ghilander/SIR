close all;
clc;
clear all;

% attiva debug:
debug = 1;

%% Setup dei parametri
if(~debug)
    b = input('Si inserisca un valore per il parametro b, i.e. rate di contagio (sugg. 0.0005): ');
    k = input('Si inserisca un valore per il parametro k, i.e. rate di guarigione (sugg. 0.01): ');
else
    b = 0.0005; %rate di contagio (1/giorni) (quante persone si ammalano al giorno)
    k = 0.01; %rate di guarigione (1/giorni) (quante persone guariscono al giorno)
end
oracle = [b,k]'; %costruisco vettore dei parametri

%% Setup dell'MCMC 
%numero di iterazioni
if(~debug)
    nIter = input("Si inserisca il numero di iterazioni dell'MCMC: ");
else
    nIter = 1000; %numero di iterazioni dell'mcmc
end

%distribuzioni a priori 
varPrior = 3;
priorB = @(b) normpdf(b,0,varPrior); 
priorK = @(k) normpdf(k,0,varPrior);

%proposta iniziale di theta
theta0 = [0,0]';

%deviazione standard della proposal
stdB = 1e-5;
stdK = 5e-4;
stdProp = [stdB, stdK]';

%% Setup dei parametri fissati e condizioni iniziali
% Parametri fissati
nSoggetti = 100; %numero di soggetti totali
epsilon = 0.005; %frazione dei soggetti inizialmente infetti

%Condizioni iniziali
S0 = nSoggetti*(1-epsilon); %99.5 % dei soggetti è inizialmente sano
I0 = nSoggetti*epsilon; %0.5 % dei soggetti è inizialmente sano
R0 = 0; %all'inizio nessuno è già guarito
X0 = [S0, I0, R0]'; %costruisco vettore delle condizioni iniziali

%% Simulazione del modello con parametri
tspan = 0:365; %simulo un totale di 365 giorni
N = length(tspan); %numero di misure
[tOracle, XOracle] = ode45('model',tspan,X0,'',oracle); %simulo il modello con i parametri scelti

%% Visualizzazione dei risultati
%Plotto le 3 equazioni di stato
figure;
hold on;
plot(tOracle,XOracle(:,1));
plot(tOracle,XOracle(:,2),'r');
plot(tOracle,XOracle(:,3),'g');
legend S I R

%% Generazione del vettore delle misure 
I = XOracle(:,2);
SDw = 5;
y = I + randn(length(XOracle(:,2)),1)*SDw; %y(t) = I(t) + w(t)  w~N(0,4)
t = tOracle;

%% Metropolis-Hastings 
bHat = zeros(nIter,1);
kHat = zeros(nIter,1);

X = theta0; %setto X alla proposta iniziale
accept = 0; %conto il numero di volte che accetto
figure;
if(debug)
    subplot(211)
    plot([1 nIter],[oracle(1) oracle(1)],'m--');
    ylabel('bHat')
    subplot(212)
    plot([1 nIter],[oracle(2) oracle(2)],'m--');
    ylabel('kHat')
end

for iter = 1:nIter
    
    %Blocco X
    [tSim, XSim] = ode45('model',tspan,X0,'',X); %simulo il modello con i parametri X
    
    %Calcolo di pi(X) (nel logaritmo per efficienza numerica)
    lX = -(N/2)*log(2*pi)-(N/2)*log((SDw^2))-0.5*sum(((XSim(:,2)-y)/SDw).^2); %calcolo log-likelihood
    piX = lX + log(priorB(X(1))) + log(priorK(X(2))); %calcolo pi(X) 
    
    %Blocco Y
    Y = X + randn(length(X),1).*stdProp;
    [tSim, XSim] = ode45('model',tspan,X0,'',Y); %simulo il modello con i parametri Y
    
    %Calcolo di pi(Y) (nel logaritmo per efficienza numerica)
    lY = -(N/2)*log(2*pi)-(N/2)*log((SDw^2))-0.5*sum(((XSim(:,2)-y)/SDw).^2); %calcolo log-likelihood
    piY = lY + log(priorB(Y(1))) + log(priorK(Y(2))); %calcolo pi(Y) 
    
    % Decido se accetto o rifiuto il nuovo candidato
    U = rand(1);
    alfa = min(1,exp(piY-piX));
    if(U<=alfa && ~isnan(exp(piY-piX)))
        X = Y;
        accept = accept + 1;
    end %if
    
    %Salvo la catena
    bHat(iter) = X(1);
    kHat(iter) = X(2);
    
    %Visualizzo come procede (SOLO DEBUG)
    if(debug)
        subplot(211)
        title(['Iterazione ' num2str(iter) ' di ' num2str(nIter)]);
        hold on
        scatter(iter,bHat(iter),'MarkerEdgeColor',[0 0 0]);
        subplot(212)
        hold on
        scatter(iter,kHat(iter),'MarkerEdgeColor',[0 0 0]);
        pause(0.001)
    end
    % ---------------------------------------------------------------------
    
end

%% Generazione delle figure richieste
%Plot I vs. y
figure;
hold on
plot(tOracle,I)
plot(t,y)
legend I y

%Plot istogrammi delle marginali
figure;
subplot(121)
hist(bHat,25)
title('bHat')
subplot(122)
hist(kHat,25)
title('kHat')
