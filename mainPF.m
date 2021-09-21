close all;
clc;
clear all;

% attiva debug:
debug = 1;

%% Setup dei parametri del modello
b = 0.0005; %rate di contagio (1/giorni) (quante persone si ammalano al giorno)
k = 0.01; %rate di guarigione (1/giorni) (quante persone guariscono al giorno)

theta = [b,k]'; %metto tutti i parametri in un vettore

%% Setup dei parametri fissati e condizioni iniziali

%Passo di campionamento del modello discreto
T = 1;

%Condizioni iniziali
S0 = 99; % 99% dei soggetti è inizialmente sano
I0 = 1; % 1% dei soggetti è inizialmente sano
R0 = 0; %all'inizio nessuno è già guarito
X0 = [S0, I0, R0]'; %costruisco vettore delle condizioni iniziali

%% Generazione del vettore delle misure
tspan = 0:365; %simulo per 1 anno
N = length(tspan); %numero di misure

X = zeros(length(tspan),3);
X(1,:) = X0;
for t = 2:length(tspan)
    X(t,:) = modelDiscrete(X(t-1,:),T,theta);
end
I = X(:,2);
SDw = 4; %errore di misura =  w~N(0,4) 
y = I + randn(length(X(:,2)),1)*SDw; %y(t ) = I(t) + w(t)

%% Visualizzazione dei risultati
%Plotto le 3 equazioni di stato
figure;
hold on;
plot(tspan,X(:,1));
plot(tspan,X(:,2),'r');
plot(tspan,X(:,3),'g');
plot(tspan,y,'b');
legend S I R y

%% Setup dei parametri del PF
if(debug)
    NPart = 1000;
else
    NPart = input('Inserire il numero di particelle: (consigliato 1000)');
end
if(debug)
    muS = 90;
    muI = 5;
    muR = 5;
else
    muS = input('Inserire mu di S (valore vero 99, consigliato 90): ');
    muI = input('Inserire mu di I (valore vero 1, consigliato 5): ');
    muR = input('Inserire mu di R (valore vero 0, consigliato 5): ');
end
if(debug)
    stdS = 0.5;
    stdI = 1;
    stdR = 0.5;
else
    stdS = input('Inserire sigma di S (valore consigliato 0.5): ');
    stdI = input('Inserire sigma di I (valore consigliato 1): ');
    stdR = input('Inserire sigma di R (valore consigliato 0.5): ');
end

%% Particle filter
particelle = zeros(length(tspan),3,NPart);
pesi = zeros(length(tspan),NPart);
Xhat = zeros(length(tspan),3);

%Setto le condizioni iniziali
Xhat(1,:) = [muS, muI, muR]';

%Std del rumore di modello
stdModello = [stdS, stdI, stdR]';

%Pesco le prime particelle
particelle(1,1,:) = S0 + randn(NPart,1)*stdModello(1);
particelle(1,2,:) = I0 + randn(NPart,1)*stdModello(2);
particelle(1,3,:) = R0 + randn(NPart,1)*stdModello(3);

for t = 2:length(tspan)

    for p = 1:NPart
        
        %Pesco una nuova particella...
        particelle(t,:,p) = modelDiscrete(particelle(t-1,:,p),T,theta);
        
        particelle(t,1,p) = particelle(t,1,p) + randn(1)*stdModello(1); %S
        particelle(t,2,p) = particelle(t,2,p) + randn(1)*stdModello(2); %I
        particelle(t,3,p) = particelle(t,3,p) + randn(1)*stdModello(3); %R
        
        %...e ne calcolo il peso
        pesi(t,p) = normpdf(y(t),particelle(t,2,p),SDw);
        
    end
    
    %Normalizzo i pesi 
    pesi(t,:) = pesi(t,:)./sum(pesi(t,:));
    
    %Stimo lo stato
    Xhat(t,1) = squeeze(particelle(t,1,:))'*pesi(t,:)'/sum(pesi(t,:));
    Xhat(t,2) = squeeze(particelle(t,2,:))'*pesi(t,:)'/sum(pesi(t,:));
    Xhat(t,3) = squeeze(particelle(t,3,:))'*pesi(t,:)'/sum(pesi(t,:));
    
    %Faccio resampling
    idx = resample(pesi(t,:)');
    particelle(t,1,:) = particelle(t,1,idx);
    particelle(t,2,:) = particelle(t,2,idx);
    particelle(t,3,:) = particelle(t,3,idx);
    
end

%% Generazione delle figure richieste
%Plot I vs. y
figure;
subplot(131)
hold on
plot(tspan,Xhat(:,1))
plot(tspan,X(:,1))
subplot(132)
hold on
plot(tspan,Xhat(:,2))
plot(tspan,X(:,2))
subplot(133)
hold on
plot(tspan,Xhat(:,3))
plot(tspan,X(:,3))
legend Xhat1 X1

figure;
plot(y)
hold on 
plot(tspan,Xhat(:,2))

legend y yhat

