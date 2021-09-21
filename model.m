function dX = model(t,X,~,param)
    
    b = param(1); %infection rate
    k = param(2); %recover rate
    
    S = X(1); %soggetti infettabili (suscettibili)
    I = X(2); %soggetti infetti
    R = X(3); %soggetti guariti (recuperati)
    
    dX(1,1) = -b*S*I;
    dX(2,1) = b*S*I - k*I;
    dX(3,1) = k*I;

end