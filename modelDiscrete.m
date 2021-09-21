function xk = modelDiscrete(xkm1,T,param)
    
    b = param(1); %infection rate
    k = param(2); %recover rate
    
    S = xkm1(1); %soggetti infettabili (suscettibili)
    I = xkm1(2); %soggetti infetti
    R = xkm1(3); %soggetti guariti (recuperati)
    
    %Eulero in avanti per discretizzare  (scelta implementativa)
    xk(1,1) = S + T *( -b*S*I );
    xk(2,1) = I + T * ( b*S*I - k*I );
    xk(3,1) = R + T * ( k*I );
   
end