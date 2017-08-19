function yf = y_forward(yi,delta_t,G,B,t_curr,t_next)
   
kv = 0.54;
rho = 2.66*10^-12;
%yf = zeros(1,9);
h = delta_t;
%syms y1 y2 y3 y4 y5 y6 y7 y8 y9;

    dy1 = -3*rho*kv*G*(yi(4)+yi(8));
    dy2 = 0;
    dy3 = G*yi(2);
    dy4 = 2*G*yi(3);
    dy5 = 3*G*yi(4);
    dy6 = B;
    dy7 = G*yi(6);
    dy8 = 2*G*yi(7);
    dy9 = 3*G*yi(8);
    
    
    y1 = yi(1) + delta_t * dy1;
    y2 = yi(2) + delta_t * dy2;
    y3 = yi(3) + delta_t * dy3;
    y4 = yi(4) + delta_t * dy4;
    y5 = yi(5) + delta_t * dy5;
    y6 = yi(6) + delta_t * dy6;
    y7 = yi(7) + delta_t * dy7;
    y8 = yi(8) + delta_t * dy8;
    y9 = yi(9) + delta_t * dy9;    
    
    
    
    options = odeset('RelTol',1e-11,'AbsTol',1e-14);
    [t,y_ode] = ode15s(@(t,y_ode)yODE(t,y_ode,G,B),[t_curr t_next],yi,options);
    
    yf = y_ode(end,:);
    %yf = [y1 y2 y3 y4 y5 y6 y7 y8 y9];

end