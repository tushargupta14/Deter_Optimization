function zf = z_backward(zi,yi,delta_t,G,T,S,t_curr,t_next)
   
    rho = 2.66*10^-12;
    kv = 0.54;
    b =1.45;
    kb = 285;
    E_b = 7517;
    %kv = 1.5;
   
    dz1 = zi(1)*(3*rho*kv*(yi(4)+yi(8))*DelG(S,T))-(zi(3)*yi(2)+2*zi(4)*yi(3)+3*zi(5)*yi(4)+zi(7)*yi(6)+2*zi(8)*yi(7)+3*zi(9)*yi(8))*DelG(S,T)-DelB(S,yi(4),yi(8),T)*zi(6);      %%%%%%%%%%%%%%%%%%%%%%%%%%
    dz2 = -zi(3)*G;
    dz3 = -2*zi(4)*G;
    dz4 =  3*zi(1)*rho*kv*G -3*zi(5)*G;
    dz5 = -zi(6)*kb*exp(-E_b/T)*S^b;
    dz6 = -zi(7)*G;
    dz7 = -2*G*zi(8);
    dz8 =  3*zi(1)*rho*kv*G-3*zi(9)*G;
    dz9 = -zi(6)*kb*exp(-E_b/T)*S^b;
    
    z1 = zi(1) - delta_t * dz1;
    z2 = zi(2) - delta_t * dz2;
    z3 = zi(3) - delta_t * dz3;
    z4 = zi(4) - delta_t * dz4;
    z5 = zi(5) - delta_t * dz5;
    z6 = zi(6) - delta_t * dz6;
    z7 = zi(7) - delta_t * dz7;
    z8 = zi(8) - delta_t * dz8;
    z9 = zi(9) - delta_t * dz9;
    
    
    tspan = [t_curr t_next];
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    [t,z] = ode15s(@(t,z)zODE(t,z,yi,G,T,S),[t_curr t_next],zi,options);
    
    zf = z(end,:);
    %zf = [z1 z2 z3 z4 z5 z6 z7 z8 z9];
end