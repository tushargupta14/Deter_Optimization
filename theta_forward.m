function thf = theta_forward(thi,delta_t,G,DelG_dT,DelB_dT,y,t_curr,t_next)
    
   rho = 2.66*10^-12;
   kv= 0.54; 
   dthe1 = -3*rho*kv*(thi(4)+thi(8))*G -3*rho*kv*(y(4)+y(8))*DelG_dT;
   dthe2 = 0;
   dthe3 = thi(2)*G+y(2)*DelG_dT;
   dthe4 = 2*y(3)*DelG_dT + 2*thi(3)*G;
   dthe5 = 3*DelG_dT*y(4)+3*thi(4)*G;
   dthe6 = DelB_dT;
   dthe7 = thi(6)*G + DelG_dT*y(6);
   dthe8 = 2*G*thi(7)+ 2*DelG_dT*y(7);
   dthe9 = 3*G*thi(8) + 3*DelG_dT*y(8);
   h = delta_t;
   thf = zeros(1,9);
  
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    [t,theta_ode] = ode15s(@(t,theta_ode)thetaODE(t,theta_ode,G,DelG_dT,DelB_dT,y),[t_curr t_next],thi,options);
        
    thf = theta_ode(end,:);
end