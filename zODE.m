function dz_dt = zODE(t,z,y,G,T,S)
    
    b =1.45;
kb = 285;
E_b = 7517;
kv = 0.54;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859;
rho = 2.66*10^-12;
   
    dz1 = z(1)*(3*rho*kv*(y(4)+y(8))*DelG(S,T))-(z(3)*y(2)+2*z(4)*y(3)+3*z(5)*y(4)+z(7)*y(6)+2*z(8)*y(7)+3*z(9)*y(8))*DelG(S,T)-DelB(S,y(4),y(8),T)*z(6);      %%%%%%%%%%%%%%%%%%%%%%%%%%
    dz2 = -z(3)*G;
    dz3 = -2*z(4)*G;
    dz4 =  3*z(1)*rho*kv*G -3*z(5)*G;
    dz5 = -z(6)*kb*exp(-E_b/T)*S^b;
    dz6 = -z(7)*G;
    dz7 = -2*G*z(8);
    dz8 =  3*z(1)*rho*kv*G-3*z(9)*G;
    dz9 = -z(6)*kb*exp(-E_b/T)*S^b;

    dz_dt = [dz1;dz2;dz3;dz4;dz5;dz6;dz7;dz8;dz9];

end
