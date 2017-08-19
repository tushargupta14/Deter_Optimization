function dtheta_dt = thetaODE(t,thi,G,DelG_dT,DelB_dT,y)

   rho = 2.66*10^-12;
   kv= 0.54; 
   
   dthe1 = -3*rho*kv*(thi(4)+thi(8))*G -3*rho*kv*(y(4)+y(8))*DelG_dT;
   dthe2 = 0;
   dthe3 = thi(2)*G+y(2)*DelG_dT;
   dthe4 = 2*y(3)*DelG_dT + 2*thi(3)*G;
   dthe5 = 3*DelG_dT*y(4)+3*thi(4)*G;
   dthe6 = DelB_dT;
   dthe7 = thi(6)*G+DelG_dT*y(6);
   dthe8 = 2*G*thi(7)+ 2*DelG_dT*y(7);
   dthe9 = 3*G*thi(8) + 3*DelG_dT*y(8);
   
   dtheta_dt = [dthe1;dthe2;dthe3;dthe4;dthe5;dthe6;dthe7;dthe8;dthe9];

end
