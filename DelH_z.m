function DelH_dz = DelH_z(G,B,rho,kv,y)
    
dh_dz1 = -3*rho*kv*G*(y(4)+y(8));
dh_dz2 = 0;
dh_dz3 = G*y(2);
dh_dz4 = 2*G*y(3);
dh_dz5 = 3*G*y(4) ;
dh_dz6 = B;
dh_dz7 = y(6)*G;
dh_dz8 = 2*G*y(7);
dh_dz9 = 3*G*y(8);

DelH_dz = [dh_dz1 dh_dz2 dh_dz3 dh_dz4 dh_dz5 dh_dz6 dh_dz7 dh_dz8 dh_dz9]; 

end