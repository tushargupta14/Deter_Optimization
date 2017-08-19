function dy_dt = yODE(t,y,G,B)
    rho = 2.66*10^-12;
    kv = 0.54;
    
    dy1 = -3*rho*kv*G*(y(4)+y(8));
    dy2 = 0;
    dy3 = G*y(2);
    dy4 = 2*G*y(3);
    dy5 = 3*G*y(4);
    dy6 = B;
    dy7 = G*y(6);
    dy8 = 2*G*y(7);
    dy9 = 3*G*y(8);
    
    dy_dt = [dy1;dy2;dy3;dy4;dy5;dy6;dy7;dy8;dy9];

    
end