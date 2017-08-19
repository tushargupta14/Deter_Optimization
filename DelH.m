function DelH_dy = DelH(G,z,y,S,T)
b =1.45;
kb = 285;
E_b = 7517;
kv = 0.54;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859;
rho = 2.66*10^-12;

B = (kb*exp(-E_b/T))*S^b*(y(5)+y(9));       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DelG_dy = DelG(S,T);

DelB_dy= DelB(S,y(4),y(8),T);          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dh_dy1 = DelG_dy*(-z(1)*3*rho*kv*(y(4)+y(8))+z(3)*y(2)+2*z(4)*y(3)+3*z(5)*y(4)+z(6)*B+z(7)*y(6)+2*z(8)*y(7)+3*z(9)*y(8)) + G * (z(6)*DelB_dy);
dh_dy2 = z(3)*G;
dh_dy3 = z(4)*2*G;
dh_dy4 = z(5)*3*G - 3*rho*kv*G*z(1);
dh_dy5 = 0 ;
dh_dy6 = z(7)*G;
dh_dy7 = 2*z(8)*G;
dh_dy8 = z(9)*3*G - 3*rho*kv*G*z(1);
dh_dy9 = 0;

DelH_dy = [dh_dy1 dh_dy2 dh_dy3 dh_dy4 dh_dy5 dh_dy6 dh_dy7 dh_dy8 dh_dy9];


end