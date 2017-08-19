function DelG_dT = DelG_T(S,theta,T,y1)

E_g = 4859;
g = 1.5;
kg = 1.44 *10^8;
Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
DCs_dT = 2.46*10^-3 -14.28*10^-6*(T-273);
exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y1 - Cs))/Cs^2;
A = kg*S^g*exp(-E_g/T)*E_g/T^2;
DelG_dT = A + kg*exp(-E_g/T)*g*(S^(g-1))*exper;

% derivative using the product rule 

end