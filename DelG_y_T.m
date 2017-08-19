function DelG_dely_delT = DelG_y_T(T,S,y,theta)


Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859; 
DCs_dT = 2.46*10^-3-14.28*10^-6*(T-273);
exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y(1) - Cs))/Cs^2;

A = kg*g*exp(-E_g/T)*S^(g-1)*(E_g/T^2)/Cs;
B = kg*g*exp(-E_g/T)*(g-1)*S^(g-2)*exper/Cs;
C = kg*g*exp(-E_g/T)*S^(g-1)*-1*DCs_dT/Cs^2;

DelG_dely_delT = A+B+C;
end 