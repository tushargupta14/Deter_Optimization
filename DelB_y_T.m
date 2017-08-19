function DelB_dely_delT = DelB_y_T(T,S,theta,y)
% All values calculated at the initial concentration 
kb = 285;
b = 1.45;
E_b = 7517;
Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
%Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
DCs_dT = 2.46*10^-3-14.28*10^-6*(T-273);
exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y(1) - Cs))/Cs^2;

A = kb*b*exp(-E_b/T)*(E_b/T^2)*S^(b-1)*(y(end,4)+y(end,8))/Cs;             %%%%%%%%%%%%%%%%%%%%%%%%%
B = kb*b*exp(-E_b/T)*(b-1)*S^(b-2)*(y(end,4)+y(end,8))*exper;               %%%%%%%%%%%%%%%%%%%
C = -kb*b*exp(-E_b/T)*S^(b-1)*DCs_dT*(y(end,4)+y(end,8))/Cs^2;            %%%%%%%%%%%%%%%%%%%%
D = kb*b*exp(-E_b/T)*S^(b-1)*(theta(end,4)+theta(end,8));                  %%%%%%%%%%%%%%%%%%%%%%

DelB_dely_delT = A+B+C+D;

end
