function DelB_dT = DelB_T(S,T,theta,y)

%display('INside function');
Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
%Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
DCs_dT = 2.46*10^-3-14.28*10^-6*(T-273);
exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y(1) - Cs))/Cs^2;


E_b = 7517;
b = 1.45;

kb = 285;
A = kb*exp(-E_b/T)*E_b/T^2*S^b*(y(5)+y(9));   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = A + kb*exp(-E_b/T)*b*S^(b-1)*exper*(y(end,5)+y(end,9));
DelB_dT = kb*exp(-E_b/T)*S^b*(theta(5)+theta(9))+B;

end