function T_new = check_constraint1(C,T,M,DH_DT)
 
%Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
%Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;

T_computed = T+(M*DH_DT);
Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_computed-273) - 7.14 * 10^-6 * (T_computed-273)^2 ;
Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T_computed-273) - 8.1*10^-6 * (T_computed-273)^2 ;

T_new =0;
if Cs < C 
    if C < Cm
        T_new = T_computed;
    else 
        % Evalute T form Cm
        p = [8.1*10^-6  -2.46*10^-3  (-7.76 * 10^-2+C)];
        disp('Inside the Cm constraint');
        r = roots(p)
        %if 27 < r(1) < 52
            %T_new = r(1)+273
            
        %end
        %if  27<r(2)< 52
            %T_new = r(2)+273
        %end
        %if 27<r(1)<52 && 27<r(2)<52 
            %T_new = min(r) +273
        
        %end
            
        if min(r) > 0 
        
            T_new = min(r)+273
        else 
            T_new = max(r)+273
        end
           
       
    end
        
    if T_new == 0
        display('Error');
    end
end
    
if C <  Cs 
    % Evaluat0e T from Cs expression
    p = [7.14*10^-6 -2.46*10^-3 (-6.29*10^-2+C)];
    disp('inside the Cs constraint');
    r = roots(p)
    if isreal(r(1)) ==0 
       C = Cs ;
       p = [7.14*10^-6 -2.46*10^-3 (-6.29*10^-2+C)];
       r = roots(p)
    end
    
    if min(r) > 0 
        
        T_new = min(r)+273
    else 
        T_new = max(r)+273
    end
    
    %Cm =  7.76 * 10^-2 + 2.46*10^-3 * (T_new-273) - 8.1*10^-6 * (T_new-273)^2 ;
    
    
    if C > Cm 
        
        p = [8.1*10^-6  -2.46*10^-3  (-7.76*10^-2+C)];
        disp('Inside the Cs Cm constraint');
        r = roots(p)
        
        t_min = min(r);
        if t_min > 0 
            T_new = t_min+273 
        else
            T_new = max(r) + 273
        end   
        %if 27 <= r(1) <= 52
            %T_new = r(1)+273
            
        %end
        %if  27<=r(2)<=52
            %T_new = r(2)+273
        %end
        %if 27<=r(1)<=52 && 27<=r(2)<=52 
            %T_new = min(r) +273
        
        %end
        
        if T_new == 0
            display('Error');
        end
    end
end
if T_new < 300 
    T_new = 300
end

if T_new > 325 
 T_new = 325

end

if T_new ==0 
 display('Here')
 T_new = T_computed;
end

end            




