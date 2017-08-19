% Tushar Gupta 13CH30023 
% Code for implementing the determinisitic approach for modelling a seeded
% batch crystallizer 

% Deterministic appproach 
%parameters   
b =1.45;
kb = 285;
E_b = 7517;
kv = 0.54;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859;
rho = 2.66*10^-12;
%%%%%%%

t0 = 0;
delta_t = 1;
batch_time = 1800;
tf = batch_time;

n = length(t0:delta_t:tf);

zf = [0 0 0 0 1 0 0 0 -1];
theta0 = [0 0 0 0 0 0 0 0 0];
fi_f = [0 0 0 0 0 0 0 0 0];
y0 = [0.1743 66.66 1.83*10^4 5.05*10^6 1.93*10^9 0.867 0 0 0];
M = -10^-7;
%-8*10^-9    
tolerance = 10^-2;
% -4*10^6
% Initialising matricess for storing values for variables at each time step
N = 7;

iteration = N;
T_vec = ones(1,length(t0:delta_t:batch_time))*323;
DH_vec = zeros(iteration,length(t0:delta_t:batch_time));
C_vec = zeros(iteration,length(t0:delta_t:batch_time));


% iteration X time steps 
iteration = 1;


%figure;hold on;
%title('Hamiltonian Derivative')

while iteration < N
   
    y_mat = zeros(length(t0:delta_t:tf),9); %1801 * 9   
    z_mat = zeros(length(t0:delta_t:tf),9);
    theta_mat = zeros(length(t0:delta_t:tf),9);
    fi_mat = zeros(length(t0:delta_t:tf),9);
    DelH_dy_mat = zeros(length(t0:delta_t:tf),9); %1801*9 
    DelH_dz_mat = zeros(length(t0:delta_t:tf),9);

    y_mat(1,:) = y0;
    z_mat(1,:) = zf;   % Backward 
    theta_mat(1,:) = theta0;
    fi_mat(1,:) = fi_f;   %Backward  
    
    %y forward integration 
    
    for t = 1:1:(n-1)
        t_curr = t;
        t_next = t+1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        %T_vec(1,t_curr)
        B = (kb*exp(-E_b/T_vec(1,t_curr)))*S^b*(y_mat(t_curr,5)+y_mat(t_curr,9));    
        y_mat(t_next,:) = y_forward(y_mat(t_curr,:),delta_t,G,B,t_curr,t_next);
        
    end
     
    
    % z backward integration     
    for t = 1:1:n-1
        t_curr = t;
        t_next = t+1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        z_mat(t_next,:) = z_backward(z_mat(t_curr,:),y_mat(t_curr,:),delta_t,G,T_vec(1,t_curr),S,t_curr,t_next)
    end
    break;
% 
    for t_curr = 1:1:n
        
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        B = (kb*exp(-E_b/T_vec(1,t_curr)))*S^b*(y_mat(t_curr,5)+y_mat(t_curr,9));                   
        DelH_dy_mat(t_curr,:) = DelH(G,z_mat(t_curr,:),y_mat(t_curr,:),S,T_vec(1,t_curr));
        DelH_dz_mat(t_curr,:) = DelH_z(G,B,rho,kv,y_mat(t_curr,:))*10^-6;
        %display('done')
        
    end
    
     T_vec(1,t_curr)
        %%% Theta Forward Integration 
     for t_curr = 1:1:(n-1) 
        
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        DelG_dT = DelG_T(S,theta_mat(t_curr,:),T_vec(1,t_curr),y_mat(t_curr,1));
        DelB_dT = DelB_T(S,T_vec(1,t_curr),theta_mat(t_curr,:),y_mat(t_curr,:));
        t_next = t_curr+1;
        theta_mat(t_next,:) = theta_forward(theta_mat(t_curr,:),delta_t,G,DelG_dT,DelB_dT,y_mat(t_curr,:),t_curr,t_next);
        
     end
     
     for t = 1:1:n-1
        t_curr = t;
        t_next = t+1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        
        fi_mat(t_next,:) = fi_backward(fi_mat(t_curr,:),delta_t,G,y_mat(t_curr,:),z_mat(t_curr,:),theta_mat(t_curr,:),T_vec(1,t_curr),t_curr,t_next);
   
     end
    %theta_mat
    %fi_mat 
    for t = 1:1:n
        sum = 0;
        for i=1:1:9
            sum = sum + DelH_dy_mat(t,i)*theta_mat(t,i) + DelH_dz_mat(i)*fi_mat(t,i);    %%$$$$$$$$$$$
            %theta_mat(t,i)
        end
        DH_vec(iteration,t)=sum;
    end
    
    for t = 1:1:n
       
        if abs(DH_vec(iteration,t)) > tolerance
            T_vec(1,t)  = check_constraint(y_mat(t,1),T_vec(1,t),M,DH_vec(iteration,t));
            %T_vec(t) = T_vec(t) + M*DH_vec(iteration,t);
        end
        
        
    end
    
    figure()
    plot(t0:delta_t:tf,DH_vec(iteration,:));
    
    iteration = iteration +1
    
 
end
figure
plot(t0:delta_t:tf,T_vec)
title('Teperature Profile');
length(T_vec)
figure

plot(t0:delta_t:tf,(y_mat(:,5)-y_mat(:,9)))
title('Objective Function');