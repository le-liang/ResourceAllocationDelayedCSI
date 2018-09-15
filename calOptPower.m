function [ Pd_opt, Pc_opt ] = calOptPower(epsi,sig2,Pc_max,Pd_max,alpha_k,alpha_mk,epsi_k,epsi_mk,h_k,h_mk,p0,gamma0)
% calOptPower: compute the optimal power allocation for a single V-UE and I-UE pair

% By Le Liang, Georgia Tech, Mar. 2, 2017

den0 = alpha_mk*(1-epsi_mk^2)*(1/p0-1)*epsi_k^2*abs(h_k)^2 - (1-epsi_k^2)*alpha_mk*epsi_mk^2*abs(h_mk)^2;
Pc0 = (1-epsi_k^2)*sig2/den0;
Pd0 = Pc0*gamma0*alpha_mk*(1-epsi_mk^2)*(1-p0)/(alpha_k*(1-epsi_k^2)*p0);

if Pd_max <= Pd0
    %% Case I
    
    B = Pd_max*alpha_k*(1-epsi_k^2);
    C = sig2+Pc_max*epsi_mk^2*alpha_mk*abs(h_mk)^2;
    D = Pc_max*alpha_mk*(1-epsi_mk^2);
    
    tmp = 1/(1-p0)*exp(epsi_k^2*abs(h_k)^2/(1-epsi_k^2));
    if exp(C*gamma0/B)*(1+D/B*gamma0) - tmp > 0
        % P_opt = (Pd_max, Pc_dmax)
        Pd_opt = Pd_max;
        
        P_left = 0;
        P_right = Pc_max;
        B = Pd_max*alpha_k*(1-epsi_k^2);
        while(abs(P_right-P_left)>epsi)
            P_mid = (P_left+P_right)/2;
            C = sig2+P_mid*epsi_mk^2*alpha_mk*abs(h_mk)^2;
            D = P_mid*alpha_mk*(1-epsi_mk^2);
            
            if exp(C*gamma0/B)*(1+D/B*gamma0)-tmp>0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
        Pc_opt = P_mid;
    else
        % P_opt = (Pd_cmax, Pc_max)
        Pc_opt = Pc_max;
        
        P_left = 0;
        P_right = Pd_max;
        C = sig2+Pc_max*alpha_mk*epsi_mk^2*abs(h_mk)^2;
        D = Pc_max*alpha_mk*(1-epsi_mk^2);
        while(abs(P_right-P_left)>epsi)
            P_mid = (P_left+P_right)/2;
            B = P_mid*alpha_k*(1-epsi_k^2);
            
            if exp(C*gamma0/B)*(1+D/B*gamma0)-tmp<0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
        Pd_opt = P_mid;
    end
    
elseif Pc_max > Pc0
    %% Case II
    
    num = (epsi_mk^2*abs(h_mk)^2)/(1-epsi_mk^2); % take log
    
    A = Pd_max*alpha_k*epsi_k^2*abs(h_k)^2;
    B = Pd_max*alpha_k*(1-epsi_k^2);
    D = Pc_max*alpha_mk*(1-epsi_mk^2);
    den1 = log(1+B/(gamma0*D));
    den2 = (A-sig2*gamma0)/(gamma0*D);
    
    if num-(den1+den2)-log(p0) > 0
        Pd_opt = Pd_max;
        
        P_left = 0;
        P_right = Pc_max;
        A = Pd_max*alpha_k*epsi_k^2*abs(h_k)^2;
        B = Pd_max*alpha_k*(1-epsi_k^2);
        while(abs(P_right-P_left)>epsi)
            P_mid = (P_left+P_right)/2;
            D = P_mid*alpha_mk*(1-epsi_mk^2);
            den1 = log(1+B/(gamma0*D));
            den2 = (A-sig2*gamma0)/(gamma0*D);
            if num-(den1+den2)-log(p0) > 0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
        Pc_opt = P_mid;
    else
        Pc_opt = Pc_max;
        P_left = 0;
        P_right = Pd_max;
        D = Pc_max*alpha_mk*(1-epsi_mk^2);
        while(abs(P_right-P_left)>epsi)
            P_mid = (P_left+P_right)/2;
            A = P_mid*alpha_k*epsi_k^2*abs(h_k)^2;
            B = P_mid*alpha_k*(1-epsi_k^2);
            
            den1 = log(1+B/(gamma0*D));
            den2 = (A-sig2*gamma0)/(gamma0*D);
            if num-(den1+den2)-log(p0)<0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
        Pd_opt = P_mid;
    end
else
    tmp = 1/(1-p0)*exp(epsi_k^2*abs(h_k)^2/(1-epsi_k^2));

    % P_opt = (Pd_cmax, Pc_max)
    Pc_opt = Pc_max;
    
    P_left = 0;
    P_right = Pd_max;
    C = sig2+Pc_max*alpha_mk*epsi_mk^2*abs(h_mk)^2;
    D = Pc_max*alpha_mk*(1-epsi_mk^2);
    while(abs(P_right-P_left)>epsi)
        P_mid = (P_left+P_right)/2;
        B = P_mid*alpha_k*(1-epsi_k^2);
        
        if exp(C*gamma0/B)*(1+D/B*gamma0)-tmp<0
            P_right = P_mid;
        else
            P_left = P_mid;
        end
    end
    Pd_opt = P_mid;
    
end

end

