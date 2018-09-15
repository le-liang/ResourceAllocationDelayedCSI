% Plot feasible region 
% Le Liang, GaTech, Nov. 2, 2016

clear; 
clc

T = 1; % reporting period in ms
fc = 2e9; % carrier frequency 2 GHz
v = 80; % vehicle speed in km/h
c = 3e8; % speed of light
epsi_k = besselj(0,2*(T*1e-3)*pi*fc*(v/3.6)/c);
epsi_mk = epsi_k;

%% everything in dB or dBm
dB_sig2 = -114; % noise power

dB_alpha_kB = -120;
dB_alpha_mB = -76;
r0 = 1; % min rate QoS for CUE

p0 = 1e-3;
dB_gamma0 = 5; % SINR_min for DUE
dB_alpha_k = -120; % -61
dB_alpha_mk = -128;

dB_Pd_max = 23;
dB_Pc_max = 23;

%% dB to linear scale %%%
gamma0 = 10^(dB_gamma0/10); % SINR_min for DUE
sig2 = 10^(dB_sig2/10);
alpha_kB = 10^(dB_alpha_kB/10);
alpha_mB = 10^(dB_alpha_mB/10);
alpha_k = 10^(dB_alpha_k/10);
alpha_mk = 10^(dB_alpha_mk/10);

h_k = 1; % fast fading component for kth D2D
h_mk = 1;
h_mB = 1;
h_kB = 1;

Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);

%% Feasible region
% P_kmin = sig2*gamma0/(alpha_k*epsi_k^2*abs(h_k)^2+alpha_k*(1-epsi_k^2)*(-log(1-p0)));
P_kmin = sig2*gamma0/(alpha_k*epsi_k^2*abs(h_k)^2);
P_k_ = P_kmin:0.001:Pd_max; 
P_k_ = P_k_(:);
P_m_ = zeros(length(P_k_),1);
P_mb_ = zeros(length(P_k_),1); % plot the case boundary
P_m2_ = zeros(length(P_k_),1);

for ind = 1 : length(P_k_)
    P_k = P_k_(ind);
    
    tmp = 1/(1-p0)*exp(epsi_k^2*abs(h_k)^2/(1-epsi_k^2));
    P_left = 0;
    P_right = Pc_max;
    
    B = P_k*alpha_k*(1-epsi_k^2);
    C = sig2+P_right*alpha_mk*epsi_mk^2*abs(h_mk)^2;
    D = P_right*alpha_mk*(1-epsi_mk^2);
    
    while(exp(C*gamma0/B)*(1+D/B*gamma0)-tmp<=0)
        P_right = 2*P_right; % exponentially find the upperbound
        
        C = sig2+P_right*alpha_mk*epsi_mk^2*abs(h_mk)^2;
        D = P_right*alpha_mk*(1-epsi_mk^2);
    end

    while(abs(P_right-P_left)>1e-6)
        P_mid = (P_left+P_right)/2;
        C = sig2+P_mid*alpha_mk*epsi_mk^2*abs(h_mk)^2;
        D = P_mid*alpha_mk*(1-epsi_mk^2);
        
        if exp(C*gamma0/B)*(1+D/B*gamma0)-tmp>0
            P_right = P_mid;
        else
            P_left = P_mid;
        end
    end
    
    P_m_(ind) = P_mid;
    P_mb_(ind) = alpha_k*epsi_k^2*abs(h_k)^2/(alpha_mk*epsi_mk^2*abs(h_mk)^2*gamma0)*P_k - sig2/(alpha_mk*epsi_mk^2*abs(h_mk)^2);
    
    %%
    if P_k < sig2*gamma0/(alpha_k*epsi_k^2*abs(h_k)^2)
        P_m2_(ind) = 0;
    else
        P_left = 0;
        P_right = Pc_max;
        num = exp((epsi_mk^2*abs(h_mk)^2)/(1-epsi_mk^2));

        A = P_k*alpha_k*epsi_k^2*abs(h_k)^2;
        B = P_k*alpha_k*(1-epsi_k^2);
        D = P_right*alpha_mk*(1-epsi_mk^2);
        
        den1 = 1+B/(gamma0*D);
        den2 = exp((A-sig2*gamma0)/(gamma0*D));
        while num/(den1*den2)-p0<0
            P_right = 2*P_right;
            D = P_right*alpha_mk*(1-epsi_mk^2);
        end

        while(abs(P_right-P_left)>1e-6)
            P_mid = (P_left+P_right)/2;
            D = P_mid*alpha_mk*(1-epsi_mk^2);

            den1 = 1+B/(gamma0*D);
            den2 = exp((A-sig2*gamma0)/(gamma0*D));
            if num/(den1*den2)-p0 > 0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
        P_m2_(ind) = P_mid;
    end
    
end

figure
% plot(10*log10(P_k),10*log10(P_m))
plot(P_k_, P_m_)
hold on
plot(P_k_, P_mb_, 'b--')
hold on
plot(P_k_, P_m2_, 'r-')
xlabel('P_k^d')
ylabel('P_m^c')
grid


%% 3D line
obj = zeros(length(P_k_), 1);

for ind = 1 : length(P_k_)
    P_kd = P_k_(ind);
    P_mc = P_m2_(ind);
    sinr = P_mc*alpha_mB*abs(h_mB)^2/(sig2+P_kd*alpha_kB*abs(h_kB)^2);
    obj(ind) = log2(1+sinr);
end
figure
% plot3(10*log10(P_k_), 10*log10(P_m2_), obj)
plot3(P_k_, P_m2_, obj)
xlabel('P_k^d (dBm)')
ylabel('P_m^c (dBm)')
zlabel('CUE rate (bps/Hz)')
grid


% %% 3D line
% P_m = zeros(length(P_k), 1);
% obj = zeros(length(P_k), 1);
% 
% for ind = 1 : length(P_k)
%     P_kd = P_k(ind);
%     
%     P_left = 0;
%     P_right = Pc_max;
%     while(abs(P_right-P_left)>1e-6)
%         P_mid = (P_left+P_right)/2;
%         B = P_kd*alpha_k*(1-epsi_k^2);
%         C = sig2+P_mid*alpha_mk*epsi_mk^2*abs(h_mk)^2;
%         D = P_mid*alpha_mk*(1-epsi_mk^2);
%         left = 1/(1-p0)*exp(epsi_k^2*abs(h_k)^2/(1-epsi_k^2));
%         right = exp(C*gamma0/B)*(1+D/B*gamma0);
%         if left >= right
%             P_left = P_mid;
%         else
%             P_right = P_mid;
%         end
%     end
%     P_m(ind) = P_left;
%     
%     P_mc = P_m(ind);
%     sinr = P_mc*alpha_mB*abs(h_mB)^2/(sig2+P_kd*alpha_kB*abs(h_kB)^2);
%     obj(ind) = log2(1+sinr);
% end
% figure
% plot3(10*log10(P_k), 10*log10(P_m), obj)
% xlabel('P_k^d (dBm)')
% ylabel('P_m^c (dBm)')
% zlabel('CUE rate (bps/Hz)')
% grid




%% 3D mesh
% P_k = 0:1:Pd_max;
% P_m = 0:1:Pc_max;
% [P_mGrid, P_kGrid] = meshgrid(P_m, P_k);
% sinr = P_mGrid.*alpha_mB.*abs(h_mB).^2./(sig2+P_kGrid.*alpha_kB.*abs(h_kB).^2);
% obj_mat = log2(1+sinr);
% figure
% grid on
% mesh(10*log10(P_kGrid), 10*log10(P_mGrid), obj_mat)
% xlabel('P_k^d (dBm)')
% ylabel('P_m^c (dBm)')
% zlabel('CUE rate (bps/Hz)')






