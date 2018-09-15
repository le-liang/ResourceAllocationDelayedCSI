% main function for D2D-enabled vehicular communications based on delayed
% CSI feedback from vehicles.
% Evaluate CDF of the V-UE's SINR.

% By Le Liang, Georgia Tech, Feb. 23, 2017

tic
clear;
clc

numCDF = 1e6;
rng(3); % control the random seed for randn, randi, rand

%% Parameters setup
infty = 2000; % used as infinity in the simulation
dB_Pd_max = 23; % max DUE transmit power in dBm
dB_Pc_max = 23; % max CUE transmit power in dBm

% cell parameter setup
fc = 2; % carrier frequency 2 GHz
radius = 500; % cell radius in meters
bsHgt = 25; % BS height in meters
disBstoHwy = 35; % BS-highway distance in meters
bsAntGain = 8; % BS antenna gain 8 dBi
bsNoiseFigure = 5; % BS noise figure 5 dB

vehHgt = 1.5; % vehicle antenna height, in meters
vehAntGain = 3; % vehicle antenna gain 3 dBi
vehNoiseFigure = 9; % vehicle noise figure 9 dB

stdV2V = 3; % shadowing std deviation
stdV2I = 8;
dB_sig2 = -114; % noise power in dBm

numLane = 6;
laneWidth = 4;
v = 100; % velocity
d_avg = 2.5*v/3.6; % average inter-vehicle distance according to TR 36.885
T = 1; % feedback period in ms

% QoS parameters for CUE and DUE
r0 = 0.5; % min rate for CUE in bps/Hz
dB_gamma0 = 5; % SINR_min for DUE in dB
p0 = 1e-4; % outage probability for DUE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% dB to linear scale conversion
sig2 = 10^(dB_sig2/10);
gamma0 = 10^(dB_gamma0/10);
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numCUE = 20;
numDUE = 20;
sinr = zeros(numCDF,1);

epsi_k = besselj(0,2*pi*(T*1e-3)*(fc*1e9)*(v/3.6)/(3e8));
epsi_mk = epsi_k;

%%
while (1)
    %% Generate traffic on the highway
    d0 = sqrt(radius^2-disBstoHwy^2);
    [genFlag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
    if genFlag == 1
        continue; % generated vehicles are not enough to perform simulation, jump to the next iteration.
    end
    
    %% random large-scale fading generation
    alpha_mB_ = zeros(1, numCUE);
    alpha_k_ = zeros(1, numDUE);
    alpha_kB_ = zeros(1, numDUE);
    alpha_mk_ = zeros(numCUE, numDUE);
    for m = 1 : numCUE
        dist_mB = sqrt(vehPos(indCUE(m),1)^2 + vehPos(indCUE(m),2)^2);
        dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, fc) + vehAntGain+bsAntGain-bsNoiseFigure;
        alpha_mB_(m) = 10^(dB_alpha_mB/10);
        
        for k = 1 : numDUE
            dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE2(k),1))^2 +  (vehPos(indCUE(m),2)-vehPos(indDUE2(k),2))^2);
            dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) + 2*vehAntGain-vehNoiseFigure;
            alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
        end
    end
    for k = 1 : numDUE
        dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 + (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
        dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) + 2*vehAntGain-vehNoiseFigure;
        alpha_k_(k) = 10^(dB_alpha_k/10);
        
        dist_kB = sqrt(vehPos(indDUE(k),1)^2 + vehPos(indDUE(k),2)^2);
        dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, fc)+ vehAntGain+bsAntGain-bsNoiseFigure;
        alpha_kB_(k) = 10^(dB_alpha_kB/10);
    end
    
    %%
    h_mB_ = (randn(numCUE,1)+1j*randn(numCUE,1))/sqrt(2);
    h_k_ = (randn(numCUE,numDUE)+1j*randn(numCUE,numDUE))/sqrt(2);
    h_kB_ = (randn(numCUE,numDUE)+1j*randn(numCUE,numDUE))/sqrt(2);
    h_mk_ = (randn(numCUE,numDUE)+1j*randn(numCUE,numDUE))/sqrt(2);
    
    Pc_opt_ = zeros(numCUE,numDUE);
    Pd_opt_ = zeros(numCUE,numDUE);
    
    %% power allocation design - single pair
    C_mk = zeros(numCUE, numDUE);
    for m = 1 : numCUE
        alpha_mB = alpha_mB_(m);
        h_mB = h_mB_(m);
        g_mB = alpha_mB*abs(h_mB)^2;
        for k = 1 : numDUE
            alpha_k = alpha_k_(k);
            alpha_kB = alpha_kB_(k);
            alpha_mk = alpha_mk_(m,k);
            
            h_k = h_k_(m,k);
            h_kB = h_kB_(m,k);
            h_mk = h_mk_(m,k);
            
            g_k = alpha_k*abs(h_k)^2;
            g_kB = alpha_kB*abs(h_kB)^2;
            g_mk = alpha_mk*abs(h_mk)^2;
            
            %% find Pd_opt, Pc_opt
            [ Pd_opt, Pc_opt ] = calOptPower(1e-6,sig2,Pc_max,Pd_max,alpha_k,alpha_mk,epsi_k,epsi_mk,h_k,h_mk,p0,gamma0);
            Pc_opt_(m,k) = Pc_opt;
            Pd_opt_(m,k) = Pd_opt;
            
            % C_mk stores the capacity result for the pair
            C_mk(m,k) = log2(1 + Pc_opt*g_mB/(sig2+Pd_opt*g_kB));
            if C_mk(m,k) < r0 % min rate for the V2I link
                C_mk(m,k) = -infty;
            end
        end
    end
    
    %% Reuse pair matching
    [assignmentSum, ~] = munkres(-C_mk);
    [sumVal_sum, ~] = sumAndMin(C_mk, assignmentSum);
    
    if sumVal_sum < 0 % infeasible problem
        continue;
    end

    m = 1; % an arbitrary pair
    k = assignmentSum(m);
    
    for cntSs = 1:numCDF
        % current channel
        e_k = sqrt(1-epsi_k^2)*(randn(1)+1j*randn(1))/sqrt(2);
        e_mk = sqrt(1-epsi_mk^2)*(randn(1)+1j*randn(1))/sqrt(2);
        
        g_k_d = alpha_k_(k)*(abs(epsi_k*h_k_(m,k))^2 + abs(e_k)^2);
        g_mk_d = alpha_mk_(m,k)*(abs(epsi_mk*h_mk_(m,k))^2 + abs(e_mk)^2);
        
        sinr(cntSs) = Pd_opt_(m,k)*g_k_d/(sig2+Pc_opt_(m,k)*g_mk_d);
    end
    break;
end


%%
figure;
% hold on
h1 = cdfplot(10*log10(sinr));
set(h1,'color','k')
set(h1,'linestyle','-')
legend('Algorithm 1')
xlabel('SINR (dB)')
ylabel('CDF')
title('')
set(gca,'YScale','log')
ylim([1e-4,1])
grid on

% save all data
% save('main_CDFvsSINR_Feb24')

toc