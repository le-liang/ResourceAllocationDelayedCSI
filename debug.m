% % test debug
% 
% m = 1;
% k = 8;
% 
% c1 = Pd_opt_(m,k)*alpha_k_(m,k)*epsi_k^2*abs(h_k_(m,k))^2;
% c2 = Pd_opt_(m,k)*alpha_k_(m,k)*(1-epsi_k^2);
% 
% c3 = sig2 + Pc_opt_(m,k)*alpha_mk_(m,k)*epsi_mk^2*abs(h_mk_(m,k))^2;
% c4 = Pc_opt_(m,k)*alpha_mk_(m,k)*(1-epsi_mk^2);
% 
% numCDF = 1e5;
% sinr = zeros(numCDF,1);
% 
% for ii = 1:numCDF
%     e1 = (randn(1)+1j*randn(1))/sqrt(2);
%     e2 = (randn(1)+1j*randn(1))/sqrt(2);
%     
%     sinr(ii) = (c1+c2*abs(e1)^2)/(c3+c4*abs(e2)^2);
% end
% 
% figure;
% h1 = cdfplot(10*log10(sinr));
% set(h1,'color','k')
% set(h1,'linestyle','-')
% legend('Algorithm 1')
% xlabel('SINR (dB)')
% ylabel('CDF')
% title('')
% set(gca,'YScale','log')
% ylim([1e-4,1])
%     
% % debug 2
% 
% m = 1;
% k = 8;
% 
% 
% tmp = 1/(1-p0)*exp(epsi_k^2*abs(h_k_(m,k))^2/(1-epsi_k^2));
% 
% B = Pd_opt_(m,k)*alpha_k_(k)*(1-epsi_k^2);
% C = sig2+Pc_opt_(m,k)*alpha_mk_(m,k)*epsi_mk^2*abs(h_mk_(m,k))^2;
% D = Pc_opt_(m,k)*alpha_mk_(m,k)*(1-epsi_mk^2);
% 
% exp(C*gamma0/B)*(1+D/B*gamma0) - tmp
% 
% Pc_opt_(m,k)


% debug
clc

num = exp((epsi_mk^2*abs(h_mk)^2)/(1-epsi_mk^2));
P_k = 80;
P_right = 29.942778937110274; % 24.940778937110274

A = P_k*alpha_k*epsi_k^2*abs(h_k)^2
B = P_k*alpha_k*(1-epsi_k^2)
D = P_right*alpha_mk*(1-epsi_mk^2)

den1 = 1 + B/(gamma0*D)
den2 = exp((A-sig2*gamma0)/(gamma0*D))

num/(den1*den2)














