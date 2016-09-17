clc
clear all;
%M4表明天线是4根
%信道误差下的
load globalM4Var_error.mat
% load individualM4Var_error.mat
load individualM4Var_error %数据名称load完为gaoM4Var_error
load individualM8Var_error
load globalM4Var_error

 er=[0:0.1:1];
% load rsnLinRobust  % linzhi IET2016 robust（针对非完全信息的窃听者）
figure(1)
% legend('Non-robust M=4','Non-robust M=8','pro-robust M=4','pro-robust M=8','tra-robust M=4','tra-robust M=8')
% plot(SNR,rsn1,SNR,rsn2,SNR,rsn3);
plot(er,gaoM4Var_error,'m*-',er,globalM4Var_error,'bs-',er,individualM8Var_error,'m*--',er,globalM8Var_error,'bs--','LineWidth',2 );
% xlabel('Channel Error $\varepsilon$','Interpreter','LaTex'),ylabel('Sum-Rate(bits)');
xlabel('Channel error ($\varepsilon^2$)','Interpreter','LaTex'),ylabel('Secrecy rate (bps/Hz)');
legend('Method in [11] with M=4','Proposed method with M=4','Method in [11] with M=8','Proposed method with M=8')
grid on;