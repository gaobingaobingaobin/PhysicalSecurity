%% 本代码区别以往要产生image，CPUtime PSNR iSNR，i五项，而是只是简单的  PSNR CPUtime，i 三项
%% 主要的，本代码的alpha=0.6，tolerance 是1e-6,1e-5和1e-4
clc
clear
close all
%% LASSO data which extracted from simulationTable1.m
load globalM8Var_snrTime
load globalM8TimeBisection
load globalM4Var_snrTimecooper
load globalM4Var_snrTime
% tempY{1}=[35,30;36,32;33,28;36,30;30,26];
tempY{1}=[globalM4Var_snrTimecooper,globalM4Var_snrTime;globalM8Var_snrTime,globalM8TimeBisection];
bar(tempY{1},1)
% xlabel('Dimension of $A$','Interpreter','latex','fontsize',12);
xlabel('Antenna number M');
ylabel('CPU time (Sec.)')
legend('DF beamforming based on bisection method','DF beamforming based on  Charnes-Cooper transformation')
set(gca,'XTickLabel',{'M=4','M=8'})
% tempY{2}=[34,31;35,32;32,29;35,31;28,25];
% for i=1:2% 
% 
%     
% subplot(1,2,i)
% bar(tempY{i},1)
% xlabel('Dimension of $A$','Interpreter','latex','fontsize',12);
% ylabel('Iteration Numbers')
% legend('PD-SADMM','PID-SADMM')
% % set(gca,'XTickLabel',{'900*3000','1050*3500','1200*4000','1350*4500','1500*5000'})
% 
% end
