clc
clear all;
N=1;
% N=1;
% M=[4 8];%中继天线数

segma=1;
% SNR=[0:2.5:20];

SNR1=2;
% er(k)=0.1;%相当于variepsilon^2

Rth=2;
% eab=(randn+i*randn);
% eae=(randn+i*randn)/sqrt(2);
% hab=sqrt(1-sqrt(er(k)))*eab+sqrt(sqrt(er(k)))*(randn+i*randn);%h_ab
% hae=sqrt(1-sqrt(er(k)))*eae+sqrt(sqrt(er(k)))*(randn+i*randn)/sqrt(2);%h_ae
load eab.mat
load hab.mat
load hab.mat
load hae.mat
for j=1:length(M)
    %     er(k)b=(randn(M(j),1)+i*randn(M(j),1));
    %     hrb=sqrt(1-sqrt(er(k)))*er(k)b+sqrt(sqrt(er(k)))*(randn(M(j),1)+i*randn(M(j),1));  %h_rb
    %     er(k)e=(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);
    %     hre=sqrt(1-sqrt(er(k)))*er(k)e+sqrt(sqrt(er(k)))*(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);%h_re
    load erb.mat
    load hrb.mat
    load ere.mat
    load hre.mat
    M=4;
    SNR=20;
    er=[0:0.1:1];%相当于variepsilon^2
%     for k=1:length(SNR)
for k=1:length(er)
% for k=5
%         snr=10.^(SNR(k)/10);
snr=10.^(SNR/10);
        snr1=10.^(SNR1/10);
        P0=snr*segma;
        Pa=P0/(M(j)+1);%第一时隙的功率
        Pr=P0-Pa;%relay的功率限制
        for tt=1:N
            %放在这儿是
            %             eab=(randn+i*randn);
            %             eae=(randn+i*randn)/sqrt(2);
            %             hab=sqrt(1-sqrt(er(k)))*eab+sqrt(sqrt(er(k)))*(randn+i*randn);%h_ab
            %             hae=sqrt(1-sqrt(er(k)))*eae+sqrt(sqrt(er(k)))*(randn+i*randn)/sqrt(2);%h_ae
            %             er(k)b=(randn(M(j),1)+i*randn(M(j),1));
            %             hrb=sqrt(1-sqrt(er(k)))*er(k)b+sqrt(sqrt(er(k)))*(randn(M(j),1)+i*randn(M(j),1));  %h_rb
            %             er(k)e=(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);
            %             hre=sqrt(1-sqrt(er(k)))*er(k)e+sqrt(sqrt(er(k)))*(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);%h_re
            V=1+2*Pa*norm(eae)^2;%文中的sigma^2+2P_A|h_ae|^2
            L=1+2*Pa*norm(eab)^2;%文中的sigma^2+2P_A|h_ab|^2
            V1=1+2*Pa*norm(hae)^2;
            L1=1+2*Pa*norm(hab)^2;
            NN=0;
            Erb=erb*erb';
            Ere=ere*ere';
            Hrb=hrb*hrb';
            Hre=hre*hre';
            
            first(j,k)=((2.^(-M(j))*exp(-1./(2*snr)))./(factorial(ceil(M(j)./2))));
            for kk=1:ceil(M(j)./2-1)
                single(kk)=((factorial(ceil(M(j)./2-1)))./(factorial(kk))).*((2./snr).^kk);
            end
            all(j,k)=sum(single);
            P(j,k)=first(j,k).*all(j,k);         %%%%%%中继处的中断概率，即未成功译码的概率
            
            
            %% subproblem for DF beamforming
            %             rmin = ( 1 + 2 * Pa * norm(hae)^2 )  /  ( 1 + 2 * Pa * norm(hae)^2 +Pr * norm(hrb)^2  );
            %             Q0=Pr/M(j) *eye( M(j), M(j) );
            %             rmax1 = 1 + 2 * Pa * norm(hae)^2 + Pr * er(k) + trace (  ( Q0 + Q0 * (Pr * eye( M(j),M(j) )-Q0  ).' * Q0 ) * hre * hre'  ); %.'为共轭转置
            %             rmax2 =  1 + 2 * Pa * norm(hab)^2 + trace ( Q0 * hrb * hrb' );
            %             rmax = rmax1 / rmax2;
            rmin =0;
            rmax=10;
            while (rmax-rmin>0.0001)
                %                 t=0.5* ( rmax + rmin );% 如果不取实数，后面CVX报错
                t=0.5* real( rmax + rmin );%
                cvx_begin sdp
                variable Phi(M(j),M(j)) symmetric;    %天线数*天线数
                variable Qx(M(j),M(j)) hermitian semidefinite;
                variable tau nonnegative;
                subject to % IET lin zhi formulation (15)
                %% 下面一个语句，也就是（16）第一个公式不符合CVX的约束规范
                %             1+2*Pa*norm(hae)^2+tau*er(k)+trace(  (Qx+Phi)*hre*hre'   ) <=t *( 1+2*Pa*norm(hab)^2 + trace(Qx * hrb * hrb') );
                % trace(  (Qx+Phi)*hre*hre'-t *( Qx * hrb * hrb' )   ) <=t-1+2*Pa*norm(hae)^2+tau*er(k);
                lefttrace= trace(  t*(Qx+Phi)*hre*hre'-Qx * hrb * hrb'    );
                rightConstraint= 1+ 2* Pa*norm(hab)^2 - t * ( 1 + 2*Pa*norm(hae)^2 + tau*er(k) );
                real(lefttrace)<=rightConstraint;
                [Phi, Qx; Qx,tau*eye( M(j),M(j) )-Qx]>=0;
                trace(Qx)<=Pr;
                %                 tau>=0;
                cvx_end
                if strfind ( cvx_status, 'Solved')
                    rmin=t;%snr针对的k=5时，t=1.56，2.34，2.74可以解。
                else
                    rmax=t;
                end
                if (strfind ( cvx_status, 'Solved')& rmax-rmin<0.1)
                    break;
                end
                    
            end
           % 接着求误差e的对偶变量lambda
            %20160902 之前一直报错，可能模型建立有问题，需要重新建立，参考huangling
            %J. Huang, and A. Lee Swindlehurst, “Robust secure transmission in MISO channels based on worst-case
            %optimization,” IEEE Trans. Signal
            %Process., vol. 60, no. 4, pp. 1696-1707, Apr. 2012.
            cvx_begin sdp
            variable eta
            variable Psi(M(j),M(j)) symmetric;    %天线数*天线数
            variable lambda
            maximize eta;
            subject to
            U11 = -Qx + lambda * eye(M(j));
            U12= Qx * hre;%目前问题是U12跟U21不对称，不知道问题所在
            U21= hre' * Qx;
            U22 = hre' * Qx * hre -lambda * er(k)  -eta;
            %                                            U_unified = [ -Qx + lambda * eye(M(j))  Qx * hre; hre.' * Qx hre' * Qx * hre -lambda * er(k)  -eta  ];
            U_unified =[U11 U12; U21 U22];
            U_unified >= 0;
            %                                         -trace( real( Qx * hre*hre.' ) )   -trace( real( Psi * hre*hre.' ) ) - lambda * er(k)>=eta;
            %                             [ Psi, Qx; Qx, lambda * eye(M(j)) -Qx]>=0;
            
            lambda>=0;
            cvx_end
            %根据对偶变量lambda求解非完全信息e
            ee0=  hre' * Qx *  (  inv(-Qx + lambda * eye(M(j)) ) );
            ee1=ee0';%因为Linzhi和huang jing的定义正好是反着的
            %             Q1=Qx;
            
            
            %                 snr_receive1(tt)=(1-P(j,k))*log2(t);         %成功译码情况下的目的节点的安全速率
            %                 snr_receive2(tt)=(1-P(j,k))*log2((L+Pr*er(k)b'*Q1*er(k)b)/(V+Pr*er(k)e'*Q1*er(k)e));       %成功译码情况下的目的节点的安全速率
            %                 snr_receive3(tt)=(1-P(j,k))*log2((L1+Pr*hrb'*Q1*hrb)/(V1+Pr*hre'*Q1*hre));     %成功译码情况下的目的节点的安全速率
            snr_receiveNumerator=(  1 + 2 * Pa * norm(hae)^2 +  hrb' * Qx * hrb );
            snr_receiveDenominator=( 1 + 2 * Pa * norm(hae)^2 + (hre + ee1)' * Qx * (hre + ee1)  );
            snr_receive1(tt)=.5 * (1-P(j,k)) * log2( snr_receiveNumerator / snr_receiveDenominator ); %林志IET2016公式（10）的前一项
            %                 if(snr_receive1(tt)==0)
            %                 NN=NN+1;
            %                 end
            % NN=NN+1;
            %% subproblem for cooper(k)ative jamming
            %完全信息计3算公式
            %                 alpha=(eye(M(j))-(hrb*hrb'/(hrb'*hrb)))*hre/norm((eye(M(j))-(hrb*hrb'/(hrb'*hrb)))*hre);%linzhi IET formulation (19)
            %                 T= alpha * alpha' * Pr; % gaobin added
            %不完全信息计算公式 %linzhi IET formulation (46)
            %实验显示，unbounded，可能scur补转换有问题
            %根据huangjing的（10）式重新构造
            cvx_begin sdp
            %             variable Psi(M(j),M(j)) symmetric;    %天线数*天线数
            variable Qz(M(j),M(j)) hermitian semidefinite;%天线数*天线数
            variable mu2 nonnegative;
            variable psi2 nonnegative;
            %             variable k nonnegative;
            maximize trace(Qz * hre * hre') - psi2 - mu2 * er(k)
            %             minimize kk * er(k)-trace( (Qz - Psi) * hre * hre' )
            subject to
            U11 = mu2 * eye( M(j) ) + Qz;
            U12 = Qz * hre;
            U21 =  hre' * Qz;
            U22 = psi2;
            U = [U11 U12; U21 U22];
            U>=0;
            trace(Qz)<=Pr;
            hrb' * Qz * hrb ==0;
            %             [ Psi, Qz; Qz, kk * eye(M(j)) + Qz ]>=0;
            %             k>=0;
            %             hrb' * Qz * hrb ==0;
            cvx_end
            % 接着求误差e的对偶变量lambda
            % 接着求误差e的对偶变量lambda,(linzhi的（23）或者huangjing的（16）)
            cvx_begin sdp
            variable eta
            variable Psi(M(j),M(j)) symmetric;    %天线数*天线数
            variable lambda
            maximize eta;
            subject to
            U11 = Qx + lambda * eye(M(j));
            U12= Qx * hre;
            U21= hre' * Qx;
            U22 = hre' * Qx * hre -lambda * er(k)  -eta;
            %                                            U_unified = [ -Qx + lambda * eye(M(j))  Qx * hre; hre.' * Qx hre' * Qx * hre -lambda * er(k)  -eta  ];
            U_unified =[U11 U12; U21 U22];
            U_unified >= 0;
            %                                         -trace( real( Qx * hre*hre.' ) )   -trace( real( Psi * hre*hre.' ) ) - lambda * er(k)>=eta;
            %                             [ Psi, Qx; Qx, lambda * eye(M(j)) -Qx]>=0;
            
            lambda>=0;
            cvx_end
            %根据对偶变量lambda求解非完全信息e
            ee02=  -hre' * Qz * ( inv(  Qz + lambda * eye(M(j)) )   );
            ee2 = ee02';
            %协同jamming 的安全速率
            snr_receive2(tt)= 0.5 * ( P(j,k) ) * log2( ( 1 + 2 * Pa * norm(hab)^2 +  hrb' * Qz * hrb )  / ...
                ( 1 + 2 * Pa * norm(hae)^2 +(hre + ee2)' * Qz * (hre + ee2) ) + (1+ (hre + ee2)' * Qz * (hre + ee2)) / ( 1+ hrb' * Qz * hrb )  ); %林志IET2016公式（10）的第二项，即协同jamming
            
            
            
            Rs(tt)=snr_receive1(tt) + snr_receive2(tt);
            %                 Rs1(tt)=snr_receive1(tt)+P(j,k)*(max(log2(L+Pr*er(k)b'*T*er(k)b)-log2(V+Pr*er(k)e'*T*er(k)e),0));   %需要修改
            %                 Rs2(tt)=snr_receive2(tt)+P(j,k)*(max(log2(L+Pr*er(k)b'*T*er(k)b)-log2(V+Pr*er(k)e'*T*er(k)e),0));   %后者为中继未译码成功时的安全速率
            %                 Rs3(tt)=snr_receive3(tt)+P(j,k)*(max(log2(L1+Pr*hrb'*T*hrb)-log2(V1+Pr*hre'*T*hre),0));   %后者为中继未译码成功时的安全速率
        end
        LinRobustM4 (j,k) =sum(Rs) / N;
        %         rsn1(j,k)=sum(Rs1)/N;
        %         rsn2(j,k)=sum(Rs2)/N;
        %         rsn3(j,k)=sum(Rs3)/N;
    end
end
gaoM4Var_error=LinRobustM4; 
save gaoM4Var_error % 不同的varepsilon产生不同的
% load rsnLinRobust  % linzhi IET2016 robust（针对非完全信息的窃听者）
figure(1)
% legend('Non-robust M=4','Non-robust M=8','pro-robust M=4','pro-robust M=8','tra-robust M=4','tra-robust M=8')
% plot(SNR,rsn1,SNR,rsn2,SNR,rsn3);
plot(er,LinRobustM4 );
xlabel('Channel Error $\varepsilon$','Interpreter','LaTex'),ylabel('Sum-Rate(bits)');
legend('LinRobustlin M=4')
grid on;