clc
clear all;
N=1;
% N=1;
% M=[4 8];%�м�������
M=4;
segma=1;
% SNR=[0:2.5:20];
SNR1=2;
er=0.1;%�൱��variepsilon^2
Rth=2;
% eab=(randn+i*randn);
% eae=(randn+i*randn)/sqrt(2);
% hab=sqrt(1-sqrt(er))*eab+sqrt(sqrt(er))*(randn+i*randn);%h_ab
% hae=sqrt(1-sqrt(er))*eae+sqrt(sqrt(er))*(randn+i*randn)/sqrt(2);%h_ae
load eab.mat
load hab.mat
load hab.mat
load hae.mat
  load erb.mat
    load hrb.mat
    load ere.mat
    load hre.mat
%     SNR=20;
for j=1:length(M)
    %     erb=(randn(M(j),1)+i*randn(M(j),1));
    %     hrb=sqrt(1-sqrt(er))*erb+sqrt(sqrt(er))*(randn(M(j),1)+i*randn(M(j),1));  %h_rb
    %     ere=(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);
    %     hre=sqrt(1-sqrt(er))*ere+sqrt(sqrt(er))*(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);%h_re
  
    for k=1:length(SNR)
% for k=5
        snr=10.^(SNR(k)/10);
% snr=10.^(SNR/10);
        snr1=10.^(SNR1/10);
        P0=snr*segma;
        Pa=P0/(M(j)+1);%��һʱ϶�Ĺ���
        Pr=P0-Pa;%relay�Ĺ�������
        for tt=1:N
            %���������
            %             eab=(randn+i*randn);
            %             eae=(randn+i*randn)/sqrt(2);
            %             hab=sqrt(1-sqrt(er))*eab+sqrt(sqrt(er))*(randn+i*randn);%h_ab
            %             hae=sqrt(1-sqrt(er))*eae+sqrt(sqrt(er))*(randn+i*randn)/sqrt(2);%h_ae
            %             erb=(randn(M(j),1)+i*randn(M(j),1));
            %             hrb=sqrt(1-sqrt(er))*erb+sqrt(sqrt(er))*(randn(M(j),1)+i*randn(M(j),1));  %h_rb
            %             ere=(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);
            %             hre=sqrt(1-sqrt(er))*ere+sqrt(sqrt(er))*(randn(M(j),1)+i*randn(M(j),1))/sqrt(2);%h_re
            V=1+2*Pa*norm(eae)^2;%���е�sigma^2+2P_A|h_ae|^2
            L=1+2*Pa*norm(eab)^2;%���е�sigma^2+2P_A|h_ab|^2
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
            P(j,k)=first(j,k).*all(j,k);         %%%%%%�м̴����жϸ��ʣ���δ�ɹ�����ĸ���
            
            
            %% subproblem for DF beamforming
            %             rmin = ( 1 + 2 * Pa * norm(hae)^2 )  /  ( 1 + 2 * Pa * norm(hae)^2 +Pr * norm(hrb)^2  );
            %             Q0=Pr/M(j) *eye( M(j), M(j) );
            %             rmax1 = 1 + 2 * Pa * norm(hae)^2 + Pr * er + trace (  ( Q0 + Q0 * (Pr * eye( M(j),M(j) )-Q0  ).' * Q0 ) * hre * hre'  ); %.'Ϊ����ת��
            %             rmax2 =  1 + 2 * Pa * norm(hab)^2 + trace ( Q0 * hrb * hrb' );
            %             rmax = rmax1 / rmax2;
%             t_start = tic;
            t_start = tic;
                cvx_begin sdp quiet
                variable GPhi(M(j),M(j)) symmetric;    %������*������
                variable Gx(M(j),M(j)) hermitian semidefinite;
                variable tau nonnegative;
                variable phi nonnegative;
                maximize phi * ( 1 + 2*Pa*norm(hab)^2 ) + trace (Gx * hrb * hrb')
                subject to % IET lin zhi formulation (15)
                %% ����һ����䣬Ҳ���ǣ�16����һ����ʽ������CVX��Լ���淶
                %             1+2*Pa*norm(hae)^2+tau*er+trace(  (Qx+Phi)*hre*hre'   ) <=t *( 1+2*Pa*norm(hab)^2 + trace(Qx * hrb8 * hrb8') );
                % trace(  (Qx+Phi)*hre*hre'-t *( Qx * hrb8 * hrb8' )   ) <=t-1+2*Pa*norm(hae)^2+tau*er;
%                 lefttrace= trace(  t*(Gx+GPhi)*hre8*hre8'-Qx * hrb8 * hrb8'    );
%                 rightConstraint= 1+ 2* Pa*norm(hab)^2 - t * ( 1 + 2*Pa*norm(hae)^2 + tau*er );
%                 real(lefttrace)<=rightConstraint;
                phi * ( 1 + 2*Pa*norm(hab)^2 ) + tau*er + trace( (Gx+GPhi)*hre*hre' )==1;
                [GPhi, Gx; Gx,tau*eye( M(j),M(j) )-Gx]>=0;
                trace(Gx)<=phi*2*Pr;
                %                 tau>=0;
                cvx_end
                globalM4Var_snrTimecooper=toc(t_start);       
            save globalM4Var_snrTimecooper
           % ���������e�Ķ�ż����lambda
            %20160902 ֮ǰһֱ��������ģ�ͽ��������⣬��Ҫ���½������ο�huangling
            %J. Huang, and A. Lee Swindlehurst, ��Robust secure transmission in MISO channels based on worst-case
            %optimization,�� IEEE Trans. Signal
            %Process., vol. 60, no. 4, pp. 1696-1707, Apr. 2012.
            cvx_begin sdp
            variable eta
            variable Psi(M(j),M(j)) symmetric;    %������*������
            variable lambda
            maximize eta;
            subject to
            U11 = -Qx + lambda * eye(M(j));
            U12= Qx * hre;%Ŀǰ������U12��U21���Գƣ���֪����������
            U21= hre' * Qx;
            U22 = hre' * Qx * hre -lambda * er  -eta;
            %                                            U_unified = [ -Qx + lambda * eye(M(j))  Qx * hre; hre.' * Qx hre' * Qx * hre -lambda * er  -eta  ];
            U_unified =[U11 U12; U21 U22];
            U_unified >= 0;
            %                                         -trace( real( Qx * hre*hre.' ) )   -trace( real( Psi * hre*hre.' ) ) - lambda * er>=eta;
            %                             [ Psi, Qx; Qx, lambda * eye(M(j)) -Qx]>=0;
            
            lambda>=0;
            cvx_end
            %���ݶ�ż����lambda������ȫ��Ϣe
            ee0=  hre' * Qx *  (  inv(-Qx + lambda * eye(M(j)) ) );
            ee1=ee0';%��ΪLinzhi��huang jing�Ķ��������Ƿ��ŵ�
            %             Q1=Qx;
            
            
            %                 snr_receive1(tt)=(1-P(j,k))*log2(t);         %�ɹ���������µ�Ŀ�Ľڵ�İ�ȫ����
            %                 snr_receive2(tt)=(1-P(j,k))*log2((L+Pr*erb'*Q1*erb)/(V+Pr*ere'*Q1*ere));       %�ɹ���������µ�Ŀ�Ľڵ�İ�ȫ����
            %                 snr_receive3(tt)=(1-P(j,k))*log2((L1+Pr*hrb'*Q1*hrb)/(V1+Pr*hre'*Q1*hre));     %�ɹ���������µ�Ŀ�Ľڵ�İ�ȫ����
            snr_receiveNumerator=(  1 + 2 * Pa * norm(hae)^2 +  hrb' * Qx * hrb );
            snr_receiveDenominator=( 1 + 2 * Pa * norm(hae)^2 + (hre + ee1)' * Qx * (hre + ee1)  );
            snr_receive1(tt)=.5 * (1-P(j,k)) * log2( snr_receiveNumerator / snr_receiveDenominator ); %��־IET2016��ʽ��10����ǰһ��
            %                 if(snr_receive1(tt)==0)
            %                 NN=NN+1;
            %                 end
            % NN=NN+1;
            %% subproblem for cooperative jamming
            %��ȫ��Ϣ��3�㹫ʽ
            %                 alpha=(eye(M(j))-(hrb*hrb'/(hrb'*hrb)))*hre/norm((eye(M(j))-(hrb*hrb'/(hrb'*hrb)))*hre);%linzhi IET formulation (19)
            %                 T= alpha * alpha' * Pr; % gaobin added
            %����ȫ��Ϣ���㹫ʽ %linzhi IET formulation (46)
            %ʵ����ʾ��unbounded������scur��ת��������
            %����huangjing�ģ�10��ʽ���¹���
            cvx_begin sdp
            %             variable Psi(M(j),M(j)) symmetric;    %������*������
            variable Qz(M(j),M(j)) hermitian semidefinite;%������*������
            variable mu2 nonnegative;
            variable psi2 nonnegative;
            %             variable k nonnegative;
            maximize trace(Qz * hre * hre') - psi2 - mu2 * er
            %             minimize kk * er-trace( (Qz - Psi) * hre * hre' )
            subject to
            U11 = mu2 * eye( M(j) ) + Qz;
            U12 = Qz * hre;
            U21 =  hre' * Qz;
            U22 = psi2;
            U = [U11 U12; U21 U22];
            U>=0;
            trace(Qz)<=2*Pr-trace(Qx);
            hrb' * Qz * hrb ==0;
            %             [ Psi, Qz; Qz, kk * eye(M(j)) + Qz ]>=0;
            %             k>=0;
            %             hrb' * Qz * hrb ==0;
            cvx_end
            % ���������e�Ķ�ż����lambda
            % ���������e�Ķ�ż����lambda,(linzhi�ģ�23������huangjing�ģ�16��)
            cvx_begin sdp
            variable eta
            variable Psi(M(j),M(j)) symmetric;    %������*������
            variable lambda
            maximize eta;
            subject to
            U11 = Qx + lambda * eye(M(j));
            U12= Qx * hre;
            U21= hre' * Qx;
            U22 = hre' * Qx * hre -lambda * er  -eta;
            %                                            U_unified = [ -Qx + lambda * eye(M(j))  Qx * hre; hre.' * Qx hre' * Qx * hre -lambda * er  -eta  ];
            U_unified =[U11 U12; U21 U22];
            U_unified >= 0;
            %                                         -trace( real( Qx * hre*hre.' ) )   -trace( real( Psi * hre*hre.' ) ) - lambda * er>=eta;
            %                             [ Psi, Qx; Qx, lambda * eye(M(j)) -Qx]>=0;
            
            lambda>=0;
            cvx_end
            %���ݶ�ż����lambda������ȫ��Ϣe
            ee02=  -hre' * Qz * ( inv(  Qz + lambda * eye(M(j)) )   );
            ee2 = ee02';
            %Эͬjamming �İ�ȫ����
            snr_receive2(tt)= 0.5 * ( P(j,k) ) * log2( ( 1 + 2 * Pa * norm(hab)^2 +  hrb' * Qz * hrb )  / ...
                ( 1 + 2 * Pa * norm(hae)^2 +(hre + ee2)' * Qz * (hre + ee2) ) + (1+ (hre + ee2)' * Qz * (hre + ee2)) / ( 1+ hrb' * Qz * hrb )  ); %��־IET2016��ʽ��10���ĵڶ����Эͬjamming
            
            
            
            Rs(tt)=snr_receive1(tt) + snr_receive2(tt);
            %                 Rs1(tt)=snr_receive1(tt)+P(j,k)*(max(log2(L+Pr*erb'*T*erb)-log2(V+Pr*ere'*T*ere),0));   %��Ҫ�޸�
            %                 Rs2(tt)=snr_receive2(tt)+P(j,k)*(max(log2(L+Pr*erb'*T*erb)-log2(V+Pr*ere'*T*ere),0));   %����Ϊ�м�δ����ɹ�ʱ�İ�ȫ����
            %                 Rs3(tt)=snr_receive3(tt)+P(j,k)*(max(log2(L1+Pr*hrb'*T*hrb)-log2(V1+Pr*hre'*T*hre),0));   %����Ϊ�м�δ����ɹ�ʱ�İ�ȫ����
        end
        PriorityRobustM4 (j,k) =sum(Rs) / N;
        %         rsn1(j,k)=sum(Rs1)/N;
        %         rsn2(j,k)=sum(Rs2)/N;
        %         rsn3(j,k)=sum(Rs3)/N;
    end
end
% save PriorityRobustM4.mat

% load rsnLinRobust  % linzhi IET2016 robust����Է���ȫ��Ϣ�������ߣ�
figure(1)
% legend('Non-robust M=4','Non-robust M=8','pro-robust M=4','pro-robust M=8','tra-robust M=4','tra-robust M=8')
% plot(SNR,rsn1,SNR,rsn2,SNR,rsn3);
plot(SNR,PriorityRobustM4 );
xlabel('System SNR��dB��'),ylabel('Sum-Rate(bits)');
legend('LinRobustPriority M=4')
grid on;