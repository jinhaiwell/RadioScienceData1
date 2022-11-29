%%%%Figure 7,8
clc;
clear;
close all;

j=sqrt(-1);
% M=5;
% R0=10000;
% for ki=1:M
%     X(ki)=abs(randn*10);
%     Y(ki)=abs(randn*10);
%     Am(ki)=abs(randn)+5;
%     phi=(120+10*randn)*pi/180;
%     dn(ki)=sqrt(R0^2+X(ki)^2-2*R0*X(ki)*cos(phi));
% end
% X=sort(X);
% Y=sort(Y);
% Am=sort(Am);
% dn=sort(dn);
Amp=[1.1 3.5 5.2 7 8.5];  
X=[1.2 2.401 4.709 6.100 9.5];
Type=[1 0.5 0 -0.5 -1];
Bw=2e9;   
f0=10e9;
fmax=f0+Bw;
n=512;
M=length(Amp);
fm=linspace(f0,fmax,n);
deltf=Bw/n;
c=3e8;
snr=1:1:30;

for kj=1:5 
% kj=3;
for sr=snr
  Xi_esprit=0;
  Ti_esprit=0;
  Ai_esprit=0;
  Xi_ssa=0;
  Ai_ssa=0;
  Api=0;
  for test=1:100 
      
  theta=0;
  Es=zeros(1,n);
  
  for ki=1:1:M
    R(ki)=Amp(ki)*exp(j*pi*Type(kj)/2)*exp(j*(-4)*pi*f0*X(ki)*cos(theta)/c);
    Z(ki)=(1+Type(kj)*deltf/f0)*exp(j*(-4)*pi*deltf*X(ki)*cos(theta)/c);
    for k=1:1:n
    Es(k)=awgn(R(ki)*(Z(ki).^k),sr)+Es(k);  
    end
  end
%%%%%%%%%%%%%%%%ESPRIT  
g=100;  
Rtemp=zeros(1,g);
for m=1:1:g     
    for i=1:1:(n-m+1) % n=512  m=100
        Rtemp(m)=Rtemp(m)+conj(Es(i))*Es(i+m-1);
    end 
    R(m)=Rtemp(m)/(n-m+1);     
end

t=g-1;        
for i=1:1:t
    for jj=1:1:t
        if (i>=jj) 
            Rxx(i,jj)=R(i-jj+1);
        end
        if (jj>i)  
            Rxx(i,jj)=conj(R(jj-i+1));
        end
    end
end

for i=1:1:t   
    for jj=1:1:t
        if (i>jj) 
            Rxy(i,jj)=R(i-jj);          
        end
        if (jj>=i)  
            Rxy(i,jj)=conj(R(jj-i+2));   
        end 
    end 
end 

[U,S,V] = svd(Rxx);
U1=U(:,1:M);
V1=V(:,1:M);
S1=S(1:M,1:M);
mm=min(diag(S1));
I =  eye(t,t);
Zg(1:1,1:t)=0;
Zg(2:t,1:t)=I(1:(t-1),1:t);
Cxy=Rxy-mm*Zg;
Y1=U1'*Cxy*V1;
Y2=S1;
Y11=inv((Y1'*Y1))*Y1';
YY=Y11*Y2;
[v,d]=eig(YY);
Pout=diag(d);
for ki=1:1:M
    Xout_esprit1(ki)=angle(Pout(ki))/(j*(-4)*pi*deltf*cos(theta)/c);
end
Xout_esprit=sort(abs(Xout_esprit1));%X=[1.2 2.401 4.709 6.100 9.5]
Xi_esprit=Xi_esprit+Xout_esprit;

for ki=1:M
%     s(:,ki)=(j*fm/f0).^Type(kj).*exp(-1j*4*pi*fm/c*Xout1(ki));
    s1(:,ki)=exp(-4*pi*fm/c*Xout_esprit1(ki));
end
Aout_esprit1=inv(s1'*s1)*s1'*Es';
Aout_esprit=sort(abs(Aout_esprit1))';%A=[1.1 3.5 5.2 7 8.5]
Ai_esprit=Ai_esprit+Aout_esprit;

hat_T=(abs(Pout)'-1)*f0/deltf; 
for ki=1:1:M
    Tout(ki)=(abs(Pout(ki))-1)*f0/deltf;%[0.5 -1 0 1 -0.5]
end
Ti_esprit=Ti_esprit+Tout;
%%%%%%%%%%%%%%%%%%%%%%%%% SSM %%%%%%%%%%
P=n/2; %
L=n-P+1;
Y0=hankel(Es(1:P),Es(P:end));
[U,S,V] = svd(Y0);
U1=U(:,1:M);
S1=S(1:M,1:M);
V1=V(1:M,:);
omga=U1*sqrt(S1);
gama=sqrt(S1)*V1;
H1=U1*S1*V1;%dominant components of the spectral norm sense
H2=omga*gama;%H1=H2
omga1=omga(2:P,:);%delete the first row
omga2=omga(1:P-1,:);%delete the last row
A=inv(omga2'*omga2)*omga2'*omga1;
% om=omga2*A;
gama1=gama(:,2:P);%delete the first column
gama2=gama(:,1:P-1);%delete the last column
% A=gama1*gama2'*pinv(gama2*gama1'); 
C=omga(1,:);
B=gama(:,1);
omgN=[];
gamN=[];
for ki=1:n
    omgan=C*A^(ki-1);
    omgN=[omgN;omgan];
end
hat_B=pinv(omgN'*omgN)*omgN'*Es';
for ki=1:n
    gaman=(A^(ki-1))*B;
    gamN=[gamN,gaman];
end
hat_C=Es*gamN'*pinv(gamN*gamN');
[vec,lamda]=eig(A);
egv=diag(lamda);
%A State Identification Method for 1-D Measurements with Gaps
fi=angle(diag(lamda));%refers to the phase of the eigenvalue lamda
Xout_ssa=-fi*c/(sqrt(-1)*4*pi*deltf);%scattering position  [1.2 2.401 4.709 6.100 9.5];
Xout=sort(abs(Xout_ssa))';
Xi_ssa=Xi_ssa+Xout;

hat_T=(abs(egv)'-1)*f0/deltf; 

s(n,M)=0;
for ki=1:M
    s(:,ki)=exp(-4*pi*fm/c*Xout_ssa(ki));
end
Aout_ssa1=pinv(s'*s)*s'*Es';
Aout_ssa=sort(abs(Aout_ssa1))';
Ai_ssa=Ai_ssa+Aout_ssa;

psi=vec;
for ki=1:M
numr(ki)=(C*psi(:,ki))*((psi(:,ki).')*hat_B);
deno(ki)=egv(ki).^(f0/deltf);
end
ai=(numr)./(deno);
[Am_hat,pos]=sort(abs(ai),'descend');
Api=Api+(Am_hat); 
%%%%%%%%%%%%%%%%%%%MUSIC%%%%%%%%%%%%%%%%
for kn=1:M
    for ki=1:n
        astar(ki,kn)=(j*(f0+(ki-1)*deltf)/f0)^Type(kj)*exp(-4*pi*j*pi*(f0+(ki-1)*deltf)*X(kn)/c);
    end
end
Xnew=astar*Amp';
Rxx1=Xnew*Xnew';
[Um,Sm,Vm] = svd(Rxx1);
Um1=Um(:,1:M);
Um2=Um(:,M+1:end);
Sm1=Sm(1:M,1:M);
Sm2=Sm(M+1:end,M+1:end);
Vm1=Vm(M,:);
Vm2=Vm(M+1,end,:);
Gn=Um2*Um2';
PwStar=astar'*astar*pinv(astar'*Gn*astar);
value=max(max(abs(PwStar)));
[row,col]=find(value==abs(PwStar));
hat_A=pinv(astar'*astar)*astar'*Xnew;
hat_A=abs(hat_A);
  end
Xi_esprit=Xi_esprit/test;
Xi_ssa=Xi_ssa/test;
Ai_esprit=Ai_esprit/test;
Ai_ssa=Ai_ssa/test;
Ti_esprit=Ti_esprit/test;
RMSE1(sr) = (sqrt(sum((Xi_esprit-X).^2)/5));
RMSE2(sr) = (sqrt(sum((Xi_ssa-X).^2)/5));
RMSE3(sr) = 10*log10(sqrt(sum((Ai_esprit-Amp).^2)/5));
RMSE4(sr) = 10*log10(sqrt(sum((Ai_ssa-Amp).^2)/5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta=0;
% Es_o(n)=0;
% Es_re_esprit(n)=0;
% Es_re_ssa(n)=0;
% for ki=1:1:M
%     R1(ki)=Amp(ki)*exp(j*pi*Type(kj)/2)*exp(j*(-4)*pi*f0*X(ki)*cos(theta)/c);
%     Z1(ki)=(1+Type(kj)*deltf/f0)*exp(j*(-4)*pi*deltf*X(ki)*cos(theta)/c);
%     R2(ki)=Ai_esprit(ki)*exp(j*pi*Type(kj)/2)*exp(j*(-4)*pi*f0*Xi_esprit(ki)*cos(theta)/c);
%     Z2(ki)=(1+Type(kj)*deltf/f0)*exp(j*(-4)*pi*deltf*Xi_esprit(ki)*cos(theta)/c);
%     R3(ki)=Ai_ssa(ki)*exp(j*pi*Type(kj)/2)*exp(j*(-4)*pi*f0*Xi_ssa(ki)*cos(theta)/c);
%     Z3(ki)=(1+Type(kj)*deltf/f0)*exp(j*(-4)*pi*deltf*Xi_ssa(ki)*cos(theta)/c);
%     for k=1:1:n
%         Es_o(k)=awgn(R1(ki)*(Z1(ki).^k),sr)+Es_o(k);  
%         Es_re_esprit(k)=awgn(R2(ki)*(Z2(ki).^k),sr)+Es_re_esprit(k);  
%         Es_re_ssa(k)=awgn(R3(ki)*(Z3(ki).^k),sr)+Es_re_ssa(k);  
%     end
% end
% Em=ifft(Es_o);
% rcs1(:,sr) = 10.0*log10(abs(Em.^2));
% Em_re_esprit=ifft(Es_re_esprit);
% rcs2(:,sr) = 10.0*log10(abs(Em_re_esprit.^2));
% Em_re_ssa=ifft(Es_re_esprit);
% rcs3(:,sr) = 10.0*log10(abs(Em_re_ssa.^2));
% figure;
% plot(fm,rcs1(:,sr),'k--');
% hold on
% plot(fm,rcs2(:,sr),'b-');
% hold on
% plot(fm,rcs3(:,sr),'r-');
% grid on
% xlabel('Frequency/Hz');
% ylabel('RCS/dBsm');
% title(['1D-GTD',' SNR=',num2str(sr),' Scattering type=',num2str(Type(kj))]);
% legend('Original','Esprit_Recovery','SSM_Recovery','Interpreter','none')%显示下划线'Interpreter','none'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ERROR Analysis
figure;
% subplot(2,1,1);
plot(snr,RMSE1,'*r-',snr,RMSE2,'*b-');
% title(['Position',' Scattering type = ',num2str(Type(kj))]);
xlabel('SNR(dB)');
ylabel('RMSE(dB)');
legend('ESPRIT','SSM');
% grid on;
figure;
% subplot(2,1,2);
plot(snr,RMSE3,'*r-',snr,RMSE4,'*b-');
% title(['Amplitude',' Scattering type = ',num2str(Type(kj))]);
xlabel('SNR(dB)');
ylabel('RMSE(dB)');
legend('ESPRIT','SSM');
% grid on;
% subplot(3,1,3);
% plot(1:30,RMSE3,'*r-');
% title('Scattering Tpye');
% xlabel('SNR/dB');
% ylabel('NaN');
% grid on;
end
