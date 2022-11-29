%%%Figure 5, Table 2
close all
clear 
clc
Type=[1 0.5 0 -0.5 -1];
Bw=2e9;   
f0=10e9;
fm=f0+Bw;
n=512;
fm=linspace(f0,fm,n);
deltf=Bw/n;
c=3e8;
snr=10:10:30;

M=5;
R0=10000;
% for typ=1:M
% for sr=snr
for ki=1:M
    X(ki)=abs(randn*10);
    Y(ki)=abs(randn*10);
    Am(ki)=abs(randn*10);
    phi=(120+10*randn)*pi/180;
    dn(ki)=sqrt(R0^2+X(ki)^2-2*R0*X(ki)*cos(phi));
end
X=sort(X)
Y=sort(Y);
Am=sort(Am)
dn=sort(dn);
dmax=max(dn);
dmin=min(dn);
% Am=[1.1 3.5 5.2 7 8.5];  
% X=[1.2 2.401 4.709 6.1 9.5];
for typ=1:M
for sr=snr
  Xi_esprit=0;
  Ti_esprit=0;
  Ai_esprit=0;
  Xi_ssa=0;
  Ti_ssa=0;
  Ai_ssa=0;
  for test=1:100 
      
  theta=0;
  Es=zeros(1,n);
  for r=1:1:M
    0
    R(r)=Am(r)*exp(j*pi*Type(typ)/2)*exp(j*(-4)*pi*f0*X(r)*cos(theta)/c);
    Z(r)=(1+Type(typ)*deltf/f0)*exp(j*(-4)*pi*deltf*X(r)*cos(theta)/c);
    for k=1:1:n
    Es(k)=R(r)*(Z(r).^k)+Es(k);  
    end
  end
%% esprit
  g=100;  
for m=1:1:g       
    Rtemp(m)=0;
end 

for m=1:1:g     
    for i=1:1:(n-m+1) 
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
Xout1=angle(Pout)/(sqrt(-1)*(-4)*pi*deltf*cos(theta)/c);
Xout=sort(abs(Xout1))';%X=[1.2 2.401 4.709 6.100 9.5]
Xi_esprit=Xi_esprit+Xout;

s_esprit=zeros(n,M);
for ki=1:M
    s_esprit(:,ki)=exp(-4*pi*fm/c*Xout1(ki));
end
Aout1_esprit=pinv(s_esprit'*s_esprit)*s_esprit'*Es';
Aout_esprit=sort(abs(Aout1_esprit))';%A=[1.1 3.5 5.2 7 8.5]
Ai_esprit=Ai_esprit+Aout_esprit;

for r=1:1:M
    Tout_esprit(r)=(abs(Pout(r))-1)*f0/deltf;
end
Ti_esprit=Ti_esprit+Tout_esprit;  %O=[0.5 -1 0 1 -0.5]
%%%%%%%%%%%%%%%
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
egv=diag(lamda);%
%A State Identification Method for 1-D Measurements with Gaps
fi=angle(diag(lamda));%refers to the phase of the eigenvalue lamda
Xout1_ssa=-fi*c/(sqrt(-1)*4*pi*deltf);%scattering position  [1.2 2.401 4.709 6.100 9.5];
Xout_ssa=sort(abs(Xout1_ssa))';
Xi_ssa=Xi_ssa+Xout_ssa;

hat_T=(abs(egv)'-1)*f0/deltf; 

s_ssa=zeros(n,M);
for ki=1:M
    s_ssa(:,ki)=exp(-4*pi*fm/c*Xout1_ssa(ki));
end
Aout1_ssa=pinv(s_ssa'*s_ssa)*s_ssa'*Es';
Aout_ssa=sort(abs(Aout1_ssa))';
Ai_ssa=Ai_ssa+Aout_ssa;

  end
Xi_esprit=Xi_esprit/test
Ai_esprit=Ai_esprit/test
Ti_esprit=Ti_esprit/test;
Xi_ssa=Xi_ssa/test
Ai_ssa=Ai_ssa/test

RMSE1(sr/10) = sqrt(sum((Xi_esprit-X).^2)/5);
RMSE2(sr/10) = sqrt(sum((Xi_ssa-X).^2)/5);
RMSE3(sr/10) = sqrt(sum((Ai_esprit-Am).^2)/5);
RMSE4(sr/10) = sqrt(sum((Ai_ssa-Am).^2)/5);
%% RCS
    theta=0;
    Es_o(n)=0;
    Es_esprit=zeros(1,n);
    Es_ssa=zeros(1,n);
    for ki=1:1:M
        R1(ki)=Am(ki)*exp(j*pi*Type(typ)/2)*exp(j*(-4)*pi*f0*X(ki)*cos(theta)/c);
        Z1(ki)=(1+Type(typ)*deltf/f0)*exp(j*(-4)*pi*deltf*X(ki)*cos(theta)/c);
        R2(ki)=Ai_esprit(ki)*exp(j*pi*Type(typ)/2)*exp(j*(-4)*pi*f0*Xi_esprit(ki)*cos(theta)/c);
        Z2(ki)=(1+Type(typ)*deltf/f0)*exp(j*(-4)*pi*deltf*Xi_esprit(ki)*cos(theta)/c);
        R3(ki)=Ai_ssa(ki)*exp(j*pi*Type(typ)/2)*exp(j*(-4)*pi*f0*Xi_ssa(ki)*cos(theta)/c);
        Z3(ki)=(1+Type(typ)*deltf/f0)*exp(j*(-4)*pi*deltf*Xi_ssa(ki)*cos(theta)/c);
        for k=1:1:n
            Es_o(k)=awgn(R1(ki)*(Z1(ki).^k),sr)+Es_o(k);  
            Es_esprit(k)=awgn(R2(ki)*(Z2(ki).^k),sr)+Es_esprit(k);  
            Es_ssa(k)=awgn(R3(ki)*(Z3(ki).^k),sr)+Es_ssa(k);  
        end
    end
    Em=ifft(Es_o);
    rcs1(:,sr) = 10.0*log10(abs(Em.^2));
    Em_esprit=ifft(Es_esprit);
    rcs2(:,sr) = 10.0*log10(abs(Em_esprit.^2));
    Em_ssa=ifft(Es_ssa);
    rcs3(:,sr) = 10.0*log10(abs(Em_ssa.^2));
    figure;
    plot(fm*1e-9,rcs1(:,sr),'r-');
    hold on
    plot(fm*1e-9,rcs2(:,sr),'b--');
    hold on 
    plot(fm*1e-9,rcs3(:,sr),'g--');
%     grid on
    xlabel('Frequency(GHz)');
    ylabel('RCS(dBsm)');
%     title(['1D-GTD',' SNR=',num2str(sr)]);
    legend('Original','ESPRIT-Recovery','SSM-Recovery')

end
%% error
figure;
% subplot(2,1,1);
plot(snr,RMSE1,'*r-');
hold on
plot(snr,RMSE2,'db-');
title('Position');
% title(['Position','Scattering type = ',num2str(Type(kj))]);
xlabel('SNR(dB)');
ylabel('RMSE(m)');
legend('ESPRIT','SSM')
% grid on;
figure;
% subplot(2,1,2);
plot(snr,10*log10(RMSE3),'*r-');
hold on
plot(snr,10*log10(RMSE4),'db-');
title('Amplitude');
% title(['Amplitude','Scattering type = ',num2str(Type(kj))]);
xlabel('SNR(dB)');
ylabel('RMSE(dB)');
legend('ESPRIT','SSM')
% grid on;
end
