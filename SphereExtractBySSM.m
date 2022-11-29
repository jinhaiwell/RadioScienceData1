%%%%%%SSM 低频可以提取爬行波，高频无法提取
%%%%%%蒙特卡罗测试随机噪声，横坐标表示各点频段0.1GHz-10GHz，纵坐标表示幅值，构建矩阵的行表示频率，列表示SNR
clc
clear 
close all
c = physconst('Lightspeed');
rad = 0.54;
az = 5.0;
el = 20.0;
j=sqrt(-1);
Am=[1 1];
%%%%%%%%%%%%%%%%%%%高频段爬行波几乎没有了
f2L=10010*1e6; 
f2H=12000*1e6;
df2=10e6;
fmH = f2L:df2:f2H;
bw2= f2H-f2L;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%以下低频段
f1L=1; 
f1H=10000e6;
df1=10e6;
fm=f1L:df1:f1H;
ln=length(fm)/20;%带宽步进数
% fq(:,:)=[fm(1:100);fm(101:200);fm(201:300);fm(301:400);fm(401:500);...
%         fm(501:600);fm(601:700);fm(701:800);fm(801:900);fm(901:1000)];
fq(:,:)=[fm(1:ln);fm(ln*1+1:ln*2);fm(ln*2+1:ln*3);fm(ln*3+1:ln*4);fm(ln*4+1:ln*5);...
        fm(ln*5+1:ln*6);fm(ln*6+1:ln*7);fm(ln*7+1:ln*8);fm(ln*8+1:ln*9);fm(ln*9+1:ln*10);...
        fm(ln*10+1:ln*11);fm(ln*11+1:ln*12);fm(ln*12+1:ln*13);fm(ln*13+1:ln*14);fm(ln*14+1:ln*15);...
        fm(ln*15+1:ln*16);fm(ln*16+1:ln*17);fm(ln*17+1:ln*18);fm(ln*18+1:ln*19);fm(ln*19+1:ln*20)];
[M,N]=size(fq);
min_th = -10 ;          % 最小方位角，最小旋转角，弧度=角度*pi/180
max_th = 10 ;           % 最大方位角，最大旋转角
phi = linspace(min_th, max_th, M);
Psi = max_th-min_th; %总角度宽度
dphi=Psi/M; %角度分辨率
winfun=hamming(N);%窗函数
snr=10:5:20;
num=0;
odrn=2; %
P=N/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sr=snr
num=num+1;
for fn=1:M
freq = fq(fn,:);%低频段
fc(fn) = (freq(1)+df1*(N/2));  % carrier frequency center
lambda = c/fc(fn);
bw= freq(end)-freq(1);
ddr1 = c/(2*bw);
dda1 = lambda/(2*Psi*pi/180);
rge1 = ddr1*(-N/2):ddr1:ddr1*(N/2);%距离向
azi1 = dda1*(-N/2):dda1:dda1*(N/2);%方位向
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if fn<=6   %%%%%低频段只能做6次MonteCarlo tests
Xi=zeros(1,odrn);
Ai=zeros(1,odrn);
Em=zeros(1,N);
for test=1:100 %MonteCarlo tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rcs,rcs_n,Ercs(fn,:),Ercs_n(fn,:),pha(fn,:),pha_n(fn,:)] = spherercs(rad,c,freq,az,phi(fn),sr);%nrcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
downrange=fftshift(abs(ifft(Ercs_n(fn,:).*winfun'))); 
downrange=downrange/max(downrange);%
%%%%%%%%%%%%%%%%%extraction peak value
IndMin=find(diff(sign(diff(downrange)))>0)+1;%
IndMax=find(diff(sign(diff(downrange)))<0)+1;% 
[vpk,pos1]=findpeaks(downrange);  %
[vpkMax,pos2]=sort(findpeaks(downrange),'descend');
vpos=pos1(pos2(1:2));%两个峰值的数组序号
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shiftnum=0;        %横坐标不平移
shiftnum=N/2-vpos(1);%横坐标平移到0点
rgx=-ddr1*(N/2-shiftnum-1):ddr1:ddr1*(N/2+shiftnum); % range vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%
vrgx(fn,:)=rgx(pos1(pos2(1:2)));      %两个峰值的原始横坐标
% spcrgx(fn)=rgx(vpos(1));              %镜面反射波的坐标
% spcval(fn)=vpos(1);                   %镜面反射波对应的横坐标
% spcvalue(num,fn)=downrange(spcval(1));%镜面反射波对应的纵坐标
% crprgx(fn)=rgx(vpos(2));              %爬行波的坐标
% crpval(fn)=vpos(2);                   %爬行波对应的横坐标
% crpvalue(num,fn)=downrange(crpval(1));%爬行波对应的纵坐标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Em=Ercs_n(fn,:);
% end
Y0=hankel(Em(1:P),Em(P:end));
[U,S,V] = svd(Y0.');
H=U*S*V';%
% omga=U*sqrt(S);
% gama=sqrt(S)*V';
U1=U(:,1:odrn);
S1=S(1:odrn,1:odrn);
V1=V(1:odrn,:);
omga=U1*sqrt(S1);
gama=sqrt(S1)*V1;
H1=U1*S1*V1;%dominant components of the spectral norm sense
% H2=omga*gama;%H1=H2
omga1=omga(2:P,:);%delete the first row
omga2=omga(1:P-1,:);%delete the last row
A=pinv(omga2'*omga2)*omga2'*omga1;
% om=omga2*A;
gama1=gama(:,2:P);%delete the first column
gama2=gama(:,1:P-1);%delete the last column
A1=gama1*gama2'*pinv(gama2*gama1'); 
C=omga(1,:);
% B=gama(:,1);
omgN=[];
gamN=[];
for ki=1:N
    omgan=C*A^(ki-1);
    omgN=[omgN;omgan];
end
B=pinv(omgN'*omgN)*omgN'*Em.';
% for ki=1:ln
%     gaman=(A^(ki-1))*B;
%     gamN=[gamN,gaman];
% end
% hat_C=Em*gamN'*pinv(gamN*gamN');
[vec,lamda]=eig(A);
egv=diag(lamda);

% A State Identification Method for 1-D Measurements with Gaps
fi=angle(diag(lamda));%refers to the phase of the eigenvalue lamda
Xout1=-fi*c/(j*4*pi*df1);%scattering position estimation [1.2 2.401 4.709 6.1 9.5]
[Xout,posx]=sort(abs(Xout1.'));%
Xi=Xi+Xout;

%%%% State–Space Spectral Estimation of Characteristic Electromagnetic Responses in Wideband Data
psi=vec;
for ki=1:odrn
numr(ki)=(C*psi(:,ki))*((psi(:,ki).')*B);
deno(ki)=egv(ki).^(f1L/df1);
end
ai=(numr)./(deno);
[Am_hat,posa]=sort(abs(ai),'descend');
Ai=Ai+abs(Am_hat); %Am_hat(1)
%%%%%%%%%%%%%%%常规方法估计,误差太大
% for ki=1:odrn
%     sig(ki,:)=exp(-4*pi*freq/c*Xout1(ki));
% end
% Aout1_ssm=pinv(sig*sig.')*sig*Em.';
% Aout=sort(abs(Aout1_ssm))';%先爬行波后镜面波
% Ai=Ai+Aout;
end
Xi=Xi/test;
Ai=Ai/test;
Xip(fn,:)=Xi;
Amp(fn,:)=Ai;
RMSE1(fn) = 10*log10(sqrt(sum(abs(Xip(fn,:)-vrgx(1,:)).^2)/2));
RMSE2(fn) = 10*log10(sqrt(sum(abs(Amp(fn,:)-Am).^2)/2));
% end
end
mag1(num,:)=10*log10(Amp(:,1)');
mag2(num,:)=10*log10(Amp(:,2)');
poshat1(num,:)=(Xip(:,1)');
poshat2(num,:)=(Xip(:,2)');
end
%%%%%%%%%%%%%%%%%%%%%镜面发射波
figure;
plot(fc*10e-10,mag1(1,:),'r--*');
hold on
plot(fc*10e-10,mag1(2,:),'g:o');
hold on
plot(fc*10e-10,mag1(3,:),'b-.x');
xlabel('Frequency(GHz)');
ylabel('Amplitude(dB)')
legend('SNR=10dB','SNR=15dB','SNR=20dB')
figure;
plot(fc*10e-10,poshat1(1,:)-mean(poshat1(1,:)),'r-*');
hold on
plot(fc*10e-10,poshat1(2,:)-mean(poshat1(2,:)),'g-*');
hold on
plot(fc*10e-10,poshat1(3,:)-mean(poshat1(3,:)),'b-*');
xlabel('Frequency(GHz)');
ylabel('Position(m)')
legend('SNR=10dB','SNR=15dB','SNR=20dB')
%%%%%%%%%%%%%%%%%%%%%%爬行波
figure;
plot(fc*10e-10,mag2(1,:),'r-*');
hold on
plot(fc*10e-10,mag2(2,:),'g-*');
hold on
plot(fc*10e-10,mag2(3,:),'b-*');
xlabel('Frequency(GHz)');
ylabel('Amplitude(dB)')
legend('SNR=10dB','SNR=15dB','SNR=20dB')
figure;
plot(fc*10e-10,(poshat2(1,:)+mean(poshat1(1,:))),'r-*');%加上均值进行平移
hold on
plot(fc*10e-10,(poshat2(2,:)+mean(poshat1(2,:))),'g-*');
hold on
plot(fc*10e-10,(poshat2(3,:)+mean(poshat1(3,:))),'b-*');
xlabel('Frequency(GHz)');
ylabel('Position(m)')
legend('SNR=10dB','SNR=15dB','SNR=20dB')
num=0;
%%%%%%%%%%%%%%%%%%提取完成后将参数带入PSF模型计算球的RCS


