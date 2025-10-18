clc
clear
close
addpath('CEEMDAN\')
fs=4;%
Ts=1/fs;
STA=1; %
[data3,fs1]=audioread('003.mp3');
% [data4,fs1]=audioread('008.mp3');

% [data5,fs1]=audioread('005.mp3');
% [data6,fs1]=audioread('006.mp3');
% data1=xlsread ('nps.xlsx');
% data1=data1(2:4001,:);
% X=data1;
% data2=data2(1:409600);
data3=data3(1501:5500);
X=data3;
L=length(X);%
t=(0:L-1)*Ts;%
%% CEEMDAN分解
Nstd = 0.2;
NR = 500;
MaxIter = 5000;
[modes,~]=ceemdan(X',0.2,500,5000);



%% 绘图
figure('Position',[100,10,300,700]);
imfn=modes;
n=size(imfn,1);
subplot(n+1,1,1);
plot(t,X); 
ylabel('Original Signal ','fontsize',12,'fontname','Times New Roman');

for n1=1:n
    subplot(n+1,1,n1+1);
    plot(t,modes(n1,:));
    ylabel(['IMF' int2str(n1)]);
end
xlabel('Time\itt/h','fontsize',12,'fontname','Times New Roman');

[a b]=size(modes);   
figure
for i=1:a
    subplot(a,1,i);
    plot(modes(i,:));
    ylabel (['IMF ' num2str(i)],'fontname','Times New Roman');
    set(gca,'Fontname','Times New Roman');
    set(gca,'XLim',[0 4000],'LineWidth',1);
    % axis([0,2^11,-inf,+inf]);
end;
 xlabel('Sample points','fontname','Times New Roman');
% u=flipud(modes);

%% 
dim = 2;   %   dim：1or2)
tau = 1;   %
for i = 1:n
	x=modes(i,:);%
    r = 0.15*std(x);  %   r：0.1*Std(data)~0.25*Std(data) )
    % Sample_Entropy(i,:) = SampleEntropy( dim, r, x, tau );
    % Sample_Entropy(i,:) = pec( x, 6, 1 );
    % Sample_Entropy(i,:) = DisEn_NCDF(x, 3, 6,'NCDF',1 );
    Sample_Entropy(i,:) = VS_Link_Dispersion_Entropy( x, 3, 6, 1,3 );

end

%% 
%matlab k-means
[idx, c] = kmeans(Sample_Entropy, 3); ,
% [idx, c] = KMeans(Sample_Entropy, 3)
Co_IMF1 = sum(modes(find(idx==2),:));  %
Co_IMF2 = sum(modes(find(idx==3),:));   %
Co_IMF3 = sum(modes(find(idx==1),:));   %

% Co2_IMF1 = sum(modes(find(U==2),:));  %
% Co2_IMF2 = sum(modes(find(U==3),:));   %
% Co2_IMF3 = sum(modes(find(U==1),:));   %


Co_IMF2=modes(4,:);
figure
subplot(3,1,1);
plot(Co_IMF1);ylabel('Co-IMF1','Fontname','Times New Roman');hold on
set(gca,'XLim',[0 4000],'LineWidth',1.5);
set(gca,'Fontname','Times New Roman');
set(gca,'xticklabel',[]);

subplot(3,1,2);
plot(Co_IMF2);ylabel('Co-IMF2','Fontname','Times New Roman');hold on
set(gca,'XLim',[0 4000],'LineWidth',1.5);
set(gca,'Fontname','Times New Roman');
set(gca,'xticklabel',[]);

subplot(3,1,3);
plot(Co_IMF3);ylabel('Co-IMF3','Fontname','Times New Roman');hold on
set(gca,'XLim',[0 4000],'LineWidth',1.5);
set(gca,'Fontname','Times New Roman');
 xlabel('Sample points','fontname','Times New Roman');

save CEEVMD_IMF_GLDE_SHIP2_2.mat Co_IMF1
save CEE_IMF_GLDE_SHIP2_2.mat Co_IMF2 Co_IMF3


%% 

maxAlpha=2500; %compactness of mode  
tau=0;%time-step of the dual ascent   
tol=1e-6; %tolerance of convergence criterion;  
stopc=4;%the type of stopping criteria 
%fs=125; % sampling frequency 
% fs=125;
T = length(Co_IMF1);% 
t = (1:T)/T;
omega_freqs = t-0.5-1/T;%
f_hat=fftshift(fft(Co_IMF1));
[svmddata,u_hat,omega]=svmd(Co_IMF1,maxAlpha,tau,tol,stopc);
Co_data = [svmddata;Co_IMF2;Co_IMF3]; %
save CEEMSVMD_Co_data_10_SHIP1_GLDE.mat Co_data

%% 
figure;
vmddata=vmddata';
[a b]=size(vmddata);   %a为模态数

for i=1:a
    subplot(a,1,i);
    plot(vmddata(i,:));
    ylabel (['SVMD-IMF ' num2str(i)],'fontname','Times New Roman');
    set(gca,'Fontname','Times New Roman');
    set(gca,'XLim',[0 4000],'LineWidth',1);
    % axis([0,2^11,-inf,+inf]);
end;
 xlabel('Sample points','fontname','Times New Roman');

imfn=vmddata;
n=size(imfn,1);
subplot(n+1,1,1);
plot(t,X); 
ylabel('Original Signal ','fontsize',12,'fontname','Times New Roman');


for n1=1:n
    subplot(n+1,1,n1+1);
    plot(t,imfn(n1,:));
    ylabel(['VMD-IMF' int2str(n1)]);
end
xlabel('时间\itt/h','fontsize',12,'fontname','宋体');

