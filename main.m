clear
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; k=1.38*10^(-23)/(1.6*10^(-19));
wyt=load('~/Desktop/temp/WYT/SiGe_energy/SiGe10.txt');%1:8 SiH	GeH	SiSi	SiGe	GeGe	Sinum	Genum	Hnum	9:11 deg	 gap	E_DFT
all_inf=[[1:size(wyt,1)]', wyt(:,[6,11,9])];%load('~/Desktop/temp/WYT/SiGe_energy/SiGe10.txt');% [index number_of_atom_A  energy  degener]
%all_inf=load('../SiGeH16/files/all_DFT.txt');
A='Ga';
B='In';
cp=30/100;% critical probability of phase diagram
dn=0.005;n_A=[0.01:dn:0.99];%%%%%%% feed ratio
frequency=[];
%frequency=load('../Al13/files/frequency.txt');
T_min=10;
T_max=2500;T=[T_min:(T_max)/100:T_max+T_min];
nm=1000;%the number of mu points
%%%%%%%%%%%%%%%%%diff mu area at diff T
mu_min=[-1.4;-2.2];mu_max=[-0.8;-0.3];
mu=zeros(size(T,2),nm+1);
for ii=1:size(mu,1)
    a=ii*[mu_min(2)-mu_min(1)]/size(mu,1)+mu_min(1);
    b=ii*[mu_max(2)-mu_max(1)]/size(mu,1)+mu_max(1);
    mu(ii,:)=[a:(b-a)/nm:b];
end
%%%%%%%%%%%same mu area at diff T
% mu_min=-2.5;mu_max=0.5; dmu=(mu_max-mu_min)/nm;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);% mu=mu(A)-mu(B)

n_mu_T=[];
[n_mu_T, n_str_max, p]=PF_fre(all_inf, T, mu, frequency, cp);%%%%%%%cal n(mu,T)
n_mu_T_plot= get_n_plot(n_mu_T,mu);%to get better n_mu_T to plot
plot_n_mu_T(n_mu_T_plot, mu, T, A, B)

% mu1= grep_mu(n_mu_T,mu,n_A,T);
% [n_mu_T1, n_str_max1, p1]=PF_fre(all_inf, T, mu1, frequency, cp);
% n_str_max1= imrotate(n_str_max1,90);
% p1= imrotate(p1,90);
% te=unique(n_str_max1);
% n_str_max_p1=zeros(size(n_str_max1,1),size(n_str_max1,2));
% for ii=1:size(n_str_max1,1)%%%%%%%adjust the color
%     for jj=1:size(n_str_max1,2)
%        n_str_max_p1(ii,jj)=find(te==n_str_max1(ii,jj));
%     end
% end  
% plot_n_str_max(n_str_max_p1, n_A, T, A, B)
% plot_n_str_max(p1, n_A, T,A, B)



