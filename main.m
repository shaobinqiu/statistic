clear
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
all_inf=load('../Al13/files/all.txt');% [index number_of_atom_A  energy  degener]
frequency=[];

cp=30/100;% critical probability of phase diagram
frequency=load('../Al13/files/sp_39.txt');
T_max=2000;T=[10:T_max/100:T_max];
nm=500;%the number of mu points
mu_min=-4.1;mu_max=-1.8;dmu=(mu_max-mu_min)/nm;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);% mu=mu(A)-mu(B)
dn=0.005;n_A=[0.01:dn:0.99];%%%%%%% feed ratio

n_mu_T=[];
[n_mu_T, n_str_max, p]=PF_fre(all_inf, T, mu, frequency, cp);
plot_n_mu_T(n_mu_T, mu, T)
mu1=zeros(size(T,2),size(n_A,2));
err=0;
for zz=1:size(n_A,2)%%%%find the function of mu(n_A, T)
Z_all_inf=[];
for ii=1:size(T,2)  
    for yy=1:size(n_mu_T,1)-1
        if n_A(1,zz)>=n_mu_T(yy,ii) && n_A(1,zz)<=n_mu_T(yy+1,ii)
            mu1(ii,zz)=min(mu(1,:))+(yy-1)*dmu;
        end
    end
    if mu1(ii,zz)==0
        err=err+1;
    end%%%could not find n_A
end
end
[n_mu_T1, n_str_max1, p1]=PF_fre(all_inf, T, mu1, frequency, cp);
n_str_max1= imrotate(n_str_max1,90);
p1= imrotate(p1,90);
te=unique(n_str_max1);
n_str_max_p1=zeros(size(n_str_max1,1),size(n_str_max1,2));
for ii=1:size(n_str_max1,1)%%%%%%%adjust the color
    for jj=1:size(n_str_max1,2)
       n_str_max_p1(ii,jj)=find(te==n_str_max1(ii,jj));
    end
end  
plot_n_str_max(n_str_max_p1, n_A, T)
% E_H=[];
% for ii=1:size(all_inf,1)
%    E_H=[E_H;all_inf(ii,3)-all_inf(1,3)/13*(13-all_inf(ii,2))-all_inf(164,3)/13*all_inf(ii,2)-all_inf(ii,2)*(-100)];
% end
% figure
% plot(E_H,'-')
% hold on 
% plot([0,164],[0,0],'r-')

