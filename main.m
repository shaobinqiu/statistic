clear
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
all_inf=load('~/Desktop/temp/LiuLu/InGaO/scell.txt');% [index number_of_atom_A  energy  degener]
%all_inf=load('../SiGeH16/files/all_DFT.txt');
A='Ga';
B='In';
frequency=[];

cp=30/100;% critical probability of phase diagram
%frequency=load('../Al13/files/frequency.txt');
%frequency=load('../SiGeH16/files/frequency.txt');
T_min=10;
T_max=5000;T=[T_min:(T_max)/200:T_max+T_min];
nm=1000;%the number of mu points
%%%%%%%%%%%%%%%%%diff mu area at diff T
mu_min=[-0.94;-2.8];mu_max=[-0.91;0.6];
mu=zeros(size(T,2),nm+1);
for ii=1:size(mu,1)
    a=ii*[mu_min(2)-mu_min(1)]/size(mu,1)+mu_min(1);
    b=ii*[mu_max(2)-mu_max(1)]/size(mu,1)+mu_max(1);
mu(ii,:)=[a:(b-a)/nm:b];
end
%%%%%%%%%%%same mu area at diff T
% mu_min=-2.5;mu_max=0.5;
% dmu=(mu_max-mu_min)/nm;
% mu=repmat([mu_min:dmu:mu_max],size(T,2),1);% mu=mu(A)-mu(B)
%%%%%%%%%%%
%mu_min=-3.7;mu_max=-2.5;dmu=(mu_max-mu_min)/nm;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);
%mu_min=-4.3;mu_max=-2.5;dmu=(mu_max-mu_min)/nm;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);
%mu_min=-1.6;mu_max=-0.6;dmu=(mu_max-mu_min)/nm;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);
dn=0.005;n_A=[0.01:dn:0.99];%%%%%%% feed ratio

n_mu_T=[];
[n_mu_T, n_str_max, p]=PF_fre(all_inf, T, mu, frequency, cp);
%%%%%%%%%%%%%%%%%%%%%%%%%to get better n_mu_T to plot
n_mu_T_plot=zeros(size(n_mu_T,1),size(n_mu_T,2));
for ii=1:size(n_mu_T,2)
    for jj=1:size(n_mu_T,1)
        for kk=1:size(mu,2)-1
             if mu(ii,jj)>=mu(end,kk) && mu(ii,jj)<=mu(end,kk+1)
                   n_mu_T_plot(kk,ii)=n_mu_T(jj,ii);
             end
        end
    end
end
n_mu_T_plot(1,:)=repmat(0,1,size(n_mu_T_plot,2));
n_mu_T_plot(end,:)=repmat(1,1,size(n_mu_T_plot,2));
    for ii=1:size(n_mu_T_plot,2)
        in=0;
        for jj=2:size(n_mu_T_plot,1)-1
            if in==1
                n_mu_T_plot(jj,ii)=1;
            end
            if n_mu_T_plot(jj,ii)>0.98 & n_mu_T_plot(jj+1,ii)==0
                in=1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
plot_n_mu_T(n_mu_T_plot, mu, T, A, B)

mu1=zeros(size(T,2),size(n_A,2));
err=0;
for zz=1:size(n_A,2)%%%%find the function of mu(n_A, T)
Z_all_inf=[];
for ii=1:size(T,2)  
    for yy=1:size(n_mu_T,1)-1
        if n_A(1,zz)>=n_mu_T(yy,ii) && n_A(1,zz)<=n_mu_T(yy+1,ii)
            mu1(ii,zz)=min(mu(ii,:))+(yy-1)*(mu(ii,2)-mu(ii,1));
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
plot_n_str_max(n_str_max_p1, n_A, T, A, B)
plot_n_str_max(p1, n_A, T,A, B)



% E_H=[];
% for ii=1:size(all_inf,1)
%    E_H=[E_H;all_inf(ii,3)-all_inf(1,3)/13*(13-all_inf(ii,2))-all_inf(164,3)/13*all_inf(ii,2)-all_inf(ii,2)*(-100)];
% end
% figure
% plot(E_H,'-')
% hold on 
% plot([0,164],[0,0],'r-')

