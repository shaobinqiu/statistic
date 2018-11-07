function [n_mu_T, n_str_max]=PF_fre(all_inf, T, mu, frequency, cp);
n_mu_T=[];
T_max=max(T);
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
n_str_max=zeros(size(delta_mu,2),size(T,2));
for zz=1:size(delta_mu,2)
    zz
    H=all(:,2)*delta_mu(1,zz)-all(:,3);
    H=H-max(H)*ones(size(H,1),1);
    H=-H;
    Z_all=[];
    for ii= 1:size(T,2)
        z_i=partitionf(H,Fre,T(ii));
        Z_all=[Z_all   z_i];
    end%%% Z = sum(g_i*z_i)        z_i from $9.7 wangzhicheng
    Z_all_norm=Z_all;
    for xx=1:size(Z_all,2)
        Z_all_norm(:,xx)=Z_all(:,xx).*degener/(Z_all(:,xx)'*degener);
    end%normalization
    [x,m]=max(Z_all_norm);
    for ww=1:size(x,2)
        if x(1,ww)>0.3%%%%structure which is more than 30% and dominant
            n_str_max(zz,ww)=m(1,ww);
        else
            n_str_max(zz,ww)=0;
        end
    end
    E_n=[];
    for ll=1:size(Z_all,2)
        en=Z_all(:,ll)'*(all(:,2).*degener)/(Z_all(:,ll)'*degener);
        E_n=[E_n en];
    end
    n_mu_T=[n_mu_T;E_n];
end
end
