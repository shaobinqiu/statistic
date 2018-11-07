function  [n_mu_T, n_str_max]=PF_nfre(all_inf, T, delta_mu, p)
n_mu_T=[];
T_max=max(T);
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
n_str_max=zeros(size(delta_mu,2),size(T,2));
for zz=1:size(delta_mu,2)
    zz
    H=all_inf(:,2)*delta_mu(1,zz)-all_inf(:,3);
    H=-(H-max(H)*ones(size(H,1),1));
    Z_all=[];
    for ii= 1:size(T,2)
        z_i=[];
        for rr=1:size(H,1)
            z=exp(-H(rr,1)/k/T(ii));
            z_i=[z_i;z];
        end
        z_i_normal=z_i.*all_inf(:,4)/(z_i'*all_inf(:,4));
        Z_all=[Z_all   z_i_normal];
    end%%% Z = sum(g_i*z_i)        z_i from $9.7 wangzhicheng
   
%     [x,m]=max(Z_all);
%     for ww=1:size(x,2)
%         if x(1,ww)>p%%%%structure which is more than p and dominant
%             n_str_max(zz,ww)=m(1,ww);
%         else
%             n_str_max(zz,ww)=0;
%         end
%     end
    
    E_n=[];
    for ll=1:size(Z_all,2)
        en=Z_all(:,ll)'*all_inf(:,2);
        E_n=[E_n en];
    end
    n_mu_T=[n_mu_T;E_n];
end
end
