%%%%%%%%%%%%%%%%%%%%%%%%%find mu(n,T) from n(mu,T)
function mu1= grep_mu(n_mu_T,mu,n_A,T)
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
end