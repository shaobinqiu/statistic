%%%%%%%%%%%%%%%%%%%%%%%%%to get better n_mu_T to plot
function n_mu_T_plot= get_n_plot(n_mu_T,mu)
zeros(size(n_mu_T,1),size(n_mu_T,2));
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
end