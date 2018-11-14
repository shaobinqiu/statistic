function plot_n_mu_T(n_mu_T, mu, T)
figure
set(gcf,'color','white');
n_mu_T= imrotate(n_mu_T,90);
image(n_mu_T,'CDataMapping','scaled')
title('n_{Si}(\Delta\mu, T)')
xlabel('\mu_{Si}-\mu_{Ge}(eV)')
ylabel('Temperature(K)')
colorbar
n_x=5;n_y=5;%%%number of axis label
xl={};yl={};
for ii=1:n_x
    xl=[xl num2str(min(mu(1,:))+(ii-1)*(max(mu(1,:))-min(mu(1,:)))/(n_x-1))]; 
end
for ii=1:n_y
    yl=[yl num2str(max(T)-(ii-1)*max(T)/(n_y-1))]; 
end
set(gca,'XTick',1:(size(n_mu_T,2)-1)/(n_x-1):size(n_mu_T,2));
set(gca,'XTicklabel',xl)
set(gca,'YTick',1:(size(n_mu_T,1)-1)/(n_y-1):size(n_mu_T,1));
set(gca,'YTicklabel',yl)
line_n=[0.01, 0.5, 0.99];%%%%%%plot the line of n=line_n(ii)
for ii=1:size(line_n,2)
    alpha=line_n(1,ii);
    line_alpha=[];
    for ii=1:size(n_mu_T,1)
        for jj=1:size(n_mu_T,2)-1
            if n_mu_T(ii,jj)<=alpha && n_mu_T(ii,jj+1)>=alpha
                line_alpha=[line_alpha;ii,jj];
            end
        end
    end
    hold on 
    plot(line_alpha(:,2),line_alpha(:,1),'r-','MarkerSize',1)
end
n_mu_T= imrotate(n_mu_T,270);
end