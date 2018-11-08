clear
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
all_inf=load('./files/all_inf_DFT.txt');% [index number_of_atom_A  energy  degener]
cp=30/100;% critical probability of phase diagram
frequency=[];
frequency=load('./files/frequency.txt');
dT=10;T_max=1000;T=[10:dT:T_max];
n_mu_T=[];
dmu=0.005;mu_min=-1.5;mu_max=-0.8;mu=repmat([mu_min:dmu:mu_max],size(T,2),1);% mu=mu(A)-mu(B)
[n_mu_T, n_str_max, p]=PF_fre(all_inf, T, mu, frequency, cp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'color','white');
n_mu_T= imrotate(n_mu_T,90);
image(n_mu_T,'CDataMapping','scaled')
title('n_{Si}(\Delta\mu, T)')
xlabel('\mu_{Si}-\mu_{Ge}(eV)')
ylabel('Temperature(K)')
colorbar
n_x=5;n_y=5;
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
for ii=1:10
    alpha=(ii-0.5)/10;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dn=0.01;n_A=[0.05:dn:0.95];
mu1=zeros(size(T,2),size(n_A,2));
err=0;
for zz=1:size(n_A,2)
Z_all_inf=[];
for ii=1:size(T,2)  
    for yy=1:size(n_mu_T,1)-1
        if n_A(1,zz)>=n_mu_T(yy,ii) && n_A(1,zz)<=n_mu_T(yy+1,ii)
            mu1(ii,zz)=min(mu(1,:))+(yy-1)*dmu;
        end
    end
    if mu1(ii,zz)==0
        err=err+1;
    end%%%could find n_A
end
end
[n_mu_T1, n_str_max1, p1]=PF_fre(all_inf, T, mu1, frequency, cp);
n_str_max1= imrotate(n_str_max1,90);
p1= imrotate(p1,90);
te=unique(n_str_max1);
n_str_max_p1=zeros(size(n_str_max1,1),size(n_str_max1,2));
for ii=1:size(n_str_max1,1)
    for jj=1:size(n_str_max1,2)
       n_str_max_p1(ii,jj)=find(te==n_str_max1(ii,jj));
    end
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure 
 set(gcf,'color','white');
 image(n_str_max_p1,'CDataMapping','scaled')
 title('predict phase(>30%)')
xlabel('n_{V}')
ylabel('Temperature(K)')
colorbar
n_x1=5;
xl1={};
for ii=1:n_x1
    xl1=[xl1 num2str(roundn((min(n_A)+(ii-1)*(max(n_A)-min(n_A))/(n_x1-1)),-2))]; 
end
set(gca,'XTick',1:(size(n_str_max_p1,2)-1)/(n_x1-1):size(n_str_max_p1,2));
set(gca,'XTicklabel',xl1)
set(gca,'YTick',1:(size(n_str_max_p1,1)-1)/(n_y-1):size(n_str_max_p1,1));
set(gca,'YTicklabel',yl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%