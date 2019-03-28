function plot_n_str_max(n_str_max_p1, n_A, T, A, B)
figure
 set(gcf,'color','white');
 image(n_str_max_p1,'CDataMapping','scaled')
 title('predict phase(>30%)')
xlabel(['n_{' A, '}'])
ylabel('Temperature(K)')
colorbar
n_x1=5;n_y=5;%%%number of axis label
xl1={};
for ii=1:n_x1
    xl1=[xl1 num2str(roundn((min(n_A)+(ii-1)*(max(n_A)-min(n_A))/(n_x1-1)),-2))]; 
end
yl={};
for ii=1:n_y
    yl=[yl num2str(max(T)-(ii-1)*(max(T)-min(T))/(n_y-1))]; 
end
set(gca,'XTick',1:(size(n_str_max_p1,2)-1)/(n_x1-1):size(n_str_max_p1,2));
set(gca,'XTicklabel',xl1)
set(gca,'YTick',1:(size(n_str_max_p1,1)-1)/(n_y-1):size(n_str_max_p1,1));
set(gca,'YTicklabel',yl)
end