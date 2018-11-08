%%%wanzcP273
function z_i=partitionf(Fi_0,Fre,T)
z_i=[];
h=6.63*10^(-34)/(1.6*10^(-19))/2/3.14; 
k=1.38*10^(-23)/(1.6*10^(-19));
for ii=1:size(Fi_0,1)
    z=exp(-Fi_0(ii,1)/k/T);
     for jj=1:size(Fre,2)
        z=z*(exp(-h*Fre(ii,jj)/(2*k*T))/(1-exp(-h*Fre(ii,jj)/(k*T))));
     end
    z_i=[z_i;z];
end
end
    