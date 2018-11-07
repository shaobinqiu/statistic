clear
all_inf=load('./files/all_inf.txt');% [index number_of_atom_A  energy  degener]
cp=30/100;% critical probability of phase diagram
frequency=[];
%frequency=load('./files/frequency.txt');
T_max=1000;dT=10;T=[10:dT:T_max];
n_mu_T=[];
dmu=0.001;mu_min=-2;mu_max=-1;mu=[mu_min:dmu:mu_max];% mu=mu(A)-mu(B)
if frequency==[]
    [n_mu_T, n_str_max]=PF_nfre(all_inf, T, mu, cp);
else
    [n_mu_T, n_str_max]=PF_fre(all_inf, T, mu, frequency, cp);
end