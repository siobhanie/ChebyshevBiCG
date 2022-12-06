function [n,b,comp_size,c,ctilde,last2,tau_list,A0,A1,A2,A3] = helm_setup(ak,ak2,N,b1,tau,f,g)

load('small_helm_bicg.mat')
n = length(b); b=b';  
btilde=randn(n,1); 
comp_size = n*(length(ak)-1); 
c=zeros(comp_size,1); c(end-n+1:end)=b;
ctilde=zeros(comp_size,1); ctilde(end-n+1:end)=btilde;

akflip=fliplr(ak); ak2flip=fliplr(ak2);

last2=[];
for i=1:length(akflip)-1
    if i==1
       	last2=[last2 A1*akflip(i+1)+A3*ak2flip(i+1)]; %prev last i==2
    elseif i==2 
        last2=[last2 A1*akflip(i+1)+A3*ak2flip(i+1)+(b1^2/2)*A2]; %prev last i==3
    elseif i==length(akflip)-3 %prev last length - 1
        last2=[last2 A1*akflip(end-2)+A3*ak2flip(end-2)-A1*akflip(end)-A3*ak2flip(end)];     
    elseif i==length(akflip)-2 %new
        last2=[last2 A1*akflip(end-1)+A3*ak2flip(end-1) + (2*tau/b1)*(A1*akflip(end)+A3*ak2flip(end))]; 
    elseif i==length(akflip)-1
        last2=[last2 A0+A1*f(tau)+A3*g(tau)+A2*tau^2];
    else
        last2=[last2 A1*akflip(i+1)+A3*ak2flip(i+1)];
    end  
end

t0 = 1; t1 = (1/b1)*tau; 
tau_list = [t0, t1]; 
for i=2:(N-2)
   tau_list = [tau_list (2/b1)*tau*tau_list(end) - tau_list(end-1)];  
end
tau_list=-tau_list(2:end); 