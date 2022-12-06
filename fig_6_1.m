%This script works with shifted Chebyshev polynomials, specifically on [-2,2]
%Note: add path to Chebfun 

clear; clc; 
randn('seed',1);

x=chebfun('x');
f = @(x) exp(-x);
a1=-2; b1=2; 
f2 = chebfun(f,[a1 b1]); 

N=length(f2); 
ak = chebpoly(f2,N); 

n=80;
 
A0=randn(n,n); A1=randn(n,n); A2=eye(n,n);   
b=randn(n,1); btilde=randn(n,1); 
comp_size = n*(length(ak)-1); 
c=zeros(comp_size,1); c(end-n+1:end)=b;
ctilde=zeros(comp_size,1); ctilde(end-n+1:end)=btilde;

mu2s=b1*[-.9,-.7,-.5,-.3,-.1,.1,.3,.5,.7,.9];  
mus=-1./(mu2s); 

eyen = speye(n); eyedn = speye(comp_size-n); 
zeros1 = sparse(n,comp_size-n); zeros2 = sparse(comp_size-n,n); 
pimat = [zeros1, eyen; eyedn, zeros2]; 

tau = 0; 
t0 = 1; t1 = (1/b1)*tau; 
tau_list = [t0, t1]; 
for i=2:(N-2)
   tau_list = [tau_list (2/b1)*tau*tau_list(end) - tau_list(end-1)];  
end
tau_list=-tau_list(2:end); 

U = zeros(comp_size,comp_size); 
for i=1:comp_size
   U(i,i)=1; 
end

for i=1:length(tau_list)
    U((i-1)*n+1:i*n,end-n+1:end) = tau_list(i)*eye(n); 
end

last2=[]; akflip=fliplr(ak); 
for i=1:length(akflip)-1
    if i==1
       	last2=[last2 A1*akflip(i+1)-b1*A2]; %linear in mu
    elseif i==2 
        last2=[last2 A1*akflip(i+1)]; %quadratic in mu
    elseif i==length(akflip)-3 %prev last length - 1
        last2=[last2 A1*akflip(end-2)-A1*akflip(end)];     
    elseif i==length(akflip)-2 %new
        last2=[last2 A1*akflip(end-1)]; 
    elseif i==length(akflip)-1
        last2=[last2 A0+A1*f(tau)-A2*tau];
    else
        last2=[last2 A1*akflip(i+1)];
    end  
end

L = zeros(comp_size,comp_size); 
for i=1:comp_size-n
    L(i,i)=1;
end
for i=2*n+1:comp_size-n
    L(i,i-2*n)=1;
end
L(end-n+1:end,:)=last2; 

%maxit=comp_size;
maxit=300; tol=1.e-10; kcount=0;  

applyPreCon = @(x) pimat*(U\(L\x));
applyPreConT = @(x) L'\(U'\(pimat'*x));

rxold=applyPreCon(c); rzold=applyPreCon(ctilde);
rhoold=1; dx=zeros(comp_size,1); dz=zeros(comp_size,1); 
mat2=zeros(comp_size); x0s=zeros(n,comp_size); 
Pd=A1*ak(1);

xi=1*ones(length(mu2s),1); xiold=1*ones(length(mu2s),1);
alphaold=1; rhoold=1;
dhatz=zeros(comp_size,length(mu2s)); dhatx=zeros(comp_size,length(mu2s)); 
xhat=zeros(comp_size,length(mu2s)); zhat=zeros(comp_size,length(mu2s)); 
xi=1*ones(length(mu2s),1); xiold=1*ones(length(mu2s),1);

for i=1:maxit
    i
    rho=(rzold')*rxold; 
    beta=-rho/rhoold;
    dx=rxold-beta*dx; dz=rzold-conj(beta)*dz; 
    vec = -applyPreCon(apply_M(ak,n,A0,A1,A2,a1,b1,dx,N,Pd)); 
    alpha=rho/((dz')*vec); 
    rx=rxold-alpha*vec;
    rz=rzold-(conj(alpha)*(-apply_MT(ak,n,A0,A1,A2,a1,b1,applyPreConT(dz),N,Pd)));
    
    for k=1:length(mus)
        mu=mus(k); 
        mu2=mu2s(k); 
        xinew(k)=(1-alpha*mu)*xi(k) + (alpha*beta/alphaold)*(xiold(k)-xi(k));
        alphahat=-alpha*(xi(k)/xinew(k));
        betahat=((xiold(k)/xi(k))^2)*beta;
        dhatx(:,k)=(1/xi(k))*rxold - betahat*dhatx(:,k);
        dhatz(:,k)=(1/conj(xi(k)))*rzold - conj(beta)*dhatz(:,k);
        xhat(:,k)=xhat(:,k)+alphahat*dhatx(:,k);
        zhat(:,k)=zhat(:,k)+conj(alphahat)*dhatz(:,k);
        xiold(k)=xi(k); xi(k)=xinew(k); 

        %Just for plotting purposes!
        res(i+1,k) = postprocess(xhat,k,b,mu2s,A0,A1,A2,f,mu);   

        if k==length(mu2s)

            res(i+1,k) = postprocess(xhat,k,b,mu2s,A0,A1,A2,f,mu); 

            %If the 'slowest' one converged 
            if res(i+1,k)<=tol
                kcount = kcount + 1; 
    
                %Check residual of the others
                for l=length(mu2s)-1:-1:1 
                    res(i+1,k) = postprocess(xhat,k,b,mu2s,A0,A1,A2,f,mu); 
                    
                    if res(i+1,k)<=tol
                        kcount = kcount + 1; 
                        
                        %If solutions have all converged, quit 
                        if kcount == length(mu2s)
                            return; 
                        end
                    end
                end
            end
        end
    end

    alphaold=alpha; 
    rxold=rx; rzold=rz; 
    rhoold=rho; 
end


res(1,1)=1; res(1,2)=1; res(1,3)=1; res(1,4)=1;
res(1,5)=1; res(1,6)=1; res(1,7)=1; res(1,8)=1;
res(1,9)=1; res(1,10)=1; 

figure(1)
semilogy([0:1:i]',res(:,1));
hold on
semilogy([0:1:i]',res(:,2));
semilogy([0:1:i]',res(:,3));
semilogy([0:1:i]',res(:,4));
semilogy([0:1:i]',res(:,5));
legendInfo{1} = ['$\mu = $' num2str(mu2s(1))]; 
legendInfo{2} = ['$\mu = $' num2str(mu2s(2))]; 
legendInfo{3} = ['$\mu = $' num2str(mu2s(3))]; 
legendInfo{4} = ['$\mu = $' num2str(mu2s(4))]; 
legendInfo{5} = ['$\mu = $' num2str(mu2s(5))]; 
legend(legendInfo,'interpreter','latex')
xlabel('its')
ylabel('rel res')

figure(2)
semilogy([0:1:i]',res(:,6));
hold on 
semilogy([0:1:i]',res(:,7));
semilogy([0:1:i]',res(:,8));
semilogy([0:1:i]',res(:,9));
semilogy([0:1:i]',res(:,10));
legendInfo{1} = ['$\mu = $' num2str(mu2s(6))]; 
legendInfo{2} = ['$\mu = $' num2str(mu2s(7))]; 
legendInfo{3} = ['$\mu = $' num2str(mu2s(8))]; 
legendInfo{4} = ['$\mu = $' num2str(mu2s(9))]; 
legendInfo{5} = ['$\mu = $' num2str(mu2s(10))]; 
legend(legendInfo,'interpreter','latex')
xlabel('its')
ylabel('rel res')

function res = postprocess(xhat,k,b,mu2s,A0,A1,A2,f,mu)
    n=length(A0); 
    x=(mu)*xhat(1:n,k);   
    A=A0+A1*f(mu2s(k))-mu2s(k)*A2;   
    res=norm(A*x-b)/norm(b);
end

function newvec = apply_M(ak,n,A0,A1,A2,a1,b1,oldvec,N,Pd_mult)
    newvec=zeros(length(oldvec),1); 
    newvec(1:n)=(1/b1)*oldvec(1:n);
    newvec(n+1:n*(length(ak)-1)) = (2/b1)*oldvec(n+1:n*(length(ak)-1)); 
    newvec(end-n+1:end)=(-2/b1)*Pd_mult*(oldvec(end-n+1:end));  
end


function newvec = apply_MT(ak,n,A0,A1,A2,a1,b1,oldvec,N,Pd_mult)
    newvec=zeros(length(oldvec),1); 
    newvec(1:n)=(1/b1)*oldvec(1:n);
    newvec(n+1:n*(length(ak)-1)) = (2/b1)*oldvec(n+1:n*(length(ak)-1)); 
    newvec(end-n+1:end)=(-2/b1)*Pd_mult'*(oldvec(end-n+1:end));  
end

