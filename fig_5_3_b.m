%This script works with shifted Chebyshev polynomials, specifically on [-a,a]
%Note: add path to Chebfun 

clear; close all; clc; 
randn('seed',1);
 
x=chebfun('x');
f = @(x) sin(x)^2;
g = @(x) cos(x)^2; 
a1=-15; b1=15;  
f2 = chebfun(f,[a1 b1]); 
g2 = chebfun(g,[a1 b1]); 

N=max(length(f2),length(g2));
ak = chebpoly(f2,N);  
ak2 = chebpoly(g2,N); 

maxit=80; 
mu2s=[12]; 
tau=11.25; 
tol = 1.e-9; 

[n,b,comp_size,c,ctilde,last2,tau_list,A0,A1,A2,A3] = helm_setup(ak,ak2,N,b1,tau,f,g); 
 
Ptau = A0+A1*f(tau)+A3*g(tau)+A2*tau^2; 
[L,U,P,Q]=lu(Ptau); 
applyPtau = @(x) Q*(U\(L\(P*x))); 
applyPtauT = @(x) P'*(L'\(U'\(Q'*x))); 

eyen = speye(n); eyedn = speye(comp_size-n); 
zeros1 = sparse(n,comp_size-n); zeros2 = sparse(comp_size-n,n); 
pimat = [zeros1, eyen; eyedn, zeros2]; 
Pd=A1*ak(1)+A3*ak2(1);
applyPreCon = @(x) pimat*(apply_Uinv(apply_Linv(x,tau,b1,N,applyPtau,last2,n),tau_list,n)); 
applyPreConT = @(x) apply_LinvT(apply_UinvT(pimat'*x,tau_list,n),tau,b1,N,applyPtauT,last2,n); 

rxold=c; rzold=ctilde;  
dx=0; dz=0;
xi=1*ones(length(mu2s),1); xiold=1*ones(length(mu2s),1);
alphaold=1; rhoold=1;
dhatz=zeros(comp_size,length(mu2s)); dhatx=zeros(comp_size,length(mu2s)); 
xhat=zeros(comp_size,length(mu2s)); zhat=zeros(comp_size,length(mu2s)); 
  
tt=tic();
for i=1:maxit
    
    i
    rho=(rzold')*rxold;
    beta=-rho/rhoold;
    dx=rxold-beta*dx;
    dz=rzold-conj(beta)*dz;
    vec=apply_M(ak,ak2,n,A0,A1,A2,a1,b1,applyPreCon(dx),N,Pd);
    alpha=rho/((dz')*vec);

    %update the seed system 
    rx=rxold-alpha*(vec);
    rz=rzold-conj(alpha)*(applyPreConT(apply_M(ak,ak2,n,A0,A1,A2,a1,b1,dz,N,Pd')));
    
    for k=1:length(mu2s)

        mu(k)=-1./(-mu2s(k)+tau); 
        xinew(k)=(1-alpha*mu(k))*xi(k) + (alpha*beta/alphaold)*(xiold(k)-xi(k));
        alphahat(k)=-alpha*(xi(k)/xinew(k)); 
        betahat(k)=((xiold(k)/xi(k))^2)*beta;
        dhatx(:,k)=(1/xi(k))*rxold - betahat(k)*dhatx(:,k);
        xhat(:,k)=xhat(:,k)+alphahat(k)*dhatx(:,k);

        %%If you are interested in the adjoint solution also
        %dhatz(:,k)=(1/conj(xi(k)))*rzold - conj(beta)*dhatz(:,k);
        %zhat(:,k)=zhat(:,k)+conj(alphahat(k))*dhatz(:,k);
        
        xiold(k)=xi(k); xi(k)=xinew(k);  
        [res,~] = postprocess(xhat(:,k),k,b,applyPreCon,mu2s,A0,A1,A2,A3,f,g,mu(k)); 
        res_mat(i+1,k)=res; 
        thist(i)=toc(tt);

        %If the 'slowest' one converged 
        if res<=tol
            return; 
        end        
        
    end
    
    alphaold=alpha; 
    rxold=rx; rzold=rz; 
    rhoold=rho;
end

semilogy(thist,res_mat(2:end),'ro')
legendInfo{1} = ['$\mu = $' num2str(mu2s(1))]; 
legend(legendInfo,'interpreter','latex')
xlabel('CPU sec')
ylabel('rel res')

function [res,xreturn] = postprocess(xhat,k,b,applyPreCon,mu2s,A0,A1,A2,A3,f,g,mu)
    A=A0+A1*f(mu2s(k))+A3*g(mu2s(k))+A2*mu2s(k)^2; 
    n=length(A0); 
    x2(:,k)=applyPreCon(xhat);
    xreturn=mu*x2(1:n,k); 
    res=norm(A*xreturn-b)/norm(b); %returns norm of the rel res
end
