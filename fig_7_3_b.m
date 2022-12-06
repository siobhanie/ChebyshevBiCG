%This script works with shifted Chebyshev polynomials, specifically on [-a,a]
%Note: add path to Chebfun 

clear; clc; close all; 
randn('seed',1);
 
x=chebfun('x');
f = @(x) sin(x)^2;
g = @(x) cos(x)^2; 
a1=-8; b1=8;  
f2 = chebfun(f,[a1 b1]); 
g2 = chebfun(g,[a1 b1]); 

N=max(length(f2),length(g2));
ak = chebpoly(f2,N);  
ak2 = chebpoly(g2,N); 

maxit=30; 
tau=6;
mu1 = 5.7;
mu2 = 5.8; 
mu3 = 5.9;
mu4 = 6.1;
mu5 = 6.2;
mu6 = 6.3;

[n,b,comp_size,c,ctilde,last2,tau_list,A0,A1,A2,A3] = helm_setup(ak,ak2,N,b1,tau,f,g); 
qold=0; r=c;
qtildeold=0; rtilde=ctilde; 

Ptau = A0+A1*f(tau)+A3*g(tau)+A2*tau^2; 
applyPtau = @(x,tol) agmg(Ptau,x,[],tol,n,1); 
applyPtauT = @(x,tol) agmg(Ptau',x,[],tol,n,1);

eyen = speye(n); eyedn = speye(comp_size-n); 
zeros1 = sparse(n,comp_size-n); zeros2 = sparse(comp_size-n,n); 
pimat = [zeros1, eyen; eyedn, zeros2]; 
Pd=A1*ak(1)+A3*ak2(1);
[L,U,P,Q2]=lu(Ptau); 
applyPtau2 = @(x,tol) Q2*(U\(L\(P*x))); 
applyPtauT2 = @(x,tol) P'*(L'\(U'\(Q2'*x))); 

Apde1=A0+A1*f(mu1)+A3*g(mu1)+(mu1^2)*A2;
Apde2=A0+A1*f(mu2)+A3*g(mu2)+(mu2^2)*A2;
Apde3=A0+A1*f(mu3)+A3*g(mu3)+(mu3^2)*A2;
Apde4=A0+A1*f(mu4)+A3*g(mu4)+(mu4^2)*A2;
Apde5=A0+A1*f(mu5)+A3*g(mu5)+(mu5^2)*A2;
Apde6=A0+A1*f(mu6)+A3*g(mu6)+(mu6^2)*A2;

applyPreCon = @(x,tol) pimat*(apply_Uinv(apply_Linv_flex(x,tau,b1,N,applyPtau,last2,n,tol),tau_list,n)); 
applyPreConT = @(x,tol) apply_LinvT_flex(apply_UinvT(pimat'*x,tau_list,n),tau,b1,N,applyPtauT,last2,n,tol); 
applyPreCon2 = @(x,tol) pimat*(apply_Uinv(apply_Linv(x,tau,b1,N,applyPtau2,last2,n),tau_list,n)); 
applyPreConT2 = @(x,tol) apply_LinvT(apply_UinvT(pimat'*x,tau_list,n),tau,b1,N,applyPtauT2,last2,n); 

tol = 1.e-14; tolvec(1) = tol; 
epsilon = 1.e-12; 
tolOuter = 1.e-9; kcount = 0; 

%flexible algorithm 
for k=1:maxit
    k
    
    beta = norm(r); 
    B(k)=beta; 
    gamma = (rtilde'*r)/beta; 
    G(k)=gamma; 
    q=r/beta; 
    qtilde=rtilde/gamma; 
    z = applyPreCon(q,tol); 
    vec0 = apply_M(ak,ak2,n,A0,A1,A2,a1,b1,z,N,Pd); 
    alpha=qtilde'*(vec0);
    A(k)=alpha; 
    
    r=vec0 - alpha*q - gamma*qold; 
    
    vec1 = applyPreConT(apply_M(ak,ak2,n,A0,A1,A2,a1,b1,qtilde,N,Pd'),tol); 
    rtilde = vec1 - alpha*qtilde - beta*qtildeold; 
    
    qold = q; qtildeold = qtilde; 
    Z(:,k)=z; 
    
    if k>1
        for i=2:k
            T(i-1,i)=G(i);
            T(i-1,i-1)=A(i-1);
            T(i,i-1)=B(i); 
        end
    end
    T(k,k)=A(k); 
    
    %flexible bicg solution 
    e1 = zeros(k,1); e1(1)=norm(c); 
    [xm1,ym1] = postprocess(k,T,Z,mu1,tau,e1); 
    [xm2,ym2] = postprocess(k,T,Z,mu2,tau,e1); 
    [xm3,ym3] = postprocess(k,T,Z,mu3,tau,e1);
    [xm4,ym4] = postprocess(k,T,Z,mu4,tau,e1);
    [xm5,ym5] = postprocess(k,T,Z,mu5,tau,e1);
    [xm6,ym6] = postprocess(k,T,Z,mu6,tau,e1);
    deltak(k) = norm(ym6(end));

    res1(k) = norm(Apde1*xm1(1:n)-b)/norm(b);

    %plotting purposes 
    res2(k) = norm(Apde2*xm2(1:n)-b)/norm(b);
    res3(k) = norm(Apde3*xm3(1:n)-b)/norm(b);
    res4(k) = norm(Apde4*xm4(1:n)-b)/norm(b);
    res5(k) = norm(Apde5*xm5(1:n)-b)/norm(b);
    res6(k) = norm(Apde6*xm6(1:n)-b)/norm(b);

    if res1(k) <= tolOuter
        kcount = kcount + 1; 
        res2(k) = norm(Apde2*xm2(1:n)-b)/norm(b);
        res3(k) = norm(Apde3*xm3(1:n)-b)/norm(b);
        res4(k) = norm(Apde4*xm4(1:n)-b)/norm(b);
        res5(k) = norm(Apde5*xm5(1:n)-b)/norm(b);
        res6(k) = norm(Apde6*xm6(1:n)-b)/norm(b);
        if (res2(k) <= tolOuter) && (res3(k) <= tolOuter) && (res4(k) <= tolOuter) && (res5(k) <= tolOuter) && (res6(k) <= tolOuter) 
            kcount = kcount + 3; 
            break; 
        end 
    end

    tol = epsilon/deltak(k); 
    tolvec(k+1) = tol; 

end

figure(1)
semilogy(tolvec)
hold on
semilogy(res1)
semilogy(res2)
semilogy(res3)
semilogy(res4)
semilogy(res5)
semilogy(res6)
xlabel('its')
ylabel('rel res')
legendInfo{1} = ['$tol_i$']; 
legendInfo{2} = ['$\mu = $' num2str(mu1)]; 
legendInfo{3} = ['$\mu = $' num2str(mu2)]; 
legendInfo{4} = ['$\mu = $' num2str(mu3)]; 
legendInfo{5} = ['$\mu = $' num2str(mu4)]; 
legendInfo{6} = ['$\mu = $' num2str(mu5)];
legendInfo{7} = ['$\mu = $' num2str(mu6)];
legend(legendInfo,'interpreter','latex')


function [x,y] = postprocess(k,T,Z,mu,tau,e1) 
    y = (eye(k) + (-mu+tau)*T)\(e1); 
    x = zeros(length(Z(:,1)),1) + Z*y;   
end
