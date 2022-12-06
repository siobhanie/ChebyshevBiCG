function newvec = apply_Linv_flex(oldvec,tau,b1,N,Ptau,last2,n,tol)

    newvec=zeros(length(oldvec),1); 
    newvec(1:n)=oldvec(1:n); 
    newvec(n+1:2*n)=oldvec(n+1:2*n) + (2*tau/b1)*newvec(1:n); 
    for i=2:(N-2)
        newvec(i*n+1:(i+1)*n) = oldvec(i*n+1:(i+1)*n) + (2*tau/b1)*newvec((i-1)*n+1:i*n) - newvec((i-2)*n+1:(i-1)*n);
    end
    summation = -last2(:,1:end-n)*newvec(1:end-n); 
    newvec(end-n+1:end)=Ptau(oldvec(end-n+1:end) + summation ,tol);
end