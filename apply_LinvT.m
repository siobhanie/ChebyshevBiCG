function newvec = apply_LinvT(oldvec,tau,b1,N,Ptau,last2,n)

    newvec=zeros(length(oldvec),1); 
    newvec(end-n+1:end) = Ptau(oldvec(end-n+1:end)); 
    newvec(end-2*n+1:end-n) = oldvec(end-2*n+1:end-n) - last2(:,end-2*n+1:end-n)*newvec(end-n+1:end);
    newvec(end-3*n+1:end-2*n) = oldvec(end-3*n+1:end-2*n) + (2*tau/b1)*newvec(end-2*n+1:end-1*n) - last2(:,end-3*n+1:end-2*n)*newvec(end-n+1:end);

    for i=3:(N-2)
       newvec(end-(i+1)*n+1:end-i*n) = oldvec(end-(i+1)*n+1:end-i*n) + (2*tau/b1)*newvec(end-(i)*n+1:end-(i-1)*n) ...
            - newvec(end-(i-1)*n+1:end-(i-2)*n) - last2(:,end-(i+1)*n+1:end-i*n)*newvec(end-n+1:end);    
    end
end
