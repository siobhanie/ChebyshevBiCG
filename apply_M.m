function newvec = apply_M(ak,ak2,n,A0,A1,A2,a1,b1,oldvec,N,Pd_mult)
    newvec=zeros(length(oldvec),1); 
    newvec(1:n)=(1/b1)*oldvec(1:n);
    newvec(n+1:n*(length(ak)-1)) = (2/b1)*oldvec(n+1:n*(length(ak)-1)); 
    newvec(end-n+1:end)=(-2/b1)*Pd_mult*(oldvec(end-n+1:end)); 
end

