function newvec = apply_UinvT(oldvec,tau_list,n) 

    newvec=zeros(length(oldvec),1);
    for i=1:length(tau_list)
        newvec((i-1)*n+1:i*n)=oldvec((i-1)*n+1:i*n);
    end
    
    for i=1:length(tau_list)
        newvec(end-n+1:end) = newvec(end-n+1:end) - oldvec((i-1)*n+1:i*n)*tau_list(i);
    end
    newvec(end-n+1:end) = newvec(end-n+1:end) + oldvec(end-n+1:end); 

end

