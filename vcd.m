function  beta = vcd ( X, Y , beta_initial ,delta,xi, lambda ,tol,maxiter)

[n,p] = size(X);
a = 3.7;
toler = 1;
iter = 0;
beta_old = beta_initial;
beta_new = beta_initial;
z0 = Y(delta==1);
N = length(z0);
I0 = (Y*ones(1,N)>=ones(n,1)*z0');
alpha = sum(repmat((1-delta).*xi,1,N).*I0)./sum(repmat((1-delta),1,N).*I0);
alpha = alpha + (alpha==0)*eps;
Rho = repmat(delta,1,N) + repmat((1-delta).*xi,1,N)./repmat(alpha,n,1);


while toler>tol && iter<maxiter
    iter = iter +1;

    expx = exp(X*beta_old)/n;
    temp1 = repmat(expx,1,N);
    temp2 = I0.*Rho./repmat(sum(I0.*Rho.*temp1),n,1).*temp1;
    nUn = delta-sum(temp2,2);
    nAn = sum(temp2,2)-sum(temp2.^2,2);
    ld = nUn';
    w = nAn';

    v = zeros(1,p);
    z = zeros(1,p);
    beta_new = zeros(p,1);
    winv = w.^(-1);
    winv(isinf(winv)) = 0;
    r_old= diag(winv)*ld';

    for j=1:p
                
        Xj = X(:,j);
        v(j) = Xj'*diag(w)*Xj/n;
        z(j) = Xj'*diag(w)*r_old/n + v(j)*beta_old(j);
    
        if abs(z(j))<= lambda*(v(j)+1)  
            beta_new(j) = Soft(z(j),lambda)/v(j);
        elseif abs(z(j))<= a*lambda*v(j)
            beta_new(j)=Soft(z(j),a*lambda/(a-1))/(v(j)-1/(a-1));            
        else
            beta_new(j) = z(j)/v(j);
        end
        
        r_new = r_old-(beta_new(j)-beta_old(j))*Xj;
        r_old = r_new;

    end
        
    beta_new = beta_new.*(abs(beta_new)>tol);
    toler = norm(beta_new-beta_old);
    beta_old = beta_new;
    
end

beta = beta_old;
    
