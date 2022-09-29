function [ nUn,nAn] = UAn(X,Y,delta,beta,xi)
    [n,p] = size(X);
    z0 = Y(delta==1);
    N = length(z0);
    I0 = (Y*ones(1,N)>=ones(n,1)*z0');
    alpha = sum(repmat((1-delta).*xi,1,N).*I0)./sum(repmat((1-delta),1,N).*I0);
    alpha = alpha + (alpha==0)*eps;
    Rho = repmat(delta,1,N) + repmat((1-delta).*xi,1,N)./repmat(alpha,n,1);
    expx = exp(X*beta);
    temp1 = expx'*(I0.*Rho);
    w = (expx*ones(1,N))./(ones(n,1)*temp1);
    temp20 = (w.*I0.*Rho)';
    temp2 = temp20 *X;
    temp3 = sum(temp20);
    nUn = -(-sum(X.*(delta*ones(1,p)))+sum(temp2))';
    nAn = -(X'*diag(temp3)*X -temp2'*temp2);
end






