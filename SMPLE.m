function [S,beta_S] = SMPLE(X,Y,delta,xi,beta_ini,k)


maxiter = 10;
u = 1;
[n,p]= size(X);
gamma = zeros(p,1);
beta_old = beta_ini;
tol = 1;
iter = 1;
while (iter<=maxiter) && (tol>1E-3)              
    [s,An] = UAn(X,Y,delta,beta_old,xi);         
    ele = (-diag(An)).^(-1);
    solveW = diag(ele);
    add = solveW * s;
    gamma = beta_old + 1/u*add;
    r = gamma.*sqrt(-diag(An));
    tmp = sort(abs(r),'descend');  
    index = find((abs(r)>=tmp(k))~=0);
    betatmp  = vcd ( X(:,index), Y , zeros(k,1) ,delta,xi, 0 ,1.0E-3,10);
    beta_new = zeros(p,1);
    beta_new(index) = betatmp;
    
    tol = norm(beta_new-beta_old);
    beta_old = beta_new;
    iter = iter + 1;


end
 
S = find(beta_old~=0);
beta_S = beta_old(S);

return