function [X,resnrm]=sym_rrgmres_iter(A,b,iterations)
%The RRGMRES algorithm for linear discrete ill-posed problems with a 
%   symmetric matrix. This version allows the user to specify the number 
%   of desired iterations.
%
%[X,resnrm]=sym_rrgmres_iter(A,b,s)
%
%X              =   computed solution to Ax=b, each column correspons to an
%                   iteration
%resnrm         =   vector containing residual norm for each iteration
%A              =   from Ax=b
%b              =   from Ax=b
%iterations     =   desired number of iterations

n=size(A,2);
s=iterations;
%Preallocate
alpha=zeros(s+1,1);
gamma=zeros(s+1,1);
X=zeros(n,s);
resnrm=zeros(s,1);
W=zeros(n,s+2);
Q=eye(s+2,s+2);
Qh=eye(s+2,s+2);
beta=zeros(s+2,1);
beta(1)=norm(b);
W(:,1)=b/norm(b);
%MATLAB indicated that sparse indexing might be slow
%{
R=sparse(zeros(s+2,s+1));
Rh=sparse(zeros(s+2,s));
T=sparse(zeros(s+2,s+1));
%}
R=zeros(s+2,s+1);
Rh=zeros(s+2,s);
T=zeros(s+2,s+1);

for k=1:s+1
    %Lanczos
    if k==1
        v=A*W(:,k);
    else
        v=A*W(:,k)-gamma(k-1)*W(:,k-1);
    end
    alpha(k)=W(:,k)'*v;
    v=v-alpha(k)*W(:,k);
    gamma(k)=norm(v);
    W(:,k+1)=v/gamma(k);

    if (k==1)
        %Givens rotation
        T(1,1)=alpha(1);
        T(2,1)=gamma(1);
        R(1:2,1)=T(1:2,1);
        mu=sqrt(R(1,1)^2+R(2,1)^2); 
        S=R(2,1)/mu;
        C=R(1,1)/mu;
        Q(1:2,:)=[C S;-S C]*Q(1:2,:);
        R(1:2,1)=[C*R(1,1)+S*R(2,1);0];
        beta(1:2)=[C S;-S C]*beta(1:2);
        continue
    else
        T(k-1,k)=gamma(k-1);
        T(k,k)=alpha(k);
        T(k+1,k)=gamma(k);
        R(k-1:k+1,k)=T(k-1:k+1,k);
        %Apply previous Givens
        if k==2
            R(k-1:k,k)=Q(k-1:k,k-1:k+1)*T(k-1:k+1,k);
        else
            R(k-2:k,k)=Q(k-2:k,k-1:k+1)*T(k-1:k+1,k);
        end
        %Next givens
        mu=sqrt(R(k,k)^2+R(k+1,k)^2);  
        S=R(k+1,k)/mu;
        C=R(k,k)/mu;
        Q(k:k+1,:)=[C S;-S C]*Q(k:k+1,:);
        R(k:k+1,k)=[C*R(k,k)+S*R(k+1,k);0];
        beta(k:k+1)=[C S;-S C]*beta(k:k+1);
    end
    %[Qh,Rh]=qr(R*Q)
    if k==2
        Rh(k-1:k,k-1)=R(k-1:k,k-1:k)*Q(k-1,k-1:k)';
    else
        Rh(k-2:k,k-1)=R(k-2:k,k-2:k)*Q(k-1,k-2:k)';
    end
    %Apply previous Givens
    if k>2
        Rh(1:k,k-1)=Qh(1:k,1:k)*Rh(1:k,k-1); 
    end
    %Next Givens
    mu=sqrt(Rh(k-1,k-1)^2+Rh(k,k-1)^2);
    S=Rh(k,k-1)/mu;
    C=Rh(k-1,k-1)/mu;
    Qh(k-1:k,:)=[C S; -S C]*Qh(k-1:k,:);
    Rh(k-1:k,k-1)=[C*Rh(k-1,k-1)+S*Rh(k,k-1);0]; 
    %compute Qh*beta but do not replace beta
    beta_res(k-1:k+1)=Qh(k-1:k+1,1:k+1)*beta(1:k+1); 
    
    %Generate X using what follows from properties of RRGMRES
    
    y(:,k-1)=W(:,1:k)*Q(k-1,1:k)';
    for i=1:k-2
        y(:,k-1)=y(:,k-1)-y(:,i)*Rh(i,k-1);
    end
    
    y(:,k-1)=y(:,k-1)/Rh(k-1,k-1);
    if k==2
        X(:,1)=y(:,1)*beta_res(1);
    else
        X(:,k-1)=X(:,k-2)+y(:,k-1)*beta_res(k-1);
    end
    
    resnrm(k-1)=norm(beta_res(k:k+1));
end    