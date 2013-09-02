function [X,resnrm]=rrgmres_iter(A,b,iterations)
%The RRGMRES algorithm for linear discrete ill-posed problems with a 
%    square nonsymmetric Matrix. This version allows the user to specify 
%    the desired number of iterations.
%
%[X,resnrm]=rrgmres_iter(A,b,iterations)
%
%X          =   computed solution to Ax=b, each column correspons to an
%               iteration
%resnrm     =   vector containing residual norm for each iteration
%A          =   from Ax=b
%b          =   from Ax=b
%iterations =   desired number of iterations

s=iterations;
n=size(A,2);
W=zeros(n,s+2);
H=zeros(s+2,s+1);
Q=eye(s+2,s+2);
R=zeros(s+2,s+1);
Rh=zeros(s+2,s);
Qh=eye(s+2,s+2);
beta=zeros(s+2,1);
resnrm=zeros(s,1);
X=zeros(n,s);
W(:,1)=b/norm(b); beta(1)=norm(b); %Begin Arnoldi
for k=1:s+1
    v=A*W(:,k); %Arnoldi
    for j=1:k
        H(j,k)=W(:,j)'*v;
        v=v-H(j,k)*W(:,j);
    end
    H(k+1,k)=norm(v);
    W(:,k+1)=v/H(k+1,k);
    
    %beta=[beta;0];
    
    %[Q,R]=qr(H);
    R(:,k)=H(:,k);
    R(1:k,k)=Q(1:k,1:k)*H(1:k,k); %Apply previous Givens
    
    mu=sqrt(R(k,k)^2+R(k+1,k)^2); %Info for next givens 
    S=R(k+1,k)/mu;
    C=R(k,k)/mu;
    Q(k:k+1,:)=[C S;-S C]*Q(k:k+1,:);
    R(k:k+1,k)=[C S;-S C]*R(k:k+1,k); %Apply new  Givens
    beta(k:k+1)=[C S;-S C]*beta(k:k+1);
    
    if k==1
        continue
    end

    Rh(1:k,k-1)=R(1:k,1:k)*Q(k-1,1:k)';
    Rh(1:k,k-1)=Qh(1:k,1:k)*Rh(1:k,k-1); %Apply previous Givens
    
    mu=sqrt(Rh(k-1,k-1)^2+Rh(k,k-1)^2); %Info for next Givens 
    S=Rh(k,k-1)/mu;
    C=Rh(k-1,k-1)/mu;
    Qh(k-1:k,:)=[C S; -S C]*Qh(k-1:k,:);
    Rh(k-1:k,k-1)=[C S;-S C]*Rh(k-1:k,k-1); %Apply new Givens
    
    beta_res(k-1:k+1,1)=Qh(k-1:k+1,1:k+1)*beta(1:k+1); %compute Qh*beta but do not replace beta
    resnrm(k-1)=norm(beta_res(k:k+1)); %residual is found in last two entries
    
    %Generate x using what follows from properties of RRGMRES
    V=W(:,1:k)*Q(1:k-1,1:k)';
    y=Rh(1:k-1,1:k-1)\beta_res(1:k-1);
    X(:,k-1)=V*y;
end