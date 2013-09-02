function [X,resnrm,iterations]=rrgmres_dp(A,b,discrepancy)
%The RRGMRES algorithm for linear discrete ill-posed problems with a 
%   square nonsymmetric matrix. This version uses the discrepancy 
%   principle to decide when to terminate the iterations.
%       
%[X,resnrm,iterations]=rrgmres_dp(A,b,discrepancy)
%
%X              =   desired solution to Ax=b, each column corresponds to an
%                   iteration
%resnrm         =   vector containing residual norm for each iteration
%iterations     =   number of iterations required to reach stopping condition
%A              =   from Ax=b
%b              =   from Ax=b
%discrepancy    =   stopping condition: resnrm < discrepancy

k=0; resnrm=2*discrepancy;
Q=1; R=[]; Rh=[]; Qh=[1,0;0,1];
W=b/norm(b); beta=norm(b); %Begin Arnoldi
X=zeros(size(A,1),1);
while (resnrm>discrepancy)&&(k<=size(A,1))
    k=k+1;
    
    v=A*W(:,k); %Arnoldi
    for j=1:k
        H(j,k)=W(:,j)'*v; 
        v=v-H(j,k)*W(:,j);
    end
    H(k+1,k)=norm(v);
    W(:,k+1)=v/H(k+1,k);
    
    beta=[beta;0];
    
    %[Q,R]=qr(H);
    R=[R,zeros(k,1);zeros(1,k)];
    R(:,k)=H(:,k);
    R(1:k,k)=Q(1:k,1:k)*H(1:k,k); %Apply previous Givens
    
    mu=sqrt(R(k,k)^2+R(k+1,k)^2); %Info for next givens 
    S=R(k+1,k)/mu;
    C=R(k,k)/mu;
    Q=[Q,zeros(k,1);zeros(1,k),1];
    Q(k:k+1,:)=[C S;-S C]*Q(k:k+1,:);
    R(k:k+1,k)=[C S;-S C]*R(k:k+1,k); %Apply new  Givens
    beta(k:k+1)=[C S;-S C]*beta(k:k+1);
    
    if k==1
        continue
    end

    %[Qh,Rh=qr(R*Q(1:k-1,1:k)');
    Rh=[Rh,zeros(k,1);zeros(1,k-1)];
    Rh(1:k,k-1)=R(1:k,:)*Q(k-1,1:k)';
    Rh(1:k,k-1)=Qh(1:k,1:k)*Rh(1:k,k-1); %Apply previous Givens
    
    mu=sqrt(Rh(k-1,k-1)^2+Rh(k,k-1)^2); %Info for next Givens 
    S=Rh(k,k-1)/mu;
    C=Rh(k-1,k-1)/mu;
    Qh=[Qh,zeros(k,1);zeros(1,k),1];
    Qh(k-1:k,:)=[C S; -S C]*Qh(k-1:k,:);
    Rh(k-1:k,k-1)=[C S;-S C]*Rh(k-1:k,k-1); %Apply new Givens
    
    beta_res=Qh*beta; %compute Qh*beta but do not replace beta
    resnrm=norm(beta_res(k:k+1)); %residual is found in last two entries
    
    %Generate X using what follows from properties of RRGMRES
    V=W*Q(1:k-1,:)';%; V=V(:,1:k-1);
    y=Rh(1:k-1,:)\beta_res(1:k-1);
    X(:,k-1)=V*y;
end
iterations=k-1; %Needed an extra iteration for Arnoldi