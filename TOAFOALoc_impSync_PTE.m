function [sol] = TOAFOALoc_impSync_PTE(rNsed, rrNsed, senPos, senVel, Qn, Qm, varargin)
% [sol] = TOAFOALoc_impSync_PTE(rNsed, rrNsed, senPos, senVel, Qn, Qm, varargin)
%
% Input:
%   rNsed:    (M x 1), TOA measurements
%   rrNsed:   (M x 1), FOA measurements
%   senPos:   (N x M), position of sensors
%   senVel    (N x M), velocity of sensors
%   Qn:       (M x M), covariance matrix of TOA measurements
%   Qm:       (M x M), covariance matrix of FOA measurements
%   varargin: covariance matrix of sensors position error, sensors velocity error
% Output:
%   sol:      (2*N + 2), estimate of source position, velocity, bias and drift.

[N,M] = size(senPos);
G11 = [-2*senPos', zeros(M,N), 2*rNsed, zeros(M,1), ones(M,1), zeros(M,1), -ones(M,1), zeros(M,1)];
b1 = rNsed.^2 - sum(senPos.^2)';
G12 = [-senVel', -senPos', rrNsed, rNsed, zeros(M,1), ones(M,1), zeros(M,1), -ones(M,1)];
b2 = rrNsed.*rNsed - diag(senVel'*senPos);
G1 = [G11;G12];
h1 = [b1;b2];

A1 = G1(:,1:2*N+2);
A2 = G1(:,2*N+3:end);
V = null(A2');
h2 = V'*h1;
G2 = V'*A1;

W2 = V'*V;
Q = blkdiag(Qn,Qm);
for itr = 1:2
    x = (G2'*W2*G2)\G2'*W2*h2;

    u_e = x(1:N);
    v_e = x(N+1:2*N);
    b_e = x(2*N+1);
    d_e = x(2*N+2);

%     y = [u_e'*u_e; u_e'*v_e; b_e^2; d_e*b_e];

    r = sqrt(sum((u_e-senPos).^2,1))';    % range, m
    rr = diag((v_e-senVel)'*(u_e-senPos))./r;  % range rate, m/s

    B1 = [2*diag(r), zeros(M);
          diag(rr),  diag(r)];
    B2 = V'*B1;
    if isempty(varargin)
        Q = blkdiag(Qn,Qm);
        
    else
        Qs = varargin{1};
        Qv = varargin{2};
    
        C1 = zeros(2*M,N*M);
        D1 = zeros(2*M,N*M);
        for j = 1:M
            C1(j,(1:N)+N*(j-1)) = -2*(senPos(:,j)-u_e)';
            C1(j+M,(1:N)+N*(j-1)) = -(senVel(:,j)-v_e)';
            D1(j+M,(1:N)+N*(j-1)) = -(senPos(:,j)-u_e)';
        end
        
    
        Q = blkdiag(Qn,Qm) + B1\C1*Qs*C1'/B1' + B1\D1*Qv*D1'/B1';
    end
    W2 = inv(B2*Q*B2');

end


for i = 1:M
    nn = norm(u_e-senPos(:,i));
    alpha = (u_e-senPos(:,i))/nn;
    beta = (v_e-senVel(:,i))/nn;

    d_x(i,:) = alpha';
    b_x(i,:) = beta' - beta'*alpha*alpha';
    b_v(i,:) = alpha';
end

K = [d_x, zeros(M,N), ones(M,1), ones(M,1);
    b_x, b_v, zeros(M,1), ones(M,1)];

rag = sqrt(sum((u_e-senPos).^2,1))';
dn = rag + b_e;
bn = diag((v_e-senVel)'*(u_e-senPos))./rag + d_e;
reconst = [dn;bn];

delta = (K'/Q*K)\K'/Q*([rNsed;rrNsed]-reconst);
sol = x(1:2*N+2) + delta;

end