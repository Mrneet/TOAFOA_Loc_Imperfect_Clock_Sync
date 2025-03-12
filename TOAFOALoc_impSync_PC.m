function [sol] = TOAFOALoc_impSync_PC(rNsed, rrNsed, senPos, senVel, Qn, Qm, varargin)
% function [sol,sol1] = TOAFOALoc_impSync_CF(rNsed, rrNsed, senPos, senVel, Qn, Qm)
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
    
    % update weighting matrix
    u_e = x(1:N);
    v_e = x(N+1:2*N);
    b_e = x(2*N+1);
    d_e = x(2*N+2);

    y = [u_e'*u_e; u_e'*v_e; b_e^2; d_e*b_e];

    r = sqrt(sum((u_e-senPos).^2,1))';    % range, m
    rr = diag((v_e-senVel)'*(u_e-senPos))./r;  % range rate, m/s

    B1 = [2*diag(r), zeros(M);
          diag(rr),  diag(r)];
    B2 = V'*B1;
    if isempty(varargin)
        Q = blkdiag(Qn,Qm);
        W2 = inv(B2*Q*B2');
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
        
    
        W2 = inv(V'*(B1'*blkdiag(Qn,Qm)*B1 + C1*Qs*C1' + D1*Qv*D1')*V);
    end
    
end

H = [
    -2*u_e', zeros(1,N), 0, 0;
    -v_e', -u_e', 0, 0;
    zeros(1,2*N), -2*b_e, 0;
    zeros(1,2*N), -d_e, -b_e
    ];
h3 = h1 - A1*x - A2*y;
G3 = -A1 + A2*H;
if isempty(varargin)
    W1 = inv(B1'*blkdiag(Qn,Qm)*B1);
else
    W1 = inv(B1'*blkdiag(Qn,Qm)*B1 + C1*Qs*C1' + D1*Qv*D1');
end
dx = (G3'*W1*G3)\G3'*W1*h3;

sol = x - dx;

end