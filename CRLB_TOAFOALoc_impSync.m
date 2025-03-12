function CRB = CRLB_TOAFOALoc_impSync(senPos, senVel, srcLoc, srcVel, Qn, Qm, varargin)
% CRB = CRLB_TOAFOALoc_impSync(senPos, senVel, srcLoc, srcVel, Qn, Qm, varargin)
%
% CRLB for TOA-FOA Localization With Clock Bias and Drift
%
% Input:
%   senPos:   (N x M), position of sensors
%   senVel    (N x M), velocity of sensors
%   srcLoc:   (N x 1), location of signal source
%   srcVel:   (N x 1), velocity of signal source
%   Qn:       (M x M), covariance matrix of TOA measurements
%   Qm:       (M x M), covariance matrix of FOA measurements
%   varargin: covariance matrix of sensors position error, sensors velocity error
% Output:
%   CRB:      (2*N + 2), estimate of CRLB for source position, velocity, bias and drift.

[N,M] = size(senPos);

if isempty(varargin)
    Q = blkdiag(Qn,Qm);
else
    Q = blkdiag(Qn,Qm);
    Qs = varargin{1};
    Qv = varargin{2};
%     Q_beta = blkdiag(Qs,Qv);


    r = sqrt(sum((srcLoc-senPos).^2,1))';    % range, m
    rr = diag((srcVel-senVel)'*(srcLoc-senPos))./r;  % range rate, m/s
    B1 = [2*diag(r), zeros(M);
          diag(rr),  diag(r)];

    C1 = zeros(2*M,N*M);
    D1 = zeros(2*M,N*M);
    for j = 1:M
        C1(j,(1:N)+N*(j-1)) = -2*(senPos(:,j)-srcLoc)';
        C1(j+M,(1:N)+N*(j-1)) = -(senVel(:,j)-srcVel)';
        D1(j+M,(1:N)+N*(j-1)) = -(senPos(:,j)-srcLoc)';
    end

    Q = blkdiag(Qn,Qm) + B1\C1*Qs*C1'/B1' + B1\D1*Qv*D1'/B1';
end
d_s = zeros(M,M*N);
b_s = zeros(M,M*N);
d_ds = zeros(M,M*N);
for i = 1:M
    alpha = (srcLoc-senPos(:,i))/norm(srcLoc-senPos(:,i));
    beta = (srcVel-senVel(:,i))/norm(srcLoc-senPos(:,i));

    d_x(i,:) = alpha';
    b_x(i,:) = beta' - beta'*alpha*alpha';
    b_v(i,:) = alpha';
    d_s(i,(i-1)*N+1:i*N) = -alpha';
    
    b_s(i,(i-1)*N+1:i*N) = -beta'+ beta'*alpha*alpha';
    
end
d_v = zeros(M,N);
d_d = ones(M,1);
d_b = zeros(M,1);
b_d = zeros(M,1);
b_b = ones(M,1);
K = [d_x, d_v, d_d, d_b;
         b_x, b_v, b_d, b_b];
CRB = inv(K'/Q*K);

	
end