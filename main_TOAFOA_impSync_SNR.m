% TOA+FOA Localization with imperfect synchronization
%
% Performance comparison when the TOA+FOA is availble by using
% - Correction From Constraints: TOAFOALoc_impSync_PC.m
% - Correction From Taylor Expansion： TOAFOALoc_impSync_PTE.m
%
% The code can reproduce the results in Figs. 1-4 of the reference paper.
%
% Number of scenarios is set to 20.
%
% Reference: Tang, B., Zou, Y., Yang, Y., Yang, X., Cong, X., & Sun, Y. (2024). 
% Asymptotically Efficient Solutions for TOA-FOA Localization 
% With Clock Bias and Drift. IEEE Sensors Journal.
%
%
%       Copyright (C) 2024
%       Sichuan University      
%       yimaosun@scu.edu.cn
%
clear; close all;
rng('default');
warning off

simMode = 'nse';

%% parameter settings
% configurations
N = 3; M = 10;
numScen = 20;

senPoss = (rand(3,M,numScen)-0.5)*800;
senVels = (rand(3,M,numScen)-0.5)*30;
srcLocs = (rand(N,numScen)-0.5)*800;
srcVels = (rand(N,numScen)-0.5)*20;

% clock bias and clock drift
clk_bias = rand*150; % m
clk_drift = rand*15; % m/s

% noise power in dB
nsePwrdB = -60:10:40;
k = 0.1;

scaling = 0.01;

% noise generating
ensembleNum = 1000;
nseRg = randn(M,ensembleNum);
nseRg = nseRg - mean(nseRg,2);
nseRR = randn(M,ensembleNum);
nseRR = nseRR - mean(nseRR,2);

totalTime = zeros(8,1);
%% Mont Carlo simulations
for is = 1:numScen
    senPos = senPoss(:,:,is);
    senVel = senVels(:,:,is);
    srcLoc = srcLocs(:,is);
    srcVel = srcVels(:,is);
    y0 = [srcLoc;srcVel;clk_bias;clk_drift];

    % true measurements
    ro = sqrt(sum((srcLoc-senPos).^2,1))';    % range, m
    rro = diag((srcVel-senVel)'*(srcLoc-senPos))./ro;  % range rate, m/s

    for in = 1:length(nsePwrdB)
        disp(['Scen: ',num2str(is),', Noise: ', num2str(nsePwrdB(in)), ' dB'])
        nsePwr = 10^(nsePwrdB(in)/10);
        Qn = nsePwr*eye(M);
        Qm = k*nsePwr*eye(M);


        CRB = CRLB_TOAFOALoc_impSync(senPos, senVel, srcLoc, srcVel, Qn, Qm);
        crlb_u_is(in,is) = trace(CRB(1:N,1:N));
        crlb_v_is(in,is) = trace(CRB(N+1:2*N,N+1:2*N));
        crlb_d_is(in,is) = CRB(2*N+1,2*N+1);
        crlb_b_is(in,is) = CRB(end,end);

       

        for im = 1:ensembleNum
            rNsed = ro + clk_bias + sqrtm(Qn)*nseRg(:,im);
            rrNsed = rro + clk_drift + sqrtm(Qm)*nseRR(:,im);

            numAlg = 1;
            % Projection
            tic
            [solCF] = TOAFOALoc_impSync_PC(rNsed, rrNsed, senPos, senVel, Qn, Qm);
            sol(:,numAlg,im) = solCF;
            totalTime(1) = totalTime(1) + toc;
            numAlg = numAlg + 1;

            % Projection with One-Step GN
            tic
            solGN = TOAFOALoc_impSync_PTE(rNsed, rrNsed, senPos, senVel, Qn, Qm);
            sol(:,numAlg,im) = solGN;
            totalTime(2) = totalTime(2) + toc;

            

        end

        for ia = 1:numAlg
            mse_u_is(in,ia,is) = mean(sum((sol(1:N,ia,:)-srcLoc).^2,1),3);
            mse_v_is(in,ia,is) = mean(sum((sol(N+1:2*N,ia,:)-srcVel).^2,1),3);
            mse_d_is(in,ia,is) = mean((sol(2*N+1,ia,:)-clk_bias).^2,3);
            mse_b_is(in,ia,is) = mean((sol(2*N+2,ia,:)-clk_drift).^2,3);

            bias_u_is(in,ia,is) = norm(mean(sol(1:N,ia,:)-srcLoc,3));
            bias_v_is(in,ia,is) = norm(mean(sol(N+1:2*N,ia,:)-srcVel,3));
            bias_d_is(in,ia,is) = abs(mean(sol(2*N+1,ia,:)-clk_bias,3));
            bias_b_is(in,ia,is) = abs(mean(sol(2*N+2,ia,:)-clk_drift,3));

            
        end
    end
end

%% result combination and plotting
% average different configures
crlb_u = mean(crlb_u_is,2);
crlb_v = mean(crlb_v_is,2);
crlb_d = mean(crlb_d_is,2);
crlb_b = mean(crlb_b_is,2);


mse_u = mean(mse_u_is,3);
mse_v = mean(mse_v_is,3);
mse_d = mean(mse_d_is,3);
mse_b = mean(mse_b_is,3);

bias_u = mean(bias_u_is,3);
bias_v = mean(bias_v_is,3);
bias_d = mean(bias_d_is,3);
bias_b = mean(bias_b_is,3);

load('plotConfigDefault.mat');
% names = {'Proj-C','Proj-TE','MLE','SDP','WLS','TSWLS','SDP-D'};
names = {'Proj-C','Proj-TE'};

switch simMode
    case 'nse'
        xlabtext = '10log(\sigma^2(m^2))';
        xdata = nsePwrdB;
        fileName = ['mat_TOAFOA_impSync_Noise_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
        yl_mse = [-50, 100; -50,80; -50, 100; -50,80];
        yl_bias = [-150,100; -180,80; -160,120; -180,100];
        mstr = 'Noise';
    case 'num'
        xlabtext = 'Num. of Sensors';
        xdata = numSens;
        fileName = ['mat_TOAFOA_impSync_Num_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
        yl_mse = [-40,10;
            -83,-40];
        yl_bias = [-85,0;
            -150,-45];
        mstr = 'Num';
end

% MSE
f1 = figure;
subplot(2,1,1);
for ia = 1:numAlg
    plot(xdata, 10*log10(mse_u(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
plot(xdata, 10*log10(crlb_u), Color=colrs{ia+1}, LineStyle="-", LineWidth=1.5, DisplayName='CRLB');

grid on;
legend('show', 'Location', 'northwest','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('10log10(MSE(u) (m^2))',FontSize=11);
ylim(yl_mse(1,:))


% figure;
subplot(2,1,2);
for ia = 1:numAlg
    plot(xdata, 10*log10(mse_v(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
plot(xdata, 10*log10(crlb_v), Color=colrs{ia+1}, LineStyle="-", LineWidth=1.5, DisplayName='CRLB');
grid on;
legend('show', 'Location', 'northwest','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('10log10(MSE(v) (m^2/s^2))',FontSize=11);
ylim(yl_mse(2,:))

f1.Position(2:4) = [300,560,600];


f2 = figure;
subplot(2,1,1);
for ia = 1:numAlg
    plot(xdata, 10*log10(mse_d(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
plot(xdata, 10*log10(crlb_d), Color=colrs{ia+1}, LineStyle="-", LineWidth=1.5, DisplayName='CRLB');
grid on;
legend('show', 'Location', 'northwest','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('10log10(MSE(b) (m^2))',FontSize=11);
ylim(yl_mse(3,:))


% figure;
subplot(2,1,2);
for ia = 1:numAlg
    plot(xdata, 10*log10(mse_b(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
plot(xdata, 10*log10(crlb_b), Color=colrs{ia+1}, LineStyle="-", LineWidth=1.5, DisplayName='CRLB');
grid on;
legend('show', 'Location', 'southeast','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('10log10(MSE(d) (m^2/s^2))',FontSize=11);
ylim(yl_mse(4,:))

f2.Position(2:4) = [300,560,600];


% bias
f3 = figure;
subplot(2,1,1);
for ia = 1:numAlg
    plot(xdata, 20*log10(bias_u(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
grid on;
legend('show', 'Location', 'southeast','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('20log10(Bias(u) (m))',FontSize=11);
ylim(yl_bias(1,:))


% figure;
subplot(2,1,2);
for ia = 1:numAlg
    plot(xdata, 20*log10(bias_v(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
grid on;
legend('show', 'Location', 'southeast','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('20log10(Bias(v) (m/s))',FontSize=11);
ylim(yl_bias(2,:))

f3.Position(2:4) = [300,560,600];


f4 = figure;
subplot(2,1,1);
for ia = 1:numAlg
    plot(xdata, 20*log10(bias_d(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
grid on;
legend('show', 'Location', 'southeast','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('20log10(Bias(b) (m))',FontSize=11);
ylim(yl_bias(3,:))


% figure;
subplot(2,1,2);
for ia = 1:numAlg
    plot(xdata, 20*log10(bias_b(:,ia)), symbs(ia), Color=colrs{ia}, LineWidth=1.5, DisplayName=names{ia}); hold on;
end
grid on;
legend('show', 'Location', 'southeast','fontsize', 11, 'numcolumns', 2);
xlabel(xlabtext,FontSize=11); ylabel('20log10(Bias(d) (m/s))',FontSize=11);
ylim(yl_bias(4,:))

f4.Position(2:4) = [300,560,600];



