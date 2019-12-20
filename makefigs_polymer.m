% Process data from Polymer sim, and make figures
%
% Created Sept 9, 2019
%
% =======================
% =======================
% TO DO
%  -- 
%

clear;

ifsave = 0;
savefile = 'polymer_hists.epsc';


% Make figures look nicer
fs = 22;  % font size
fsleg = 16;  % font size for legend
ms = 12;  % marker size
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',fs)
%cc = get(gca,'colororder');


% Parameters
n = 6;    % number of spheres
dim = 3;  % dimension of spheres

nignoreS = 1e3;  % ignore the first this many data points (Strat)
nignoreB = 1e3;  % ignore the first this many data points (BD)


% --------------------------------
%     Stratification data
% --------------------------------

%datafile1 = 'PolymerStrat_n6d3_k1p25_a_data.txt'; kap=1.25;  % 10^7 points
%datafile2 = 'PolymerStrat_n6d3_k1p25_b_data.txt'; kap=1.25;  % 10^7 points
% datafile1 = 'PolymerStrat_n6d3_k1p25_c_data.txt';  kap=2;  % 10^7 points
% datafile2 = 'PolymerStrat_n6d3_k1p25_d_data.txt';  kap=2;  % 10^7 points
% data1 = load(datafile1);
% data2 = load(datafile2);
% data = [data1;data2];

%datafile = 'Data/PolymerStrat_s1_data.txt'; kap = 2; % to generate lines
datafile = 'Data/PolymerStrat_s2_data.txt'; kap = 2.2885; % to compare error bars
dataS = load(datafile);

nbin = 8;  % number of bins for calculating mean


% Compute dimension statistics
neqns = dataS(nignoreS:end,1);
neuclid = dim*(dim+1)/2;  % number of euclidean motions (rotations + translations)
dint = dim*n - neuclid - neqns;  % instrinsic dimension of each state

dimlist = min(dint):max(dint);  % all the different dimensions observed
nd = length(dimlist);

% Estimate error bars
dstats = zeros(1,nbin);  % holds the means
npts = nbin*floor(length(neqns) / nbin);  % total number of points to use for statistics
for dd=1:nd
    state = (dint==dd-1);
    for jj=1:nbin
        dstats(dd,jj) = mean(state((jj-1)*npts/nbin+1:jj*npts/nbin));
    end
end
dmean = mean(dstats,2);  % mean time in each dimension
dsig = std(dstats,1,2)/sqrt(nbin);  % 1-sigma error bar for this estimate of mean time

% Display dimension statistics
disp('Statistics:');
for ii=1:nd
    txt = sprintf('dim %d: %.5f +- %.5f. Relative error = %.4f.',...
        dimlist(ii), dmean(ii), dsig(ii),dsig(ii)/dmean(ii));
    disp(txt);
end
dsig0=dsig;
%disp('Full statistics:')
%disp(dstats);
%disp('Ratios/kappa:');
%ratios = dmean(1:end-1)./dmean(2:end)/kap;
%disp(ratios);

% Estimate geometric partition functions at each dimension, setting Z0=1
zg = ones(nd,1);
for ii=2:nd
    zg(ii) = zg(ii-1) * kap * dmean(ii) / dmean(ii-1);
end
disp('Geometric partition functions:');
disp(zg);

% Infer dimension statistics for all values of kappa
nkap = 200;
kapmax = 10;
kk = linspace(0,kapmax,nkap);
iflogscale = 0;

% construct total partition function
eqnlist = dim*n - dimlist - neuclid;
Ztot = zeros(size(kk));
for dd=1:nd
    Ztot = Ztot + kk.^eqnlist(dd).*zg(dd);
end


%return;


% --------------------------------
%     Plot
% --------------------------------
h1 = figure(1);
clf
if(iflogscale)
    set(gca,'yscale','log');
end
hold on
legtxt = cell(nd,1);
for dd=1:nd  % loop through intrinsic dimensions
    pp = kk.^eqnlist(dd)*zg(dd) ./ Ztot;
    if(dd < nd) pp(1) = 0; 
    else pp(1) = 1; end
    legtxt{dd} = [num2str(eqnlist(dd)), ' bonds'];
    if(~iflogscale) plot(kk,pp); end
    if(iflogscale) semilogy(kk,pp); end
end
ax = gca;   % restsart colour order
%ax.ColorOrderIndex = 1;  
%plot(kap,dmean,'kx');
xlabel('$\kappa$','interpreter','latex');
ylabel('probability');
legend(legtxt,'Location','Best','FontSize',fsleg);



% --------------------------------
%     Brownian dynamics data
% --------------------------------

Elist = [3.447689, 4.24483, 5.02205, 5.78599, 6.54042];  % list of E values for each experiment
% #3 is the one we compare for error bars
%Elist = [5.78599, 5.02205];  % list of E values for each experiment
nk = length(Elist);
kaplist = NaN(nk,1);

for jk=3%1:nk   % loop through experiments
    
    bdfile = ['Data/PolymerBD_b',num2str(jk),'_data.txt'];
    %bdfile = ['/Users/mirandaholmes-cerfon/Dropbox/Work/Sampling/BrownianDynamics/Tests/PolymerBD_expt',num2str(jk),'_data.txt'];
    
    data = load(bdfile);
    neqns = data(nignoreB:end,1);
    dint = dim*n - neuclid - neqns;  % instrinsic dimension of each state
    
    % Work out actual kappa
    E = Elist(jk);
    rho = 60;
    morse = @(r) E*(1-exp(-rho*(r-1))).^2 - E;
    zmorse = @(r) exp(-morse(r));
    kaplist(jk) = integral(zmorse,0.9,1+2.5/60);
    
    
    % Estimate error bars
    dstats = zeros(1,nbin);  % holds the means
    npts = nbin*floor(length(neqns) / nbin);  % total number of points to use for statistics
    for dd=1:nd
        state = (dint==dd-1);
        for jj=1:nbin
            dstats(dd,jj) = mean(state((jj-1)*npts/nbin+1:jj*npts/nbin));
        end
    end
    dmean = mean(dstats,2);  % mean time in each dimension
    dsig = std(dstats,1,2)/sqrt(nbin);  % 1-sigma error bar for this estimate of mean time
    
    % Display dimension statistics
    disp(['Brownian Dynamics, Experiment ',num2str(jk)]);
    disp('Statistics:');
    for ii=1:nd
        txt = sprintf('dim %d: %.5f +- %.5f. Relative error = %.4f.',...
            dimlist(ii), dmean(ii), dsig(ii),dsig(ii)/dmean(ii));
        disp(txt);
    end
    
    
    figure(1)
    ax.ColorOrderIndex = 1;  
    plot(kaplist(jk),dmean,'x','MarkerSize',ms,'HandleVisibility','off'); 
    %errorbar(kaplist(jk)*ones(size(dmean)),dmean,2*dsig,'kx','HandleVisibility','off'); 
    % must set colour order by hand, for error bar plot
    set(gca,'ylim',[0,1]);
end


if(ifsave)
    saveas(h1,savefile);
end


% Return to default values
set(0,'DefaultLineLineWidth','remove')
set(0,'DefaultAxesLineWidth','remove')
set(0,'DefaultAxesFontSize','remove')



