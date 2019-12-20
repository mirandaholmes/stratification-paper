% Make polymer-wall data figures

% created november 8, 2019

clear;

whichfig = 1;   % 1 = timing data, 0 = others
ifsave = 1;

savefile1 = 'polywall_timing.epsc';
savefile2 = 'polywall_nwall.epsc';
savefile3 = 'polywall_rg.epsc';

istart = 1e3;  % start processing data at this point
sig = sqrt(5)/4;


% Make figures look nicer
fs = 22;  % font size
fsleg = 16;  % font size for legend
ms = 12;  % marker size
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',fs)


if(whichfig == 1)
    nn = [10,20,40,80];
    timing = [1530,7554,46632,286406]/60;
    
    h1 = figure(1);
    loglog(nn,timing,'^-','MarkerSize',ms);
    xlabel('N');
    ylabel('minutes');
    p = polyfit(log(nn),log(timing),1);
    f = polyval(p,log(nn));
    hold on
    plot(nn,exp(f)*1.5,'k--');
    hold off
    legend({'timing';['slope = ',num2str(round(p(1)*100)/100)]},...
        'Location','NorthWest');
    
    if(ifsave)
        saveas(h1,savefile1);
    end
    
    return;
end

% ---------  Setup  --------- %

kap0list = [0.0894,0.2,0.4472,1,2.2361,5,11.1803];


% construct kappa for plotting curves
nkap = 100;
kap = logspace(-2,2,nkap);



% ---------  Load data  --------- %

n10wall = NaN(7,1e6);
n10rg = NaN(7,1e6);
n20wall = NaN(7,1e6);
n20rg = NaN(7,1e6);
n40wall = NaN(7,1e6);
n40rg = NaN(7,1e6);

% load n=10 data
filelist = {'PolymerWall_n10_k2_data.txt';...
    'PolymerWall_n10_k3_data.txt';...
    'PolymerWall_n10_k4_data.txt';...
    'PolymerWall_n10_k1_data.txt';...
    'PolymerWall_n10_k5_data.txt';...
    'PolymerWall_n10_k6_data.txt';...
    'PolymerWall_n10_k7_data.txt'};
tic
for ii=1:7
    datafile = ['Data/',filelist{ii}];
    data = load(datafile);
    n10wall(ii,:) = data(:,1)';
    n10rg(ii,:) = data(:,2)';
end
toc
n10wall = n10wall(:,istart:end);
n10rg = n10rg(:,istart:end);
nd = size(n10rg,2);


% load n=20 data
filelist = {'PolymerWall_n20_k4_data.txt';...
    'PolymerWall_n20_k2_data.txt';...
    'PolymerWall_n20_k5_data.txt';...
    'PolymerWall_n20_k1_data.txt';...
    'PolymerWall_n20_k6_data.txt';...
    'PolymerWall_n20_k3_data.txt';...
    'PolymerWall_n20_k7_data.txt'};
tic
for ii=1:7
    datafile = ['Data/',filelist{ii}];
    data = load(datafile);
    n20wall(ii,:) = data(:,1)';
    n20rg(ii,:) = data(:,2)';
end
toc
n20wall = n20wall(:,istart:end);
n20rg = n20rg(:,istart:end);


% load n=40 data
filelist = {'PolymerWall_n40_k4_data.txt';...
    'PolymerWall_n40_k2_data.txt';...
    'PolymerWall_n40_k5_data.txt';...
    'PolymerWall_n40_k1_data.txt';...
    'PolymerWall_n40_k6_data.txt';...
    'PolymerWall_n40_k3_data.txt';...
    'PolymerWall_n40_k7_data.txt'};
tic
for ii=1:7
    datafile = ['Data/',filelist{ii}];
    data = load(datafile);
    n40wall(ii,:) = data(:,1)';
    n40rg(ii,:) = data(:,2)';
end
toc
n40wall = n40wall(:,istart:end);
n40rg = n40rg(:,istart:end);


% load n=20 bending data
filelist = {'PolymerWall_kbend_n20_k2_data.txt';...
    'PolymerWall_kbend_n20_k3_data.txt';...
    'PolymerWall_kbend_n20_k4_data.txt';...
    'PolymerWall_kbend_n20_k1_data.txt';...
    'PolymerWall_kbend_n20_k5_data.txt';...
    'PolymerWall_kbend_n20_k6_data.txt';...
    'PolymerWall_kbend_n20_k7_data.txt'};
tic
for ii=1:7
    datafile = ['Data/',filelist{ii}];
    data = load(datafile);
    n20kbwall(ii,:) = data(:,1)';
    n20kbrg(ii,:) = data(:,2)';
end
toc
n20kbwall = n20kbwall(:,istart:end);
n20kbrg = n20kbrg(:,istart:end);



% ---------  Compute averages vs kappa  --------- %
tic
n10Ewall = NaN(size(kap));
n10Erg = NaN(size(kap));
for jk = 1:nkap
    weights = (kap(jk)./repmat(kap0list',[1,nd])).^n10wall;
    tempwall = sum(weights.*n10wall,2) ./ sum(weights,2);
    temprg = sum(weights.*n10rg,2) ./ sum(weights,2);
    kweights = exp(-(log(kap0list)-log(kap(jk))).^2./2./sig^2);
    n10Ewall(jk) = (tempwall')*(kweights')/sum(kweights);
    n10Erg(jk) = (temprg')*(kweights')/sum(kweights);
end
toc

tic
n20Ewall = NaN(size(kap));
n20Erg = NaN(size(kap));
for jk = 1:nkap
    weights = (kap(jk)./repmat(kap0list',[1,nd])).^n20wall;
    tempwall = sum(weights.*n20wall,2) ./ sum(weights,2);
    temprg = sum(weights.*n20rg,2) ./ sum(weights,2);
    kweights = exp(-(log(kap0list)-log(kap(jk))).^2./2./sig^2);
    n20Ewall(jk) = (tempwall')*(kweights')/sum(kweights);
    n20Erg(jk) = (temprg')*(kweights')/sum(kweights);
end
toc

tic
n40Ewall = NaN(size(kap));
n40Erg = NaN(size(kap));
for jk = 1:nkap
    weights = (kap(jk)./repmat(kap0list',[1,nd])).^n40wall;
    tempwall = sum(weights.*n40wall,2) ./ sum(weights,2);
    temprg = sum(weights.*n40rg,2) ./ sum(weights,2);
    kweights = exp(-(log(kap0list)-log(kap(jk))).^2./2./sig^2);
    n40Ewall(jk) = (tempwall')*(kweights')/sum(kweights);
    n40Erg(jk) = (temprg')*(kweights')/sum(kweights);
end
toc


tic
n20kbEwall = NaN(size(kap));
n20kbErg = NaN(size(kap));
for jk = 1:nkap
    weights = (kap(jk)./repmat(kap0list',[1,nd])).^n20kbwall;
    tempwall = sum(weights.*n20kbwall,2) ./ sum(weights,2);
    temprg = sum(weights.*n20kbrg,2) ./ sum(weights,2);
    kweights = exp(-(log(kap0list)-log(kap(jk))).^2./2./sig^2);
    n20kbEwall(jk) = (tempwall')*(kweights')/sum(kweights);
    n20kbErg(jk) = (temprg')*(kweights')/sum(kweights);
end
toc


% ---------  Plot  --------- %

% NWall
h2 = figure(2);
clf
hold on
semilogx(kap,n10Ewall/10,'-');
semilogx(kap,n20Ewall/20,'-');
semilogx(kap,n40Ewall/40,'-');
semilogx(kap,n20kbEwall/20,'-');

ax = gca;
ax.ColorOrderIndex = 1;
semilogx(kap0list,mean(n10wall/10,2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n20wall/20,2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n40wall/40,2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n20kbwall/20,2),'d','MarkerSize',ms);
set(gca,'xscale','log');
set(gca,'ylim',[0 1]);
xlabel('\kappa');
ylabel('fraction on wall');
hold off
legend({'N=10';'N=20';'N=40';'N=20, $k_{\rm bend}{=}2$'},'Location','NorthWest','interpreter','latex');


% Rg
h3 = figure(3);
clf
hold on
semilogx(kap,n10Erg/sqrt(10-1),'-');
semilogx(kap,n20Erg/sqrt(20-1),'-');
semilogx(kap,n40Erg/sqrt(40-1),'-');
semilogx(kap,n20kbErg/sqrt(20-1)/2,'--');

ax = gca;
ax.ColorOrderIndex = 1;
semilogx(kap0list,mean(n10rg/sqrt(10-1),2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n20rg/sqrt(20-1),2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n40rg/sqrt(40-1),2),'x','MarkerSize',ms);
semilogx(kap0list,mean(n20kbrg/sqrt(20-1)/2,2),'d','MarkerSize',ms);
text(kap(end),n20kbErg(end)/sqrt(20-1)/2,'$\times 2$','interpreter','latex','FontSize',fs,'color',[0.4940    0.1840    0.5560]);
set(gca,'xscale','log');
ylim([0.72 1.05]);

xlabel('\kappa');
ylabel('$R_g/\sqrt{N-1}$','interpreter','latex');   %ylabel('radius of gyration');
hold off
legend({'N=10';'N=20';'N=40';'N=20, $k_{\rm bend}{=}2$'},'Location','SouthEast','interpreter','latex');



if(ifsave)
    saveas(h2,savefile2);
    saveas(h3,savefile3);
end


% % Return to default values
% set(0,'DefaultLineLineWidth','remove')
% set(0,'DefaultAxesLineWidth','remove')
% set(0,'DefaultAxesFontSize','remove')
