% make trimer figures for stratification paper

% created aug 15, 2019

clear;

ifsave = 0;
savefile = 'trimer_ptriangle.epsc';
savefile2 = 'trimer_anglehist.epsc';

z2t = 4*pi/9;
z3t = 4*sqrt(3)/9;


% Make figures look nicer
fs = 22;  % font size
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',fs)
cc = get(gca,'colororder');


% Plot angle distribution, k=1
loadfile = 'Data/trimer_d2_k1_data2.txt';  % data2 has 10x more points
data = load(loadfile);
state = data(:,3);  % 0 (polymer), 1 (triangle)
angle = data(:,4);

th = acos(angle(state==0)); 
nb = 20;  % number of bins for histogram
edges = linspace(pi/3,pi,nb+1);
h3 = figure(3);
clf
ht = histogram(th,edges,'normalization','pdf');
hold on
%plot([pi/3,pi],[1,1]/(2*pi/3),'k--');
hold off
set(gca,'box','off');
xlabel('$\theta$','interpreter','latex');
ylabel('$p(\theta)$','interpreter','latex');

if(ifsave)
    saveas(h3,savefile2);
end


h2=figure(2);
plot([pi/3,pi],[1,1]/(2*pi/3),'k--');
hold on
plot((edges(1:end-1)+edges(2:end))/2, ht.Values,'x-');
%ht = histogram(th,edges,'normalization','pdf');
hold off
set(gca,'box','off');
xlim([pi/3,pi]);
ylim([0 0.6]);
xlabel('$\theta$','interpreter','latex');
ylabel('$p(\theta)$','interpreter','latex');


%return;

% Get data for P(Triangle)

nbin = 10;  % number of bins for calculating mean
klist = [1,2,4,8,16];  % empirical values of k
ptri = NaN(size(klist));
sigs = NaN(size(klist));

ii=1;
for k=klist
    
    loadfile = ['Data/trimer_d2_k',num2str(k),'_data.txt'];
    data = load(loadfile);
    state = data(:,3);  % 0 (polymer), 1 (triangle)
    angle = data(:,4);
    npts = length(state);
    
    ptrilist = zeros(1,nbin);  % holds the means
    for jj=1:nbin
        ptrilist(jj) = mean(state((jj-1)*npts/nbin+1:jj*npts/nbin));
    end
    
    ptri(ii) = mean(ptrilist);
    sigs(ii) = std(ptrilist)/sqrt(nbin);
    

    
    disp(['k = ',num2str(k)]);
    disp(['mean = ',num2str(ptri(ii))]);
    disp(['true = ',num2str(k.^3*z3t ./(k.^3*z3t + k.^2*z2t))])
    disp(['std  = ',num2str(sigs(ii))])
    disp(ptrilist)
    ii = ii+1;
end


% Plot P(Triangle)
h1=figure(1);
clf
k = 0.1:0.1:20;
ptrue = k.^3*z3t ./(k.^3*z3t + k.^2*z2t);
plot([0,k],[0,ptrue],'k--');
hold on
plot(klist,ptri,'x','Linewidth',2,'color',cc(2,:),'MarkerSize',10);
%errorbar(klist,ptri,2*sigs,'x','Linewidth',2,'color',cc(2,:),'MarkerSize',8);
hold off
set(gca,'box','off');
xlabel('$\kappa$','interpreter','latex');
ylabel('P(triangle)','interpreter','latex');
legend({'Analytical';'Empirical'},'Location','southeast');


if(ifsave)
    saveas(h1,savefile);
end


% Return to default values
set(0,'DefaultLineLineWidth','remove')
set(0,'DefaultAxesLineWidth','remove')
set(0,'DefaultAxesFontSize','remove')



