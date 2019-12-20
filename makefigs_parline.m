% make par/line figures for stratification paper

% created august 12, 2019

clear;

whichfig = 2;
% 1 = par/line, no point
% 2 = examples of points
% 3 = print statistics from big run


ifsave = 0;    % save figure

% Make figures look nicer
fs = 22;  % font size
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',fs)


if(whichfig == 1)
    
    savefile = 'parline.epsc';
    
    h1 = figure(1);
    cc = get(gca,'colororder');
    
    xt = linspace(-sqrt(2),sqrt(2),200);
    yt = xt.^2;
    yt2 = 2*ones(size(xt));
    fill([xt,fliplr(xt)],[yt,fliplr(yt2)],cc(1,:));
    hold on
    plot(xt,yt,'-',xt,yt2,'-',[-sqrt(2) sqrt(2)],[2 2],'.','Linewidth',8,'MarkerSize',72);
    hold off
    box off
    axis off
    
    axis([-1.5 1.5 -0.15 2.2]);
    axis equal
    
    fs1 = 48;
    fs2 = 48;
    
    % add labels for lines
    %t1 = text(0,1.85,'y=2','color',cc(3,:),'FontSize',fs1,'HorizontalAlignment','center');
    %t2 = text(0.05,0.3,'y=x^2','color',cc(2,:),'FontSize',fs1,'HorizontalAlignment','center');%,'Rotation',0);
    t1 = text(0,2.03,'y=2','color',cc(3,:),'FontSize',fs1,'HorizontalAlignment','center','VerticalAlignment','bottom');
    t2 = text(0.06,-0.00,'y=x^2','color',cc(2,:),'FontSize',fs1,'HorizontalAlignment','center','VerticalAlignment','top');%,'Rotation',0);
    t3 = text(-sqrt(2),2.03,{'$(-\sqrt{2},2)$'},'color',cc(4,:),'interpreter','latex','FontSize',fs2,'HorizontalAlignment','center','VerticalAlignment','bottom');
    t4 = text(sqrt(2),2.03,{'$(\sqrt{2},2)$'},'color',cc(4,:),'interpreter','latex','FontSize',fs2,'HorizontalAlignment','center','VerticalAlignment','bottom');
    
    
    if(ifsave)
        saveas(h1,savefile);
    end
end

if(whichfig == 2)
    
    savefile = 'parline_points.epsc';
    
    numpl = 500;  % number of points to plot
    
    pts = load('Data/parline1_pts.txt');
    labels = load('Data/parline1_labels.txt');
    pts = pts(1:numpl,:);
    labels = labels(1:numpl,:);
    
    % Constants labeling type of function: equation or inequality
    cEq = 1; cIn = 2;
    
    iint = find(ismember(labels,[cIn,cIn],'rows')); % interior
    ipar = find(ismember(labels,[cEq,cIn],'rows')); % parabola
    ilin = find(ismember(labels,[cIn,cEq],'rows')); % line
    icor = find(ismember(labels,[cEq,cEq],'rows'));   % corners
    
    npts = size(pts,1);
    
    pts_int = pts(iint,:);
    pts_par = pts(ipar,:);
    pts_lin = pts(ilin,:);
    pts_cor = pts(icor,:);
    
    x_int = pts_int(:,1);
    y_int = pts_int(:,2);
    x_par = pts_par(:,1);
    y_par = pts_par(:,2);
    x_lin = pts_lin(:,1);
    y_lin = pts_lin(:,2);
    x_cor = pts_cor(:,1);
    y_cor = pts_cor(:,2);
    
    % Calculate ratios in each state
    npar = length(x_par);
    nint = length(x_int);
    ncor = length(x_cor);
    nlin = length(x_lin);
    
    xt = linspace(-sqrt(2),sqrt(2),200);
    yt = xt.^2;
    yt2 = 2*ones(size(xt));
    
    h1 = figure(2);
    clf
    %plot(xt,yt,'k:',xt,yt2,'k:','Linewidth',1);  % theory
    hold on
    %plot(x_int,y_int,'*',x_par,y_par,'*',x_lin,y_lin,'*',x_cor,y_cor,'*');
    ms = 10;  % marker size
    small = 1e-2 * 6;
    %xrand = small*normrnd(0,1,size(x_cor));
    %yrand = small*normrnd(0,1,size(x_cor));
    xrand = small*cos(1:ncor)';
    yrand = small*sin(1:ncor)';
    plot(x_int,y_int,'.',x_par,y_par,'.',x_lin,y_lin,'.',x_cor+xrand,y_cor+yrand,'.','MarkerSize',ms);
    hold off
    box off
    axis off
    axis([-1.5 1.5 -0.15 2.2]);
    axis equal
    
    if(ifsave)
        saveas(h1,savefile);
    end
end

if(whichfig == 3)
    
    savefile = 'parline_hists.epsc';
    r = 0.5;  % for saving
    %r=1;   % for testing
    pos = [25 10 50*r 15*r];
    
    % Make figures look nicer
    fs = 28;  % font size
    set(0,'DefaultAxesFontSize',fs)
    
    nbox = 30;  % number of boxes in histogram
    ms = 12;  % marker size
    
    pts = load('Data/parline1_pts.txt');
    labels = load('Data/parline1_labels.txt');
    
    % Constants labeling type of function: equation or inequality
    cEq = 1; cIn = 2;
    
    iint = find(ismember(labels,[cIn,cIn],'rows')); % interior
    ipar = find(ismember(labels,[cEq,cIn],'rows')); % parabola
    ilin = find(ismember(labels,[cIn,cEq],'rows')); % line
    icor = find(ismember(labels,[cEq,cEq],'rows'));   % corners
    
    npts = size(pts,1);
    
    pts_int = pts(iint,:);
    pts_par = pts(ipar,:);
    pts_lin = pts(ilin,:);
    pts_cor = pts(icor,:);
    
    x_int = pts_int(:,1);
    y_int = pts_int(:,2);
    x_par = pts_par(:,1);
    y_par = pts_par(:,2);
    x_lin = pts_lin(:,1);
    y_lin = pts_lin(:,2);
    x_cor = pts_cor(:,1);
    y_cor = pts_cor(:,2);
    
    % Calculate ratios in each state
    npar = length(x_par);
    nint = length(x_int);
    ncor = length(x_cor);
    nlin = length(x_lin);
     % -------  Display counts in each state  -------
    tcor = 2 * (ncor > 0);
    tlin = 2*sqrt(2) * (nlin > 0);
    tpar = (3*sqrt(2) + 1/2*asinh(2*sqrt(2))) * (npar > 0);
    tint = 8*sqrt(2)/3 * (nint > 0);
    ttot = tcor + tlin + tpar + tint;
    display(['Interior: theory = ',num2str(tint/ttot),', sampler = ',num2str(nint/npts), ', nint = ',num2str(nint)]);
    display(['Parabola: theory = ',num2str(tpar/ttot),', sampler = ',num2str(npar/npts), ', npar = ',num2str(npar)]);
    display(['Line    : theory = ',num2str(tlin/ttot),', sampler = ',num2str(nlin/npts), ', nlin = ',num2str(nlin)]);
    display(['Corners : theory = ',num2str(tcor/ttot),', sampler = ',num2str(ncor/npts), ', ncor = ',num2str(ncor)]);

    
    % ------   plot histogram   ------ 
    h1 = figure(3);
    clf
    set(h1,'Units','centimeters');
    set(h1,'Position',pos);
    
    
    % ----- Interior histogram -----
    % empirical histogram
    tt = linspace(-sqrt(2),sqrt(2),nbox);
    dt = tt(2)-tt(1);
    ph = histcounts(x_int,tt);
    ph = ph/nint/dt; % normalize so sum is 1
    
    % theoretical histogram
    tt1 = tt(1:end-1) + dt/2;
    ptheory = 2-tt1.^2;
    ptheory = ptheory/sum(ptheory)/dt;
    
    % plot histogram
    ha1 = subplot(1,3,1);
    plot(tt1,ptheory,'-',tt1,ph,'x-','MarkerSize',ms);
    legend({'$\propto 2-x^2$';'empirical'},'Interpreter','latex','Location','SouthEast');
    title('Interior ($M_1$)','Interpreter','latex');
    xlabel('x');
    ylabel('p(x)')
    ylim([0 0.65]);
    xlim([-sqrt(2),sqrt(2)]);
    
    % ----- Parabola histogram -----
    % empirical histogram
    tt = linspace(-sqrt(2),sqrt(2),nbox);
    dt = tt(2)-tt(1);
    ph = histcounts(x_par,tt);
    ph = ph/npar/dt; % normalize so sum is 1
    
    % theoretical histogram
    tt1 = tt(1:end-1) + dt/2;
    ptheory = sqrt(1+4*tt1.^2);
    ptheory = ptheory/sum(ptheory)/dt;
    

    ha2 = subplot(1,3,2);
    plot(tt1,ptheory,'-',tt1,ph,'x-','MarkerSize',ms);
    legend({'$\propto\sqrt{1+4x^2}$';'empirical'},'Interpreter','latex','Location','SouthEast');
    title('Parabola ($M_2$)','Interpreter','latex');
    xlabel('x');
    ylabel('p(x)')
    xlim([-sqrt(2),sqrt(2)]);
    ylim([0 0.6]);
    
    
    % ----- Line histogram -----
    % empirical histogram
    tt = linspace(-sqrt(2),sqrt(2),nbox);
    dt = tt(2)-tt(1);
    ph = histcounts(x_lin,tt);
    ph = ph/nlin/dt; % normalize so sum is 1
    
    % theoretical histogram
    tt1 = tt(1:end-1) + dt/2;
    ptheory = ones(size(tt1));
    ptheory = ptheory/sum(ptheory)/dt;
    
    % plot histogram
    ha3 = subplot(1,3,3);
    plot(tt1,ptheory,'-',tt1,ph,'x-','MarkerSize',ms);
    legend({'$\propto 1$';'empirical'},'Interpreter','latex','Location','SouthEast');
    title('Line ($M_3$)','Interpreter','latex');
    xlabel('x');
    ylabel('p(x)')
    ylim([0 0.5]);
    xlim([-sqrt(2),sqrt(2)]);
    
    % ----- Corner histogram -----
    ph1 = sum(abs(x_cor + sqrt(2)) < 1e-5) / ncor;
    ph2 = sum(abs(x_cor - sqrt(2)) < 1e-5) / ncor;
    disp(['[Left,Right]: [',num2str(ph1),',',num2str(ph2),...
        '], out of ',num2str(ncor),' pts in the corners.']);
    
    
    % ------ Move subplots around to remove space bewteen them ------
    
    set(ha1,'position',[0.1 .15 .22 .8])
    set(ha2,'position',[.37 .15 .22 .8])
    set(ha3,'position',[.64 .15 .22 .8])
    
    
    % Save, if desired
    if(ifsave)
        set(h1,'PaperPosition',pos);
        saveas(h1,savefile);
    end
    
end


% h1 = figure(1);
% clf
% set(h1,'Units','centimeters');
% set(h1,'Position',pos);
% 
% set(h1,'PaperPosition',pos);
%     saveas(h1,sfile);



% Return to default values
set(0,'DefaultLineLineWidth','remove')
set(0,'DefaultAxesLineWidth','remove')
set(0,'DefaultAxesFontSize','remove')
