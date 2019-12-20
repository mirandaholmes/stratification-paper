% Make figure showing probability of octahedron at different choices of 
% sticky parameters kAA, kAB, kBB

% Created November 9, 2019

%clear;

ifsave = 1;
kbb0 = 0.0;

if(kbb0 == 0)
    savefile1 = 'octahedronAB_b0.epsc';  % use this for kbb0 = 0
    savefile2 = 'rigidAB_b0.epsc';
else
    savefile1 = 'octahedronAB_b1.epsc';  % use this for kbb0 = 0.1
    savefile2 = 'rigidAB_b1.epsc';
end



% kappa parameters
nk = 20;
kaa0 = logspace(log(0.01)/log(10),log(1e3)/log(10),nk);
kab0 = logspace(log(0.01)/log(10),log(1e3)/log(10),nk);
%kaa0 = [0.05,0.1,0.5,1,2,5,10,20,40];
%kab0 = [0.05,0.1,0.5,1,2,5,10,20,40,80,200,500,1e3];
%kbb0 = 0.05;%0.01;%[0.05,0.1,0.5,1,2,5,10,20,40];


labelfile = 'Data/PolymerStrat_s1_labels.txt'; 
kap = 2; 

disp('Loading data');
tic
%labels = load(labelfile);
toc
nd = size(labels,1);

n = 6;    % number of spheres
dim = 3;  % dimension of spheres



% --------------------------------
%     Parameters
% --------------------------------

% Particle type
A = 1; B = 2;
ptype = [B,A,A,A,A,B]; % B - A - A - A - A - B 

  
% edges are (c-indexing, starts at 0)
edges = [[0 0 0 0 0 1 1 1 1 2 2 2 3 3 4];
[ 1 2 3 4 5 2 3 4 5 3 4 5 4 5 5]] + 1;

% Edge Codes
tf = 1;   % fixed edge
taa = 2;  % A-A
tab = 3;  % A-B
tbb = 4;  % B-B

% Particular edge types
edgetype = NaN(1,15);
for ip=1:15
    ir = edges(1,ip);
    ic = edges(2,ip);
    if(ic == ir+1)
        edgetype(ip) = tf;
    elseif(ptype(ir) == A && ptype(ic) == A)
        edgetype(ip) = taa;
    elseif(ptype(ir) == B && ptype(ic) == B)
        edgetype(ip) = tbb;
    else
        edgetype(ip) = tab;
    end
end

% B - A - A - A - A - B 
%edgetype = [tf,tab,tab,tab,tbb, tf,taa,taa,tab, tf,taa,tab, tf,tab, tf];
  

% Label constants
cEq = 1;  % equation
cIn = 0;  % inequality



% --------------------------------
%     Find indices of octahedron
% --------------------------------

i0 = find(sum(labels==cEq,2)==12);  % indices of rigid clusters
i1 = find(sum(labels(:,[1,2,3,4,5]),2)==5);  % particle 1 touches everyone
i2 = find(sum(labels(:,[6,7,8,9,1]),2)==5);  % particle 2 touches everyone
i3 = find(sum(labels(:,[10,11,12,2,6]),2)==5);  % particle 3 touches everyone
i4 = find(sum(labels(:,[3,7,10,13,14]),2)==5);  % particle 4 touches everyone
i5 = find(sum(labels(:,[4,8,11,13,15]),2)==5);  % particle 5 touches everyone
i6 = find(sum(labels(:,[5,9,12,14,15]),2)==5);  % particle 6 touches everyone

itouchall = sort(unique([i1;i2;i3;i4;i5;i6]));
ioct = setdiff(i0,itouchall);
ipoly = setdiff(i0,ioct);


% --------------------------------
%     Compute types of edges
% --------------------------------
% 
edgetypeM = repmat(edgetype,[nd,1]);
naa = sum((labels==cEq).*(edgetypeM==taa),2);
nab = sum((labels==cEq).*(edgetypeM==tab),2);
nbb = sum((labels==cEq).*(edgetypeM==tbb),2);



% -----------------------------------------
%     Compute probability of octahedron
% -----------------------------------------

poct = NaN(length(kab0),length(kaa0),length(kbb0));
ppoly = NaN(length(kab0),length(kaa0),length(kbb0));
poct0 = NaN(length(kab0),length(kaa0),length(kbb0));

tic
for kk = 1:length(kbb0)
    
    kbb = kbb0(kk);
    kbbexp = ((kbb/kap).^nbb);
    
    for ik = 1:length(kab0)
        
        kab = kab0(ik);
        kabexp = ((kab/kap).^nab);
               
        for jk = 1:length(kaa0)
            
            kaa = kaa0(jk);
            kaaexp = ((kaa/kap).^naa);
            
            % Compute partition function weights
            zweights = kaaexp.*kabexp.*kbbexp;
            %zweights = (kaa.^naa).*(kab.^nbb).*(kbb.^nbb);
            
            % Compute probability of the octahedron
            poct(ik,jk,kk) = sum(zweights(ioct))/sum(zweights);
            ppoly(ik,jk,kk) = sum(zweights(ipoly))/sum(zweights);
            poct0(ik,jk,kk) = poct(ik,jk,kk)/(poct(ik,jk,kk) + ppoly(ik,jk,kk));
            
        end  
    end
    toc
end


disp('Isotropic interactions, kap=2:');
disp(['P(oct) = ',num2str(length(ioct)/nd)]);
disp(['P(poly) = ',num2str(length(ipoly)/nd)]);
disp(['P(oct | rigid) = ',num2str(length(ioct)/length(i0))]);

%return;


% --------------------------------
%     Plot
% --------------------------------

[xx,yy] = meshgrid(kaa0,kab0);


% Make figures look nicer
fs = 22;  % font size
fsleg = 16;  % font size for legend
ms = 12;  % marker size
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',fs)
%cc = get(gca,'colororder');


h1 = figure(1);
clf

surf(xx,yy,poct0);
xlabel('$\kappa_{AA}$','interpreter','latex');
ylabel('$\kappa_{AB}$','interpreter','latex');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[-inf,inf]);
set(gca,'ylim',[-inf,inf]);
set(gca,'zlim',[0 1]);
title(['$\kappa_{BB} = $',num2str(kbb0)],'interpreter','latex');
zlabel('Prob(octahedron | rigid)');
shading interp



h2 = figure(2);
clf
surf(xx,yy,poct+ppoly);
xlabel('$\kappa_{AA}$','interpreter','latex');
ylabel('$\kappa_{AB}$','interpreter','latex');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[-inf,inf]);
set(gca,'ylim',[-inf,inf]);
set(gca,'zlim',[0 1]);
title(['$\kappa_{BB} = $',num2str(kbb0)],'interpreter','latex');
zlabel('Prob(rigid)');
shading interp


if(ifsave)
    saveas(h1,savefile1);
    saveas(h2,savefile2);
end


% Return to default values
set(0,'DefaultLineLineWidth','remove')
set(0,'DefaultAxesLineWidth','remove')
set(0,'DefaultAxesFontSize','remove')





