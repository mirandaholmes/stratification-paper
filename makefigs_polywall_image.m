% plot polymer on a wall

% created nov 8 2019

clear;

whichfig = 1;    % 1 = no bending stiffness, 2 = bending stiffness
ifsave = 1;

if(whichfig == 1)
    savefile = 'polywall_free_image.epsc';
    ptfile = 'Data/PolymerWall_n20_k1_pts.txt';

    ix = 192;   % which one to plot
    
%     xmin = -4;
%     xmax = 3;
%     ymin = -5;
%     ymax = 2;

    xmin = -6;
    xmax = 5;
    ymin = -7;
    ymax = 4;
    
    zmin = -0.5;
    zmax = 1;
    ax = [xmin xmax ymin ymax zmin zmax];
    az = -25;
    el = 25;
    zth = 0;
end

if(whichfig == 2)
    savefile = 'polywall_bend_image.epsc';
    ptfile = 'Data/PolymerWall_kbend_n20_k1_pts.txt';

    ix = 200;   % which one to plot
    
    xmin = -3;
    xmax = 8;
    ymin = -3;
    ymax = 8;
    zmin = -0.5;
    zmax = 1;
    ax = [xmin xmax ymin ymax zmin zmax];
    az = -25;
    el = 25;
    zth = -0.6;
end


cwall = [0,0,1];  % wall color
cfree = [0,1,0];  % free space color


xlist = load(ptfile);
nx = size(xlist,1);   % total number of points
dim = 3;
n = size(xlist,2) / dim;  % number of spheres


x = xlist(ix,:);
x(1:3:end) = x(1:3:end) - x(1);
x(2:3:end) = x(2:3:end) - x(2);
z = x(3:3:end);
iwall = (z<1e-9);

scolr = ones(n,3);

for ii=1:n
    if(iwall(ii) == 1)  % on wall
        scolr(ii,:) = cwall;
    else
        scolr(ii,:) = cfree;
    end
end


% surface
[xx,yy] = meshgrid(xmin:0.5:xmax,ymin:0.5:ymax);
zz = zeros(size(xx));

% adjacency matrix
a = zeros(n);
for k=1:n-1
    a(k,k+1) = 1;
    a(k+1,k) = 1;
end

% Options for plotting
opts = struct('srad',0.18,'fig',1,'lcolr',0.7*[1,1,1],'lrad',0.05,...
    'scolr',scolr,'salph',0.5,'ifgrid',0,'iftext',0,'bondeps',1e-5,...
    'ax',ax,'iftight',0,'a',a,'pos',[24,10,20,12],'az',az,'el',el,'zth',zth);


h1 = figure(1);
clf
plotcluster3d(x,opts);
hold on
surf(xx,yy,zz,'FaceAlpha',0.5,'EdgeColor','none','FaceColor',[0.2 0.2 0.8]);
hold off
%axis equal
set(gca,'visible','off')
drawnow

if(ifsave)
    saveas(h1,savefile);
end
        