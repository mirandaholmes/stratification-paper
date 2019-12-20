% Plot an image of the polymer

% Created September 11, 2019

clear;

ifsave = 1;
savefile = 'polymer_image.epsc';

datafile = 'Data/PolymerStrat_s3_pts.txt';
xlist = load(datafile);
nx = size(xlist,1);   % total number of points to animate

n = 6;
dim = 3;

pos = [24,10,40,12];

swidth = 0.3;
sheight = 1.2;
sdiff = 0.2;
sbot = -0.1;
psub = [0, sbot, swidth, sheight;
    1.05*sdiff, sbot, swidth, sheight;
    1.95*sdiff, sbot, swidth, sheight;
    2.8*sdiff, sbot, swidth, sheight;
    3.6*sdiff, sbot, swidth, sheight];
        
twidth = 0.1;
theight = 0.1;
tbot = 0.8;
tdiff = 0.2;
t0 = 0.09;
tsub = [t0, tbot, twidth, theight;
    t0+1.1*tdiff, tbot, twidth, theight;
    t0+2*tdiff, tbot, twidth, theight;
    t0+2.85*tdiff, tbot, twidth, theight;
    t0+3.7*tdiff, tbot, twidth, theight];

fs = 18;

h1 = figure(1);
clf
set(h1,'Units','centimeters');
set(h1,'Position',pos);




Rx = @(th) ([1, 0, 0, ; ...
    0, cos(th), -sin(th) ;...
    0, sin(th), cos(th) ]);
Ry = @(th) ([cos(th), 0, sin(th);...
    0, 1, 0;...
    -sin(th), 0, cos(th)]);
Rz = @(th) ([cos(th), -sin(th), 0; ...
    sin(th), cos(th), 0; ...
    0, 0, 1]);

xthlist = [0,0.8,0,-0.1,0];
ythlist = [0,0,0,0,0];
zthlist = [-0.6,-0.5,0,0.4,-0.2];

ik=1;
for k=[2, 65, 160, 500, 850] %1:nx
    
    x = center(xlist(k,:),dim);
    
    % Options for plotting
    %ff = 1;
    ff = ik;  % which figure
    xth = xthlist(ik);
    yth = ythlist(ik);
    zth = zthlist(ik);
    
    
    % -------------------------
    % This is all copy-pasted from plotcluster3d.m
    
    ax = 2.5*[-1 1 -1 1 -1 1];
    bondeps = 1e-5;
    iftight = 0;
    iflines = 1;
    iftext = 1;
    tfs = 24;  % text font size
    tshift = 0.12;  % text shift
    srad = 0.5;  % radius of spheres
    salph = 0.5;  % transparency for spheres, 0 for transparent (0.65 orig)
    lalph = 1;  % transparency for lines
    lightcolr = 0.5*[1 1 1];  % colour of light
    lightpos = [-1 0.25 1]; %[0 0.25 1];  % light position
    ambstrength = 0.5;  % intensity of ambient component of light reflected from object
    specstrength = 0.8;  % intensity of specular component of reflected light
    diffstrength = 1;  % intensity of diffuse component of reflected light
    specexp = 2;    % specular exponent (large = small light spots)
    lcolr = 0.7*[1,1,1];  % line colour (cylinder)
    lrad = 0.05;  % radius of lines (cylinder)  % increased from 0.03 feb 27 2014
    
    scolr=[1,0,0;...
        0,1,0;...
        0,0,1;...
        1,1,0;...
        1,0,1;...
        0,1,1;...
        0.25,0.25,0;...
        0,0.25,0.25;...
        0.25,0,0.25;...
        0,0,0.5;...
        0,0.5,0;...
        0.5,0,0;...
        ];  % default color scheme for spheres

    
    % get adjacency matrix
    r = zeros(n);
    for ii=1:n
        for jj=ii+1:n
            r(ii,jj) = norm(x(3*ii-2:3*ii)-x(3*jj-2:3*jj));
            r(jj,ii) = r(ii,jj);
        end
    end
    a0 = +( (r + 2*(1+bondeps)*eye(n))< 1 + bondeps);
    
    nb=sum(sum(triu(a0)));
    
    % Set up spheres
    [sx,sy,sz]= sphere(40);  %-- sphere facets
    X = Rz(zth)*Ry(yth)*Rx(xth)*reshape(x,3,n);
    
    
    pax = axes('Position',psub(ik,:));
    %subplot('Position',psub(ik,:));
    %subplot(1,5,ik)
    
    annotation('textbox', tsub(ik,:), 'String',...
        {[num2str(nb),' bonds']}, 'FontSize',fs,...
        'FitBoxToText', 'on','HorizontalAlignment','center','LineStyle','none');
    %text({[num2str(k), ' steps'];[num2str(nb),' bonds']},'FontSize',fs);
    %title([num2str(k), ' steps'],'FontSize',fs);
    
    %del = 0;
    %axC = axes(pax,'Position',[0+del 0+del 1-2*del 1-2*del]);  % figure goes from 0 to 1, takes full window
    
    
    
    hold on
    % plot spheres
    for jn=1:n
        if(size(scolr,1)==1)
            scolt = scolr;
        else
            scolt = scolr(jn,:);
        end
        surf(sx*srad+X(1,jn), sy*srad+X(2,jn), sz*srad+X(3,jn),...
            'LineStyle','none',...
            'FaceColor',scolt,...
            'FaceAlpha',salph,...
            'DiffuseStrength',diffstrength,...
            'AmbientStrength',ambstrength,...
            'SpecularStrength',specstrength,...
            'SpecularExponent',specexp);
    end
    
    % plot lines
    if(iflines)
        ind = find(triu(a0));
        siz = size(a0);
        for jb=1:length(ind)  % plot non-broken bonds
            [ir,ic] = ind2sub(siz,ind(jb));
            v1 = [X(1,ir),X(2,ir),X(3,ir)];
            v2 = [X(1,ic),X(2,ic),X(3,ic)];
            [xc,yc,zc] = cylinder2P(lrad,20,v1,v2);
            if(size(lcolr,1)==1)
                lcolt = lcolr;
            else
                lcolt = lcolr(jb,:);
            end
            surf(xc,yc,zc,...
                'LineStyle','none',...
                'FaceColor',lcolt,...
                'FaceAlpha',lalph,...
                'DiffuseStrength',diffstrength,...
                'AmbientStrength',ambstrength,...
                'SpecularStrength',specstrength,...
                'SpecularExponent',specexp);
            % Previous code
            %     line([X(1,ir) X(1,ic)],[X(2,ir) X(2,ic)],[X(3,ir) X(3,ic)],...
            %         'Color',lcolr,'LineWidth',lw);
        end
    end
    hold off
    axis(ax);
    daspect([1,1,1]);
    light('Position',lightpos,'Style','infinit','Color',lightcolr);
    lighting phong
    axis off
    %axes(axC);
    % sphere labels
    if(iftext)
        for jn=1:n
            text(X(1,jn)-tshift,X(2,jn)-tshift,X(3,jn)+tshift,num2str(jn),'FontSize',tfs,...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    
    %--------------------

    drawnow
    ik=ik+1;
end



if(ifsave)
    saveas(h1,savefile);
end



% Function to put center of mass at origin

% Created April 25, 2012

function xc = center(x,varargin)

if(nargin > 1)
    d = varargin{1};
else
    d = 3;
end

xc = x;
for jj=1:d
    xc(jj:d:end) = xc(jj:d:end) - mean(xc(jj:d:end));
end

end


%  CYLINDER:  A function to draw a N-sided cylinder based on the
%             generator curve in the vector R.
%
%  Usage:      [X, Y, Z] = cylinder(R, N)
%
%  Arguments:  R - The vector of radii used to define the radius of
%                  the different segments of the cylinder.
%              N - The number of points around the circumference.
%
%  Returns:    X - The x-coordinates of each facet in the cylinder.
%              Y - The y-coordinates of each facet in the cylinder.
%              Z - The z-coordinates of each facet in the cylinder.
%
%  Author:     Luigi Barone
%  Date:       9 September 2001
%  Modified:   Per Sundqvist July 2004

function [X, Y, Z] = cylinder2P(R, N,r1,r2)

% The parametric surface will consist of a series of N-sided
% polygons with successive radii given by the array R.
% Z increases in equal sized steps from 0 to 1.

% Set up an array of angles for the polygon.
theta = linspace(0,2*pi,N);

m = length(R);                 % Number of radius values
% supplied.

if m == 1                      % Only one radius value supplied.
    R = [R; R];                % Add a duplicate radius to make
    m = 2;                     % a cylinder.
end


X = zeros(m, N);             % Preallocate memory.
Y = zeros(m, N);
Z = zeros(m, N);

v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
%cylinder axis described by: r(t)=r1+v*t for 0<t<1
R2=rand(1,3);              %linear independent vector (of v)
x2=v-R2/(R2*v');    %orthogonal vector to v
x2=x2/sqrt(x2*x2');     %orthonormal vector to v
x3=cross(v,x2);     %vector orthonormal to v and x2
x3=x3/sqrt(x3*x3');

r1x=r1(1);r1y=r1(2);r1z=r1(3);
r2x=r2(1);r2y=r2(2);r2z=r2(3);
vx=v(1);vy=v(2);vz=v(3);
x2x=x2(1);x2y=x2(2);x2z=x2(3);
x3x=x3(1);x3y=x3(2);x3z=x3(3);

time=linspace(0,1,m);
for j = 1 : m
    t=time(j);
    X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x;
    Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y;
    Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
end

%surf(X, Y, Z);

end


