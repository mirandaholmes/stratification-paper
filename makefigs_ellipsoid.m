% plot a 3-dimensional ellipsoid, and intersecting planes

% created nov 13, 2019

clear;

ifsave = 1;
savefile = 'ellipsoid.epsc';


a1 = 3;
a2 = 2;
a3 = 1;

c = 1.2;
xmin = -c*a1;
xmax = c*a1;
ymin = -c*a2;
ymax = c*a2;
zmin = -c*a3;
zmax = c*a3+0.2;

pos = [24,10,18,12];

salph = 0.6;  % transparency for spheres, 0 for transparent (0.65 orig)
lalph = 1;  % transparency for lines
lightcolr = 0.5*[1 1 1];  % colour of light
lightpos = [-1 0.25 1]; %[0 0.25 1];  % light position
ambstrength = 0.5;  % intensity of ambient component of light reflected from object
specstrength = 0.8;  % intensity of specular component of reflected light
diffstrength = 1;  % intensity of diffuse component of reflected light
specexp = 2;    % specular exponent (large = small light spots)


h1=figure(1);
clf
set(h1,'Units','centimeters');
set(h1,'Position',pos);
hold on

[x,y,z] = ellipsoid(0,0,0,a1,a2,a3,30);
surf(x,y,z,...
        'LineStyle','none',...
        'FaceColor',0.5*[1 1 1],...
        'FaceAlpha',salph,...
        'DiffuseStrength',diffstrength,...
        'AmbientStrength',ambstrength,...
        'SpecularStrength',specstrength,...
        'SpecularExponent',specexp);

[xx,yy] = meshgrid(xmin:0.5:xmax,ymin:0.5:ymax);
zz = zeros(size(xx));
surf(xx,yy,zz,'FaceAlpha',0.4,'EdgeColor','none','FaceColor',0.5*[1 0 0]);

[x2,z2] = meshgrid(xmin:0.5:xmax,zmin:0.5:zmax);
y2 = zeros(size(x2));
surf(x2,y2,z2,'FaceAlpha',0.4,'EdgeColor','none','FaceColor',0.5*[0 0 1]);


hold off
set(gca,'visible','off')

daspect([1,1,1]);
view(3);
view([-39,22]);
light('Position',lightpos,'Style','infinit','Color',lightcolr);
lighting phong
drawnow


if(ifsave)
    saveas(h1,savefile);
end



