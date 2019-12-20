% play around with ellipsoid data

% created nov 13, 2019

clear;
whichcase = 1;

n=10;

data = load(['Data/ellipsoid_whichcase',num2str(whichcase),'_data.txt']);

te = data(:,2);  % in ellipsoid
t0 = data(:,3);  % in 0-dimensional
eqns = data(:,1);  % # of equations
nd = length(t0);

pfac = 2*2*2*2*3*3*3*1*1*1;
vfac = 2.550164;

nbin = 10;
vest = NaN(nbin,1);

for jb=1:nbin
    if(whichcase == 0)
        vest(jb) = mean(te((jb-1)*nd/nbin+1:jb*nd/nbin))/mean(t0((jb-1)*nd/nbin+1:jb*nd/nbin)) * 2  * exp(n*0.94)/exp(1*0.94);
    end
    if(whichcase == 1)
        vest(jb) = mean(te((jb-1)*nd/nbin+1:jb*nd/nbin))/mean(t0((jb-1)*nd/nbin+1:jb*nd/nbin)) * vfac * pfac;
    end
        
end

disp('estimates: ');
disp(vest');
disp(['mean = ',num2str(mean(vest))]);
disp(['std = ',num2str(std(vest)/sqrt(nbin))]);



return;

figure(1)
subplot(1,2,1)
hh = histogram(eqns,[min(eqns):max(eqns)+1]);

subplot(1,2,2)
plot([min(eqns):max(eqns)],1./hh.Values/sum(1./hh.Values));
set(gca,'yscale','log');

 
 