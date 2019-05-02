clear;
clc;
close;

IntCon = [1,5];
opts = optimoptions('ga');
opts.UseParallel = true();
opts.Display = "iter";
opts.EliteCount = 5;
oldbest = 9e9;
bestval = 100001;
bestcount = 2;
tic
RunOptimization([1,55,55,10000,2])

% while bestval >  1000
% [xbest,bestval] = ga(@RunOptimization,5,[],[],[],[],[2,10,45,10000,2],[8,120,90,10000,7],[],IntCon,opts)
% if bestval < oldbest(bestcount-1)
% opts.InitialPopulationMatrix = xbest;
% beststate(bestcount-1,:) = xbest;
% oldbest(bestcount)=bestval;
% bestcount = bestcount+1;
% save('orbit params','beststate')
% end
% end
function cost = RunOptimization(x)
n_planes = ceil(x(1));
d_RAAN = x(2);
i = x(3);
a = x(4);
nsats= ceil(x(5));
fsep = 360/(nsats);
outputs = struct([]);
covCone = 80;
% xgrid = 0:1:360;
% ygrid = 90:-1:-90;
% [~,ygrid]=meshgrid(xgrid,ygrid);
counter = 1;
for planes = 1:n_planes
    RAAN = (planes-1)*d_RAAN;
    orbitAlt = a;
    inc = i;
    for sats = 1:nsats
        f = fsep*(sats-1);
        [outputs(counter).alfa,outputs(counter).delta]=GroundTrack(orbitAlt,RAAN,inc,f);
        counter = counter + 1;
    end
    
end

coverage = zeros(length(outputs(1).alfa),1);
cost  = 0;
for timestep = 1:length(outputs(1).alfa)
    fun = @(lam,phi) 0;
    for sat = 1:length(outputs)
        alfa = outputs(sat).alfa(timestep);
        delta = outputs(sat).delta(timestep);
%         testVal = @(lam,phi) -Vincenty(alfa,delta,lam,phi)+covCone
        fun = @(lam,phi) fun(lam,phi) + heaviside(-Vincenty(alfa,delta,lam,phi)+covCone);
    end
    fun = @(lam,phi) (abs(fun(lam,phi)));
    coverage(timestep)=simp2D(fun,0,360,-90,90,180,90);
    cost = trapz(coverage)/360/180+cost;
    if any(timestep == [1])
        h=figure(timestep)
        plotcov(fun)
    end
end
cost = cost*nsats*n_planes;
end

function plotcov(fun)
h=figure(1)
colormap parula
x = 0:.5:360;
y= -90:.5:90;
[x,y] = meshgrid(x,y);
z=fun(x,y);
surf(x,y,z,'lines','none')
colorbar
axis([0,360,-90,90])
set(gca,'fontsize',18)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figure','-dpng','-r0')
end