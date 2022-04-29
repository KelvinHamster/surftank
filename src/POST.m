
sim = sim.update_characteristics();
PTS = 40;

xmin = min(sim.boundary.x);
xmax = max(sim.boundary.x);
ymin = min(sim.boundary.z);
ymax = max(sim.boundary.z);


xmin = 10; xmax = 30;
ymin = -0.9;
ymax = 0.2;

x = linspace(xmin,xmax,PTS*2 + 1);
y = linspace(ymin,ymax,PTS*2 + 1);

x = x(2:2:PTS*2);
y = y(2:2:PTS*2);

s = zeros(length(y),length(x));
s_interior = zeros(length(y),length(x));
for j = 1:length(y)-1
    for i=1:length(x)-1
        s_interior(j,i) = 1*sim.is_point_inside(x(i),y(j));
        if s_interior(j,i)
            s(j,i) = sim.interior_eval(x(i),y(j));
        end
        %fprintf('%d ',i)
    end
    fprintf('%d ',j)
end
fprintf('\n')

figure(1)
splot = surface(x,y,s,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');
splot.AlphaData = s_interior;
hold on

sim.plot_full();
xlim([xmin, xmax])
ylim([ymin, ymax])
