function h = mc_scatter(x,y,varargin)

bad = isnan(x) | isnan(y);
x(bad) = [];
y(bad) = [];

size=25;
color=[0 0 1];
Xlabel = 'X';
Ylabel = 'Y';
Title = 'Title';
cmap = [];
edge = [0.5 0.5 0.5];
alpha = repmat(1,numel(x),1);

p = inputParser;

addParameter(p,'size',[]);
addParameter(p,'color',[]);
addParameter(p,'Xlabel',[]);
addParameter(p,'Ylabel',[]);
addParameter(p,'Title',[]);
addParameter(p,'cmap',[]);
addParameter(p,'MarkerEdgeColor',[]);
addParameter(p,'AlphaData',[]);

parse(p,varargin{:});

if (~isempty(p.Results.size))
    size = p.Results.size;
end
if (~isempty(p.Results.color))
    color = p.Results.color;
end
if (~isempty(p.Results.Xlabel))
    Xlabel = p.Results.Xlabel;
end
if (~isempty(p.Results.Ylabel))
    Ylabel = p.Results.Ylabel;
end
if (~isempty(p.Results.Title))
    Title = p.Results.Title;
end
if (~isempty(p.Results.cmap))
    cmap = p.Results.cmap;
end
if (~isempty(p.Results.MarkerEdgeColor))
    edge = p.Results.MarkerEdgeColor;
end
if (~isempty(p.Results.AlphaData))
    alpha = p.Results.AlphaData;
end

h = figure;
scatter(x,y,size,color,'filled','MarkerEdgeColor',edge,'LineWidth',0.25,'MarkerFaceAlpha','flat','AlphaData',alpha);

if (~isempty(cmap))
    %[cmap,labels] = mc_cmap_gordon();
    colormap(cmap);
    %colorbar('Ticks',[1:16],'TickLabels',labels);
end

hold on;
h = refline(0,0);
h.Color = 'red';
h.LineWidth = 2;
title(Title,'FontSize',30,'FontWeight','bold');
xlabel(Xlabel,'FontSize',24,'FontWeight','bold');
ylabel(Ylabel,'FontSize',24,'FontWeight','bold');
fit = polyfit(x,y,1);
p = polyval(fit,x);
imin = find(x==nanmin(x));
imax = find(x==nanmax(x));

plot([nanmin(x) nanmax(x)],[p(imin(1)) p(imax(1))],'k--','LineWidth',2);
set(gcf,'color','w');
set(gcf,'position',[10 10 1024 576]);

