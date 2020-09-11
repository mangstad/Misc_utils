function QC_RSFC_plot(px,py,Group,graph,titlesuffix)
scatter(px,py,graph.size,'fill');refline(0,0);
if graph.iffit
    p=polyfit(px,py,1);
    fitval=p(1)*px+p(2);
    hold on;
    plot(px,fitval,'r');
end 
if ~isfield(graph,'axisize')
    graph.axisize=11;
end
if ~isfield(graph,'labelsize')
    graph.labelsize=11;
end
set(gca,'FontSize',graph.axisize)
xlabel('distance','FontSize',graph.labelsize);ylabel('qc-rsfc correlation','FontSize',graph.labelsize);ylim([-0.8,0.8]);
title(strcat(graph.titleprefix,{' '},Group,titlesuffix));


end

