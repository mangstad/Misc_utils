function results = mc_QC_RSFC(featuremat,meanFD,coords)

[x,y] = size(coords);
if (x==3)
    coords = coords';
end

nROI = size(coords,1);
distanceCutoff = 150;

distanceMat=zeros(nROI,nROI);
for i=1:(nROI-1)
    for j=(i+1):nROI
        distanceMat(i,j)=sqrt((coords(i,1)-coords(j,1))^2+(coords(i,2)-coords(j,2))^2+(coords(i,3)-coords(j,3))^2);
    end
end
distance=mc_flatten_upper_triangle(distanceMat);
distanceMask = distance > distanceCutoff;

distance(distanceMask) = [];

c = corr(featuremat,meanFD);
c(distanceMask) = [];

fprintf(1,'%d connections out of %d total were removed by distanceCutoff (%6.4f%%)\n',sum(distanceMask),size(distanceMask,2),sum(distanceMask)/size(distanceMask,2));

mc_scatter(distance',c,'Xlabel',{'Distance (mm)'},'Ylabel',{'qc-rsfc correlation'},'Title',{'QC-RSFC';''});