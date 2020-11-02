function [output] = mc_summarize_rois(data,rois,varargin)
%this function will take a matrix with unique values per ROI
%and will return the mean value within each ROI, indexed according 
%to the ROI numbers
%data must have subjects/timepoints in the first dimension

scale = 1;
if (nargin>2)
    scale = varargin{1};
end

s = size(data);
N = s(1);
P = prod(s(2:end));

if (numel(rois)~=P)
    error('Size of data does not match size of ROI matrix');
end

data = reshape(data,N,P);
rois = reshape(rois,P,1);

rois(rois==0) = max(rois)+1;

u = unique(rois);
mapping = [u [1:numel(u)]'];

mask = mc_Vector2Mask(rois,scale);

mean = data*mask;

output = zeros(size(mean,1),max(u));
output(:,mapping(:,1)) = mean;
