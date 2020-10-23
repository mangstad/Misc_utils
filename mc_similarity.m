function [dc,jc] = mc_similarity(X,Y)
%takes vectors or matrices X and Y. Both should be either n x 1 or n x p
%will return the similarity coefficient for each set (size either 1 or p)
%if one is a matrix and the other a vector it will calculate the vector
%against each column of the matrix.
%if both are matrices they must be the same size and it will calculate for
%each pair of columns 1, 2, ..., p
%the first output is the dice, the second is the jaccard

%coerce to non-logical
X = X + 1 - 1;
Y = Y + 1 - 1;

inter = X'*Y;
total = sum(X) + sum(Y);
union = total - inter;

jc = inter./union;
dc = 2*inter./total;
%dc = 2*jc./(1+jc);

if (min(size(jc))>1)
    jc = diag(jc);
    dc = diag(dc);
end

jc = jc(:);
dc = dc(:);
