function o = orth_poly(x,degree)

if (size(x,2)>1 & size(x,1)==1)
    x = x';
end

ans = [];
n = numel(x);
for k = 0:degree
    v = x.^k;
    cmps = zeros(n,1);
    if (k>0) 
        for j = 0:(k-1)
            cmps = cmps + v'*ans(:,j+1).*ans(:,j+1);
        end
    end
    p = v - cmps;
    p = p/sum(p.^2).^0.5;
    ans = [ans p];
end
o = ans(:,2:end);

    