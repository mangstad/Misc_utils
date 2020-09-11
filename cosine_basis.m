function cdrift = cosine_basis(high,low, frametimes)
    len_tim = numel(frametimes);
    n_times = 0:(len_tim-1);
    dt = frametimes(2)-frametimes(1);
    sr = 1/dt;
    nyq = sr/2;
    nfct = sqrt(2/len_tim);
    
    order = max(floor(2*len_tim*nyq*dt),1);
    cdrift = zeros(len_time,order);

    for k = 1:(order-1)
        cdrift(:,k) = nfct*cos((pi/len_tim)*(n_times+0.5)*k);
    end
    cdrift(:,order) = 1;
    
    order2 = max(floor(2*len_tim*low*dt),1);
    order3 = max(floor(2*len_tim*high*dt),1);

    cdrift(:,order2:order3) = [];
    
    
    
    