function NRMSE = calc_NRMSE(gold,rec)
% (c) E. Shimron 2019

gold = abs(gold);
rec = abs(rec);

RMSE = sqrt(sum((rec(:)-gold(:)).^2 ) / (size(gold,1)*size(gold,2)) );

range = max(gold(:)) - min(gold(:)) ;

if range<=0
    disp('problem in NRMSE calculation: range cannot be negative')
    return
end

NRMSE = RMSE/range;  

end