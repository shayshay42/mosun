function [out] = PK_v26(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time, fvalidation) 

if fvalidation == 0
    out = PK(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time);
elseif fvalidation == 1
    out = PK_validation(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time);
elseif fvalidation == 2
    out = PK_validation_2(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time);
elseif fvalidation == 5
    out = PK_mean(TDBc_ugperkg, Vc, 0, VPid, ticker, end_time);    
else
    out = 0;
end

end

