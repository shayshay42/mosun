function [out] = PK_mean(TDBc_ugperkg, Vc, PKflag, case_no, ticker, end_time) 

% persistent data
if ( PKflag == 1 )
    out = (TDBc_ugperkg/Vc>1e-5)*TDBc_ugperkg/Vc;      
elseif ( PKflag == 0 )
    
    vIDvec = zeros(1, 14); % need for acceleration
    if case_no == 103 || case_no == 113
            vIDvec = [2001, 2002, 2501, 2502,   0,   0,  0,  0,  0,  0,  0,  0, 0,  0];
    elseif case_no == 104 || case_no == 114
            vIDvec = [3001, 3002, 3501, 3502, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0];
    elseif case_no == 105 || case_no == 115
            vIDvec = [4001, 4002, 4003, 4004, 4501, 4502, 4503,20013,20023,20033,20043, 0, 0,  0];
%             vIDvec = [3001, 3002, 3501, 3502];
    elseif case_no == 106
            vIDvec = [10012,10022,10032,10015,10025,10035,20015,20025,20035,30015,30025,30035, 0,  0];
%             vIDvec = [3001, 3002, 3501, 3502];
    end

%     pk_ary = zeros(1, length(vIDvec));
%     for vidpos = 1 : length(vIDvec)
%         pk_ary(vidpos) = PK(10, 1, 0, vIDvec(vidpos), ticker, end_time);
%     end
%     out = mean(pk_ary);

    pk_sum = 0;
    pk_mul = 1;
    cnt = 0;
    for vidpos = 1 : length(vIDvec)
        if vIDvec(vidpos) ~= 0
            switch case_no
                case {103, 104, 105, 106, ...
                        113, 114, 115}
                    pk_sum = pk_sum + PK(10, 1, 0, vIDvec(vidpos), ticker, end_time);
                    pk_mul = pk_mul * PK(10, 1, 0, vIDvec(vidpos), ticker, end_time);
            end
            cnt = cnt + 1;
        end
    end
    % out = pk_sum/cnt;
    out = pk_mul^(1/cnt);
elseif ( PKflag == 2 )
    out = PK(10, 1, 0, case_no, ticker, end_time);
else
    out = 0;
end

end

