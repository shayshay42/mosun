function [ h_plot, G_cnt, Obj, val, ob_ary_tot] = objective_run(case_no, Pvec, param, expmodel, data, datavec, outvec, dose_names, color_str)

s_med_col_ary =  {[192 080 077]/255, [247 150 070]/255, [075 172 198]/255, [155 187 089]/255, [128, 100, 162]/ 255};

switch case_no
    case {10, 11, 12, 13}
        color_str = s_med_col_ary{1};
    case {20, 21, 22}
        color_str = s_med_col_ary{4};
    case {30}
        color_str = s_med_col_ary{4};
end

h_plot = 1;
val = [];
ob_ary_tot = [];

IDvec = unique(data.ID);

switch case_no
    case {10}
        G = {[1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504]};
        G_str = 'Fit - Placebo'; title_str = '  Control - qwx4  ';    
    case {11}
        G = {[2001, 2002, 2501, 2502]};
        G_str = 'Fit - Low'; title_str = '  0.01 mg/kg - qwx4  ';
    case {12}
        G = {[3001, 3002, 3501, 3502]};
        G_str = 'Fit - Med';  title_str = '  0.1 mg/kg - qwx4  ';   
    case {13}
        G = {[4001, 4002, 4003, 4004, 4501, 4502, 4503, 20013, 20023, 20033, 20043]};
        G_str = 'Fit - High'; title_str = '  1.0 mg/kg - qwx4  ';
        
    case {20}
        G = {[3001, 3002, 3003]};
        G_str = 'Valid - Very Low'; title_str = '  0.001 mg/kg - single  ';   
    case {21}
        G = {[2001, 2002, 2003, 2004]};
        G_str = 'Valid - Low'; title_str = '  0.01 mg/kg - single  ';     
    case {22}
        G = {[1001, 1002, 1003, 1004]};
        G_str = 'Valid - Med'; title_str = '  0.1 mg/kg - single  ';
    case {30}
        G = {[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20]}; 
        G_str = 'Valid - Blin';  title_str = '  15 \mug/m^{\fontsize{12}2}/day - 4 weeks  '; 
end

Obj = 0;

G_cnt = length(G);
for G_no = 1 : length(G)
    IDvec = G{G_no};
   
    ind = find(ismember(data.ID, IDvec));
    Tdata_tmp = sort(unique(data.Time(ind)));
    SimTime = max(Tdata_tmp)+7;
    if case_no == 30
        SimTime = max(Tdata_tmp);
    end
    outtime = sort(unique([(0:0.1:SimTime)'; Tdata_tmp]));
    outtime = outtime(outtime >= 0);
    
    Xpksim_ary = zeros(length(outtime), length(IDvec));
    Xsim_ary = zeros(length(outtime), length(outvec)*length(IDvec));
    
    Xdata_ary = {};
    Tdata_ary = {};
    em = expmodel;
    em.SimulationOptions.OutputTimes = outtime;
    em.SimulationOptions.StopTime = SimTime;
    em.SimulationOptions.MaximumWallClock = 30;
    
    for idpos = 1 : length(IDvec)

        if (case_no ~= 30)
            [Tdata, Xpk] = getIDdata(data, IDvec(idpos), {'PK'});
            Tdata = Tdata(isnan(Xpk)==0);
            Xpk   = Xpk(isnan(Xpk)==0);
            [Tdata, dpos] = unique(Tdata);
            Xpk = Xpk(dpos);
            figure(1);
            semilogy(Tdata, Xpk, '.', 'Markersize', 35, 'color', [0.7, 0.7, 0.7], 'linewidth', 2);
            hold on;
            Xpksim = zeros(1, length(outtime));
            for ii = 1:length(outtime)
                switch case_no

                    case {10, 11, 12, 13}
                        Xpksim(ii) = PK(10, 1, 0, IDvec(idpos), outtime(ii), SimTime);
                    case {20, 21, 22}
                        Xpksim(ii) = PK_validation(10, 1, 0, IDvec(idpos), outtime(ii), SimTime);
                end
            end
            Xpksim_ary(:, idpos) = Xpksim';
        end        
        
        [Tdata, Xdata] = getIDdata(data, IDvec(idpos), datavec);
        
        pos = find(data.ID == IDvec(idpos));
        init_Bpbo = nanmean(data.Bcells(pos(data.Time(pos)==0)))*1000;
        init_Tpbo = nanmean(data.CD8Tcells(pos(data.Time(pos)==0)))*1000;
        init_frac = nanmean(data.fracCD69(pos(data.Time(pos)==0)));
                
        switch case_no
            case {10, 11, 12, 13, 20, 21, 22}
                param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 0, IDvec(idpos), SimTime];  
            case {30}
                param_values = [init_Bpbo, init_Bpbo, init_Tpbo, init_Tpbo, 1]; 
        end
        
        if strcmp(unique(data.schedule(pos)),'4qw')==1
            switch case_no
                case {10, 11, 12, 13}
                    if ismember(IDvec(idpos), [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                        fixinj_dose = getdose(em, dose_names{1});
                    else
                        fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
                    end
                case {20, 21, 22}
                    fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{3})];
                otherwise
                    fixinj_dose = getdose(em, dose_names{1});
            end
        elseif strcmp(unique(data.schedule(pos)),'single dose')==1
            switch case_no
                case {10, 11, 12, 13}
                    if ismember(IDvec(idpos), [1001, 1002, 1003, 1004, 1501, 1502, 1503, 1504])
                        fixinj_dose = getdose(em, dose_names{2});
                    else
                        fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
                    end
                case {20, 21, 22}
                    fixinj_dose = [getdose(em, dose_names{2}), getdose(em, dose_names{4})];
                otherwise
                    fixinj_dose = getdose(em, dose_names{2});
            end
        elseif strcmp(unique(data.schedule(pos)),'cont_infusion')==1
            fixinj_dose = [getdose(em, dose_names{1}), getdose(em, dose_names{2}), getdose(em, dose_names{3})];
        end
                
        simData = simulate(em, param_values, fixinj_dose);
        [Tsim, Xsim] = selectbyname(simData, outvec);
        
        switch case_no
            case {10, 11, 12, 13}
                Xsim(:,1:2)=Xsim(:,1:2)/1000;
                Xsim(:,3) = (Xsim(:,3)+init_frac) * 100;                
                Xsim(:,6) = (Xsim(:,6)*4+2); % IL6 levels, baseline = 2 is added                
            case {20, 21, 22, 30}
                Xsim(:,1:2)=Xsim(:,1:2)/1000;
                Xsim(:,3) = (Xsim(:,3)+init_frac) * 100;                               
        end
        Xdata(:,3)=Xdata(:,3)*100;
        
        for ii = 1 : size(Xdata, 2)
            figure(ii+1);
            T1 = Tdata(isfinite(Xdata(:, ii)));
            X1 = Xdata(isfinite(Xdata(:, ii)), ii);
            plot(T1, X1, '.', 'Markersize', 35, 'color', [0.7, 0.7, 0.7], 'linewidth', 2);
            hold on
        end
        
        Xdata_ary{idpos} = Xdata;
        Tdata_ary{idpos} = Tdata;
        Xsim_ary(:, [1:length(outvec)]+(idpos-1)*length(outvec)) = Xsim;
    end
    
    if (case_no ~= 30)
        figure(1);
        semilogy(Tsim, mean(Xpksim_ary,2), 'color', color_str, 'linewidth', 4);
        tmp_time = [];
        for tmp_i = 1 : length(Tdata_ary)
            tmp_time = [tmp_time; Tdata_ary{tmp_i}];
        end
        tmp_time = unique(tmp_time);
        [~, PK_ind] = ismember(tmp_time, Tsim);
        PK_mean = mean(Xpksim_ary,2);
        PK_std = std(Xpksim_ary,[], 2);
        h_plot(G_no) = errorbar(Tsim(PK_ind), PK_mean(PK_ind), PK_std(PK_ind), 'color', 'k', 'linewidth', 2, 'linestyle', 'none');
        xlim([0, SimTime]);
        set(gcf, 'Position', [1 1 500 300]);
        set(gcf, 'PaperPositionMode', 'auto');
        set(gca, 'Fontsize', 20, 'Fontweight', 'Bold');
        ylim([10^-5, 10^2]);
        set(gca, 'YTick', [10^-4, 10^-2, 10^0, 10^2]);
        set(gca, 'YTickLabel', {'10^{-4}', '10^{-2}', '10^{0}', '10^{2}'});
        box off
        grid minor
        grid minor
        grid on
        xlabel('Time (days)', 'Fontsize', 20, 'Fontweight', 'Normal');
        ylabel('TDB Conc [\mug/ml]', 'Fontsize', 20, 'Fontweight', 'Normal');
        print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//PK - %s.tiff', G_str));
        title(title_str, 'Fontsize', 24, 'Fontweight', 'Normal', 'Backgroundcolor', [0.85, 0.85, 0.85])
        print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//PK with title - %s.tiff', G_str));

    end
    
    for ii = 1 : size(Xsim, 2)
        figure(ii+1);
        Xdata_tot_tmp = {};
        for jj = 1 : length(Tdata_tmp)
            Xdata_tot_tmp{jj} = [];
        end
        for kk = 1 : length(IDvec)
            Xdata_tmp = Xdata_ary{kk};
            for jj = 1 : length(Tdata_tmp)
                [d3, d2] = ismember(Tdata_ary{kk}, Tdata_tmp(jj));
                Xdata_tot_tmp{jj} = [Xdata_tot_tmp{jj}, Xdata_tmp(find(d3 == 1), ii)'];
            end
        end
        mean_tmp = cellfun(@(x) nanmean(x(:)), Xdata_tot_tmp);
        std_tmp = cellfun(@(x) nanstd(x(:)), Xdata_tot_tmp);
        t_tmp = Tdata_tmp(isfinite(mean_tmp));
        std_tmp = std_tmp(isfinite(mean_tmp));
        mean_tmp = mean_tmp(isfinite(mean_tmp));
        
        if (case_no == 30)
            if ii == 3
                std_tmp = [8.26, 16.02];
            end
        end
        errorbar(t_tmp, mean_tmp, std_tmp, 'linestyle', 'none', 'color', 'k', 'linewidth', 2);
        xlim([0, SimTime]);
        if ii == 2
            if (case_no == 30)
                ylim([0, 1500]);
            else
                ylim([0, 9000]);
            end
        end
        if ii == 1
            if (case_no == 30)
                ylim([0, 600]);
            else
                ylim([0, 2000]);
            end
        end
    end
    
    for ii = 1 : size(Xsim, 2)
        figure(ii+1);
        for jj = 1 : length(IDvec)
            Tmax_tmp = max(Tdata_ary{jj});
            find(Tsim > Tmax_tmp);
        end
        
        mean_tmp = nanmean(Xsim_ary(:, ii+length(outvec)*[0:length(IDvec)-1]),2);
        if ii == 6
            Tsim_tmp = Tsim;
            mean_tmp_tmp = mean_tmp;
            ind_tmp = find(0 < Tsim_tmp & Tsim_tmp < 0.2);
            Tsim_tmp(ind_tmp) = [];
            mean_tmp_tmp(ind_tmp) = [];
            ind_tmp = find(0.2 < Tsim_tmp & Tsim_tmp < 0.5);
            Tsim_tmp(ind_tmp) = [];
            mean_tmp_tmp(ind_tmp) = [];
            plot(Tsim_tmp, mean_tmp_tmp, '-', 'color', color_str, 'linewidth', 4);
        else
            plot(Tsim, mean_tmp, '-', 'color', color_str, 'linewidth', 4);
        end
        hold on
        set(gca, 'Fontsize', 20, 'Fontweight', 'Bold');
        xlabel('Time (days)', 'Fontsize', 20, 'Fontweight', 'Normal');
        set(gcf, 'Position', [1 1 500 300]);
        set(gcf, 'PaperPositionMode', 'auto');
        grid minor
        grid minor
        grid on
        box off
        switch ii
            case 1
                ylabel('B cell count [/\mul]', 'Fontsize', 20, 'Fontweight', 'Normal');
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//B cells - %s.tiff', G_str))
            case 2
                ylabel('CD8+ T cell count [/\mul]', 'Fontsize', 20, 'Fontweight', 'Normal');
                title(title_str, 'Fontsize', 24, 'Fontweight', 'Normal', 'Backgroundcolor', [0.85, 0.85, 0.85])
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//T cells - %s.tiff', G_str));
            case 3
                ylabel('CD69+ CD8+ T cells [%]', 'Fontsize', 20, 'Fontweight', 'Normal');
                ylim([0, 100]);
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//Act Frac - %s.tiff', G_str));
            case 4
                switch case_no
                    case {10, 11, 12, 13}
                        ylabel('Spleen B/T ratio', 'Fontsize', 20, 'Fontweight', 'Normal');
                        ylim([0, 2]);
                        print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//Spleen B_T - %s.tiff', G_str));
                        title(title_str, 'Fontsize', 24, 'Fontweight', 'Normal', 'Backgroundcolor', [0.85, 0.85, 0.85])
                        print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//Spleen B_T with title - %s.tiff', G_str));
                end
                
            case 5
                ylabel('LNs B/T ratio', 'Fontsize', 20, 'Fontweight', 'Normal');
                ylim([0, 2]);
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//LNs B_T - %s.tiff', G_str));
            case 6
                ylabel('IL6 levels (pg/mL)', 'Fontsize', 20, 'Fontweight', 'Normal');
                set(gca, 'YScale', 'Log');
                xlim([0, 29]);
                ylim([10^0, 10^4]);
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//IL6 Levels - %s.tiff', G_str));
                title(title_str, 'Fontsize', 24, 'Fontweight', 'Normal', 'Backgroundcolor', [0.85, 0.85, 0.85])
                print(gcf, '-loose', '-dtiff', '-r150', sprintf('Figures//IL6 Levels with title - %s.tiff', G_str));
        end
        
        
    end
    
end


