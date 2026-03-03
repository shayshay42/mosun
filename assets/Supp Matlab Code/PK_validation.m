function [out] = PK_validation(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time) 

% persistent data

if ( PKflag == 1 )
    out = (TDBc_ugperkg/Vc>1e-5)*TDBc_ugperkg/Vc;      
else

    data.AnimalID = [1001	1001	1001	1001	1001	1001	1001	1001	1001	1001	1001	1001	1002	1002	1002	1002	1002	1002	1002	1002	1002	1002	1002	1002	1003	1003	1003	1003	1003	1003	1003	1003	1003	1003	1003	1003	1004	1004	1004	1004	1004	1004	1004	1004	1004	1004	1004	1004	2001	2001	2001	2001	2001	2001	2001	2001	2001	2001	2001	2001	2002	2002	2002	2002	2002	2002	2002	2002	2002	2002	2002	2002	2003	2003	2003	2003	2003	2003	2003	2003	2003	2003	2003	2003	2004	2004	2004	2004	2004	2004	2004	2004	2004	2004	2004	2004	3001	3001	3001	3001	3001	3001	3001	3001	3001	3001	3001	3001	3002	3002	3002	3002	3002	3002	3002	3002	3002	3002	3002	3002	3003	3003	3003	3003	3003	3003	3003	3003	3003	3003	3003	3003	3004	3004	3004	3004	3004	3004	3004	3004	3004	3004	3004	3004]';
    data.Time = [0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29	0	0.0035	0.083	0.25	1	2	4	7	10	14	21	29]';
    data.PK = [NaN	2230.83	1209.39	1003.75	493.41	299.18	182.89	135.48	75.14	NaN	NaN	NaN	NaN	2072.39	1261.22	926.72	454.26	231.39	177.70	147.42	61.99	NaN	NaN	NaN	NaN	1953.68	1543.82	842.05	468.12	302.23	204.55	143.34	93.20	NaN	NaN	NaN	NaN	1793.37	1214.90	781.56	406.61	224.53	156.67	106.55	10.70	NaN	NaN	NaN	NaN	159.19	106.50	91.04	54.43	48.07	23.22	17.42	9.63	NaN	NaN	NaN	NaN	144.74	121.67	69.21	38.82	32.67	20.35	12.79	NaN	NaN	NaN	NaN	NaN	79.20	76.01	61.59	53.77	37.96	24.95	16.36	11.75	NaN	NaN	NaN	NaN	128.29	110.97	83.32	45.18	27.36	19.91	12.52	6.76	NaN	NaN	NaN	NaN	10.95	9.08	7.26	5.13	3.09	NaN	NaN	NaN	NaN	NaN	NaN	NaN	8.52	7.34	6.44	5.19	3.14	NaN	NaN	NaN	NaN	NaN	NaN	NaN	12.96	8.43	7.45	4.30	3.11	NaN	NaN	NaN	NaN	NaN	NaN	NaN	12.26	9.46	8.44	6.11	4.53	3.57	NaN	NaN	NaN	NaN	NaN]'/1000;

    cpos=data.AnimalID==VPid;
    Tvec=data.Time(cpos);
    PKvec=data.PK(cpos);
    Tvec=Tvec(isnan(PKvec)==0);
    PKvec=PKvec(isnan(PKvec)==0);
    [TvecR,dpos]=unique(Tvec);
    PKvecR=PKvec(dpos);
   
    % Iraj: updated the last time point
    PKvecR=max(1e-6,PKvecR);
    TvecR=[-1; TvecR];
    PKvecR=[PKvecR(1); PKvecR];    
    PKvecR=log10(PKvecR);

    m = (PKvecR(end)-PKvecR(end-1))/(TvecR(end)-TvecR(end-1));
    m = min(m, -0.02);
    TvecR=[TvecR; end_time];
    % PKvecR=[PKvecR; max(-6, PKvecR(end)+m*(end_time-TvecR(end-1)))];    
    % Sep 19, 2016, I changed this because didn't want to force the PK to
    % be high for 0.001 mg/kg
    PKvecR=[PKvecR; max(-9, PKvecR(end)+m*(end_time-TvecR(end-1)))];   
    
    
    out = 10^interp1(TvecR, PKvecR, ticker);
end

end

