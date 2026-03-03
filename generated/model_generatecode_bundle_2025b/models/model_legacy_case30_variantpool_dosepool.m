function m1 = model_legacy_case30_variantpool_dosepool()
% Created in MATLAB R2025b

%% Create the model.
m1 = sbiomodel("TDB");
m1.Notes = "";
m1.Tag = "";

%% Create compartments.
c1 = addcompartment(m1, "unnamed");
c1.Constant = true;
c1.Value = 1;
c1.Units = "";
c1.BoundaryCondition = false;
c1.Notes = "";
c1.Tag = "";

%% Create species.
s1 = addspecies(c1, "actTtiss");
s1.Constant = false;
s1.Value = 0;
s1.Units = "";
s1.BoundaryCondition = false;
s1.Notes = "";
s1.Tag = "";

s2 = addspecies(c1, "Btiss");
s2.Constant = false;
s2.Value = 0;
s2.Units = "";
s2.BoundaryCondition = false;
s2.Notes = "";
s2.Tag = "";

s3 = addspecies(c1, "TDBc_ugperkg");
s3.Constant = false;
s3.Value = 0;
s3.Units = "";
s3.BoundaryCondition = false;
s3.Notes = "";
s3.Tag = "";

s4 = addspecies(c1, "actTpb");
s4.Constant = false;
s4.Value = 0;
s4.Units = "";
s4.BoundaryCondition = false;
s4.Notes = "";
s4.Tag = "";

s5 = addspecies(c1, "drugBtisskill");
s5.Constant = false;
s5.Value = 0;
s5.Units = "";
s5.BoundaryCondition = false;
s5.Notes = "";
s5.Tag = "";

s6 = addspecies(c1, "drugTtissact");
s6.Constant = false;
s6.Value = 0;
s6.Units = "";
s6.BoundaryCondition = false;
s6.Notes = "";
s6.Tag = "";

s7 = addspecies(c1, "restTtiss");
s7.Constant = false;
s7.Value = 0;
s7.Units = "";
s7.BoundaryCondition = false;
s7.Notes = "";
s7.Tag = "";

s8 = addspecies(c1, "TDBp_ugperkg");
s8.Constant = false;
s8.Value = 0;
s8.Units = "";
s8.BoundaryCondition = false;
s8.Notes = "";
s8.Tag = "";

s9 = addspecies(c1, "restTpb");
s9.Constant = false;
s9.Value = 2000000;
s9.Units = "";
s9.BoundaryCondition = false;
s9.Notes = "";
s9.Tag = "";

s10 = addspecies(c1, "Bpb");
s10.Constant = false;
s10.Value = 0;
s10.Units = "";
s10.BoundaryCondition = false;
s10.Notes = "";
s10.Tag = "";

s11 = addspecies(c1, "TDBc_ugperml");
s11.Constant = false;
s11.Value = 0;
s11.Units = "";
s11.BoundaryCondition = false;
s11.Notes = "";
s11.Tag = "";

s12 = addspecies(c1, "TDBt_ugperml");
s12.Constant = false;
s12.Value = 0;
s12.Units = "";
s12.BoundaryCondition = false;
s12.Notes = "";
s12.Tag = "";

s13 = addspecies(c1, "drugTpbact");
s13.Constant = false;
s13.Value = 0;
s13.Units = "";
s13.BoundaryCondition = false;
s13.Notes = "";
s13.Tag = "";

s14 = addspecies(c1, "drugBpbkill");
s14.Constant = false;
s14.Value = 0;
s14.Units = "";
s14.BoundaryCondition = false;
s14.Notes = "";
s14.Tag = "";

s15 = addspecies(c1, "restTtiss_perml");
s15.Constant = false;
s15.Value = 0;
s15.Units = "";
s15.BoundaryCondition = false;
s15.Notes = "";
s15.Tag = "";

s16 = addspecies(c1, "actTtiss_perml");
s16.Constant = false;
s16.Value = 0;
s16.Units = "";
s16.BoundaryCondition = false;
s16.Notes = "";
s16.Tag = "";

s17 = addspecies(c1, "Btiss_perml");
s17.Constant = false;
s17.Value = 0;
s17.Units = "";
s17.BoundaryCondition = false;
s17.Notes = "";
s17.Tag = "";

s18 = addspecies(c1, "restTpb_perml");
s18.Constant = false;
s18.Value = 0;
s18.Units = "";
s18.BoundaryCondition = false;
s18.Notes = "";
s18.Tag = "";

s19 = addspecies(c1, "actTpb_perml");
s19.Constant = false;
s19.Value = 0;
s19.Units = "";
s19.BoundaryCondition = false;
s19.Notes = "";
s19.Tag = "";

s20 = addspecies(c1, "Bpb_perml");
s20.Constant = false;
s20.Value = 0;
s20.Units = "";
s20.BoundaryCondition = false;
s20.Notes = "";
s20.Tag = "";

s21 = addspecies(c1, "BTrratio_pb");
s21.Constant = false;
s21.Value = 0;
s21.Units = "";
s21.BoundaryCondition = false;
s21.Notes = "";
s21.Tag = "";

s22 = addspecies(c1, "BTrratio_tiss");
s22.Constant = false;
s22.Value = 0;
s22.Units = "";
s22.BoundaryCondition = false;
s22.Notes = "";
s22.Tag = "";

s23 = addspecies(c1, "TaBratio_pb");
s23.Constant = false;
s23.Value = 0;
s23.Units = "";
s23.BoundaryCondition = false;
s23.Notes = "";
s23.Tag = "";

s24 = addspecies(c1, "TaBratio_tiss");
s24.Constant = false;
s24.Value = 0;
s24.Units = "";
s24.BoundaryCondition = false;
s24.Notes = "";
s24.Tag = "";

s25 = addspecies(c1, "Tafraction_pb");
s25.Constant = false;
s25.Value = 0;
s25.Units = "";
s25.BoundaryCondition = false;
s25.Notes = "";
s25.Tag = "";

s26 = addspecies(c1, "totTpb_perml");
s26.Constant = false;
s26.Value = 0;
s26.Units = "";
s26.BoundaryCondition = false;
s26.Notes = "";
s26.Tag = "";

s27 = addspecies(c1, "totTtiss_perml");
s27.Constant = false;
s27.Value = 0;
s27.Units = "";
s27.BoundaryCondition = false;
s27.Notes = "";
s27.Tag = "";

s28 = addspecies(c1, "act0Ttiss_perml");
s28.Constant = false;
s28.Value = 0;
s28.Units = "";
s28.BoundaryCondition = false;
s28.Notes = "";
s28.Tag = "";

s29 = addspecies(c1, "act0Tpb_perml");
s29.Constant = false;
s29.Value = 0;
s29.Units = "";
s29.BoundaryCondition = false;
s29.Notes = "";
s29.Tag = "";

s30 = addspecies(c1, "Tafraction_tiss");
s30.Constant = false;
s30.Value = 0;
s30.Units = "";
s30.BoundaryCondition = false;
s30.Notes = "";
s30.Tag = "";

s31 = addspecies(c1, "BAFF");
s31.Constant = false;
s31.Value = 0.1;
s31.Units = "";
s31.BoundaryCondition = false;
s31.Notes = "";
s31.Tag = "";

s32 = addspecies(c1, "Baffconsumption");
s32.Constant = false;
s32.Value = 0;
s32.Units = "";
s32.BoundaryCondition = false;
s32.Notes = "";
s32.Tag = "";

s33 = addspecies(c1, "act0Tpb");
s33.Constant = false;
s33.Value = 0;
s33.Units = "";
s33.BoundaryCondition = false;
s33.Notes = "";
s33.Tag = "";

s34 = addspecies(c1, "act0Ttiss");
s34.Constant = false;
s34.Value = 0;
s34.Units = "";
s34.BoundaryCondition = false;
s34.Notes = "";
s34.Tag = "";

s35 = addspecies(c1, "injection_effect");
s35.Constant = false;
s35.Value = 0;
s35.Units = "";
s35.BoundaryCondition = false;
s35.Notes = "";
s35.Tag = "";

s36 = addspecies(c1, "Btiss2");
s36.Constant = false;
s36.Value = 0;
s36.Units = "";
s36.BoundaryCondition = false;
s36.Notes = "";
s36.Tag = "";

s37 = addspecies(c1, "act0Tpb_1");
s37.Constant = false;
s37.Value = 0;
s37.Units = "";
s37.BoundaryCondition = false;
s37.Notes = "";
s37.Tag = "";

s38 = addspecies(c1, "Bpb_1");
s38.Constant = false;
s38.Value = 0;
s38.Units = "";
s38.BoundaryCondition = false;
s38.Notes = "";
s38.Tag = "";

s39 = addspecies(c1, "restTpb_1");
s39.Constant = false;
s39.Value = 2000000;
s39.Units = "";
s39.BoundaryCondition = false;
s39.Notes = "";
s39.Tag = "";

s40 = addspecies(c1, "actTpb_1");
s40.Constant = false;
s40.Value = 0;
s40.Units = "";
s40.BoundaryCondition = false;
s40.Notes = "";
s40.Tag = "";

s41 = addspecies(c1, "act0Ttiss2");
s41.Constant = false;
s41.Value = 0;
s41.Units = "";
s41.BoundaryCondition = false;
s41.Notes = "";
s41.Tag = "";

s42 = addspecies(c1, "restTtiss2");
s42.Constant = false;
s42.Value = 0;
s42.Units = "";
s42.BoundaryCondition = false;
s42.Notes = "";
s42.Tag = "";

s43 = addspecies(c1, "TDBt2_ugperml");
s43.Constant = false;
s43.Value = 0;
s43.Units = "";
s43.BoundaryCondition = false;
s43.Notes = "";
s43.Tag = "";

s44 = addspecies(c1, "drugTtissact2");
s44.Constant = false;
s44.Value = 0;
s44.Units = "";
s44.BoundaryCondition = false;
s44.Notes = "";
s44.Tag = "";

s45 = addspecies(c1, "drugBtisskill2");
s45.Constant = false;
s45.Value = 0;
s45.Units = "";
s45.BoundaryCondition = false;
s45.Notes = "";
s45.Tag = "";

s46 = addspecies(c1, "actTtiss2");
s46.Constant = false;
s46.Value = 0;
s46.Units = "";
s46.BoundaryCondition = false;
s46.Notes = "";
s46.Tag = "";

s47 = addspecies(c1, "Tafraction_tiss2");
s47.Constant = false;
s47.Value = 0;
s47.Units = "";
s47.BoundaryCondition = false;
s47.Notes = "";
s47.Tag = "";

s48 = addspecies(c1, "act0Ttiss2_perml");
s48.Constant = false;
s48.Value = 0;
s48.Units = "";
s48.BoundaryCondition = false;
s48.Notes = "";
s48.Tag = "";

s49 = addspecies(c1, "totTtiss2_perml");
s49.Constant = false;
s49.Value = 0;
s49.Units = "";
s49.BoundaryCondition = false;
s49.Notes = "";
s49.Tag = "";

s50 = addspecies(c1, "TaBratio_tiss2");
s50.Constant = false;
s50.Value = 0;
s50.Units = "";
s50.BoundaryCondition = false;
s50.Notes = "";
s50.Tag = "";

s51 = addspecies(c1, "BTrratio_tiss2");
s51.Constant = false;
s51.Value = 0;
s51.Units = "";
s51.BoundaryCondition = false;
s51.Notes = "";
s51.Tag = "";

s52 = addspecies(c1, "Btiss2_perml");
s52.Constant = false;
s52.Value = 0;
s52.Units = "";
s52.BoundaryCondition = false;
s52.Notes = "";
s52.Tag = "";

s53 = addspecies(c1, "actTtiss2_perml");
s53.Constant = false;
s53.Value = 0;
s53.Units = "";
s53.BoundaryCondition = false;
s53.Notes = "";
s53.Tag = "";

s54 = addspecies(c1, "restTtiss2_perml");
s54.Constant = false;
s54.Value = 0;
s54.Units = "";
s54.BoundaryCondition = false;
s54.Notes = "";
s54.Tag = "";

s55 = addspecies(c1, "RTXc_ugperkg");
s55.Constant = false;
s55.Value = 0;
s55.Units = "";
s55.BoundaryCondition = false;
s55.Notes = "";
s55.Tag = "";

s56 = addspecies(c1, "RTXp_ugperkg");
s56.Constant = false;
s56.Value = 0;
s56.Units = "";
s56.BoundaryCondition = false;
s56.Notes = "";
s56.Tag = "";

s57 = addspecies(c1, "RTXc_ugperml");
s57.Constant = false;
s57.Value = 0;
s57.Units = "";
s57.BoundaryCondition = false;
s57.Notes = "";
s57.Tag = "";

s58 = addspecies(c1, "RTXt_ugperml");
s58.Constant = false;
s58.Value = 0;
s58.Units = "";
s58.BoundaryCondition = false;
s58.Notes = "";
s58.Tag = "";

s59 = addspecies(c1, "RTXt2_ugperml");
s59.Constant = false;
s59.Value = 0;
s59.Units = "";
s59.BoundaryCondition = false;
s59.Notes = "";
s59.Tag = "";

s60 = addspecies(c1, "drug_effect");
s60.Constant = false;
s60.Value = 0;
s60.Units = "";
s60.BoundaryCondition = false;
s60.Notes = "";
s60.Tag = "";

s61 = addspecies(c1, "B1920tiss3");
s61.Constant = false;
s61.Value = 0;
s61.Units = "";
s61.BoundaryCondition = false;
s61.Notes = "";
s61.Tag = "";

s62 = addspecies(c1, "B19no20tiss3");
s62.Constant = false;
s62.Value = 0;
s62.Units = "";
s62.BoundaryCondition = false;
s62.Notes = "";
s62.Tag = "";

s63 = addspecies(c1, "restTtiss3");
s63.Constant = false;
s63.Value = 0;
s63.Units = "";
s63.BoundaryCondition = false;
s63.Notes = "";
s63.Tag = "";

s64 = addspecies(c1, "act0Ttiss3");
s64.Constant = false;
s64.Value = 0;
s64.Units = "";
s64.BoundaryCondition = false;
s64.Notes = "";
s64.Tag = "";

s65 = addspecies(c1, "actTtiss3");
s65.Constant = false;
s65.Value = 0;
s65.Units = "";
s65.BoundaryCondition = false;
s65.Notes = "";
s65.Tag = "";

s66 = addspecies(c1, "drugTtissact3");
s66.Constant = false;
s66.Value = 0;
s66.Units = "";
s66.BoundaryCondition = false;
s66.Notes = "";
s66.Tag = "";

s67 = addspecies(c1, "TDBt3_ugperml");
s67.Constant = false;
s67.Value = 0;
s67.Units = "";
s67.BoundaryCondition = false;
s67.Notes = "";
s67.Tag = "";

s68 = addspecies(c1, "B19Trratio_tiss3");
s68.Constant = false;
s68.Value = 0;
s68.Units = "";
s68.BoundaryCondition = false;
s68.Notes = "";
s68.Tag = "";

s69 = addspecies(c1, "B19tiss3");
s69.Constant = false;
s69.Value = 0;
s69.Units = "";
s69.BoundaryCondition = false;
s69.Notes = "";
s69.Tag = "";

s70 = addspecies(c1, "restTtiss3_perml");
s70.Constant = false;
s70.Value = 0;
s70.Units = "";
s70.BoundaryCondition = false;
s70.Notes = "";
s70.Tag = "";

s71 = addspecies(c1, "act0Ttiss3_perml");
s71.Constant = false;
s71.Value = 0;
s71.Units = "";
s71.BoundaryCondition = false;
s71.Notes = "";
s71.Tag = "";

s72 = addspecies(c1, "actTtiss3_perml");
s72.Constant = false;
s72.Value = 0;
s72.Units = "";
s72.BoundaryCondition = false;
s72.Notes = "";
s72.Tag = "";

s73 = addspecies(c1, "B19tiss3_perml");
s73.Constant = false;
s73.Value = 0;
s73.Units = "";
s73.BoundaryCondition = false;
s73.Notes = "";
s73.Tag = "";

s74 = addspecies(c1, "Tafraction_tiss3");
s74.Constant = false;
s74.Value = 0;
s74.Units = "";
s74.BoundaryCondition = false;
s74.Notes = "";
s74.Tag = "";

s75 = addspecies(c1, "totTtiss3_perml");
s75.Constant = false;
s75.Value = 0;
s75.Units = "";
s75.BoundaryCondition = false;
s75.Notes = "";
s75.Tag = "";

s76 = addspecies(c1, "TaB19ratio_tiss3");
s76.Constant = false;
s76.Value = 0;
s76.Units = "";
s76.BoundaryCondition = false;
s76.Notes = "";
s76.Tag = "";

s77 = addspecies(c1, "B19TtotRatio_tiss3");
s77.Constant = false;
s77.Value = 0;
s77.Units = "";
s77.BoundaryCondition = false;
s77.Notes = "";
s77.Tag = "";

s78 = addspecies(c1, "RTXt3_ugperml");
s78.Constant = false;
s78.Value = 0;
s78.Units = "";
s78.BoundaryCondition = false;
s78.Notes = "";
s78.Tag = "";

s79 = addspecies(c1, "drugBtisskill3");
s79.Constant = false;
s79.Value = 0;
s79.Units = "";
s79.BoundaryCondition = false;
s79.Notes = "";
s79.Tag = "";

s80 = addspecies(c1, "TaB1920ratio_tiss3");
s80.Constant = false;
s80.Value = 0;
s80.Units = "";
s80.BoundaryCondition = false;
s80.Notes = "";
s80.Tag = "";

s81 = addspecies(c1, "B1920tiss3_perml");
s81.Constant = false;
s81.Value = 0;
s81.Units = "";
s81.BoundaryCondition = false;
s81.Notes = "";
s81.Tag = "";

s82 = addspecies(c1, "B1920Trratio_tiss3");
s82.Constant = false;
s82.Value = 0;
s82.Units = "";
s82.BoundaryCondition = false;
s82.Notes = "";
s82.Tag = "";

s83 = addspecies(c1, "B19no20tiss3_perml");
s83.Constant = false;
s83.Value = 0;
s83.Units = "";
s83.BoundaryCondition = false;
s83.Notes = "";
s83.Tag = "";

s84 = addspecies(c1, "Blinc_ug");
s84.Constant = false;
s84.Value = 0;
s84.Units = "";
s84.BoundaryCondition = false;
s84.Notes = "";
s84.Tag = "";

s85 = addspecies(c1, "Blinc_ngperml");
s85.Constant = false;
s85.Value = 0;
s85.Units = "";
s85.BoundaryCondition = false;
s85.Notes = "";
s85.Tag = "";

s86 = addspecies(c1, "BlinTpbact");
s86.Constant = false;
s86.Value = 0;
s86.Units = "";
s86.BoundaryCondition = false;
s86.Notes = "";
s86.Tag = "";

s87 = addspecies(c1, "BlinBpbkill");
s87.Constant = false;
s87.Value = 0;
s87.Units = "";
s87.BoundaryCondition = false;
s87.Notes = "";
s87.Tag = "";

s88 = addspecies(c1, "BlinTtissact");
s88.Constant = false;
s88.Value = 0;
s88.Units = "";
s88.BoundaryCondition = false;
s88.Notes = "";
s88.Tag = "";

s89 = addspecies(c1, "BlinBtisskill");
s89.Constant = false;
s89.Value = 0;
s89.Units = "";
s89.BoundaryCondition = false;
s89.Notes = "";
s89.Tag = "";

s90 = addspecies(c1, "Blint_ngperml");
s90.Constant = false;
s90.Value = 0;
s90.Units = "";
s90.BoundaryCondition = false;
s90.Notes = "";
s90.Tag = "";

s91 = addspecies(c1, "Blint2_ngperml");
s91.Constant = false;
s91.Value = 0;
s91.Units = "";
s91.BoundaryCondition = false;
s91.Notes = "";
s91.Tag = "";

s92 = addspecies(c1, "BlinBtisskill2");
s92.Constant = false;
s92.Value = 0;
s92.Units = "";
s92.BoundaryCondition = false;
s92.Notes = "";
s92.Tag = "";

s93 = addspecies(c1, "BlinTtissact2");
s93.Constant = false;
s93.Value = 0;
s93.Units = "";
s93.BoundaryCondition = false;
s93.Notes = "";
s93.Tag = "";

s94 = addspecies(c1, "BlinBtisskill3");
s94.Constant = false;
s94.Value = 0;
s94.Units = "";
s94.BoundaryCondition = false;
s94.Notes = "";
s94.Tag = "";

s95 = addspecies(c1, "Blint3_ngperml");
s95.Constant = false;
s95.Value = 0;
s95.Units = "";
s95.BoundaryCondition = false;
s95.Notes = "";
s95.Tag = "";

s96 = addspecies(c1, "BlinTtissact3");
s96.Constant = false;
s96.Value = 0;
s96.Units = "";
s96.BoundaryCondition = false;
s96.Notes = "";
s96.Tag = "";

s97 = addspecies(c1, "Bpb_norm");
s97.Constant = false;
s97.Value = 0;
s97.Units = "";
s97.BoundaryCondition = false;
s97.Notes = "";
s97.Tag = "";

s98 = addspecies(c1, "Tafraction_pb_init");
s98.Constant = false;
s98.Value = 0;
s98.Units = "";
s98.BoundaryCondition = false;
s98.Notes = "";
s98.Tag = "";

s99 = addspecies(c1, "totTtiss");
s99.Constant = false;
s99.Value = 0;
s99.Units = "";
s99.BoundaryCondition = false;
s99.Notes = "";
s99.Tag = "";

s100 = addspecies(c1, "totTtiss2");
s100.Constant = false;
s100.Value = 0;
s100.Units = "";
s100.BoundaryCondition = false;
s100.Notes = "";
s100.Tag = "";

s101 = addspecies(c1, "totTtiss3");
s101.Constant = false;
s101.Value = 0;
s101.Units = "";
s101.BoundaryCondition = false;
s101.Notes = "";
s101.Tag = "";

s102 = addspecies(c1, "restTtumor");
s102.Constant = false;
s102.Value = 0;
s102.Units = "";
s102.BoundaryCondition = false;
s102.Notes = "";
s102.Tag = "";

s103 = addspecies(c1, "actTtumor");
s103.Constant = false;
s103.Value = 0;
s103.Units = "";
s103.BoundaryCondition = false;
s103.Notes = "";
s103.Tag = "";

s104 = addspecies(c1, "Btumor");
s104.Constant = false;
s104.Value = 0;
s104.Units = "";
s104.BoundaryCondition = false;
s104.Notes = "";
s104.Tag = "";

s105 = addspecies(c1, "act0Ttumor");
s105.Constant = false;
s105.Value = 0;
s105.Units = "";
s105.BoundaryCondition = false;
s105.Notes = "";
s105.Tag = "";

s106 = addspecies(c1, "drugTtumoract");
s106.Constant = false;
s106.Value = 0;
s106.Units = "";
s106.BoundaryCondition = false;
s106.Notes = "";
s106.Tag = "";

s107 = addspecies(c1, "drugBtumorkill");
s107.Constant = false;
s107.Value = 0;
s107.Units = "";
s107.BoundaryCondition = false;
s107.Notes = "";
s107.Tag = "";

s108 = addspecies(c1, "BlinTtumoract");
s108.Constant = false;
s108.Value = 0;
s108.Units = "";
s108.BoundaryCondition = false;
s108.Notes = "";
s108.Tag = "";

s109 = addspecies(c1, "BlinBtumorkill");
s109.Constant = false;
s109.Value = 0;
s109.Units = "";
s109.BoundaryCondition = false;
s109.Notes = "";
s109.Tag = "";

s110 = addspecies(c1, "TDBtumor_ugperml");
s110.Constant = false;
s110.Value = 0;
s110.Units = "";
s110.BoundaryCondition = false;
s110.Notes = "";
s110.Tag = "";

s111 = addspecies(c1, "restTtumor_perml");
s111.Constant = false;
s111.Value = 0;
s111.Units = "";
s111.BoundaryCondition = false;
s111.Notes = "";
s111.Tag = "";

s112 = addspecies(c1, "act0Ttumor_perml");
s112.Constant = false;
s112.Value = 0;
s112.Units = "";
s112.BoundaryCondition = false;
s112.Notes = "";
s112.Tag = "";

s113 = addspecies(c1, "actTtumor_perml");
s113.Constant = false;
s113.Value = 0;
s113.Units = "";
s113.BoundaryCondition = false;
s113.Notes = "";
s113.Tag = "";

s114 = addspecies(c1, "Btumor_perml");
s114.Constant = false;
s114.Value = 0;
s114.Units = "";
s114.BoundaryCondition = false;
s114.Notes = "";
s114.Tag = "";

s115 = addspecies(c1, "Tafraction_tumor");
s115.Constant = false;
s115.Value = 0;
s115.Units = "";
s115.BoundaryCondition = false;
s115.Notes = "";
s115.Tag = "";

s116 = addspecies(c1, "BTrratio_tumor");
s116.Constant = false;
s116.Value = 0;
s116.Units = "";
s116.BoundaryCondition = false;
s116.Notes = "";
s116.Tag = "";

s117 = addspecies(c1, "TaBratio_tumor");
s117.Constant = false;
s117.Value = 0;
s117.Units = "";
s117.BoundaryCondition = false;
s117.Notes = "";
s117.Tag = "";

s118 = addspecies(c1, "totTtumor");
s118.Constant = false;
s118.Value = 0;
s118.Units = "";
s118.BoundaryCondition = false;
s118.Notes = "";
s118.Tag = "";

s119 = addspecies(c1, "totTtumor_perml");
s119.Constant = false;
s119.Value = 0;
s119.Units = "";
s119.BoundaryCondition = false;
s119.Notes = "";
s119.Tag = "";

s120 = addspecies(c1, "Blintumor_ngperml");
s120.Constant = false;
s120.Value = 0;
s120.Units = "";
s120.BoundaryCondition = false;
s120.Notes = "";
s120.Tag = "";

s121 = addspecies(c1, "RTXtumor_ugperml");
s121.Constant = false;
s121.Value = 0;
s121.Units = "";
s121.BoundaryCondition = false;
s121.Notes = "";
s121.Tag = "";

s122 = addspecies(c1, "BTtotRatio_tiss");
s122.Constant = false;
s122.Value = 0;
s122.Units = "";
s122.BoundaryCondition = false;
s122.Notes = "";
s122.Tag = "";

s123 = addspecies(c1, "BTtotRatio_tiss2");
s123.Constant = false;
s123.Value = 0;
s123.Units = "";
s123.BoundaryCondition = false;
s123.Notes = "";
s123.Tag = "";

s124 = addspecies(c1, "TDBsc_ugperkg");
s124.Constant = false;
s124.Value = 0;
s124.Units = "";
s124.BoundaryCondition = false;
s124.Notes = "";
s124.Tag = "";

s125 = addspecies(c1, "TDBc_ugperml_AUC");
s125.Constant = false;
s125.Value = 0;
s125.Units = "";
s125.BoundaryCondition = false;
s125.Notes = "";
s125.Tag = "";

s126 = addspecies(c1, "IL6pb");
s126.Constant = false;
s126.Value = 0;
s126.Units = "";
s126.BoundaryCondition = false;
s126.Notes = "";
s126.Tag = "";

s127 = addspecies(c1, "IL6tiss");
s127.Constant = false;
s127.Value = 0;
s127.Units = "";
s127.BoundaryCondition = false;
s127.Notes = "";
s127.Tag = "";

s128 = addspecies(c1, "IL6tiss2");
s128.Constant = false;
s128.Value = 0;
s128.Units = "";
s128.BoundaryCondition = false;
s128.Notes = "";
s128.Tag = "";

s129 = addspecies(c1, "IL6tiss3");
s129.Constant = false;
s129.Value = 0;
s129.Units = "";
s129.BoundaryCondition = false;
s129.Notes = "";
s129.Tag = "";

s130 = addspecies(c1, "IL6tumor");
s130.Constant = false;
s130.Value = 0;
s130.Units = "";
s130.BoundaryCondition = false;
s130.Notes = "";
s130.Tag = "";

s131 = addspecies(c1, "IL6combo");
s131.Constant = false;
s131.Value = 0;
s131.Units = "";
s131.BoundaryCondition = false;
s131.Notes = "";
s131.Tag = "";

%% Create reactions.
r1 = addreaction(m1, "Btiss -> 2 Btiss");
r1.ReactionRate = "kBprolif*KBp*Bpbref_perml*Vtissue*(max(0, 1 - Btiss/(KBp*Bpbref_perml*Vtissue)))^1";
r1.Active = true;
r1.Name = "Btiss_prolf";
r1.Notes = "";
r1.Tag = "";
klaw1 = addkineticlaw(r1, 'Unknown');
klaw1.ParameterVariableNames = strings(0,1);
klaw1.SpeciesVariableNames = strings(0,1);
klaw1.Name = "";
klaw1.Notes = "";
klaw1.Tag = "";

r2 = addreaction(m1, "restTtiss + drugTtissact + BlinTtissact <-> actTtiss + drugTtissact + BlinTtissact");
r2.ReactionRate = "kTact*((drugTtissact+BlinTtissact)*restTtiss-fTadeact*actTtiss)";
r2.Active = true;
r2.Name = "Ttiss_activn";
r2.Notes = "";
r2.Tag = "";
klaw2 = addkineticlaw(r2, 'Unknown');
klaw2.ParameterVariableNames = strings(0,1);
klaw2.SpeciesVariableNames = strings(0,1);
klaw2.Name = "";
klaw2.Notes = "";
klaw2.Tag = "";

r3 = addreaction(m1, "Btiss -> null");
r3.ReactionRate = "0*fBprolif*kBprolif*Btiss*kBapop_cll";
r3.Active = true;
r3.Name = "Btiss_apop";
r3.Notes = "";
r3.Tag = "";
klaw3 = addkineticlaw(r3, 'Unknown');
klaw3.ParameterVariableNames = strings(0,1);
klaw3.SpeciesVariableNames = strings(0,1);
klaw3.Name = "";
klaw3.Notes = "";
klaw3.Tag = "";

r4 = addreaction(m1, "TDBc_ugperkg <-> null");
r4.ReactionRate = "(Cl_tdb+Vm_tdb/(Km_tdb+TDBc_ugperkg/Vc_tdb)/BW)/Vc_tdb*TDBc_ugperkg";
r4.Active = true;
r4.Name = "`";
r4.Notes = "";
r4.Tag = "";
klaw4 = addkineticlaw(r4, 'Unknown');
klaw4.ParameterVariableNames = strings(0,1);
klaw4.SpeciesVariableNames = strings(0,1);
klaw4.Name = "";
klaw4.Notes = "";
klaw4.Tag = "";

r5 = addreaction(m1, "TDBc_ugperkg <-> TDBp_ugperkg");
r5.ReactionRate = "Cld_tdb*(TDBc_ugperkg/Vc_tdb-TDBp_ugperkg/Vp_tdb)";
r5.Active = true;
r5.Name = "TDBdistr";
r5.Notes = "";
r5.Tag = "";
klaw5 = addkineticlaw(r5, 'Unknown');
klaw5.ParameterVariableNames = strings(0,1);
klaw5.SpeciesVariableNames = strings(0,1);
klaw5.Name = "";
klaw5.Notes = "";
klaw5.Tag = "";

r6 = addreaction(m1, "restTtiss + Btiss + drugTtissact + TDBt_ugperml -> restTtiss + Btiss + drugTtissact + TDBt_ugperml + restTtiss");
r6.ReactionRate = "0";
r6.Active = true;
r6.Name = "UnnamedReaction_6";
r6.Notes = "";
r6.Tag = "";

r7 = addreaction(m1, "actTtiss + Btiss + drugBtisskill + TDBt_ugperml -> actTtiss + Btiss + drugBtisskill + TDBt_ugperml");
r7.ReactionRate = "0";
r7.Active = true;
r7.Name = "UnnamedReaction_7";
r7.Notes = "";
r7.Tag = "";

r8 = addreaction(m1, "restTpb + injection_effect + drug_effect <-> restTtiss + injection_effect + drug_effect");
r8.ReactionRate = "tissue1on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp-restTtiss/Vtissue)";
r8.Active = true;
r8.Name = "restTtraffic";
r8.Notes = "";
r8.Tag = "";
klaw6 = addkineticlaw(r8, 'Unknown');
klaw6.ParameterVariableNames = strings(0,1);
klaw6.SpeciesVariableNames = strings(0,1);
klaw6.Name = "";
klaw6.Notes = "";
klaw6.Tag = "";

r9 = addreaction(m1, "totTpb_perml -> restTpb + totTpb_perml");
r9.ReactionRate = "kTgen*Vpb*Trpbref_perml*(fTgenbl+max(0,1-totTpb_perml/Trpbref_perml)^2)";
r9.Active = true;
r9.Name = "Tgen";
r9.Notes = "";
r9.Tag = "";
klaw7 = addkineticlaw(r9, 'Unknown');
klaw7.ParameterVariableNames = strings(0,1);
klaw7.SpeciesVariableNames = strings(0,1);
klaw7.Name = "";
klaw7.Notes = "";
klaw7.Tag = "";

r10 = addreaction(m1, "Bpb -> Btiss");
r10.ReactionRate = "tissue1on*fBexit*kTaexit*Vpb*max(0, Bpb/Vpb*KBp-Btiss/Vtissue)";
r10.Active = true;
r10.Name = "Btraffic";
r10.Notes = "";
r10.Tag = "";
klaw8 = addkineticlaw(r10, 'Unknown');
klaw8.ParameterVariableNames = strings(0,1);
klaw8.SpeciesVariableNames = strings(0,1);
klaw8.Name = "";
klaw8.Notes = "";
klaw8.Tag = "";

r11 = addreaction(m1, "restTpb + drugTpbact + BlinTpbact <-> actTpb + drugTpbact + BlinTpbact");
r11.ReactionRate = "kTact*((drugTpbact+BlinTpbact)*restTpb-fTadeact*actTpb)";
r11.Active = true;
r11.Name = "Tpb_activn";
r11.Notes = "";
r11.Tag = "";
klaw9 = addkineticlaw(r11, 'Unknown');
klaw9.ParameterVariableNames = strings(0,1);
klaw9.SpeciesVariableNames = strings(0,1);
klaw9.Name = "";
klaw9.Notes = "";
klaw9.Tag = "";

r12 = addreaction(m1, "Bpb -> 2 Bpb");
r12.ReactionRate = "tissue3on*fBtissue3_v1*kBapop*Bpb*kBapop_cll";
r12.Active = true;
r12.Name = "Bpb_prolif";
r12.Notes = "";
r12.Tag = "";

r13 = addreaction(m1, "Bpb -> null");
r13.ReactionRate = "kBapop*Bpb*kBapop_cll";
r13.Active = true;
r13.Name = "Bpb_apop";
r13.Notes = "";
r13.Tag = "";

r14 = addreaction(m1, "TDBc_ugperml + drugTpbact + Bpb -> TDBc_ugperml + drugTpbact + Bpb");
r14.ReactionRate = "0";
r14.Active = true;
r14.Name = "UnnamedReaction_14";
r14.Notes = "";
r14.Tag = "";

r15 = addreaction(m1, "TDBc_ugperml + actTpb + drugBpbkill + Bpb -> TDBc_ugperml + actTpb + drugBpbkill + Bpb");
r15.ReactionRate = "0";
r15.Active = true;
r15.Name = "UnnamedReaction_15";
r15.Notes = "";
r15.Tag = "";

r16 = addreaction(m1, "actTpb + drugBpbkill + Bpb + RTXc_ugperml + BlinBpbkill -> actTpb + drugBpbkill + RTXc_ugperml + BlinBpbkill");
r16.ReactionRate = "kBkill*(drugBpbkill + BlinBpbkill + RTXc_ugperml/(Kmkill_rtx/1000+RTXc_ugperml))*Bpb";
r16.Active = true;
r16.Name = "Bpb_kill";
r16.Notes = "";
r16.Tag = "";
klaw10 = addkineticlaw(r16, 'Unknown');
klaw10.ParameterVariableNames = strings(0,1);
klaw10.SpeciesVariableNames = strings(0,1);
klaw10.Name = "";
klaw10.Notes = "";
klaw10.Tag = "";

r17 = addreaction(m1, "actTtiss + drugBtisskill + Btiss + RTXt_ugperml + BlinBtisskill -> actTtiss + drugBtisskill + RTXt_ugperml + BlinBtisskill");
r17.ReactionRate = "fBkill*kBkill*(drugBtisskill + BlinBtisskill + RTXt_ugperml/(Kmkill_rtx/1000+RTXt_ugperml))*Btiss";
r17.Active = true;
r17.Name = "Btiss_kill";
r17.Notes = "";
r17.Tag = "";

r18 = addreaction(m1, "restTpb -> null");
r18.ReactionRate = "kTgen*fTgenbl*restTpb + fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb)";
r18.Active = true;
r18.Name = "restTpb_apop";
r18.Notes = "";
r18.Tag = "";
klaw11 = addkineticlaw(r18, 'Unknown');
klaw11.ParameterVariableNames = strings(0,1);
klaw11.SpeciesVariableNames = strings(0,1);
klaw11.Name = "";
klaw11.Notes = "";
klaw11.Tag = "";

r19 = addreaction(m1, "actTpb + injection_effect + drug_effect <-> actTtiss + injection_effect + drug_effect");
r19.ReactionRate = "tissue1on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp*fTap-actTtiss/Vtissue)";
r19.Active = true;
r19.Name = "actTtraffic";
r19.Notes = "";
r19.Tag = "";

r20 = addreaction(m1, "actTpb -> null");
r20.ReactionRate = "kTaapop*(actTpb + fAICD*actTpb^2/(Vpb*Trpbref_perml))";
r20.Active = true;
r20.Name = "actTpb_apop";
r20.Notes = "";
r20.Tag = "";
klaw12 = addkineticlaw(r20, 'Unknown');
klaw12.ParameterVariableNames = strings(0,1);
klaw12.SpeciesVariableNames = strings(0,1);
klaw12.Name = "";
klaw12.Notes = "";
klaw12.Tag = "";

r21 = addreaction(m1, "BAFF + Btiss + Bpb + Btiss2 -> Btiss + Bpb + Btiss2");
r21.ReactionRate = "log(2)/(thBAFF/24/60)*(fBAFFo*BAFF/BAFFo+(1-fBAFFo)*(Btiss+ Bpb + Btiss2)/(Bpbref_perml*(Vpb+KBp*Vtissue + KBp2*Vtissue2)))";
r21.Active = true;
r21.Name = "BAFFcons";
r21.Notes = "";
r21.Tag = "";
klaw13 = addkineticlaw(r21, 'Unknown');
klaw13.ParameterVariableNames = strings(0,1);
klaw13.SpeciesVariableNames = strings(0,1);
klaw13.Name = "";
klaw13.Notes = "";
klaw13.Tag = "";

r22 = addreaction(m1, "null -> BAFF");
r22.ReactionRate = "log(2)/(thBAFF/24/60) * (fBAFFo+(1-fBAFFo)*Bpbo_perml/Bpbref_perml)";
r22.Active = true;
r22.Name = "BAFFprod";
r22.Notes = "";
r22.Tag = "";
klaw14 = addkineticlaw(r22, 'Unknown');
klaw14.ParameterVariableNames = strings(0,1);
klaw14.SpeciesVariableNames = strings(0,1);
klaw14.Name = "";
klaw14.Notes = "";
klaw14.Tag = "";

r23 = addreaction(m1, "act0Tpb -> restTpb");
r23.ReactionRate = "act0on*fTa0deact*act0Tpb";
r23.Active = true;
r23.Name = "act0Tpb_deact";
r23.Notes = "";
r23.Tag = "";
klaw15 = addkineticlaw(r23, 'Unknown');
klaw15.ParameterVariableNames = strings(0,1);
klaw15.SpeciesVariableNames = strings(0,1);
klaw15.Name = "";
klaw15.Notes = "";
klaw15.Tag = "";

r24 = addreaction(m1, "act0Tpb + drugTpbact + BlinTpbact <-> actTpb + drugTpbact + BlinTpbact");
r24.ReactionRate = "act0on*kTact*((drugTpbact+BlinTpbact)*act0Tpb-fTadeact*actTpb)";
r24.Active = true;
r24.Name = "act0Tpb_activn";
r24.Notes = "";
r24.Tag = "";

r25 = addreaction(m1, "act0Tpb -> null");
r25.ReactionRate = "fTa0apop*kTaapop*act0Tpb";
r25.Active = true;
r25.Name = "act0Tpb_apop";
r25.Notes = "";
r25.Tag = "";
klaw16 = addkineticlaw(r25, 'Unknown');
klaw16.ParameterVariableNames = strings(0,1);
klaw16.SpeciesVariableNames = strings(0,1);
klaw16.Name = "";
klaw16.Notes = "";
klaw16.Tag = "";

r26 = addreaction(m1, "act0Ttiss -> null");
r26.ReactionRate = "fTa0apop*kTaapop*act0Ttiss";
r26.Active = true;
r26.Name = "act0Ttiss_apop";
r26.Notes = "";
r26.Tag = "";
klaw17 = addkineticlaw(r26, 'Unknown');
klaw17.ParameterVariableNames = strings(0,1);
klaw17.SpeciesVariableNames = strings(0,1);
klaw17.Name = "";
klaw17.Notes = "";
klaw17.Tag = "";

r27 = addreaction(m1, "act0Ttiss + drugTtissact + BlinTtissact <-> actTtiss + drugTtissact + BlinTtissact");
r27.ReactionRate = "act0on*kTact*((drugTtissact+BlinTtissact)*act0Ttiss-fTadeact*actTtiss)";
r27.Active = true;
r27.Name = "act0Ttiss_activn";
r27.Notes = "";
r27.Tag = "";

r28 = addreaction(m1, "act0Ttiss -> restTtiss");
r28.ReactionRate = "act0on*fTa0deact*act0Ttiss";
r28.Active = true;
r28.Name = "act0Ttiss_deact";
r28.Notes = "";
r28.Tag = "";
klaw18 = addkineticlaw(r28, 'Unknown');
klaw18.ParameterVariableNames = strings(0,1);
klaw18.SpeciesVariableNames = strings(0,1);
klaw18.Name = "";
klaw18.Notes = "";
klaw18.Tag = "";

r29 = addreaction(m1, "act0Tpb + injection_effect + drug_effect <-> act0Ttiss + injection_effect + drug_effect");
r29.ReactionRate = "tissue1on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp-act0Ttiss/Vtissue)";
r29.Active = true;
r29.Name = "act0Ttraffic";
r29.Notes = "";
r29.Tag = "";
klaw19 = addkineticlaw(r29, 'Unknown');
klaw19.ParameterVariableNames = strings(0,1);
klaw19.SpeciesVariableNames = strings(0,1);
klaw19.Name = "";
klaw19.Notes = "";
klaw19.Tag = "";

r30 = addreaction(m1, "actTpb -> actTpb + act0Tpb");
r30.ReactionRate = "act0on*fTaprolif*kTprolif*actTpb";
r30.Active = true;
r30.Name = "actTpb_prolif";
r30.Notes = "";
r30.Tag = "";

r31 = addreaction(m1, "actTtiss -> actTtiss + act0Ttiss");
r31.ReactionRate = "act0on*fTaprolif*kTprolif*actTtiss";
r31.Active = true;
r31.Name = "actTtiss_prolif";
r31.Notes = "";
r31.Tag = "";

r32 = addreaction(m1, "injection_effect -> null");
r32.ReactionRate = "log(2)/tinjhalf*injection_effect";
r32.Active = true;
r32.Name = "UnnamedReaction_32";
r32.Notes = "";
r32.Tag = "";
klaw20 = addkineticlaw(r32, 'Unknown');
klaw20.ParameterVariableNames = strings(0,1);
klaw20.SpeciesVariableNames = strings(0,1);
klaw20.Name = "";
klaw20.Notes = "";
klaw20.Tag = "";

r33 = addreaction(m1, "actTtiss2 -> actTtiss2 + act0Ttiss2");
r33.ReactionRate = "act0on*fTaprolif*kTprolif*actTtiss2";
r33.Active = true;
r33.Name = "actTtiss_prolif_2";
r33.Notes = "";
r33.Tag = "";

r34 = addreaction(m1, "act0Tpb + injection_effect + drug_effect <-> act0Ttiss2 + injection_effect + drug_effect");
r34.ReactionRate = "tissue2on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp2-act0Ttiss2/Vtissue2)";
r34.Active = true;
r34.Name = "act0Ttraffic_2";
r34.Notes = "";
r34.Tag = "";
klaw21 = addkineticlaw(r34, 'Unknown');
klaw21.ParameterVariableNames = strings(0,1);
klaw21.SpeciesVariableNames = strings(0,1);
klaw21.Name = "";
klaw21.Notes = "";
klaw21.Tag = "";

r35 = addreaction(m1, "act0Ttiss2 -> restTtiss2");
r35.ReactionRate = "act0on*fTa0deact*act0Ttiss2";
r35.Active = true;
r35.Name = "act0Ttiss_deact_2";
r35.Notes = "";
r35.Tag = "";
klaw22 = addkineticlaw(r35, 'Unknown');
klaw22.ParameterVariableNames = strings(0,1);
klaw22.SpeciesVariableNames = strings(0,1);
klaw22.Name = "";
klaw22.Notes = "";
klaw22.Tag = "";

r36 = addreaction(m1, "act0Ttiss2 + drugTtissact2 + BlinTtissact2 <-> actTtiss2 + drugTtissact2 + BlinTtissact2");
r36.ReactionRate = "act0on*kTact*((drugTtissact2+BlinTtissact2)*act0Ttiss2-fTadeact*actTtiss2)";
r36.Active = true;
r36.Name = "act0Ttiss_activn_2";
r36.Notes = "";
r36.Tag = "";

r37 = addreaction(m1, "act0Ttiss2 -> null");
r37.ReactionRate = "fTa0apop*kTaapop*act0Ttiss2";
r37.Active = true;
r37.Name = "act0Ttiss_apop_2";
r37.Notes = "";
r37.Tag = "";
klaw23 = addkineticlaw(r37, 'Unknown');
klaw23.ParameterVariableNames = strings(0,1);
klaw23.SpeciesVariableNames = strings(0,1);
klaw23.Name = "";
klaw23.Notes = "";
klaw23.Tag = "";

r38 = addreaction(m1, "actTpb + injection_effect + drug_effect <-> actTtiss2 + injection_effect + drug_effect");
r38.ReactionRate = "tissue2on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp2*fTap-actTtiss2/Vtissue2)";
r38.Active = true;
r38.Name = "actTtraffic_2";
r38.Notes = "";
r38.Tag = "";

r39 = addreaction(m1, "actTtiss2 + drugBtisskill2 + Btiss2 + RTXt2_ugperml + BlinBtisskill2 -> actTtiss2 + drugBtisskill2 + RTXt2_ugperml + BlinBtisskill2");
r39.ReactionRate = "fBkill*kBkill*(drugBtisskill2+ BlinBtisskill2 + RTXt2_ugperml/(Kmkill_rtx/1000+RTXt2_ugperml))*Btiss2";
r39.Active = true;
r39.Name = "Btiss_kill_2";
r39.Notes = "";
r39.Tag = "";

r40 = addreaction(m1, "Bpb -> Btiss2");
r40.ReactionRate = "tissue2on*fBexit*kTaexit*Vpb*max(0, Bpb/Vpb*KBp2-Btiss2/Vtissue2)";
r40.Active = true;
r40.Name = "Btraffic_2";
r40.Notes = "";
r40.Tag = "";
klaw24 = addkineticlaw(r40, 'Unknown');
klaw24.ParameterVariableNames = strings(0,1);
klaw24.SpeciesVariableNames = strings(0,1);
klaw24.Name = "";
klaw24.Notes = "";
klaw24.Tag = "";

r41 = addreaction(m1, "restTpb + injection_effect + drug_effect <-> restTtiss2 + injection_effect + drug_effect");
r41.ReactionRate = "tissue2on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp2-restTtiss2/Vtissue2)";
r41.Active = true;
r41.Name = "restTtraffic_2";
r41.Notes = "";
r41.Tag = "";
klaw25 = addkineticlaw(r41, 'Unknown');
klaw25.ParameterVariableNames = strings(0,1);
klaw25.SpeciesVariableNames = strings(0,1);
klaw25.Name = "";
klaw25.Notes = "";
klaw25.Tag = "";

r42 = addreaction(m1, "actTtiss2 + Btiss2 + drugBtisskill2 + TDBt2_ugperml -> actTtiss2 + Btiss2 + drugBtisskill2 + TDBt2_ugperml");
r42.ReactionRate = "0";
r42.Active = true;
r42.Name = "UnnamedReaction_42";
r42.Notes = "";
r42.Tag = "";

r43 = addreaction(m1, "restTtiss2 + Btiss2 + drugTtissact2 + TDBt2_ugperml -> restTtiss2 + Btiss2 + drugTtissact2 + TDBt2_ugperml");
r43.ReactionRate = "0";
r43.Active = true;
r43.Name = "UnnamedReaction_43";
r43.Notes = "";
r43.Tag = "";

r44 = addreaction(m1, "Btiss2 -> null");
r44.ReactionRate = "0*fBprolif*kBprolif*Btiss2*kBapop_cll";
r44.Active = true;
r44.Name = "Btiss_apop_2";
r44.Notes = "";
r44.Tag = "";
klaw26 = addkineticlaw(r44, 'Unknown');
klaw26.ParameterVariableNames = strings(0,1);
klaw26.SpeciesVariableNames = strings(0,1);
klaw26.Name = "";
klaw26.Notes = "";
klaw26.Tag = "";

r45 = addreaction(m1, "restTtiss2 + drugTtissact2 + BlinTtissact2 <-> actTtiss2 + drugTtissact2 + BlinTtissact2");
r45.ReactionRate = "kTact*((drugTtissact2+BlinTtissact2)*restTtiss2-fTadeact*actTtiss2)";
r45.Active = true;
r45.Name = "Ttiss_activn_2";
r45.Notes = "";
r45.Tag = "";
klaw27 = addkineticlaw(r45, 'Unknown');
klaw27.ParameterVariableNames = strings(0,1);
klaw27.SpeciesVariableNames = strings(0,1);
klaw27.Name = "";
klaw27.Notes = "";
klaw27.Tag = "";

r46 = addreaction(m1, "Btiss2 -> 2 Btiss2");
r46.ReactionRate = "kBprolif*KBp2*Bpbref_perml*Vtissue2*(max(0, 1 - Btiss2/(KBp2*Bpbref_perml*Vtissue2)))^1";
r46.Active = true;
r46.Name = "Btiss_prolf_2";
r46.Notes = "";
r46.Tag = "";
klaw28 = addkineticlaw(r46, 'Unknown');
klaw28.ParameterVariableNames = strings(0,1);
klaw28.SpeciesVariableNames = strings(0,1);
klaw28.Name = "";
klaw28.Notes = "";
klaw28.Tag = "";

r47 = addreaction(m1, "RTXc_ugperkg <-> RTXp_ugperkg");
r47.ReactionRate = "Cld_rtx*(RTXc_ugperkg/Vc_rtx-RTXp_ugperkg/Vp_rtx)";
r47.Active = true;
r47.Name = "RTXdistr";
r47.Notes = "";
r47.Tag = "";
klaw29 = addkineticlaw(r47, 'Unknown');
klaw29.ParameterVariableNames = strings(0,1);
klaw29.SpeciesVariableNames = strings(0,1);
klaw29.Name = "";
klaw29.Notes = "";
klaw29.Tag = "";

r48 = addreaction(m1, "RTXc_ugperkg <-> null");
r48.ReactionRate = "(Cl_rtx+Vm_rtx/(Km_rtx+RTXc_ugperkg/Vc_rtx)/BW)/Vc_rtx*RTXc_ugperkg";
r48.Active = true;
r48.Name = "RTXclear";
r48.Notes = "";
r48.Tag = "";
klaw30 = addkineticlaw(r48, 'Unknown');
klaw30.ParameterVariableNames = strings(0,1);
klaw30.SpeciesVariableNames = strings(0,1);
klaw30.Name = "";
klaw30.Notes = "";
klaw30.Tag = "";

r49 = addreaction(m1, "drug_effect -> null");
r49.ReactionRate = "log(2)/tinjhalf*drug_effect";
r49.Active = true;
r49.Name = "UnnamedReaction_49";
r49.Notes = "";
r49.Tag = "";

r50 = addreaction(m1, "null -> B19no20tiss3");
r50.ReactionRate = "Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3*(kBapop+(kBapop/kBmat_kBapop_ratio))* (1+5*(max(0, 1 - (Bpb+Btiss+Btiss2)/(Bpbref_perml*Vpb+KBp*Bpbo_perml*Vtissue+KBp2*Bpbo_perml*Vtissue2)))^1)";
r50.Active = true;
r50.Name = "B19no20tiss3_gen";
r50.Notes = "";
r50.Tag = "";

r51 = addreaction(m1, "B19no20tiss3 -> B1920tiss3");
r51.ReactionRate = "(kBapop/kBmat_kBapop_ratio)*B19no20tiss3";
r51.Active = true;
r51.Name = "B19n020tiss3_mat";
r51.Notes = "";
r51.Tag = "";

r52 = addreaction(m1, "B1920tiss3 -> Bpb");
r52.ReactionRate = "tissue3on*fBtissue3_v1*kBtiss3exit*Vpb*max(0,B1920tiss3/Vtissue3-Bpb/Vpb*KBp3)";
r52.Active = true;
r52.Name = "Btraffic_3_v1";
r52.Notes = "";
r52.Tag = "";

r53 = addreaction(m1, "B19no20tiss3 -> null");
r53.ReactionRate = "kBapop*B19no20tiss3*kBapop_cll";
r53.Active = true;
r53.Name = "B19no20tiss3_apop";
r53.Notes = "";
r53.Tag = "";

r54 = addreaction(m1, "B1920tiss3 -> null");
r54.ReactionRate = "kBapop*B1920tiss3*kBapop_cll";
r54.Active = true;
r54.Name = "B1920tiss3_apop";
r54.Notes = "";
r54.Tag = "";

r55 = addreaction(m1, "restTpb + injection_effect + drug_effect <-> restTtiss3 + drug_effect + injection_effect");
r55.ReactionRate = "tissue3on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrp3-restTtiss3/Vtissue3)";
r55.Active = true;
r55.Name = "restTtraffic_3";
r55.Notes = "";
r55.Tag = "";

r56 = addreaction(m1, "actTtiss3 + RTXt3_ugperml + drugBtisskill3 + B1920tiss3 + BlinBtisskill3 -> actTtiss3 + RTXt3_ugperml + drugBtisskill3 + BlinBtisskill3");
r56.ReactionRate = "fBkill*kBkill*(drugBtisskill3+ BlinBtisskill3 + RTXt3_ugperml/(Kmkill_rtx/1000+RTXt3_ugperml))*B1920tiss3";
r56.Active = true;
r56.Name = "B1920tiss3_kill";
r56.Notes = "";
r56.Tag = "";

r57 = addreaction(m1, "act0Tpb + injection_effect + drug_effect <-> act0Ttiss3 + injection_effect + drug_effect");
r57.ReactionRate = "tissue3on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrp3-act0Ttiss3/Vtissue3)";
r57.Active = true;
r57.Name = "act0Ttraffic_3";
r57.Notes = "";
r57.Tag = "";

r58 = addreaction(m1, "actTpb + injection_effect + drug_effect <-> actTtiss3 + drug_effect + injection_effect");
r58.ReactionRate = "tissue3on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrp3*fTap-actTtiss3/Vtissue3)";
r58.Active = true;
r58.Name = "actTtraffic_3";
r58.Notes = "";
r58.Tag = "";

r59 = addreaction(m1, "act0Ttiss3 -> null");
r59.ReactionRate = "fTa0apop*kTaapop*act0Ttiss3";
r59.Active = true;
r59.Name = "act0Ttiss_apop_3";
r59.Notes = "";
r59.Tag = "";

r60 = addreaction(m1, "act0Ttiss3 + drugTtissact3 + BlinTtissact3 <-> actTtiss3 + drugTtissact3 + BlinTtissact3");
r60.ReactionRate = "act0on*kTact*((drugTtissact3+BlinTtissact3)*act0Ttiss3-fTadeact*actTtiss3)";
r60.Active = true;
r60.Name = "act0Ttiss_activn_3";
r60.Notes = "";
r60.Tag = "";

r61 = addreaction(m1, "act0Ttiss3 -> restTtiss3");
r61.ReactionRate = "act0on*fTa0deact*act0Ttiss3";
r61.Active = true;
r61.Name = "act0Ttiss_deact_3";
r61.Notes = "";
r61.Tag = "";

r62 = addreaction(m1, "restTtiss3 + drugTtissact3 + BlinTtissact3 <-> actTtiss3 + drugTtissact3 + BlinTtissact3");
r62.ReactionRate = "kTact*((drugTtissact3+BlinTtissact3)*restTtiss3-fTadeact*actTtiss3)";
r62.Active = true;
r62.Name = "Ttiss_activn_3";
r62.Notes = "";
r62.Tag = "";

r63 = addreaction(m1, "actTtiss3 -> act0Ttiss3 + actTtiss3");
r63.ReactionRate = "act0on*fTaprolif*kTprolif*actTtiss3";
r63.Active = true;
r63.Name = "actTtiss_prolif_3";
r63.Notes = "";
r63.Tag = "";

r64 = addreaction(m1, "drugTtissact3 + TDBt3_ugperml + restTtiss3 + B1920tiss3 -> drugTtissact3 + TDBt3_ugperml + restTtiss3 + B1920tiss3");
r64.ReactionRate = "0";
r64.Active = true;
r64.Name = "UnnamedReaction_64";
r64.Notes = "";
r64.Tag = "";

r65 = addreaction(m1, "B19no20tiss3 + actTtiss3 + BlinBtisskill3 -> actTtiss3 + BlinBtisskill3");
r65.ReactionRate = "fBkill*kBkill*BlinBtisskill3*B19no20tiss3";
r65.Active = true;
r65.Name = "Btiss3_19_kill";
r65.Notes = "";
r65.Tag = "";

r66 = addreaction(m1, "drugBtisskill3 + TDBt3_ugperml + actTtiss3 + B1920tiss3 -> drugBtisskill3 + TDBt3_ugperml + actTtiss3 + B1920tiss3");
r66.ReactionRate = "0";
r66.Active = true;
r66.Name = "UnnamedReaction_66";
r66.Notes = "";
r66.Tag = "";

r67 = addreaction(m1, "B19no20tiss3 + B1920tiss3 + B19tiss3 -> B1920tiss3 + B19no20tiss3 + B19tiss3");
r67.ReactionRate = "0";
r67.Active = true;
r67.Name = "UnnamedReaction_67";
r67.Notes = "";
r67.Tag = "";

r68 = addreaction(m1, "restTtiss -> null");
r68.ReactionRate = "fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue*KTrp)*max(0,restTtiss/(Trpbref_perml*Vtissue*KTrp)-1)";
r68.Active = true;
r68.Name = "restTtiss_apop_v24";
r68.Notes = "";
r68.Tag = "";

r69 = addreaction(m1, "restTtiss2 -> null");
r69.ReactionRate = "fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue2*KTrp2)*max(0,restTtiss2/(Trpbref_perml*Vtissue2*KTrp2)-1)";
r69.Active = true;
r69.Name = "restTtiss_apop_v24_2";
r69.Notes = "";
r69.Tag = "";

r70 = addreaction(m1, "actTtiss -> null");
r70.ReactionRate = "fapop_v24*kTaapop*(actTtiss+fAICD*actTtiss^2/(KTrp*Vtissue*Trpbref_perml))";
r70.Active = true;
r70.Name = "actTtiss_apop_v24";
r70.Notes = "";
r70.Tag = "";

r71 = addreaction(m1, "actTtiss2 -> null");
r71.ReactionRate = "fapop_v24*kTaapop*(actTtiss2+fAICD*actTtiss2^2/(Vtissue2*KTrp2*Trpbref_perml))";
r71.Active = true;
r71.Name = "actTtiss_apop_v24_2";
r71.Notes = "";
r71.Tag = "";

r72 = addreaction(m1, "restTtiss3 -> null");
r72.ReactionRate = "fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtissue3*KTrp3)*max(0,restTtiss3/(Trpbref_perml*Vtissue3*KTrp3)-1)";
r72.Active = true;
r72.Name = "restTtiss_apop_v24_3";
r72.Notes = "";
r72.Tag = "";

r73 = addreaction(m1, "actTtiss3 -> null");
r73.ReactionRate = "fapop_v24*kTaapop*(actTtiss3+fAICD*actTtiss3^2/(Vtissue3*KTrp3*Trpbref_perml))";
r73.Active = true;
r73.Name = "actTtiss_apop_v24_3";
r73.Notes = "";
r73.Tag = "";

r74 = addreaction(m1, "null -> Bpb");
r74.ReactionRate = "tissue3on*fBtissue3_v1*kBprolif*kBapop_cll*Bpbref_perml*Vpb*(max(0, 1 -Bpb/(Bpbref_perml*Vpb)))^1";
r74.Active = true;
r74.Name = "Bpb_prolif_test";
r74.Notes = "";
r74.Tag = "";

r75 = addreaction(m1, "Blinc_ug <-> null");
r75.ReactionRate = "Cl_blin*Blinc_ug/Vz_blin";
r75.Active = true;
r75.Name = "Blinclear";
r75.Notes = "";
r75.Tag = "";
klaw31 = addkineticlaw(r75, 'Unknown');
klaw31.ParameterVariableNames = strings(0,1);
klaw31.SpeciesVariableNames = strings(0,1);
klaw31.Name = "";
klaw31.Notes = "";
klaw31.Tag = "";

r76 = addreaction(m1, "Bpb + BlinBpbkill + actTpb + Blinc_ngperml -> Bpb + BlinBpbkill + actTpb + Blinc_ngperml");
r76.ReactionRate = "0";
r76.Active = true;
r76.Name = "UnnamedReaction_76";
r76.Notes = "";
r76.Tag = "";

r77 = addreaction(m1, "Blinc_ngperml + BlinTpbact + Bpb + restTpb -> Blinc_ngperml + BlinTpbact + Bpb + restTpb");
r77.ReactionRate = "0";
r77.Active = true;
r77.Name = "UnnamedReaction_77";
r77.Notes = "";
r77.Tag = "";

r78 = addreaction(m1, "Btiss + Blint_ngperml + actTtiss + BlinBtisskill -> BlinBtisskill + Blint_ngperml + actTtiss + Btiss");
r78.ReactionRate = "0";
r78.Active = true;
r78.Name = "UnnamedReaction_78";
r78.Notes = "";
r78.Tag = "";

r79 = addreaction(m1, "Btiss + restTtiss + Blint_ngperml + BlinTtissact -> Blint_ngperml + Btiss + BlinTtissact + restTtiss");
r79.ReactionRate = "0";
r79.Active = true;
r79.Name = "UnnamedReaction_79";
r79.Notes = "";
r79.Tag = "";

r80 = addreaction(m1, "Blint2_ngperml + BlinTtissact2 + Btiss2 + restTtiss2 -> Blint2_ngperml + BlinTtissact2 + Btiss2 + restTtiss2");
r80.ReactionRate = "0";
r80.Active = true;
r80.Name = "UnnamedReaction_80";
r80.Notes = "";
r80.Tag = "";

r81 = addreaction(m1, "Blint2_ngperml + BlinBtisskill2 + Btiss2 + actTtiss2 -> BlinBtisskill2 + Blint2_ngperml + Btiss2 + actTtiss2");
r81.ReactionRate = "0";
r81.Active = true;
r81.Name = "UnnamedReaction_81";
r81.Notes = "";
r81.Tag = "";

r82 = addreaction(m1, "BlinBtisskill3 + Blint3_ngperml + actTtiss3 + B19tiss3 -> BlinBtisskill3 + Blint3_ngperml + actTtiss3 + B19tiss3");
r82.ReactionRate = "0";
r82.Active = true;
r82.Name = "UnnamedReaction_82";
r82.Notes = "";
r82.Tag = "";

r83 = addreaction(m1, "BlinTtissact3 + Blint3_ngperml + restTtiss3 + B19tiss3 -> BlinTtissact3 + Blint3_ngperml + restTtiss3 + B19tiss3");
r83.ReactionRate = "0";
r83.Active = true;
r83.Name = "UnnamedReaction_83";
r83.Notes = "";
r83.Tag = "";

r84 = addreaction(m1, "null -> B1920tiss3");
r84.ReactionRate = "kBprolif*KBp3*Bpbref_perml*Vtissue3*(max(0, 1 - B1920tiss3/(KBp3*Bpbref_perml*Vtissue3)))^1";
r84.Active = true;
r84.Name = "B1920tiss3_prolf";
r84.Notes = "";
r84.Tag = "";

r85 = addreaction(m1, "null -> B19no20tiss3");
r85.ReactionRate = "kBprolif*KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3*(max(0, 1 - B19no20tiss3/(KBp3*Bpbref_perml*B19no20_B1920_ratio*Vtissue3)))^1";
r85.Active = true;
r85.Name = "B19no20tiss3_prolif";
r85.Notes = "";
r85.Tag = "";

r86 = addreaction(m1, "act0Ttumor -> restTtumor");
r86.ReactionRate = "act0on*fTa0deact*act0Ttumor";
r86.Active = true;
r86.Name = "act0Ttumor_deact";
r86.Notes = "";
r86.Tag = "";
klaw32 = addkineticlaw(r86, 'Unknown');
klaw32.ParameterVariableNames = strings(0,1);
klaw32.SpeciesVariableNames = strings(0,1);
klaw32.Name = "";
klaw32.Notes = "";
klaw32.Tag = "";

r87 = addreaction(m1, "act0Ttumor + drugTtumoract + BlinTtumoract <-> actTtumor + drugTtumoract + BlinTtumoract");
r87.ReactionRate = "act0on*kTact*((drugTtumoract+BlinTtumoract)*act0Ttumor-fTadeact*actTtumor)";
r87.Active = true;
r87.Name = "act0Ttumor_activn";
r87.Notes = "";
r87.Tag = "";
klaw33 = addkineticlaw(r87, 'Unknown');
klaw33.ParameterVariableNames = strings(0,1);
klaw33.SpeciesVariableNames = strings(0,1);
klaw33.Name = "";
klaw33.Notes = "";
klaw33.Tag = "";

r88 = addreaction(m1, "actTtumor -> act0Ttumor + actTtumor");
r88.ReactionRate = "act0on*fTaprolif*kTprolif*actTtumor";
r88.Active = true;
r88.Name = "actTtumor_prolif";
r88.Notes = "";
r88.Tag = "";
klaw34 = addkineticlaw(r88, 'Unknown');
klaw34.ParameterVariableNames = strings(0,1);
klaw34.SpeciesVariableNames = strings(0,1);
klaw34.Name = "";
klaw34.Notes = "";
klaw34.Tag = "";

r89 = addreaction(m1, "restTtumor + drugTtumoract + BlinTtumoract <-> actTtumor + drugTtumoract + BlinTtumoract");
r89.ReactionRate = "kTact*((drugTtumoract+BlinTtumoract)*restTtumor-fTadeact*actTtumor)";
r89.Active = true;
r89.Name = "Ttumor_activn";
r89.Notes = "";
r89.Tag = "";
klaw35 = addkineticlaw(r89, 'Unknown');
klaw35.ParameterVariableNames = strings(0,1);
klaw35.SpeciesVariableNames = strings(0,1);
klaw35.Name = "";
klaw35.Notes = "";
klaw35.Tag = "";

r90 = addreaction(m1, "BlinTtumoract + Blintumor_ngperml + restTtumor + Btumor -> BlinTtumoract + Blintumor_ngperml + Btumor + restTtumor");
r90.ReactionRate = "0";
r90.Active = true;
r90.Name = "UnnamedReaction_90";
r90.Notes = "";
r90.Tag = "";
klaw36 = addkineticlaw(r90, 'Unknown');
klaw36.ParameterVariableNames = strings(0,1);
klaw36.SpeciesVariableNames = strings(0,1);
klaw36.Name = "";
klaw36.Notes = "";
klaw36.Tag = "";

r91 = addreaction(m1, "Blintumor_ngperml + Btumor + BlinBtumorkill + actTtumor -> Blintumor_ngperml + BlinBtumorkill + actTtumor + Btumor");
r91.ReactionRate = "0";
r91.Active = true;
r91.Name = "UnnamedReaction_91";
r91.Notes = "";
r91.Tag = "";
klaw37 = addkineticlaw(r91, 'Unknown');
klaw37.ParameterVariableNames = strings(0,1);
klaw37.SpeciesVariableNames = strings(0,1);
klaw37.Name = "";
klaw37.Notes = "";
klaw37.Tag = "";

r92 = addreaction(m1, "drugTtumoract + Btumor + restTtumor + TDBtumor_ugperml -> drugTtumoract + Btumor + TDBtumor_ugperml + restTtumor");
r92.ReactionRate = "0";
r92.Active = true;
r92.Name = "UnnamedReaction_92";
r92.Notes = "";
r92.Tag = "";
klaw38 = addkineticlaw(r92, 'Unknown');
klaw38.ParameterVariableNames = strings(0,1);
klaw38.SpeciesVariableNames = strings(0,1);
klaw38.Name = "";
klaw38.Notes = "";
klaw38.Tag = "";

r93 = addreaction(m1, "TDBtumor_ugperml + Btumor + actTtumor + drugBtumorkill -> TDBtumor_ugperml + Btumor + drugBtumorkill + actTtumor");
r93.ReactionRate = "0";
r93.Active = true;
r93.Name = "UnnamedReaction_93";
r93.Notes = "";
r93.Tag = "";
klaw39 = addkineticlaw(r93, 'Unknown');
klaw39.ParameterVariableNames = strings(0,1);
klaw39.SpeciesVariableNames = strings(0,1);
klaw39.Name = "";
klaw39.Notes = "";
klaw39.Tag = "";

r94 = addreaction(m1, "Btumor -> 2 Btumor");
r94.ReactionRate = "kBtumorprolif*Btumor + 0*kBprolif*KBptumor*Bpbref_perml*Vtumor*(max(0, 1 - Btumor/(KBptumor*Bpbref_perml*Vtumor)))^1";
r94.Active = true;
r94.Name = "Btumor_prolf";
r94.Notes = "";
r94.Tag = "";
klaw40 = addkineticlaw(r94, 'Unknown');
klaw40.ParameterVariableNames = strings(0,1);
klaw40.SpeciesVariableNames = strings(0,1);
klaw40.Name = "";
klaw40.Notes = "";
klaw40.Tag = "";

r95 = addreaction(m1, "Btumor -> null");
r95.ReactionRate = "0*fBprolif*kBprolif*Btumor";
r95.Active = true;
r95.Name = "Btumor_apop";
r95.Notes = "";
r95.Tag = "";
klaw41 = addkineticlaw(r95, 'Unknown');
klaw41.ParameterVariableNames = strings(0,1);
klaw41.SpeciesVariableNames = strings(0,1);
klaw41.Name = "";
klaw41.Notes = "";
klaw41.Tag = "";

r96 = addreaction(m1, "Btumor + actTtumor + BlinBtumorkill + drugBtumorkill + RTXtumor_ugperml -> actTtumor + BlinBtumorkill + drugBtumorkill + RTXtumor_ugperml");
r96.ReactionRate = "fBkill*kBkill*(drugBtumorkill+ BlinBtumorkill + RTXtumor_ugperml/(Kmkill_rtx/1000+RTXtumor_ugperml))*Btumor";
r96.Active = true;
r96.Name = "Btumor_kill";
r96.Notes = "";
r96.Tag = "";
klaw42 = addkineticlaw(r96, 'Unknown');
klaw42.ParameterVariableNames = strings(0,1);
klaw42.SpeciesVariableNames = strings(0,1);
klaw42.Name = "";
klaw42.Notes = "";
klaw42.Tag = "";

r97 = addreaction(m1, "actTtumor -> null");
r97.ReactionRate = "fapop_v24*kTaapop*(actTtumor+fAICD*actTtumor^2/(Vtumor*KTrptumor*Trpbref_perml))";
r97.Active = true;
r97.Name = "actTtumor_apop_v24";
r97.Notes = "";
r97.Tag = "";
klaw43 = addkineticlaw(r97, 'Unknown');
klaw43.ParameterVariableNames = strings(0,1);
klaw43.SpeciesVariableNames = strings(0,1);
klaw43.Name = "";
klaw43.Notes = "";
klaw43.Tag = "";

r98 = addreaction(m1, "Bpb -> Btumor");
r98.ReactionRate = "Bcell_tumor_trafficking_on*tumor_on*fBexit*kTaexit*Vpb*max(0, Bpb/Vpb*KBptumor-Btumor/Vtumor)";
r98.Active = true;
r98.Name = "Btraffic_tumor";
r98.Notes = "";
r98.Tag = "";
klaw44 = addkineticlaw(r98, 'Unknown');
klaw44.ParameterVariableNames = strings(0,1);
klaw44.SpeciesVariableNames = strings(0,1);
klaw44.Name = "";
klaw44.Notes = "";
klaw44.Tag = "";

r99 = addreaction(m1, "restTpb + drug_effect + injection_effect <-> restTtumor + drug_effect + injection_effect");
r99.ReactionRate = "tumor_on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*restTpb/Vpb*KTrptumor-restTtumor/Vtumor)";
r99.Active = true;
r99.Name = "restTtraffic_tumor";
r99.Notes = "";
r99.Tag = "";
klaw45 = addkineticlaw(r99, 'Unknown');
klaw45.ParameterVariableNames = strings(0,1);
klaw45.SpeciesVariableNames = strings(0,1);
klaw45.Name = "";
klaw45.Notes = "";
klaw45.Tag = "";

r100 = addreaction(m1, "actTpb + drug_effect + injection_effect <-> actTtumor + injection_effect + drug_effect");
r100.ReactionRate = "tumor_on*kTaexit*Vpb*((1+fa0*(finj*injection_effect+fdrug*drug_effect))*actTpb/Vpb*KTrptumor*fTap-actTtumor/Vtumor)";
r100.Active = true;
r100.Name = "actTtraffic_tumor";
r100.Notes = "";
r100.Tag = "";
klaw46 = addkineticlaw(r100, 'Unknown');
klaw46.ParameterVariableNames = strings(0,1);
klaw46.SpeciesVariableNames = strings(0,1);
klaw46.Name = "";
klaw46.Notes = "";
klaw46.Tag = "";

r101 = addreaction(m1, "act0Tpb + drug_effect + injection_effect <-> act0Ttumor + injection_effect + drug_effect");
r101.ReactionRate = "tumor_on*act0on*kTrexit*Vpb*((1+finj*injection_effect+fdrug*drug_effect)*act0Tpb/Vpb*KTrptumor-act0Ttumor/Vtumor)";
r101.Active = true;
r101.Name = "act0Ttraffic_tumor";
r101.Notes = "";
r101.Tag = "";
klaw47 = addkineticlaw(r101, 'Unknown');
klaw47.ParameterVariableNames = strings(0,1);
klaw47.SpeciesVariableNames = strings(0,1);
klaw47.Name = "";
klaw47.Notes = "";
klaw47.Tag = "";

r102 = addreaction(m1, "act0Ttumor -> null");
r102.ReactionRate = "fTa0apop*kTaapop*act0Ttumor";
r102.Active = true;
r102.Name = "act0Ttumor_apop";
r102.Notes = "";
r102.Tag = "";
klaw48 = addkineticlaw(r102, 'Unknown');
klaw48.ParameterVariableNames = strings(0,1);
klaw48.SpeciesVariableNames = strings(0,1);
klaw48.Name = "";
klaw48.Notes = "";
klaw48.Tag = "";

r103 = addreaction(m1, "restTtumor -> null");
r103.ReactionRate = "fapop_v24*fTrapop*kTaapop*(Trpbref_perml*Vtumor*KTrptumor)*max(0,restTtumor/(Trpbref_perml*Vtumor*KTrptumor)-1)";
r103.Active = true;
r103.Name = "restTtumor_apop_v24";
r103.Notes = "";
r103.Tag = "";
klaw49 = addkineticlaw(r103, 'Unknown');
klaw49.ParameterVariableNames = strings(0,1);
klaw49.SpeciesVariableNames = strings(0,1);
klaw49.Name = "";
klaw49.Notes = "";
klaw49.Tag = "";

r104 = addreaction(m1, "TDBsc_ugperkg -> TDBc_ugperkg");
r104.ReactionRate = "kabs_TDB*fbio_TDB*TDBsc_ugperkg";
r104.Active = true;
r104.Name = "TDB_abs";
r104.Notes = "";
r104.Tag = "";
klaw50 = addkineticlaw(r104, 'Unknown');
klaw50.ParameterVariableNames = strings(0,1);
klaw50.SpeciesVariableNames = strings(0,1);
klaw50.Name = "";
klaw50.Notes = "";
klaw50.Tag = "";

r105 = addreaction(m1, "TDBsc_ugperkg -> null");
r105.ReactionRate = "kabs_TDB*(1-fbio_TDB)*TDBsc_ugperkg";
r105.Active = true;
r105.Name = "TDBSC_loss";
r105.Notes = "";
r105.Tag = "";
klaw51 = addkineticlaw(r105, 'Unknown');
klaw51.ParameterVariableNames = strings(0,1);
klaw51.SpeciesVariableNames = strings(0,1);
klaw51.Name = "";
klaw51.Notes = "";
klaw51.Tag = "";

r106 = addreaction(m1, "null -> TDBc_ugperml_AUC");
r106.ReactionRate = "TDBc_ugperml";
r106.Active = true;
r106.Name = "UnnamedReaction_106";
r106.Notes = "";
r106.Tag = "";
klaw52 = addkineticlaw(r106, 'Unknown');
klaw52.ParameterVariableNames = strings(0,1);
klaw52.SpeciesVariableNames = strings(0,1);
klaw52.Name = "";
klaw52.Notes = "";
klaw52.Tag = "";

r107 = addreaction(m1, "null -> IL6pb");
r107.ReactionRate = "kIL6prod*actTpb/Vpb*(Bpb_perml/Bpbref_perml)*((TDBc_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBc_ugperml*1000)^ndrugactT)+Blinc_ngperml^ndrugactT_blin/(Blinc_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))";
r107.Active = true;
r107.Name = "UnnamedReaction_107";
r107.Notes = "";
r107.Tag = "";
klaw53 = addkineticlaw(r107, 'Unknown');
klaw53.ParameterVariableNames = strings(0,1);
klaw53.SpeciesVariableNames = strings(0,1);
klaw53.Name = "";
klaw53.Notes = "";
klaw53.Tag = "";

r108 = addreaction(m1, "IL6pb -> null");
r108.ReactionRate = "log(2)/(thalfIL6/60/24)*IL6pb";
r108.Active = true;
r108.Name = "UnnamedReaction_108";
r108.Notes = "";
r108.Tag = "";
klaw54 = addkineticlaw(r108, 'Unknown');
klaw54.ParameterVariableNames = strings(0,1);
klaw54.SpeciesVariableNames = strings(0,1);
klaw54.Name = "";
klaw54.Notes = "";
klaw54.Tag = "";

r109 = addreaction(m1, "null -> IL6tiss");
r109.ReactionRate = "kIL6prod*actTtiss/Vtissue*(Btiss_perml/(Bpbref_perml*KBp))*((TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT)+Blint_ngperml^ndrugactT_blin/(Blint_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))";
r109.Active = true;
r109.Name = "UnnamedReaction_109";
r109.Notes = "";
r109.Tag = "";
klaw55 = addkineticlaw(r109, 'Unknown');
klaw55.ParameterVariableNames = strings(0,1);
klaw55.SpeciesVariableNames = strings(0,1);
klaw55.Name = "";
klaw55.Notes = "";
klaw55.Tag = "";

r110 = addreaction(m1, "IL6tiss -> null");
r110.ReactionRate = "log(2)/(thalfIL6/60/24)*IL6tiss";
r110.Active = true;
r110.Name = "UnnamedReaction_110";
r110.Notes = "";
r110.Tag = "";
klaw56 = addkineticlaw(r110, 'Unknown');
klaw56.ParameterVariableNames = strings(0,1);
klaw56.SpeciesVariableNames = strings(0,1);
klaw56.Name = "";
klaw56.Notes = "";
klaw56.Tag = "";

r111 = addreaction(m1, "null -> IL6tiss2");
r111.ReactionRate = "kIL6prod*actTtiss2/Vtissue2*(Btiss2_perml/(Bpbref_perml*KBp2))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT)+Blint2_ngperml^ndrugactT_blin/(Blint2_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))";
r111.Active = true;
r111.Name = "UnnamedReaction_111";
r111.Notes = "";
r111.Tag = "";
klaw57 = addkineticlaw(r111, 'Unknown');
klaw57.ParameterVariableNames = strings(0,1);
klaw57.SpeciesVariableNames = strings(0,1);
klaw57.Name = "";
klaw57.Notes = "";
klaw57.Tag = "";

r112 = addreaction(m1, "IL6tiss2 -> null");
r112.ReactionRate = "log(2)/(thalfIL6/60/24)*IL6tiss2";
r112.Active = true;
r112.Name = "UnnamedReaction_112";
r112.Notes = "";
r112.Tag = "";
klaw58 = addkineticlaw(r112, 'Unknown');
klaw58.ParameterVariableNames = strings(0,1);
klaw58.SpeciesVariableNames = strings(0,1);
klaw58.Name = "";
klaw58.Notes = "";
klaw58.Tag = "";

r113 = addreaction(m1, "null -> IL6tiss3");
r113.ReactionRate = "kIL6prod*actTtiss3/Vtissue3*((B1920tiss3_perml/(Bpbref_perml*KBp3))*((TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT))+(B19tiss3_perml/(Bpbref_perml*KBp3*(1+B19no20_B1920_ratio)))*(Blint3_ngperml^ndrugactT_blin/(Blint3_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin)))";
r113.Active = true;
r113.Name = "UnnamedReaction_113";
r113.Notes = "";
r113.Tag = "";
klaw59 = addkineticlaw(r113, 'Unknown');
klaw59.ParameterVariableNames = strings(0,1);
klaw59.SpeciesVariableNames = strings(0,1);
klaw59.Name = "";
klaw59.Notes = "";
klaw59.Tag = "";

r114 = addreaction(m1, "IL6tiss3 -> null");
r114.ReactionRate = "log(2)/(thalfIL6/60/24)*IL6tiss3";
r114.Active = true;
r114.Name = "UnnamedReaction_114";
r114.Notes = "";
r114.Tag = "";
klaw60 = addkineticlaw(r114, 'Unknown');
klaw60.ParameterVariableNames = strings(0,1);
klaw60.SpeciesVariableNames = strings(0,1);
klaw60.Name = "";
klaw60.Notes = "";
klaw60.Tag = "";

r115 = addreaction(m1, "null -> IL6tumor");
r115.ReactionRate = "kIL6prod*actTtumor/Vtumor*(Btumor_perml/(Bpbref_perml*KBptumor))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT)+Blintumor_ngperml^ndrugactT_blin/(Blintumor_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))";
r115.Active = true;
r115.Name = "UnnamedReaction_115";
r115.Notes = "";
r115.Tag = "";
klaw61 = addkineticlaw(r115, 'Unknown');
klaw61.ParameterVariableNames = strings(0,1);
klaw61.SpeciesVariableNames = strings(0,1);
klaw61.Name = "";
klaw61.Notes = "";
klaw61.Tag = "";

r116 = addreaction(m1, "IL6tumor -> null");
r116.ReactionRate = "log(2)/(thalfIL6/60/24)*IL6tumor";
r116.Active = true;
r116.Name = "UnnamedReaction_116";
r116.Notes = "";
r116.Tag = "";
klaw62 = addkineticlaw(r116, 'Unknown');
klaw62.ParameterVariableNames = strings(0,1);
klaw62.SpeciesVariableNames = strings(0,1);
klaw62.Name = "";
klaw62.Notes = "";
klaw62.Tag = "";

%% Create parameters.
p1 = addparameter(m1, "VmB");
p1.Constant = true;
p1.Value = 0.95;
p1.Units = "";
p1.BoundaryCondition = false;
p1.Notes = "";
p1.Tag = "";

p2 = addparameter(m1, "KmTB_kill");
p2.Constant = true;
p2.Value = 0.75;
p2.Units = "";
p2.BoundaryCondition = false;
p2.Notes = "";
p2.Tag = "";

p3 = addparameter(m1, "KdrugB");
p3.Constant = true;
p3.Value = 10;
p3.Units = "nanogram/milliliter";
p3.BoundaryCondition = false;
p3.Notes = "";
p3.Tag = "";

p4 = addparameter(m1, "kBapop");
p4.Constant = true;
p4.Value = 0.03;
p4.Units = "1/day";
p4.BoundaryCondition = false;
p4.Notes = "";
p4.Tag = "";

p5 = addparameter(m1, "kBprolif");
p5.Constant = true;
p5.Value = 0.7;
p5.Units = "1/day";
p5.BoundaryCondition = false;
p5.Notes = "";
p5.Tag = "";

p6 = addparameter(m1, "kTprolif");
p6.Constant = true;
p6.Value = 0.7;
p6.Units = "1/day";
p6.BoundaryCondition = false;
p6.Notes = "";
p6.Tag = "";

p7 = addparameter(m1, "kBkill");
p7.Constant = true;
p7.Value = 1;
p7.Units = "";
p7.BoundaryCondition = false;
p7.Notes = "";
p7.Tag = "";

p8 = addparameter(m1, "KdrugactT");
p8.Constant = true;
p8.Value = 50;
p8.Units = "nanogram/milliliter";
p8.BoundaryCondition = false;
p8.Notes = "";
p8.Tag = "";

p9 = addparameter(m1, "VmT");
p9.Constant = true;
p9.Value = 0.9;
p9.Units = "";
p9.BoundaryCondition = false;
p9.Notes = "";
p9.Tag = "";

p10 = addparameter(m1, "KmBT_act");
p10.Constant = true;
p10.Value = 0.1;
p10.Units = "";
p10.BoundaryCondition = false;
p10.Notes = "";
p10.Tag = "";

p11 = addparameter(m1, "kTaexit");
p11.Constant = true;
p11.Value = 1;
p11.Units = "";
p11.BoundaryCondition = false;
p11.Notes = "";
p11.Tag = "";

p12 = addparameter(m1, "kTact");
p12.Constant = true;
p12.Value = 1;
p12.Units = "";
p12.BoundaryCondition = false;
p12.Notes = "";
p12.Tag = "";

p13 = addparameter(m1, "fTadeact");
p13.Constant = true;
p13.Value = 0;
p13.Units = "";
p13.BoundaryCondition = false;
p13.Notes = "";
p13.Tag = "";

p14 = addparameter(m1, "fTap");
p14.Constant = true;
p14.Value = 1;
p14.Units = "";
p14.BoundaryCondition = false;
p14.Notes = "";
p14.Tag = "";

p15 = addparameter(m1, "Cl_tdb");
p15.Constant = true;
p15.Value = 3;
p15.Units = "milliliter/day/kilogram";
p15.BoundaryCondition = false;
p15.Notes = "";
p15.Tag = "";

p16 = addparameter(m1, "Cld_tdb");
p16.Constant = true;
p16.Value = 25;
p16.Units = "milliliter/day/kilogram";
p16.BoundaryCondition = false;
p16.Notes = "";
p16.Tag = "";

p17 = addparameter(m1, "Vc_tdb");
p17.Constant = true;
p17.Value = 35;
p17.Units = "milliliter/kilogram";
p17.BoundaryCondition = false;
p17.Notes = "";
p17.Tag = "";

p18 = addparameter(m1, "Vp_tdb");
p18.Constant = true;
p18.Value = 110;
p18.Units = "milliliter/kilogram";
p18.BoundaryCondition = false;
p18.Notes = "";
p18.Tag = "";

p19 = addparameter(m1, "kTgen");
p19.Constant = true;
p19.Value = 1;
p19.Units = "1/day";
p19.BoundaryCondition = false;
p19.Notes = "";
p19.Tag = "";

p20 = addparameter(m1, "kTaapop");
p20.Constant = true;
p20.Value = 0.5;
p20.Units = "";
p20.BoundaryCondition = false;
p20.Notes = "";
p20.Tag = "";

p21 = addparameter(m1, "KTrp");
p21.Constant = true;
p21.Value = 1;
p21.Units = "";
p21.BoundaryCondition = false;
p21.Notes = "";
p21.Tag = "";

p22 = addparameter(m1, "Trpbo_perml");
p22.Constant = true;
p22.Value = 2000000;
p22.Units = "1/milliliter";
p22.BoundaryCondition = false;
p22.Notes = "";
p22.Tag = "";

p23 = addparameter(m1, "fTaprolif");
p23.Constant = true;
p23.Value = 2;
p23.Units = "";
p23.BoundaryCondition = false;
p23.Notes = "";
p23.Tag = "";

p24 = addparameter(m1, "fBprolif");
p24.Constant = true;
p24.Value = 0;
p24.Units = "";
p24.BoundaryCondition = false;
p24.Notes = "";
p24.Tag = "";

p25 = addparameter(m1, "Bpbo_perml");
p25.Constant = true;
p25.Value = 1000000;
p25.Units = "1/milliliter";
p25.BoundaryCondition = false;
p25.Notes = "";
p25.Tag = "";

p26 = addparameter(m1, "fBexit");
p26.Constant = true;
p26.Value = 1;
p26.Units = "";
p26.BoundaryCondition = false;
p26.Notes = "";
p26.Tag = "";

p27 = addparameter(m1, "KBp");
p27.Constant = true;
p27.Value = 1;
p27.Units = "";
p27.BoundaryCondition = false;
p27.Notes = "";
p27.Tag = "";

p28 = addparameter(m1, "Kp");
p28.Constant = true;
p28.Value = 0.14;
p28.Units = "";
p28.BoundaryCondition = false;
p28.Notes = "";
p28.Tag = "";

p29 = addparameter(m1, "Vpb");
p29.Constant = true;
p29.Value = 40;
p29.Units = "milliliter";
p29.BoundaryCondition = false;
p29.Notes = "";
p29.Tag = "";

p30 = addparameter(m1, "Vtissue");
p30.Constant = true;
p30.Value = 40;
p30.Units = "milliliter";
p30.BoundaryCondition = false;
p30.Notes = "";
p30.Tag = "";

p31 = addparameter(m1, "fTrapop");
p31.Constant = true;
p31.Value = 0.1;
p31.Units = "1/day";
p31.BoundaryCondition = false;
p31.Notes = "";
p31.Tag = "";

p32 = addparameter(m1, "Trpbref_perml");
p32.Constant = true;
p32.Value = 2000000;
p32.Units = "1/milliliter";
p32.BoundaryCondition = false;
p32.Notes = "";
p32.Tag = "";

p33 = addparameter(m1, "nkill");
p33.Constant = true;
p33.Value = 2;
p33.Units = "";
p33.BoundaryCondition = false;
p33.Notes = "";
p33.Tag = "";

p34 = addparameter(m1, "KTrp2");
p34.Constant = true;
p34.Value = 0.02;
p34.Units = "";
p34.BoundaryCondition = false;
p34.Notes = "";
p34.Tag = "";

p35 = addparameter(m1, "Vm_tdb");
p35.Constant = true;
p35.Value = 212;
p35.Units = "microgram/day";
p35.BoundaryCondition = false;
p35.Notes = "";
p35.Tag = "";

p36 = addparameter(m1, "Km_tdb");
p36.Constant = true;
p36.Value = 4.74;
p36.Units = "microgram/milliliter";
p36.BoundaryCondition = false;
p36.Notes = "";
p36.Tag = "";

p37 = addparameter(m1, "BW");
p37.Constant = true;
p37.Value = 3.4;
p37.Units = "kilogram";
p37.BoundaryCondition = false;
p37.Notes = "";
p37.Tag = "";

p38 = addparameter(m1, "act0on");
p38.Constant = true;
p38.Value = 1;
p38.Units = "";
p38.BoundaryCondition = false;
p38.Notes = "";
p38.Tag = "";

p39 = addparameter(m1, "kIL6prod");
p39.Constant = true;
p39.Value = 1;
p39.Units = "";
p39.BoundaryCondition = false;
p39.Notes = "";
p39.Tag = "";

p40 = addparameter(m1, "fa0");
p40.Constant = true;
p40.Value = 0.5;
p40.Units = "";
p40.BoundaryCondition = false;
p40.Notes = "";
p40.Tag = "";

p41 = addparameter(m1, "fAICD");
p41.Constant = true;
p41.Value = 1;
p41.Units = "";
p41.BoundaryCondition = false;
p41.Notes = "";
p41.Tag = "";

p42 = addparameter(m1, "kBAFFprod");
p42.Constant = true;
p42.Value = 1;
p42.Units = "";
p42.BoundaryCondition = false;
p42.Notes = "";
p42.Tag = "";

p43 = addparameter(m1, "Bpbref_perml");
p43.Constant = true;
p43.Value = 1000000;
p43.Units = "";
p43.BoundaryCondition = false;
p43.Notes = "";
p43.Tag = "";

p44 = addparameter(m1, "thBAFF");
p44.Constant = true;
p44.Value = 30;
p44.Units = "minute";
p44.BoundaryCondition = false;
p44.Notes = "";
p44.Tag = "";

p45 = addparameter(m1, "thalfIL6");
p45.Constant = true;
p45.Value = 20;
p45.Units = "minute";
p45.BoundaryCondition = false;
p45.Notes = "";
p45.Tag = "";

p46 = addparameter(m1, "fTa0deact");
p46.Constant = true;
p46.Value = 1;
p46.Units = "";
p46.BoundaryCondition = false;
p46.Notes = "";
p46.Tag = "";

p47 = addparameter(m1, "fTa0apop");
p47.Constant = true;
p47.Value = 1;
p47.Units = "";
p47.BoundaryCondition = false;
p47.Notes = "";
p47.Tag = "";

p48 = addparameter(m1, "BAFFo");
p48.Constant = true;
p48.Value = 1;
p48.Units = "";
p48.BoundaryCondition = false;
p48.Notes = "";
p48.Tag = "";

p49 = addparameter(m1, "fBAFFo");
p49.Constant = true;
p49.Value = 0.1;
p49.Units = "";
p49.BoundaryCondition = false;
p49.Notes = "";
p49.Tag = "";

p50 = addparameter(m1, "tinjhalf");
p50.Constant = true;
p50.Value = 1;
p50.Units = "day";
p50.BoundaryCondition = false;
p50.Notes = "";
p50.Tag = "";

p51 = addparameter(m1, "finj");
p51.Constant = true;
p51.Value = 0.1;
p51.Units = "";
p51.BoundaryCondition = false;
p51.Notes = "";
p51.Tag = "";

p52 = addparameter(m1, "KBp2");
p52.Constant = true;
p52.Value = 1;
p52.Units = "";
p52.BoundaryCondition = false;
p52.Notes = "";
p52.Tag = "";

p53 = addparameter(m1, "Vtissue2");
p53.Constant = true;
p53.Value = 1;
p53.Units = "";
p53.BoundaryCondition = false;
p53.Notes = "";
p53.Tag = "";

p54 = addparameter(m1, "Kp2");
p54.Constant = true;
p54.Value = 0.07000000000000001;
p54.Units = "";
p54.BoundaryCondition = false;
p54.Notes = "";
p54.Tag = "";

p55 = addparameter(m1, "tissue2on");
p55.Constant = true;
p55.Value = 0;
p55.Units = "";
p55.BoundaryCondition = false;
p55.Notes = "";
p55.Tag = "";

p56 = addparameter(m1, "Vc_rtx");
p56.Constant = true;
p56.Value = 28;
p56.Units = "";
p56.BoundaryCondition = false;
p56.Notes = "";
p56.Tag = "";

p57 = addparameter(m1, "Vp_rtx");
p57.Constant = true;
p57.Value = 9;
p57.Units = "";
p57.BoundaryCondition = false;
p57.Notes = "";
p57.Tag = "";

p58 = addparameter(m1, "Cl_rtx");
p58.Constant = true;
p58.Value = 8;
p58.Units = "";
p58.BoundaryCondition = false;
p58.Notes = "";
p58.Tag = "";

p59 = addparameter(m1, "Cld_rtx");
p59.Constant = true;
p59.Value = 23;
p59.Units = "";
p59.BoundaryCondition = false;
p59.Notes = "";
p59.Tag = "";

p60 = addparameter(m1, "Vm_rtx");
p60.Constant = true;
p60.Value = 0;
p60.Units = "";
p60.BoundaryCondition = false;
p60.Notes = "";
p60.Tag = "";

p61 = addparameter(m1, "Km_rtx");
p61.Constant = true;
p61.Value = 1;
p61.Units = "";
p61.BoundaryCondition = false;
p61.Notes = "";
p61.Tag = "";

p62 = addparameter(m1, "Kmkill_rtx");
p62.Constant = true;
p62.Value = 25;
p62.Units = "";
p62.BoundaryCondition = false;
p62.Notes = "";
p62.Tag = "";

p63 = addparameter(m1, "fTgenbl");
p63.Constant = true;
p63.Value = 0;
p63.Units = "";
p63.BoundaryCondition = false;
p63.Notes = "";
p63.Tag = "";

p64 = addparameter(m1, "depleteTpb");
p64.Constant = true;
p64.Value = 0;
p64.Units = "";
p64.BoundaryCondition = false;
p64.Notes = "";
p64.Tag = "";

p65 = addparameter(m1, "depleteBpb");
p65.Constant = true;
p65.Value = 0;
p65.Units = "";
p65.BoundaryCondition = false;
p65.Notes = "";
p65.Tag = "";

p66 = addparameter(m1, "kTrexit");
p66.Constant = true;
p66.Value = 1;
p66.Units = "";
p66.BoundaryCondition = false;
p66.Notes = "";
p66.Tag = "";

p67 = addparameter(m1, "PKflag");
p67.Constant = true;
p67.Value = 1;
p67.Units = "";
p67.BoundaryCondition = false;
p67.Notes = "";
p67.Tag = "";

p68 = addparameter(m1, "VPid");
p68.Constant = true;
p68.Value = 1;
p68.Units = "";
p68.BoundaryCondition = false;
p68.Notes = "";
p68.Tag = "";

p69 = addparameter(m1, "end_time");
p69.Constant = true;
p69.Value = 1;
p69.Units = "";
p69.BoundaryCondition = false;
p69.Notes = "";
p69.Tag = "";

p70 = addparameter(m1, "fdrug");
p70.Constant = true;
p70.Value = 1;
p70.Units = "";
p70.BoundaryCondition = false;
p70.Notes = "";
p70.Tag = "";

p71 = addparameter(m1, "Vtissue3");
p71.Constant = true;
p71.Value = 1;
p71.Units = "";
p71.BoundaryCondition = false;
p71.Notes = "";
p71.Tag = "";

p72 = addparameter(m1, "KTrp3");
p72.Constant = true;
p72.Value = 500;
p72.Units = "";
p72.BoundaryCondition = false;
p72.Notes = "";
p72.Tag = "";

p73 = addparameter(m1, "tissue3on");
p73.Constant = true;
p73.Value = 0;
p73.Units = "";
p73.BoundaryCondition = false;
p73.Notes = "";
p73.Tag = "";

p74 = addparameter(m1, "Kp3");
p74.Constant = true;
p74.Value = 0.07000000000000001;
p74.Units = "";
p74.BoundaryCondition = false;
p74.Notes = "";
p74.Tag = "";

p75 = addparameter(m1, "kBtiss3exit");
p75.Constant = true;
p75.Value = 0.03;
p75.Units = "";
p75.BoundaryCondition = false;
p75.Notes = "";
p75.Tag = "";

p76 = addparameter(m1, "KBp3");
p76.Constant = true;
p76.Value = 600;
p76.Units = "";
p76.BoundaryCondition = false;
p76.Notes = "";
p76.Tag = "";

p77 = addparameter(m1, "B19no20_B1920_ratio");
p77.Constant = true;
p77.Value = 0.25;
p77.Units = "";
p77.BoundaryCondition = false;
p77.Notes = "";
p77.Tag = "";

p78 = addparameter(m1, "fvalidation");
p78.Constant = true;
p78.Value = 0;
p78.Units = "";
p78.BoundaryCondition = false;
p78.Notes = "";
p78.Tag = "";

p79 = addparameter(m1, "fapop_v24");
p79.Constant = true;
p79.Value = 0;
p79.Units = "";
p79.BoundaryCondition = false;
p79.Notes = "";
p79.Tag = "";

p80 = addparameter(m1, "fBtissue3_v1");
p80.Constant = true;
p80.Value = 1;
p80.Units = "";
p80.BoundaryCondition = false;
p80.Notes = "";
p80.Tag = "";

p81 = addparameter(m1, "kBmat_kBapop_ratio");
p81.Constant = true;
p81.Value = 0.25;
p81.Units = "";
p81.BoundaryCondition = false;
p81.Notes = "";
p81.Tag = "";

p82 = addparameter(m1, "ndrugactT");
p82.Constant = true;
p82.Value = 1;
p82.Units = "";
p82.BoundaryCondition = false;
p82.Notes = "";
p82.Tag = "";

p83 = addparameter(m1, "Cl_blin");
p83.Constant = true;
p83.Value = 22300;
p83.Units = "";
p83.BoundaryCondition = false;
p83.Notes = "";
p83.Tag = "";

p84 = addparameter(m1, "Vz_blin");
p84.Constant = true;
p84.Value = 1610;
p84.Units = "";
p84.BoundaryCondition = false;
p84.Notes = "";
p84.Tag = "";

p85 = addparameter(m1, "KdrugactT_blin");
p85.Constant = true;
p85.Value = 0.1;
p85.Units = "nanogram/milliliter";
p85.BoundaryCondition = false;
p85.Notes = "";
p85.Tag = "";

p86 = addparameter(m1, "ndrugactT_blin");
p86.Constant = true;
p86.Value = 1;
p86.Units = "";
p86.BoundaryCondition = false;
p86.Notes = "";
p86.Tag = "";

p87 = addparameter(m1, "KmTB_kill_blin");
p87.Constant = true;
p87.Value = 0.75;
p87.Units = "";
p87.BoundaryCondition = false;
p87.Notes = "";
p87.Tag = "";

p88 = addparameter(m1, "KdrugB_blin");
p88.Constant = true;
p88.Value = 0.015;
p88.Units = "nanogram/milliliter";
p88.BoundaryCondition = false;
p88.Notes = "";
p88.Tag = "";

p89 = addparameter(m1, "nkill_blin");
p89.Constant = true;
p89.Value = 1.027;
p89.Units = "";
p89.BoundaryCondition = false;
p89.Notes = "";
p89.Tag = "";

p90 = addparameter(m1, "KmBT_act_blin");
p90.Constant = true;
p90.Value = 0.716;
p90.Units = "";
p90.BoundaryCondition = false;
p90.Notes = "";
p90.Tag = "";

p91 = addparameter(m1, "kBapop_cll");
p91.Constant = true;
p91.Value = 1;
p91.Units = "";
p91.BoundaryCondition = false;
p91.Notes = "";
p91.Tag = "";

p92 = addparameter(m1, "kBgen_cll");
p92.Constant = true;
p92.Value = 1;
p92.Units = "";
p92.BoundaryCondition = false;
p92.Notes = "";
p92.Tag = "";

p93 = addparameter(m1, "fBkill");
p93.Constant = true;
p93.Value = 1;
p93.Units = "";
p93.BoundaryCondition = false;
p93.Notes = "";
p93.Tag = "";

p94 = addparameter(m1, "fTact");
p94.Constant = true;
p94.Value = 1;
p94.Units = "";
p94.BoundaryCondition = false;
p94.Notes = "";
p94.Tag = "";

p95 = addparameter(m1, "fKmTB_kill");
p95.Constant = true;
p95.Value = 1;
p95.Units = "";
p95.BoundaryCondition = false;
p95.Notes = "";
p95.Tag = "";

p96 = addparameter(m1, "Kptumor");
p96.Constant = true;
p96.Value = 1;
p96.Units = "";
p96.BoundaryCondition = false;
p96.Notes = "";
p96.Tag = "";

p97 = addparameter(m1, "Vtumor");
p97.Constant = true;
p97.Value = 1;
p97.Units = "";
p97.BoundaryCondition = false;
p97.Notes = "";
p97.Tag = "";

p98 = addparameter(m1, "KBptumor");
p98.Constant = true;
p98.Value = 1;
p98.Units = "";
p98.BoundaryCondition = false;
p98.Notes = "";
p98.Tag = "";

p99 = addparameter(m1, "KTrptumor");
p99.Constant = true;
p99.Value = 1;
p99.Units = "";
p99.BoundaryCondition = false;
p99.Notes = "";
p99.Tag = "";

p100 = addparameter(m1, "tumor_on");
p100.Constant = true;
p100.Value = 0;
p100.Units = "";
p100.BoundaryCondition = false;
p100.Notes = "";
p100.Tag = "";

p101 = addparameter(m1, "IL6_tiss_contribution");
p101.Constant = true;
p101.Value = 0.02;
p101.Units = "";
p101.BoundaryCondition = false;
p101.Notes = "";
p101.Tag = "";

p102 = addparameter(m1, "kBtumorprolif");
p102.Constant = true;
p102.Value = 0.025;
p102.Units = "";
p102.BoundaryCondition = false;
p102.Notes = "";
p102.Tag = "";

p103 = addparameter(m1, "S");
p103.Constant = true;
p103.Value = 1;
p103.Units = "";
p103.BoundaryCondition = false;
p103.Notes = "";
p103.Tag = "";

p104 = addparameter(m1, "tissue1on");
p104.Constant = true;
p104.Value = 1;
p104.Units = "";
p104.BoundaryCondition = false;
p104.Notes = "";
p104.Tag = "";

p105 = addparameter(m1, "kabs_TDB");
p105.Constant = true;
p105.Value = 1.4;
p105.Units = "";
p105.BoundaryCondition = false;
p105.Notes = "";
p105.Tag = "";

p106 = addparameter(m1, "fbio_TDB");
p106.Constant = true;
p106.Value = 0.6;
p106.Units = "";
p106.BoundaryCondition = false;
p106.Notes = "";
p106.Tag = "";

p107 = addparameter(m1, "Bcell_tumor_trafficking_on");
p107.Constant = true;
p107.Value = 1;
p107.Units = "";
p107.BoundaryCondition = false;
p107.Notes = "";
p107.Tag = "";

%% Create rules.
rule1 = addrule(m1, "restTpb = (1-depleteTpb)*Trpbo_perml*Vpb");
rule1.RuleType = "initialAssignment";
rule1.Active = true;
rule1.Name = "Rule_1";
rule1.Notes = "";
rule1.Tag = "";

rule2 = addrule(m1, "restTtiss = KTrp*Trpbo_perml*Vtissue");
rule2.RuleType = "initialAssignment";
rule2.Active = true;
rule2.Name = "Rule_2";
rule2.Notes = "";
rule2.Tag = "";

rule3 = addrule(m1, "restTtiss2 = KTrp2*Trpbo_perml*Vtissue2");
rule3.RuleType = "initialAssignment";
rule3.Active = true;
rule3.Name = "Rule_3";
rule3.Notes = "";
rule3.Tag = "";

rule4 = addrule(m1, "restTtiss3 = KTrp3*Trpbo_perml*Vtissue3");
rule4.RuleType = "initialAssignment";
rule4.Active = true;
rule4.Name = "Rule_4";
rule4.Notes = "";
rule4.Tag = "";

rule5 = addrule(m1, "restTtumor = KTrptumor*Trpbo_perml*Vtumor");
rule5.RuleType = "initialAssignment";
rule5.Active = true;
rule5.Name = "Rule_5";
rule5.Notes = "";
rule5.Tag = "";

rule6 = addrule(m1, "Bpb = (1-depleteBpb)*Bpbo_perml*Vpb");
rule6.RuleType = "initialAssignment";
rule6.Active = true;
rule6.Name = "Rule_6";
rule6.Notes = "";
rule6.Tag = "";

rule7 = addrule(m1, "Btiss = KBp*Bpbo_perml*Vtissue");
rule7.RuleType = "initialAssignment";
rule7.Active = true;
rule7.Name = "Rule_7";
rule7.Notes = "";
rule7.Tag = "";

rule8 = addrule(m1, "Btiss2 = KBp2*Bpbo_perml*Vtissue2");
rule8.RuleType = "initialAssignment";
rule8.Active = true;
rule8.Name = "Rule_8";
rule8.Notes = "";
rule8.Tag = "";

rule9 = addrule(m1, "B1920tiss3 = Bpbref_perml*KBp3*Vtissue3");
rule9.RuleType = "initialAssignment";
rule9.Active = true;
rule9.Name = "Rule_9";
rule9.Notes = "";
rule9.Tag = "";

rule10 = addrule(m1, "B19no20tiss3 = Bpbref_perml*KBp3*B19no20_B1920_ratio*Vtissue3");
rule10.RuleType = "initialAssignment";
rule10.Active = true;
rule10.Name = "Rule_10";
rule10.Notes = "";
rule10.Tag = "";

rule11 = addrule(m1, "Btumor = KBptumor*Bpbo_perml*Vtumor");
rule11.RuleType = "initialAssignment";
rule11.Active = true;
rule11.Name = "Rule_11";
rule11.Notes = "";
rule11.Tag = "";

rule12 = addrule(m1, "BAFF = BAFFo");
rule12.RuleType = "initialAssignment";
rule12.Active = true;
rule12.Name = "Rule_12";
rule12.Notes = "";
rule12.Tag = "";

rule13 = addrule(m1, "TDBc_ugperml = PK_v26(TDBc_ugperkg, Vc_tdb, PKflag, VPid, time, end_time, fvalidation)");
rule13.RuleType = "repeatedAssignment";
rule13.Active = true;
rule13.Name = "Rule_13";
rule13.Notes = "";
rule13.Tag = "";

rule14 = addrule(m1, "TDBt_ugperml = Kp*TDBc_ugperml");
rule14.RuleType = "repeatedAssignment";
rule14.Active = true;
rule14.Name = "Rule_14";
rule14.Notes = "";
rule14.Tag = "";

rule15 = addrule(m1, "TDBt2_ugperml = Kp2*TDBc_ugperml");
rule15.RuleType = "repeatedAssignment";
rule15.Active = true;
rule15.Name = "Rule_15";
rule15.Notes = "";
rule15.Tag = "";

rule16 = addrule(m1, "TDBt3_ugperml = tissue3on*Kp3*TDBc_ugperml");
rule16.RuleType = "repeatedAssignment";
rule16.Active = true;
rule16.Name = "Rule_16";
rule16.Notes = "";
rule16.Tag = "";

rule17 = addrule(m1, "TDBtumor_ugperml = tumor_on*Kptumor*TDBc_ugperml");
rule17.RuleType = "repeatedAssignment";
rule17.Active = true;
rule17.Name = "Rule_17";
rule17.Notes = "";
rule17.Tag = "";

rule18 = addrule(m1, "RTXc_ugperml = (RTXc_ugperkg/Vc_rtx>1e-5)*RTXc_ugperkg/Vc_rtx");
rule18.RuleType = "repeatedAssignment";
rule18.Active = true;
rule18.Name = "Rule_18";
rule18.Notes = "";
rule18.Tag = "";

rule19 = addrule(m1, "RTXt_ugperml = RTXc_ugperml*Kp");
rule19.RuleType = "repeatedAssignment";
rule19.Active = true;
rule19.Name = "Rule_19";
rule19.Notes = "";
rule19.Tag = "";

rule20 = addrule(m1, "RTXt2_ugperml = RTXc_ugperml*Kp2");
rule20.RuleType = "repeatedAssignment";
rule20.Active = true;
rule20.Name = "Rule_20";
rule20.Notes = "";
rule20.Tag = "";

rule21 = addrule(m1, "RTXt3_ugperml = tissue3on*RTXc_ugperml*Kp3");
rule21.RuleType = "repeatedAssignment";
rule21.Active = true;
rule21.Name = "Rule_21";
rule21.Notes = "";
rule21.Tag = "";

rule22 = addrule(m1, "RTXtumor_ugperml = tumor_on*RTXc_ugperml*Kptumor");
rule22.RuleType = "repeatedAssignment";
rule22.Active = true;
rule22.Name = "Rule_22";
rule22.Notes = "";
rule22.Tag = "";

rule23 = addrule(m1, "Blinc_ngperml = Blinc_ug/Vz_blin*1000");
rule23.RuleType = "repeatedAssignment";
rule23.Active = true;
rule23.Name = "Rule_23";
rule23.Notes = "";
rule23.Tag = "";

rule24 = addrule(m1, "Blint_ngperml = Kp*Blinc_ngperml");
rule24.RuleType = "repeatedAssignment";
rule24.Active = true;
rule24.Name = "Rule_24";
rule24.Notes = "";
rule24.Tag = "";

rule25 = addrule(m1, "Blint2_ngperml = Kp2*Blinc_ngperml");
rule25.RuleType = "repeatedAssignment";
rule25.Active = true;
rule25.Name = "Rule_25";
rule25.Notes = "";
rule25.Tag = "";

rule26 = addrule(m1, "Blint3_ngperml = tissue3on*Kp3*Blinc_ngperml");
rule26.RuleType = "repeatedAssignment";
rule26.Active = true;
rule26.Name = "Rule_26";
rule26.Notes = "";
rule26.Tag = "";

rule27 = addrule(m1, "Blintumor_ngperml = tumor_on*Kptumor*Blinc_ngperml");
rule27.RuleType = "repeatedAssignment";
rule27.Active = true;
rule27.Name = "Rule_27";
rule27.Notes = "";
rule27.Tag = "";

rule28 = addrule(m1, "restTpb_perml = restTpb/Vpb");
rule28.RuleType = "repeatedAssignment";
rule28.Active = true;
rule28.Name = "Rule_28";
rule28.Notes = "";
rule28.Tag = "";

rule29 = addrule(m1, "restTtiss_perml = restTtiss/Vtissue");
rule29.RuleType = "repeatedAssignment";
rule29.Active = true;
rule29.Name = "Rule_29";
rule29.Notes = "";
rule29.Tag = "";

rule30 = addrule(m1, "restTtiss2_perml = restTtiss2/Vtissue2");
rule30.RuleType = "repeatedAssignment";
rule30.Active = true;
rule30.Name = "Rule_30";
rule30.Notes = "";
rule30.Tag = "";

rule31 = addrule(m1, "restTtiss3_perml = restTtiss3/Vtissue3");
rule31.RuleType = "repeatedAssignment";
rule31.Active = true;
rule31.Name = "Rule_31";
rule31.Notes = "";
rule31.Tag = "";

rule32 = addrule(m1, "restTtumor_perml = restTtumor/Vtumor");
rule32.RuleType = "repeatedAssignment";
rule32.Active = true;
rule32.Name = "Rule_32";
rule32.Notes = "";
rule32.Tag = "";

rule33 = addrule(m1, "act0Tpb_perml = act0Tpb/Vpb");
rule33.RuleType = "repeatedAssignment";
rule33.Active = true;
rule33.Name = "Rule_33";
rule33.Notes = "";
rule33.Tag = "";

rule34 = addrule(m1, "act0Ttiss_perml = act0Ttiss/Vtissue");
rule34.RuleType = "repeatedAssignment";
rule34.Active = true;
rule34.Name = "Rule_34";
rule34.Notes = "";
rule34.Tag = "";

rule35 = addrule(m1, "act0Ttiss2_perml = act0Ttiss2/Vtissue2");
rule35.RuleType = "repeatedAssignment";
rule35.Active = true;
rule35.Name = "Rule_35";
rule35.Notes = "";
rule35.Tag = "";

rule36 = addrule(m1, "act0Ttiss3_perml = act0Ttiss3/Vtissue3");
rule36.RuleType = "repeatedAssignment";
rule36.Active = true;
rule36.Name = "Rule_36";
rule36.Notes = "";
rule36.Tag = "";

rule37 = addrule(m1, "act0Ttumor_perml = act0Ttumor/Vtumor");
rule37.RuleType = "repeatedAssignment";
rule37.Active = true;
rule37.Name = "Rule_37";
rule37.Notes = "";
rule37.Tag = "";

rule38 = addrule(m1, "actTpb_perml = actTpb/Vpb");
rule38.RuleType = "repeatedAssignment";
rule38.Active = true;
rule38.Name = "Rule_38";
rule38.Notes = "";
rule38.Tag = "";

rule39 = addrule(m1, "actTtiss_perml = actTtiss/Vtissue");
rule39.RuleType = "repeatedAssignment";
rule39.Active = true;
rule39.Name = "Rule_39";
rule39.Notes = "";
rule39.Tag = "";

rule40 = addrule(m1, "actTtiss2_perml = actTtiss2/Vtissue2");
rule40.RuleType = "repeatedAssignment";
rule40.Active = true;
rule40.Name = "Rule_40";
rule40.Notes = "";
rule40.Tag = "";

rule41 = addrule(m1, "actTtiss3_perml = actTtiss3/Vtissue3");
rule41.RuleType = "repeatedAssignment";
rule41.Active = true;
rule41.Name = "Rule_41";
rule41.Notes = "";
rule41.Tag = "";

rule42 = addrule(m1, "actTtumor_perml = actTtumor/Vtumor");
rule42.RuleType = "repeatedAssignment";
rule42.Active = true;
rule42.Name = "Rule_42";
rule42.Notes = "";
rule42.Tag = "";

rule43 = addrule(m1, "B19tiss3 = B19no20tiss3+B1920tiss3");
rule43.RuleType = "repeatedAssignment";
rule43.Active = true;
rule43.Name = "Rule_43";
rule43.Notes = "";
rule43.Tag = "";

rule44 = addrule(m1, "Bpb_perml = Bpb/Vpb");
rule44.RuleType = "repeatedAssignment";
rule44.Active = true;
rule44.Name = "Rule_44";
rule44.Notes = "";
rule44.Tag = "";

rule45 = addrule(m1, "Btiss_perml = Btiss/Vtissue");
rule45.RuleType = "repeatedAssignment";
rule45.Active = true;
rule45.Name = "Rule_45";
rule45.Notes = "";
rule45.Tag = "";

rule46 = addrule(m1, "Btiss2_perml = Btiss2/Vtissue2");
rule46.RuleType = "repeatedAssignment";
rule46.Active = true;
rule46.Name = "Rule_46";
rule46.Notes = "";
rule46.Tag = "";

rule47 = addrule(m1, "B19tiss3_perml = B19tiss3/Vtissue3");
rule47.RuleType = "repeatedAssignment";
rule47.Active = true;
rule47.Name = "Rule_47";
rule47.Notes = "";
rule47.Tag = "";

rule48 = addrule(m1, "B1920tiss3_perml = B1920tiss3/Vtissue3");
rule48.RuleType = "repeatedAssignment";
rule48.Active = true;
rule48.Name = "Rule_48";
rule48.Notes = "";
rule48.Tag = "";

rule49 = addrule(m1, "B19no20tiss3_perml = B19no20tiss3/Vtissue3");
rule49.RuleType = "repeatedAssignment";
rule49.Active = true;
rule49.Name = "Rule_49";
rule49.Notes = "";
rule49.Tag = "";

rule50 = addrule(m1, "Btumor_perml = Btumor/Vtumor");
rule50.RuleType = "repeatedAssignment";
rule50.Active = true;
rule50.Name = "Rule_50";
rule50.Notes = "";
rule50.Tag = "";

rule51 = addrule(m1, "BTrratio_pb = Bpb/max(restTpb+act0Tpb,1)");
rule51.RuleType = "repeatedAssignment";
rule51.Active = true;
rule51.Name = "Rule_51";
rule51.Notes = "";
rule51.Tag = "";

rule52 = addrule(m1, "BTrratio_tiss = Btiss/max(restTtiss+act0Ttiss,1)");
rule52.RuleType = "repeatedAssignment";
rule52.Active = true;
rule52.Name = "Rule_52";
rule52.Notes = "";
rule52.Tag = "";

rule53 = addrule(m1, "BTrratio_tiss2 = Btiss2/max(restTtiss2+act0Ttiss2,1)");
rule53.RuleType = "repeatedAssignment";
rule53.Active = true;
rule53.Name = "Rule_53";
rule53.Notes = "";
rule53.Tag = "";

rule54 = addrule(m1, "B19Trratio_tiss3 = B19tiss3/max(restTtiss3+act0Ttiss3,1)");
rule54.RuleType = "repeatedAssignment";
rule54.Active = true;
rule54.Name = "Rule_54";
rule54.Notes = "";
rule54.Tag = "";

rule55 = addrule(m1, "B1920Trratio_tiss3 = B1920tiss3/max(restTtiss3+act0Ttiss3,1)");
rule55.RuleType = "repeatedAssignment";
rule55.Active = true;
rule55.Name = "Rule_55";
rule55.Notes = "";
rule55.Tag = "";

rule56 = addrule(m1, "BTrratio_tumor = Btumor/max(restTtumor+act0Ttumor,1)");
rule56.RuleType = "repeatedAssignment";
rule56.Active = true;
rule56.Name = "Rule_56";
rule56.Notes = "";
rule56.Tag = "";

rule57 = addrule(m1, "TaBratio_pb = actTpb/max(Bpb,1)");
rule57.RuleType = "repeatedAssignment";
rule57.Active = true;
rule57.Name = "Rule_57";
rule57.Notes = "";
rule57.Tag = "";

rule58 = addrule(m1, "TaBratio_tiss = actTtiss/max(Btiss,1)");
rule58.RuleType = "repeatedAssignment";
rule58.Active = true;
rule58.Name = "Rule_58";
rule58.Notes = "";
rule58.Tag = "";

rule59 = addrule(m1, "TaBratio_tiss2 = actTtiss2/max(Btiss2,1)");
rule59.RuleType = "repeatedAssignment";
rule59.Active = true;
rule59.Name = "Rule_59";
rule59.Notes = "";
rule59.Tag = "";

rule60 = addrule(m1, "TaB19ratio_tiss3 = actTtiss3/max(B19tiss3,1)");
rule60.RuleType = "repeatedAssignment";
rule60.Active = true;
rule60.Name = "Rule_60";
rule60.Notes = "";
rule60.Tag = "";

rule61 = addrule(m1, "TaB1920ratio_tiss3 = actTtiss3/max(B1920tiss3,1)");
rule61.RuleType = "repeatedAssignment";
rule61.Active = true;
rule61.Name = "Rule_61";
rule61.Notes = "";
rule61.Tag = "";

rule62 = addrule(m1, "TaBratio_tumor = actTtumor/max(Btumor,1)");
rule62.RuleType = "repeatedAssignment";
rule62.Active = true;
rule62.Name = "Rule_62";
rule62.Notes = "";
rule62.Tag = "";

rule63 = addrule(m1, "Tafraction_pb = actTpb_perml/max(totTpb_perml,1)");
rule63.RuleType = "repeatedAssignment";
rule63.Active = true;
rule63.Name = "Rule_63";
rule63.Notes = "";
rule63.Tag = "";

rule64 = addrule(m1, "Tafraction_tiss = actTtiss_perml/max(totTtiss_perml,1)");
rule64.RuleType = "repeatedAssignment";
rule64.Active = true;
rule64.Name = "Rule_64";
rule64.Notes = "";
rule64.Tag = "";

rule65 = addrule(m1, "Tafraction_tiss2 = actTtiss2_perml/max(totTtiss2_perml,1)");
rule65.RuleType = "repeatedAssignment";
rule65.Active = true;
rule65.Name = "Rule_65";
rule65.Notes = "";
rule65.Tag = "";

rule66 = addrule(m1, "Tafraction_tiss3 = actTtiss3_perml/max(totTtiss3_perml,1)");
rule66.RuleType = "repeatedAssignment";
rule66.Active = true;
rule66.Name = "Rule_66";
rule66.Notes = "";
rule66.Tag = "";

rule67 = addrule(m1, "Tafraction_tumor = actTtumor_perml/max(totTtumor_perml,1)");
rule67.RuleType = "repeatedAssignment";
rule67.Active = true;
rule67.Name = "Rule_67";
rule67.Notes = "";
rule67.Tag = "";

rule68 = addrule(m1, "totTpb_perml = restTpb_perml+act0Tpb_perml + actTpb_perml");
rule68.RuleType = "repeatedAssignment";
rule68.Active = true;
rule68.Name = "Rule_68";
rule68.Notes = "";
rule68.Tag = "";

rule69 = addrule(m1, "totTtiss_perml = restTtiss_perml + act0Ttiss_perml + actTtiss_perml");
rule69.RuleType = "repeatedAssignment";
rule69.Active = true;
rule69.Name = "Rule_69";
rule69.Notes = "";
rule69.Tag = "";

rule70 = addrule(m1, "totTtiss2_perml = restTtiss2_perml + act0Ttiss2_perml + actTtiss2_perml");
rule70.RuleType = "repeatedAssignment";
rule70.Active = true;
rule70.Name = "Rule_70";
rule70.Notes = "";
rule70.Tag = "";

rule71 = addrule(m1, "totTtiss3_perml = restTtiss3_perml + act0Ttiss3_perml + actTtiss3_perml");
rule71.RuleType = "repeatedAssignment";
rule71.Active = true;
rule71.Name = "Rule_71";
rule71.Notes = "";
rule71.Tag = "";

rule72 = addrule(m1, "totTtumor_perml = restTtumor_perml + act0Ttumor_perml + actTtumor_perml");
rule72.RuleType = "repeatedAssignment";
rule72.Active = true;
rule72.Name = "Rule_72";
rule72.Notes = "";
rule72.Tag = "";

rule73 = addrule(m1, "drugTpbact = VmT*(BTrratio_pb^S/(KmBT_act^S+BTrratio_pb^S))*((TDBc_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBc_ugperml*1000)^ndrugactT))");
rule73.RuleType = "repeatedAssignment";
rule73.Active = true;
rule73.Name = "Rule_73";
rule73.Notes = "";
rule73.Tag = "";

rule74 = addrule(m1, "drugTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act^S+BTrratio_tiss^S))*((TDBt_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt_ugperml*1000)^ndrugactT))");
rule74.RuleType = "repeatedAssignment";
rule74.Active = true;
rule74.Name = "Rule_74";
rule74.Notes = "";
rule74.Tag = "";

rule75 = addrule(m1, "drugTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act^S+BTrratio_tiss2^S))*((TDBt2_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt2_ugperml*1000)^ndrugactT))");
rule75.RuleType = "repeatedAssignment";
rule75.Active = true;
rule75.Name = "Rule_75";
rule75.Notes = "";
rule75.Tag = "";

rule76 = addrule(m1, "drugTtissact3 = fTact*VmT*(B1920Trratio_tiss3^S/(KmBT_act^S+B1920Trratio_tiss3^S))*((TDBt3_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBt3_ugperml*1000)^ndrugactT))");
rule76.RuleType = "repeatedAssignment";
rule76.Active = true;
rule76.Name = "Rule_76";
rule76.Notes = "";
rule76.Tag = "";

rule77 = addrule(m1, "drugTtumoract = fTact*VmT*(max(BTrratio_tumor, 0)^S/(KmBT_act^S+max(BTrratio_tumor, 0)^S))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT))");
rule77.RuleType = "repeatedAssignment";
rule77.Active = true;
rule77.Name = "Rule_77";
rule77.Notes = "";
rule77.Tag = "";

rule78 = addrule(m1, "drugBpbkill = VmB*max(0, TaBratio_pb)^nkill/(max(0, KmTB_kill)^nkill+max(0, TaBratio_pb)^nkill)*((TDBc_ugperml*1000)/(KdrugB+(TDBc_ugperml*1000)))");
rule78.RuleType = "repeatedAssignment";
rule78.Active = true;
rule78.Name = "Rule_78";
rule78.Notes = "";
rule78.Tag = "";

rule79 = addrule(m1, "drugBtisskill = VmB*max(0, TaBratio_tiss)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss)^nkill)*((TDBt_ugperml*1000)/(KdrugB+(TDBt_ugperml*1000)))");
rule79.RuleType = "repeatedAssignment";
rule79.Active = true;
rule79.Name = "Rule_79";
rule79.Notes = "";
rule79.Tag = "";

rule80 = addrule(m1, "drugBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tiss2)^nkill)*((TDBt2_ugperml*1000)/(KdrugB+(TDBt2_ugperml*1000)))");
rule80.RuleType = "repeatedAssignment";
rule80.Active = true;
rule80.Name = "Rule_80";
rule80.Notes = "";
rule80.Tag = "";

rule81 = addrule(m1, "drugBtisskill3 = VmB*max(0,TaB1920ratio_tiss3)^nkill/(max(0,KmTB_kill*fKmTB_kill)^nkill+max(0,TaB1920ratio_tiss3)^nkill)*((TDBt3_ugperml*1000)/(KdrugB+(TDBt3_ugperml*1000)))");
rule81.RuleType = "repeatedAssignment";
rule81.Active = true;
rule81.Name = "Rule_81";
rule81.Notes = "";
rule81.Tag = "";

rule82 = addrule(m1, "drugBtumorkill = VmB*max(0, TaBratio_tumor)^nkill/(max(0, KmTB_kill*fKmTB_kill)^nkill+max(0, TaBratio_tumor)^nkill)*((TDBtumor_ugperml*1000)/(KdrugB+(TDBtumor_ugperml*1000)))");
rule82.RuleType = "repeatedAssignment";
rule82.Active = true;
rule82.Name = "Rule_82";
rule82.Notes = "";
rule82.Tag = "";

rule83 = addrule(m1, "BlinTpbact = VmT*(BTrratio_pb^S/(KmBT_act_blin^S+BTrratio_pb^S))*(Blinc_ngperml^ndrugactT_blin/(Blinc_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))");
rule83.RuleType = "repeatedAssignment";
rule83.Active = true;
rule83.Name = "Rule_83";
rule83.Notes = "";
rule83.Tag = "";

rule84 = addrule(m1, "BlinTtissact = fTact*VmT*(BTrratio_tiss^S/(KmBT_act_blin^S+BTrratio_tiss^S))*(Blint_ngperml^ndrugactT_blin/(Blint_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))");
rule84.RuleType = "repeatedAssignment";
rule84.Active = true;
rule84.Name = "Rule_84";
rule84.Notes = "";
rule84.Tag = "";

rule85 = addrule(m1, "BlinTtissact2 = fTact*VmT*(BTrratio_tiss2^S/(KmBT_act_blin^S+BTrratio_tiss2^S))*(Blint2_ngperml^ndrugactT_blin/(Blint2_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))");
rule85.RuleType = "repeatedAssignment";
rule85.Active = true;
rule85.Name = "Rule_85";
rule85.Notes = "";
rule85.Tag = "";

rule86 = addrule(m1, "BlinTtissact3 = fTact*VmT*(B19Trratio_tiss3^S/(KmBT_act_blin^S+B19Trratio_tiss3^S))*(Blint3_ngperml^ndrugactT_blin/(Blint3_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))");
rule86.RuleType = "repeatedAssignment";
rule86.Active = true;
rule86.Name = "Rule_86";
rule86.Notes = "";
rule86.Tag = "";

rule87 = addrule(m1, "BlinTtumoract = fTact*VmT*(max(BTrratio_tumor,0)^S/(KmBT_act_blin^S+max(BTrratio_tumor, 0)^S))*(Blintumor_ngperml^ndrugactT_blin/(Blintumor_ngperml^ndrugactT_blin+KdrugactT_blin^ndrugactT_blin))");
rule87.RuleType = "repeatedAssignment";
rule87.Active = true;
rule87.Name = "Rule_87";
rule87.Notes = "";
rule87.Tag = "";

rule88 = addrule(m1, "BlinBpbkill = VmB*max(0, TaBratio_pb)^nkill_blin/(max(0, KmTB_kill_blin)^nkill_blin+max(0, TaBratio_pb)^nkill_blin)*(Blinc_ngperml/(KdrugB_blin+Blinc_ngperml))");
rule88.RuleType = "repeatedAssignment";
rule88.Active = true;
rule88.Name = "Rule_88";
rule88.Notes = "";
rule88.Tag = "";

rule89 = addrule(m1, "BlinBtisskill = VmB*max(0, TaBratio_tiss)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss)^nkill_blin)*(Blint_ngperml/(KdrugB_blin+Blint_ngperml))");
rule89.RuleType = "repeatedAssignment";
rule89.Active = true;
rule89.Name = "Rule_89";
rule89.Notes = "";
rule89.Tag = "";

rule90 = addrule(m1, "BlinBtisskill2 = VmB*max(0, TaBratio_tiss2)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tiss2)^nkill_blin)*(Blint2_ngperml/(KdrugB_blin+Blint2_ngperml))");
rule90.RuleType = "repeatedAssignment";
rule90.Active = true;
rule90.Name = "Rule_90";
rule90.Notes = "";
rule90.Tag = "";

rule91 = addrule(m1, "BlinBtisskill3 = VmB*max(0,TaB19ratio_tiss3)^nkill_blin/(max(0,KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0,TaB19ratio_tiss3)^nkill_blin)*(Blint3_ngperml/(KdrugB_blin+Blint3_ngperml))");
rule91.RuleType = "repeatedAssignment";
rule91.Active = true;
rule91.Name = "Rule_91";
rule91.Notes = "";
rule91.Tag = "";

rule92 = addrule(m1, "BlinBtumorkill = VmB*max(0, TaBratio_tumor)^nkill_blin/(max(0, KmTB_kill_blin*fKmTB_kill)^nkill_blin+max(0, TaBratio_tumor)^nkill_blin)*(Blintumor_ngperml/(KdrugB_blin+Blintumor_ngperml))");
rule92.RuleType = "repeatedAssignment";
rule92.Active = true;
rule92.Name = "Rule_92";
rule92.Notes = "";
rule92.Tag = "";

rule93 = addrule(m1, "Bpb_norm = Bpb_perml/Bpbref_perml");
rule93.RuleType = "repeatedAssignment";
rule93.Active = true;
rule93.Name = "Rule_93";
rule93.Notes = "";
rule93.Tag = "";

rule94 = addrule(m1, "totTtiss = restTtiss + act0Ttiss + actTtiss");
rule94.RuleType = "repeatedAssignment";
rule94.Active = true;
rule94.Name = "Rule_94";
rule94.Notes = "";
rule94.Tag = "";

rule95 = addrule(m1, "totTtiss2 = restTtiss2 + act0Ttiss2 + actTtiss2");
rule95.RuleType = "repeatedAssignment";
rule95.Active = true;
rule95.Name = "Rule_95";
rule95.Notes = "";
rule95.Tag = "";

rule96 = addrule(m1, "totTtiss3 = restTtiss3 + act0Ttiss3 + actTtiss3");
rule96.RuleType = "repeatedAssignment";
rule96.Active = true;
rule96.Name = "Rule_96";
rule96.Notes = "";
rule96.Tag = "";

rule97 = addrule(m1, "totTtumor = restTtumor + act0Ttumor + actTtumor");
rule97.RuleType = "repeatedAssignment";
rule97.Active = true;
rule97.Name = "Rule_97";
rule97.Notes = "";
rule97.Tag = "";

rule98 = addrule(m1, "Baffconsumption = log(2)/thBAFF*BAFF*(Btiss+Bpb + Btiss2)/(Bpbref_perml*(Vpb+KBp*Vtissue+KBp2*Vtissue2))");
rule98.RuleType = "repeatedAssignment";
rule98.Active = true;
rule98.Name = "Rule_98";
rule98.Notes = "";
rule98.Tag = "";

rule99 = addrule(m1, "BTtotRatio_tiss = Btiss_perml/totTtiss_perml");
rule99.RuleType = "repeatedAssignment";
rule99.Active = true;
rule99.Name = "rule_2";
rule99.Notes = "";
rule99.Tag = "";

rule100 = addrule(m1, "BTtotRatio_tiss2 = Btiss2_perml/totTtiss2_perml");
rule100.RuleType = "repeatedAssignment";
rule100.Active = true;
rule100.Name = "rule_3";
rule100.Notes = "";
rule100.Tag = "";

rule101 = addrule(m1, "unnamed.IL6combo = IL6pb+IL6_tiss_contribution*(IL6tiss*Vtissue+IL6tiss2*Vtissue2+IL6tiss3*Vtissue3+IL6tumor*Vtumor)/Vpb");
rule101.RuleType = "repeatedAssignment";
rule101.Active = true;
rule101.Name = "rule_4";
rule101.Notes = "";
rule101.Tag = "";

%% Create and configure configsets.
cs1 = getconfigset(m1, "default");
setactiveconfigset(m1, cs1);
cs1.Notes = "";
cs1.SolverType = "ode15s";
cs1.StopTime = 40;
cs1.MaximumNumberOfLogs = Inf;
cs1.MaximumWallClock = Inf;
cs1.TimeUnits = "day";
cs1.AmountUnits = "<automatic>";
cs1.MassUnits = "<automatic>";
cs1.CompileOptions.DefaultSpeciesDimension = "concentration";
cs1.CompileOptions.DimensionalAnalysis = false;
cs1.CompileOptions.UnitConversion = false;
cs1.SolverOptions.AbsoluteTolerance = 1e-06;
cs1.SolverOptions.AbsoluteToleranceScaling = true;
cs1.SolverOptions.MaxStep = zeros(0,1);
cs1.SolverOptions.OutputTimes = zeros(0,1);
cs1.SolverOptions.RelativeTolerance = 0.001;
cs1.SolverOptions.SensitivityAnalysis = false;
cs1.SolverOptions.AbsoluteToleranceStepSize = zeros(0,1);
cs1.SensitivityAnalysisOptions.Normalization = "None";
cs1.RunTimeOptions.StatesToLog = ["actTtiss";"Btiss";"TDBc_ugperkg";"actTpb";"drugBtisskill";"drugTtissact";"restTtiss";"TDBp_ugperkg";"restTpb";"Bpb";"TDBc_ugperml";"TDBt_ugperml";"drugTpbact";"drugBpbkill";"restTtiss_perml";"actTtiss_perml";"Btiss_perml";"restTpb_perml";"actTpb_perml";"Bpb_perml";"BTrratio_pb";"BTrratio_tiss";"TaBratio_pb";"TaBratio_tiss";"Tafraction_pb";"totTpb_perml";"totTtiss_perml";"act0Ttiss_perml";"act0Tpb_perml";"Tafraction_tiss";"BAFF";"Baffconsumption";"act0Tpb";"act0Ttiss";"injection_effect";"Btiss2";"act0Tpb_1";"Bpb_1";"restTpb_1";"actTpb_1";"act0Ttiss2";"restTtiss2";"TDBt2_ugperml";"drugTtissact2";"drugBtisskill2";"actTtiss2";"Tafraction_tiss2";"act0Ttiss2_perml";"totTtiss2_perml";"TaBratio_tiss2";"BTrratio_tiss2";"Btiss2_perml";"actTtiss2_perml";"restTtiss2_perml";"RTXc_ugperkg";"RTXp_ugperkg";"RTXc_ugperml";"RTXt_ugperml";"RTXt2_ugperml";"drug_effect";"B1920tiss3";"B19no20tiss3";"restTtiss3";"act0Ttiss3";"actTtiss3";"drugTtissact3";"TDBt3_ugperml";"B19Trratio_tiss3";"B19tiss3";"restTtiss3_perml";"act0Ttiss3_perml";"actTtiss3_perml";"B19tiss3_perml";"Tafraction_tiss3";"totTtiss3_perml";"TaB19ratio_tiss3";"B19TtotRatio_tiss3";"RTXt3_ugperml";"drugBtisskill3";"TaB1920ratio_tiss3";"B1920tiss3_perml";"B1920Trratio_tiss3";"B19no20tiss3_perml";"Blinc_ug";"Blinc_ngperml";"BlinTpbact";"BlinBpbkill";"BlinTtissact";"BlinBtisskill";"Blint_ngperml";"Blint2_ngperml";"BlinBtisskill2";"BlinTtissact2";"BlinBtisskill3";"Blint3_ngperml";"BlinTtissact3";"Bpb_norm";"Tafraction_pb_init";"totTtiss";"totTtiss2";"totTtiss3";"restTtumor";"actTtumor";"Btumor";"act0Ttumor";"drugTtumoract";"drugBtumorkill";"BlinTtumoract";"BlinBtumorkill";"TDBtumor_ugperml";"restTtumor_perml";"act0Ttumor_perml";"actTtumor_perml";"Btumor_perml";"Tafraction_tumor";"BTrratio_tumor";"TaBratio_tumor";"totTtumor";"totTtumor_perml";"Blintumor_ngperml";"RTXtumor_ugperml";"BTtotRatio_tiss";"BTtotRatio_tiss2";"TDBsc_ugperkg";"TDBc_ugperml_AUC";"IL6pb";"IL6tiss";"IL6tiss2";"IL6tiss3";"IL6tumor";"IL6combo"];

%% Create variants.
v1 = addvariant(m1, "act0 unopt");
v1.addcontent({'parameter', 'fTaprolif', 'Value', 0.65});
v1.addcontent({'parameter', 'kTaapop', 'Value', 0.05});
v1.addcontent({'parameter', 'fTrapop', 'Value', 5000});
v1.addcontent({'parameter', 'fTa0apop', 'Value', 1});
v1.addcontent({'parameter', 'fa0', 'Value', 0});
v1.addcontent({'parameter', 'fAICD', 'Value', 1});
v1.addcontent({'parameter', 'kTact', 'Value', 4});
v1.addcontent({'parameter', 'fTadeact', 'Value', 0.01});
v1.addcontent({'parameter', 'fTa0deact', 'Value', 0.1});
v1.addcontent({'parameter', 'kBkill', 'Value', 10000});
v1.addcontent({'parameter', 'kTaexit', 'Value', 100});
v1.addcontent({'parameter', 'finj', 'Value', 0.1});
v1.addcontent({'parameter', 'fBexit', 'Value', 0});
v1.addcontent({'parameter', 'kTgen', 'Value', 0});
v1.addcontent({'parameter', 'KdrugactT', 'Value', 150});
v1.addcontent({'parameter', 'KdrugB', 'Value', 10});
v1.addcontent({'parameter', 'KmBT_act', 'Value', 1});
v1.addcontent({'parameter', 'nkill', 'Value', 2});
v1.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v1.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v1.addcontent({'parameter', 'kBapop', 'Value', 0.015});
v1.addcontent({'parameter', 'kBprolif', 'Value', 0.7});
v1.addcontent({'parameter', 'fBprolif', 'Value', 0});
v1.addcontent({'parameter', 'act0on', 'Value', 1});
v1.Active = false;
v1.Name = "act0 unopt";
v1.Notes = "";
v1.Tag = "";

v2 = addvariant(m1, "combo tissues 1");
v2.addcontent({'parameter', 'Vpb', 'Value', 380});
v2.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v2.addcontent({'parameter', 'Trpbref_perml', 'Value', 2000000});
v2.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v2.addcontent({'parameter', 'Bpbref_perml', 'Value', 1000000});
v2.addcontent({'parameter', 'tissue1on', 'Value', 1});
v2.addcontent({'parameter', 'Vtissue', 'Value', 7});
v2.addcontent({'parameter', 'KTrp', 'Value', 500});
v2.addcontent({'parameter', 'fTap', 'Value', 11.1});
v2.addcontent({'parameter', 'KBp', 'Value', 900});
v2.addcontent({'parameter', 'tissue2on', 'Value', 1});
v2.addcontent({'parameter', 'Vtissue2', 'Value', 25});
v2.addcontent({'parameter', 'KTrp2', 'Value', 500});
v2.addcontent({'parameter', 'KBp2', 'Value', 600});
v2.Active = false;
v2.Name = "combo tissues 1";
v2.Notes = "";
v2.Tag = "";

v3 = addvariant(m1, "act optimizn 1");
v3.addcontent({'parameter', 'fTaprolif', 'Value', 0.63});
v3.addcontent({'parameter', 'kTaapop', 'Value', 0.05});
v3.addcontent({'parameter', 'fTrapop', 'Value', 3.8});
v3.addcontent({'parameter', 'fTa0apop', 'Value', 1.3});
v3.addcontent({'parameter', 'fa0', 'Value', 0.26});
v3.addcontent({'parameter', 'fAICD', 'Value', 1});
v3.addcontent({'parameter', 'kTact', 'Value', 9.23});
v3.addcontent({'parameter', 'fTadeact', 'Value', 0.1});
v3.addcontent({'parameter', 'fTa0deact', 'Value', 0.02});
v3.addcontent({'parameter', 'kBkill', 'Value', 915});
v3.addcontent({'parameter', 'kTaexit', 'Value', 6577});
v3.addcontent({'parameter', 'fBexit', 'Value', 0});
v3.addcontent({'parameter', 'kTgen', 'Value', 1});
v3.addcontent({'parameter', 'KdrugactT', 'Value', 122});
v3.addcontent({'parameter', 'KdrugB', 'Value', 0.3});
v3.addcontent({'parameter', 'KmBT_act', 'Value', 0.43});
v3.addcontent({'parameter', 'nkill', 'Value', 1});
v3.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v3.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v3.addcontent({'parameter', 'kBapop', 'Value', 0.015});
v3.addcontent({'parameter', 'kBprolif', 'Value', 0.7});
v3.addcontent({'parameter', 'fBprolif', 'Value', 0});
v3.addcontent({'parameter', 'act0on', 'Value', 1});
v3.addcontent({'parameter', 'finj', 'Value', 0.66});
v3.Active = false;
v3.Name = "act optimizn 1";
v3.Notes = "";
v3.Tag = "";

v4 = addvariant(m1, "act optimizn 2a");
v4.addcontent({'parameter', 'fTaprolif', 'Value', 9.5});
v4.addcontent({'parameter', 'kTaapop', 'Value', 0.4});
v4.addcontent({'parameter', 'fTrapop', 'Value', 33});
v4.addcontent({'parameter', 'fTa0apop', 'Value', 3.2});
v4.addcontent({'parameter', 'fa0', 'Value', 0.54});
v4.addcontent({'parameter', 'fAICD', 'Value', 1});
v4.addcontent({'parameter', 'kTact', 'Value', 7.6});
v4.addcontent({'parameter', 'fTadeact', 'Value', 0.001});
v4.addcontent({'parameter', 'fTa0deact', 'Value', 1});
v4.addcontent({'parameter', 'kBkill', 'Value', 470});
v4.addcontent({'parameter', 'kTaexit', 'Value', 1});
v4.addcontent({'parameter', 'fBexit', 'Value', 0});
v4.addcontent({'parameter', 'kTgen', 'Value', 1});
v4.addcontent({'parameter', 'KdrugactT', 'Value', 40});
v4.addcontent({'parameter', 'KdrugB', 'Value', 7.5});
v4.addcontent({'parameter', 'KmBT_act', 'Value', 5.7});
v4.addcontent({'parameter', 'nkill', 'Value', 1});
v4.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v4.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v4.addcontent({'parameter', 'kBapop', 'Value', 0.015});
v4.addcontent({'parameter', 'kBprolif', 'Value', 0.7});
v4.addcontent({'parameter', 'fBprolif', 'Value', 0});
v4.addcontent({'parameter', 'act0on', 'Value', 1});
v4.addcontent({'parameter', 'finj', 'Value', 0.2});
v4.Active = false;
v4.Name = "act optimizn 2a";
v4.Notes = "";
v4.Tag = "";

v5 = addvariant(m1, "combo tissues 2a");
v5.addcontent({'parameter', 'Vpb', 'Value', 380});
v5.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v5.addcontent({'parameter', 'Trpbref_perml', 'Value', 2000000});
v5.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v5.addcontent({'parameter', 'Bpbref_perml', 'Value', 1000000});
v5.addcontent({'parameter', 'tissue1on', 'Value', 1});
v5.addcontent({'parameter', 'Vtissue', 'Value', 7});
v5.addcontent({'parameter', 'KTrp', 'Value', 500});
v5.addcontent({'parameter', 'fTap', 'Value', 2.4});
v5.addcontent({'parameter', 'KBp', 'Value', 900});
v5.addcontent({'parameter', 'tissue2on', 'Value', 1});
v5.addcontent({'parameter', 'Vtissue2', 'Value', 25});
v5.addcontent({'parameter', 'KTrp2', 'Value', 500});
v5.addcontent({'parameter', 'KBp2', 'Value', 600});
v5.Active = true;
v5.Name = "combo tissues 2a";
v5.Notes = "";
v5.Tag = "";

v6 = addvariant(m1, "act optimizn 2b");
v6.addcontent({'parameter', 'fTaprolif', 'Value', 3.76});
v6.addcontent({'parameter', 'kTaapop', 'Value', 0.3});
v6.addcontent({'parameter', 'fTrapop', 'Value', 446});
v6.addcontent({'parameter', 'fTa0apop', 'Value', 0.3});
v6.addcontent({'parameter', 'fa0', 'Value', 0.24});
v6.addcontent({'parameter', 'fAICD', 'Value', 1});
v6.addcontent({'parameter', 'kTact', 'Value', 1.9});
v6.addcontent({'parameter', 'fTadeact', 'Value', 0.004});
v6.addcontent({'parameter', 'fTa0deact', 'Value', 0.012});
v6.addcontent({'parameter', 'kBkill', 'Value', 661});
v6.addcontent({'parameter', 'kTaexit', 'Value', 158});
v6.addcontent({'parameter', 'fBexit', 'Value', 0});
v6.addcontent({'parameter', 'kTgen', 'Value', 1});
v6.addcontent({'parameter', 'KdrugactT', 'Value', 163});
v6.addcontent({'parameter', 'KdrugB', 'Value', 0.4});
v6.addcontent({'parameter', 'KmBT_act', 'Value', 0.025});
v6.addcontent({'parameter', 'nkill', 'Value', 1.44});
v6.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v6.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v6.addcontent({'parameter', 'kBapop', 'Value', 0.015});
v6.addcontent({'parameter', 'kBprolif', 'Value', 0.7});
v6.addcontent({'parameter', 'fBprolif', 'Value', 0});
v6.addcontent({'parameter', 'act0on', 'Value', 1});
v6.addcontent({'parameter', 'finj', 'Value', 0.44});
v6.Active = false;
v6.Name = "act optimizn 2b";
v6.Notes = "";
v6.Tag = "";

v7 = addvariant(m1, "combo tissues 2b");
v7.addcontent({'parameter', 'Vpb', 'Value', 380});
v7.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v7.addcontent({'parameter', 'Trpbref_perml', 'Value', 2000000});
v7.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v7.addcontent({'parameter', 'Bpbref_perml', 'Value', 1000000});
v7.addcontent({'parameter', 'tissue1on', 'Value', 1});
v7.addcontent({'parameter', 'Vtissue', 'Value', 7});
v7.addcontent({'parameter', 'KTrp', 'Value', 500});
v7.addcontent({'parameter', 'fTap', 'Value', 2.6});
v7.addcontent({'parameter', 'KBp', 'Value', 900});
v7.addcontent({'parameter', 'tissue2on', 'Value', 1});
v7.addcontent({'parameter', 'Vtissue2', 'Value', 25});
v7.addcontent({'parameter', 'KTrp2', 'Value', 500});
v7.addcontent({'parameter', 'KBp2', 'Value', 600});
v7.Active = false;
v7.Name = "combo tissues 2b";
v7.Notes = "";
v7.Tag = "";

v8 = addvariant(m1, "initialTdepleted");
v8.addcontent({'parameter', 'depleteTpb', 'Value', 1});
v8.Active = false;
v8.Name = "initialTdepleted";
v8.Notes = "";
v8.Tag = "";

v9 = addvariant(m1, "Tgeneration");
v9.addcontent({'parameter', 'fTgenbl', 'Value', 0});
v9.addcontent({'parameter', 'kTgen', 'Value', 1});
v9.Active = true;
v9.Name = "Tgeneration";
v9.Notes = "";
v9.Tag = "";

v10 = addvariant(m1, "Tgeneration 2b");
v10.addcontent({'parameter', 'fTgenbl', 'Value', 1});
v10.addcontent({'parameter', 'kTgen', 'Value', 1});
v10.Active = false;
v10.Name = "Tgeneration 2b";
v10.Notes = "";
v10.Tag = "";

v11 = addvariant(m1, "act optimizn 2a_modified");
v11.addcontent({'parameter', 'fTaprolif', 'Value', 10});
v11.addcontent({'parameter', 'kTaapop', 'Value', 0.4});
v11.addcontent({'parameter', 'fTrapop', 'Value', 24});
v11.addcontent({'parameter', 'fTa0apop', 'Value', 2.4});
v11.addcontent({'parameter', 'fa0', 'Value', 0.017});
v11.addcontent({'parameter', 'fAICD', 'Value', 1});
v11.addcontent({'parameter', 'kTact', 'Value', 7.6});
v11.addcontent({'parameter', 'fTadeact', 'Value', 0.001});
v11.addcontent({'parameter', 'fTa0deact', 'Value', 1});
v11.addcontent({'parameter', 'kBkill', 'Value', 470});
v11.addcontent({'parameter', 'kTaexit', 'Value', 1});
v11.addcontent({'parameter', 'fBexit', 'Value', 0});
v11.addcontent({'parameter', 'kTgen', 'Value', 1});
v11.addcontent({'parameter', 'KdrugactT', 'Value', 100});
v11.addcontent({'parameter', 'KdrugB', 'Value', 7.5});
v11.addcontent({'parameter', 'KmBT_act', 'Value', 5.7});
v11.addcontent({'parameter', 'nkill', 'Value', 1});
v11.addcontent({'parameter', 'Trpbo_perml', 'Value', 2000000});
v11.addcontent({'parameter', 'Bpbo_perml', 'Value', 1000000});
v11.addcontent({'parameter', 'kBapop', 'Value', 0.03});
v11.addcontent({'parameter', 'kBprolif', 'Value', 0.7});
v11.addcontent({'parameter', 'fBprolif', 'Value', 0});
v11.addcontent({'parameter', 'act0on', 'Value', 1});
v11.addcontent({'parameter', 'fdrug', 'Value', 5.26});
v11.Active = true;
v11.Name = "act optimizn 2a_modified";
v11.Notes = "";
v11.Tag = "";

v12 = addvariant(m1, "initialBdepleted");
v12.addcontent({'parameter', 'depleteBpb', 'Value', 1});
v12.Active = false;
v12.Name = "initialBdepleted";
v12.Notes = "";
v12.Tag = "";

v13 = addvariant(m1, "finj optimization");
v13.addcontent({'parameter', 'kTrexit', 'Value', 0.03});
v13.addcontent({'parameter', 'finj', 'Value', 1.2});
v13.Active = true;
v13.Name = "finj optimization";
v13.Notes = "";
v13.Tag = "";

v14 = addvariant(m1, "tissue3");
v14.addcontent({'parameter', 'KTrp3', 'Value', 50});
v14.addcontent({'parameter', 'KBp3', 'Value', 60});
v14.addcontent({'parameter', 'tissue3on', 'Value', 1});
v14.addcontent({'parameter', 'Vtissue3', 'Value', 50});
v14.addcontent({'parameter', 'kBtiss3exit', 'Value', 0.041});
v14.addcontent({'parameter', 'fBtissue3_v1', 'Value', 1});
v14.addcontent({'parameter', 'kBmat_kBapop_ratio', 'Value', 0.25});
v14.Active = true;
v14.Name = "tissue3";
v14.Notes = "";
v14.Tag = "";

v15 = addvariant(m1, "validation");
v15.addcontent({'parameter', 'fvalidation', 'Value', 1});
v15.Active = false;
v15.Name = "validation";
v15.Notes = "";
v15.Tag = "";

v16 = addvariant(m1, "act optimizn v24");
v16.addcontent({'parameter', 'fTaprolif', 'Value', 9.805999999999999});
v16.addcontent({'parameter', 'kTaapop', 'Value', 0.639});
v16.addcontent({'parameter', 'fTrapop', 'Value', 0.249});
v16.addcontent({'parameter', 'fTa0apop', 'Value', 0.313});
v16.addcontent({'parameter', 'fa0', 'Value', 0.002});
v16.addcontent({'parameter', 'fAICD', 'Value', 22.027});
v16.addcontent({'parameter', 'kTact', 'Value', 8.278});
v16.addcontent({'parameter', 'fTadeact', 'Value', 0.002});
v16.addcontent({'parameter', 'fTa0deact', 'Value', 1});
v16.addcontent({'parameter', 'kBkill', 'Value', 510.411});
v16.addcontent({'parameter', 'kTaexit', 'Value', 1.506});
v16.addcontent({'parameter', 'KdrugactT', 'Value', 600});
v16.addcontent({'parameter', 'KdrugB', 'Value', 0.228});
v16.addcontent({'parameter', 'nkill', 'Value', 1.002});
v16.addcontent({'parameter', 'KmBT_act', 'Value', 3.681});
v16.addcontent({'parameter', 'fTap', 'Value', 1.52});
v16.addcontent({'parameter', 'kBapop', 'Value', 0.182});
v16.addcontent({'parameter', 'fdrug', 'Value', 6.999});
v16.addcontent({'parameter', 'kTrexit', 'Value', 0.056});
v16.addcontent({'parameter', 'finj', 'Value', 0.585});
v16.addcontent({'parameter', 'fapop_v24', 'Value', 1});
v16.addcontent({'parameter', 'ndrugactT', 'Value', 0.7});
v16.addcontent({'parameter', 'fBkill', 'Value', 0.001959});
v16.Active = false;
v16.Name = "act optimizn v24";
v16.Notes = "";
v16.Tag = "";

v17 = addvariant(m1, "tissue3_v2");
v17.addcontent({'parameter', 'KTrp3', 'Value', 50});
v17.addcontent({'parameter', 'KBp3', 'Value', 60});
v17.addcontent({'parameter', 'tissue3on', 'Value', 1});
v17.addcontent({'parameter', 'Vtissue3', 'Value', 50});
v17.addcontent({'parameter', 'kBtiss3exit', 'Value', 1});
v17.addcontent({'parameter', 'fBtissue3_v1', 'Value', 0});
v17.addcontent({'parameter', 'kBmat_kBapop_ratio', 'Value', 0.2218935});
v17.addcontent({'parameter', 'kBapop', 'Value', 0.005});
v17.Active = false;
v17.Name = "tissue3_v2";
v17.Notes = "";
v17.Tag = "";

v18 = addvariant(m1, "B cell migration");
v18.addcontent({'parameter', 'fBexit', 'Value', 1});
v18.Active = false;
v18.Name = "B cell migration";
v18.Notes = "";
v18.Tag = "";

v19 = addvariant(m1, "B cell proliferation");
v19.addcontent({'parameter', 'kBprolif', 'Value', 0.02});
v19.Active = false;
v19.Name = "B cell proliferation";
v19.Notes = "";
v19.Tag = "";

v20 = addvariant(m1, "human parameters");
v20.addcontent({'parameter', 'Vpb', 'Value', 5000});
v20.addcontent({'parameter', 'Vtissue', 'Value', 210});
v20.addcontent({'parameter', 'Vtissue2', 'Value', 400});
v20.addcontent({'parameter', 'Vtissue3', 'Value', 500});
v20.addcontent({'parameter', 'Trpbref_perml', 'Value', 500000});
v20.addcontent({'parameter', 'Trpbo_perml', 'Value', 500000});
v20.addcontent({'parameter', 'Bpbref_perml', 'Value', 500000});
v20.addcontent({'parameter', 'Bpbo_perml', 'Value', 500000});
v20.addcontent({'parameter', 'KTrp', 'Value', 200});
v20.addcontent({'parameter', 'KBp', 'Value', 333});
v20.addcontent({'parameter', 'KTrp2', 'Value', 190});
v20.addcontent({'parameter', 'KBp2', 'Value', 190});
v20.addcontent({'parameter', 'KTrp3', 'Value', 60});
v20.addcontent({'parameter', 'KBp3', 'Value', 80});
v20.addcontent({'parameter', 'Cl_tdb', 'Value', 1.5});
v20.Active = true;
v20.Name = "human parameters";
v20.Notes = "";
v20.Tag = "";

v21 = addvariant(m1, "human ALL parameters");
v21.addcontent({'parameter', 'Vpb', 'Value', 5000});
v21.addcontent({'parameter', 'Vtissue', 'Value', 210});
v21.addcontent({'parameter', 'Vtissue2', 'Value', 400});
v21.addcontent({'parameter', 'Vtissue3', 'Value', 500});
v21.addcontent({'parameter', 'Trpbref_perml', 'Value', 200000});
v21.addcontent({'parameter', 'Trpbo_perml', 'Value', 200000});
v21.addcontent({'parameter', 'Bpbref_perml', 'Value', 70000});
v21.addcontent({'parameter', 'Bpbo_perml', 'Value', 70000});
v21.addcontent({'parameter', 'KTrp', 'Value', 200});
v21.addcontent({'parameter', 'KBp', 'Value', 333});
v21.addcontent({'parameter', 'KTrp2', 'Value', 190});
v21.addcontent({'parameter', 'KBp2', 'Value', 190});
v21.addcontent({'parameter', 'KTrp3', 'Value', 60});
v21.addcontent({'parameter', 'KBp3', 'Value', 80});
v21.addcontent({'parameter', 'Cl_tdb', 'Value', 1.5});
v21.Active = false;
v21.Name = "human ALL parameters";
v21.Notes = "";
v21.Tag = "";

v22 = addvariant(m1, "Blin parameters");
v22.addcontent({'parameter', 'KdrugactT_blin', 'Value', 0.1});
v22.addcontent({'parameter', 'ndrugactT_blin', 'Value', 1});
v22.addcontent({'parameter', 'KmBT_act_blin', 'Value', 0.716});
v22.addcontent({'parameter', 'KdrugB_blin', 'Value', 0.15});
v22.addcontent({'parameter', 'nkill_blin', 'Value', 1.027});
v22.addcontent({'parameter', 'KmTB_kill_blin', 'Value', 0.75});
v22.Active = true;
v22.Name = "Blin parameters";
v22.Notes = "";
v22.Tag = "";

v23 = addvariant(m1, "human CLL parameters");
v23.addcontent({'parameter', 'Vpb', 'Value', 5000});
v23.addcontent({'parameter', 'Vtissue', 'Value', 210});
v23.addcontent({'parameter', 'Vtissue2', 'Value', 400});
v23.addcontent({'parameter', 'Vtissue3', 'Value', 500});
v23.addcontent({'parameter', 'Trpbref_perml', 'Value', 500000});
v23.addcontent({'parameter', 'Trpbo_perml', 'Value', 500000});
v23.addcontent({'parameter', 'Bpbref_perml', 'Value', 500000});
v23.addcontent({'parameter', 'Bpbo_perml', 'Value', 500000});
v23.addcontent({'parameter', 'KTrp', 'Value', 200});
v23.addcontent({'parameter', 'KBp', 'Value', 333});
v23.addcontent({'parameter', 'KTrp2', 'Value', 190});
v23.addcontent({'parameter', 'KBp2', 'Value', 190});
v23.addcontent({'parameter', 'KTrp3', 'Value', 60});
v23.addcontent({'parameter', 'KBp3', 'Value', 80});
v23.addcontent({'parameter', 'Cl_tdb', 'Value', 1.5});
v23.addcontent({'parameter', 'kBapop_cll', 'Value', 0.1215});
v23.Active = false;
v23.Name = "human CLL parameters";
v23.Notes = "";
v23.Tag = "";

v24 = addvariant(m1, "act optimizn v22_2");
v24.addcontent({'parameter', 'fTaprolif', 'Value', 2.106});
v24.addcontent({'parameter', 'kTaapop', 'Value', 0.061});
v24.addcontent({'parameter', 'fTrapop', 'Value', 0.2});
v24.addcontent({'parameter', 'fTa0apop', 'Value', 2.019});
v24.addcontent({'parameter', 'fa0', 'Value', 0.002});
v24.addcontent({'parameter', 'fAICD', 'Value', 1.497});
v24.addcontent({'parameter', 'kTact', 'Value', 9.827999999999999});
v24.addcontent({'parameter', 'fTadeact', 'Value', 0.005});
v24.addcontent({'parameter', 'fTa0deact', 'Value', 0.001});
v24.addcontent({'parameter', 'kBkill', 'Value', 275.189});
v24.addcontent({'parameter', 'kTaexit', 'Value', 0.123});
v24.addcontent({'parameter', 'KdrugactT', 'Value', 130.083});
v24.addcontent({'parameter', 'KdrugB', 'Value', 1.302});
v24.addcontent({'parameter', 'nkill', 'Value', 1.027});
v24.addcontent({'parameter', 'KmBT_act', 'Value', 0.716});
v24.addcontent({'parameter', 'fTap', 'Value', 3.773});
v24.addcontent({'parameter', 'kBapop', 'Value', 0.02});
v24.addcontent({'parameter', 'fdrug', 'Value', 6.3839});
v24.addcontent({'parameter', 'kTrexit', 'Value', 0.045});
v24.addcontent({'parameter', 'finj', 'Value', 1});
v24.addcontent({'parameter', 'kBprolif', 'Value', 0.05});
v24.addcontent({'parameter', 'fBexit', 'Value', 0.333});
v24.addcontent({'parameter', 'fBkill', 'Value', 0.1});
v24.addcontent({'parameter', 'fTact', 'Value', 0.25});
v24.addcontent({'parameter', 'fKmTB_kill', 'Value', 10});
v24.addcontent({'parameter', 'S', 'Value', 1.4});
v24.addcontent({'parameter', 'ndrugactT', 'Value', 0.8});
v24.addcontent({'parameter', 'fapop_v24', 'Value', 1});
v24.addcontent({'parameter', 'IL6_tiss_contribution', 'Value', 0.003});
v24.Active = true;
v24.Name = "act optimizn v22_2";
v24.Notes = "";
v24.Tag = "";

v25 = addvariant(m1, "human DLBCL parameters");
v25.addcontent({'parameter', 'Trpbref_perml', 'Value', 500000});
v25.addcontent({'parameter', 'Trpbo_perml', 'Value', 500000});
v25.addcontent({'parameter', 'Bpbref_perml', 'Value', 250000});
v25.addcontent({'parameter', 'Bpbo_perml', 'Value', 250000});
v25.Active = false;
v25.Name = "human DLBCL parameters";
v25.Notes = "";
v25.Tag = "";

v26 = addvariant(m1, "Cyno PK Lin (GUI Fit - all values normalized)");
v26.addcontent({'parameter', 'Cl_tdb', 'Value', 8.5});
v26.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v26.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v26.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v26.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v26.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v26.Active = false;
v26.Name = "Cyno PK Lin (GUI Fit - all values normalized)";
v26.Notes = "";
v26.Tag = "";

v27 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized)");
v27.addcontent({'parameter', 'Cl_tdb', 'Value', 5.4});
v27.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v27.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v27.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v27.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v27.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v27.addcontent({'parameter', 'kabs_TDB', 'Value', 0});
v27.addcontent({'parameter', 'fbio_TDB', 'Value', 0});
v27.Active = false;
v27.Name = "Human PK Lin (GUI Fit - all values normalized)";
v27.Notes = "";
v27.Tag = "";

v28 = addvariant(m1, "Tumor compartment - DLBCL");
v28.addcontent({'parameter', 'tumor_on', 'Value', 1});
v28.addcontent({'parameter', 'Vtumor', 'Value', 100});
v28.addcontent({'parameter', 'KTrptumor', 'Value', 100});
v28.addcontent({'parameter', 'KBptumor', 'Value', 1600});
v28.addcontent({'parameter', 'Kptumor', 'Value', 0.05});
v28.addcontent({'parameter', 'kBtumorprolif', 'Value', 0.025});
v28.addcontent({'parameter', 'IL6_tiss_contribution', 'Value', 0.0004});
v28.Active = false;
v28.Name = "Tumor compartment - DLBCL";
v28.Notes = "";
v28.Tag = "";

v29 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized) - SC - Ave Abs - Ave fbio");
v29.addcontent({'parameter', 'Cl_tdb', 'Value', 5.4});
v29.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v29.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v29.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v29.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v29.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v29.addcontent({'parameter', 'kabs_TDB', 'Value', 0.35});
v29.addcontent({'parameter', 'fbio_TDB', 'Value', 0.5});
v29.Active = false;
v29.Name = "Human PK Lin (GUI Fit - all values normalized) - SC - Ave Abs - Ave fbio";
v29.Notes = "";
v29.Tag = "";

v30 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized) - SC - Fast Abs - High fbio");
v30.addcontent({'parameter', 'Cl_tdb', 'Value', 5.4});
v30.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v30.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v30.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v30.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v30.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v30.addcontent({'parameter', 'kabs_TDB', 'Value', 0.5});
v30.addcontent({'parameter', 'fbio_TDB', 'Value', 0.8});
v30.Active = false;
v30.Name = "Human PK Lin (GUI Fit - all values normalized) - SC - Fast Abs - High fbio";
v30.Notes = "";
v30.Tag = "";

v31 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized) - SC - Slow Abs - Low fbio");
v31.addcontent({'parameter', 'Cl_tdb', 'Value', 5.4});
v31.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v31.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v31.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v31.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v31.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v31.addcontent({'parameter', 'kabs_TDB', 'Value', 0.2});
v31.addcontent({'parameter', 'fbio_TDB', 'Value', 0.2});
v31.Active = false;
v31.Name = "Human PK Lin (GUI Fit - all values normalized) - SC - Slow Abs - Low fbio";
v31.Notes = "";
v31.Tag = "";

v32 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized) - SC - Fast Abs - Ave fbio");
v32.addcontent({'parameter', 'Cl_tdb', 'Value', 5.4});
v32.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v32.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v32.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v32.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v32.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v32.addcontent({'parameter', 'kabs_TDB', 'Value', 0.5});
v32.addcontent({'parameter', 'fbio_TDB', 'Value', 0.5});
v32.Active = false;
v32.Name = "Human PK Lin (GUI Fit - all values normalized) - SC - Fast Abs - Ave fbio";
v32.Notes = "";
v32.Tag = "";

v33 = addvariant(m1, "Human PK Lin (GUI Fit - all values normalized) - SC - Slow Abs - Ave fbio");
v33.addcontent({'parameter', 'Bcell_tumor_trafficking_on', 'Value', 5.4});
v33.addcontent({'parameter', 'Cld_tdb', 'Value', 24.17});
v33.addcontent({'parameter', 'Vc_tdb', 'Value', 36.78});
v33.addcontent({'parameter', 'Vp_tdb', 'Value', 173.47});
v33.addcontent({'parameter', 'Vm_tdb', 'Value', 0});
v33.addcontent({'parameter', 'Km_tdb', 'Value', 5.72});
v33.addcontent({'parameter', 'kabs_TDB', 'Value', 0.2});
v33.addcontent({'parameter', 'fbio_TDB', 'Value', 0.5});
v33.Active = false;
v33.Name = "Human PK Lin (GUI Fit - all values normalized) - SC - Slow Abs - Ave fbio";
v33.Notes = "";
v33.Tag = "";

v34 = addvariant(m1, "Tumor B cell off");
v34.addcontent({'parameter', 'Bcell_tumor_trafficking_on', 'Value', 0});
v34.Active = false;
v34.Name = "Tumor B cell off";
v34.Notes = "";
v34.Tag = "";

v35 = addvariant(m1, "chimp - blin");
v35.addcontent({'parameter', 'IL6_tiss_contribution', 'Value', 0.0004});
v35.Active = false;
v35.Name = "chimp - blin";
v35.Notes = "";
v35.Tag = "";

%% Create doses.
dose1 = adddose(m1, "1mg/kg 4qw", "repeat");
dose1.Amount = 1000;
dose1.Interval = 7;
dose1.Rate = 0;
dose1.RepeatCount = 3;
dose1.StartTime = 0;
dose1.Active = false;
dose1.AmountUnits = "";
dose1.DurationParameterName = "";
dose1.EventMode = "stop";
dose1.LagParameterName = "";
dose1.RateUnits = "";
dose1.TargetName = "TDBc_ugperkg";
dose1.TimeUnits = "";
dose1.Name = "1mg/kg 4qw";
dose1.Notes = "";
dose1.Tag = "";

dose2 = adddose(m1, "0.1 mg/kg 4qw", "repeat");
dose2.Amount = 100;
dose2.Interval = 7;
dose2.Rate = 0;
dose2.RepeatCount = 3;
dose2.StartTime = 0;
dose2.Active = false;
dose2.AmountUnits = "";
dose2.DurationParameterName = "";
dose2.EventMode = "stop";
dose2.LagParameterName = "";
dose2.RateUnits = "";
dose2.TargetName = "TDBc_ugperkg";
dose2.TimeUnits = "";
dose2.Name = "0.1 mg/kg 4qw";
dose2.Notes = "";
dose2.Tag = "";

dose3 = adddose(m1, "1mg/kg single", "repeat");
dose3.Amount = 1000;
dose3.Interval = 365;
dose3.Rate = 0;
dose3.RepeatCount = 0;
dose3.StartTime = 0;
dose3.Active = false;
dose3.AmountUnits = "";
dose3.DurationParameterName = "";
dose3.EventMode = "stop";
dose3.LagParameterName = "";
dose3.RateUnits = "";
dose3.TargetName = "TDBc_ugperkg";
dose3.TimeUnits = "";
dose3.Name = "1mg/kg single";
dose3.Notes = "";
dose3.Tag = "";

dose4 = adddose(m1, "0.1 mg/kg single", "repeat");
dose4.Amount = 100;
dose4.Interval = 7;
dose4.Rate = 0;
dose4.RepeatCount = 0;
dose4.StartTime = 0;
dose4.Active = false;
dose4.AmountUnits = "";
dose4.DurationParameterName = "";
dose4.EventMode = "stop";
dose4.LagParameterName = "";
dose4.RateUnits = "";
dose4.TargetName = "TDBc_ugperkg";
dose4.TimeUnits = "";
dose4.Name = "0.1 mg/kg single";
dose4.Notes = "";
dose4.Tag = "";

dose5 = adddose(m1, "0.01 mg/kg 4qw", "repeat");
dose5.Amount = 10;
dose5.Interval = 7;
dose5.Rate = 0;
dose5.RepeatCount = 3;
dose5.StartTime = 0;
dose5.Active = false;
dose5.AmountUnits = "";
dose5.DurationParameterName = "";
dose5.EventMode = "stop";
dose5.LagParameterName = "";
dose5.RateUnits = "";
dose5.TargetName = "TDBc_ugperkg";
dose5.TimeUnits = "";
dose5.Name = "0.01 mg/kg 4qw";
dose5.Notes = "";
dose5.Tag = "";

dose6 = adddose(m1, "0.01 mg/kg 2qw", "repeat");
dose6.Amount = 10;
dose6.Interval = 7;
dose6.Rate = 0;
dose6.RepeatCount = 1;
dose6.StartTime = 0;
dose6.Active = false;
dose6.AmountUnits = "";
dose6.DurationParameterName = "";
dose6.EventMode = "stop";
dose6.LagParameterName = "";
dose6.RateUnits = "";
dose6.TargetName = "TDBc_ugperkg";
dose6.TimeUnits = "";
dose6.Name = "0.01 mg/kg 2qw";
dose6.Notes = "";
dose6.Tag = "";

dose7 = adddose(m1, "4qw inj", "repeat");
dose7.Amount = 1;
dose7.Interval = 7;
dose7.Rate = 0;
dose7.RepeatCount = 3;
dose7.StartTime = 0;
dose7.Active = false;
dose7.AmountUnits = "";
dose7.DurationParameterName = "";
dose7.EventMode = "stop";
dose7.LagParameterName = "";
dose7.RateUnits = "";
dose7.TargetName = "injection_effect";
dose7.TimeUnits = "";
dose7.Name = "4qw inj";
dose7.Notes = "";
dose7.Tag = "";

dose8 = adddose(m1, "single inj", "schedule");
dose8.Amount = 1;
dose8.Rate = 0;
dose8.Time = 0;
dose8.Active = true;
dose8.AmountUnits = "";
dose8.DurationParameterName = "";
dose8.EventMode = "stop";
dose8.LagParameterName = "";
dose8.RateUnits = "";
dose8.TargetName = "injection_effect";
dose8.TimeUnits = "";
dose8.Name = "single inj";
dose8.Notes = "";
dose8.Tag = "";

dose9 = adddose(m1, "2qw inj", "schedule");
dose9.Amount = [1;1];
dose9.Rate = [0;0];
dose9.Time = [0;7];
dose9.Active = false;
dose9.AmountUnits = "";
dose9.DurationParameterName = "";
dose9.EventMode = "stop";
dose9.LagParameterName = "";
dose9.RateUnits = "";
dose9.TargetName = "injection_effect";
dose9.TimeUnits = "";
dose9.Name = "2qw inj";
dose9.Notes = "";
dose9.Tag = "";

dose10 = adddose(m1, "0.01 mg/kg single", "repeat");
dose10.Amount = 10;
dose10.Interval = 7;
dose10.Rate = 0;
dose10.RepeatCount = 0;
dose10.StartTime = 0;
dose10.Active = false;
dose10.AmountUnits = "";
dose10.DurationParameterName = "";
dose10.EventMode = "stop";
dose10.LagParameterName = "";
dose10.RateUnits = "";
dose10.TargetName = "TDBc_ugperkg";
dose10.TimeUnits = "";
dose10.Name = "0.01 mg/kg single";
dose10.Notes = "";
dose10.Tag = "";

dose11 = adddose(m1, "0.001 mg/kg single", "repeat");
dose11.Amount = 1;
dose11.Interval = 7;
dose11.Rate = 0;
dose11.RepeatCount = 0;
dose11.StartTime = 0;
dose11.Active = false;
dose11.AmountUnits = "";
dose11.DurationParameterName = "";
dose11.EventMode = "stop";
dose11.LagParameterName = "";
dose11.RateUnits = "";
dose11.TargetName = "TDBc_ugperkg";
dose11.TimeUnits = "";
dose11.Name = "0.001 mg/kg single";
dose11.Notes = "";
dose11.Tag = "";

dose12 = adddose(m1, "0.003 mg/kg single", "repeat");
dose12.Amount = 3;
dose12.Interval = 7;
dose12.Rate = 0;
dose12.RepeatCount = 0;
dose12.StartTime = 0;
dose12.Active = false;
dose12.AmountUnits = "";
dose12.DurationParameterName = "";
dose12.EventMode = "stop";
dose12.LagParameterName = "";
dose12.RateUnits = "";
dose12.TargetName = "TDBc_ugperkg";
dose12.TimeUnits = "";
dose12.Name = "0.003 mg/kg single";
dose12.Notes = "";
dose12.Tag = "";

dose13 = adddose(m1, "RTX 10 mg/kg 2qw", "repeat");
dose13.Amount = 10000;
dose13.Interval = 7;
dose13.Rate = 0;
dose13.RepeatCount = 1;
dose13.StartTime = 0;
dose13.Active = false;
dose13.AmountUnits = "";
dose13.DurationParameterName = "";
dose13.EventMode = "stop";
dose13.LagParameterName = "";
dose13.RateUnits = "";
dose13.TargetName = "RTXc_ugperkg";
dose13.TimeUnits = "";
dose13.Name = "RTX 10 mg/kg 2qw";
dose13.Notes = "";
dose13.Tag = "";

dose14 = adddose(m1, "RTX 10 mg/kg single dose", "repeat");
dose14.Amount = 10000;
dose14.Interval = 7;
dose14.Rate = 0;
dose14.RepeatCount = 0;
dose14.StartTime = 0;
dose14.Active = false;
dose14.AmountUnits = "";
dose14.DurationParameterName = "";
dose14.EventMode = "stop";
dose14.LagParameterName = "";
dose14.RateUnits = "";
dose14.TargetName = "RTXc_ugperkg";
dose14.TimeUnits = "";
dose14.Name = "RTX 10 mg/kg single dose";
dose14.Notes = "";
dose14.Tag = "";

dose15 = adddose(m1, "RTX 0.2 mg/kg 2qw", "repeat");
dose15.Amount = 200;
dose15.Interval = 7;
dose15.Rate = 0;
dose15.RepeatCount = 1;
dose15.StartTime = 0;
dose15.Active = false;
dose15.AmountUnits = "";
dose15.DurationParameterName = "";
dose15.EventMode = "stop";
dose15.LagParameterName = "";
dose15.RateUnits = "";
dose15.TargetName = "RTXc_ugperkg";
dose15.TimeUnits = "";
dose15.Name = "RTX 0.2 mg/kg 2qw";
dose15.Notes = "";
dose15.Tag = "";

dose16 = adddose(m1, "RTX 2 mg/kg 2qw", "repeat");
dose16.Amount = 2000;
dose16.Interval = 7;
dose16.Rate = 0;
dose16.RepeatCount = 1;
dose16.StartTime = 0;
dose16.Active = false;
dose16.AmountUnits = "";
dose16.DurationParameterName = "";
dose16.EventMode = "stop";
dose16.LagParameterName = "";
dose16.RateUnits = "";
dose16.TargetName = "RTXc_ugperkg";
dose16.TimeUnits = "";
dose16.Name = "RTX 2 mg/kg 2qw";
dose16.Notes = "";
dose16.Tag = "";

dose17 = adddose(m1, "RTX 0.05 mg/kg 2qw", "repeat");
dose17.Amount = 50;
dose17.Interval = 7;
dose17.Rate = 0;
dose17.RepeatCount = 1;
dose17.StartTime = 0;
dose17.Active = false;
dose17.AmountUnits = "";
dose17.DurationParameterName = "";
dose17.EventMode = "stop";
dose17.LagParameterName = "";
dose17.RateUnits = "";
dose17.TargetName = "RTXc_ugperkg";
dose17.TimeUnits = "";
dose17.Name = "RTX 0.05 mg/kg 2qw";
dose17.Notes = "";
dose17.Tag = "";

dose18 = adddose(m1, "4qw inj drug", "repeat");
dose18.Amount = 1;
dose18.Interval = 7;
dose18.Rate = 0;
dose18.RepeatCount = 3;
dose18.StartTime = 0;
dose18.Active = false;
dose18.AmountUnits = "";
dose18.DurationParameterName = "";
dose18.EventMode = "stop";
dose18.LagParameterName = "";
dose18.RateUnits = "";
dose18.TargetName = "drug_effect";
dose18.TimeUnits = "";
dose18.Name = "4qw inj drug";
dose18.Notes = "";
dose18.Tag = "";

dose19 = adddose(m1, "single inj drug", "schedule");
dose19.Amount = 1;
dose19.Rate = 0;
dose19.Time = 0;
dose19.Active = true;
dose19.AmountUnits = "";
dose19.DurationParameterName = "";
dose19.EventMode = "stop";
dose19.LagParameterName = "";
dose19.RateUnits = "";
dose19.TargetName = "drug_effect";
dose19.TimeUnits = "";
dose19.Name = "single inj drug";
dose19.Notes = "";
dose19.Tag = "";

dose20 = adddose(m1, "Blin 15 IV 4w", "schedule");
dose20.Amount = 420;
dose20.Rate = 15;
dose20.Time = 0;
dose20.Active = true;
dose20.AmountUnits = "";
dose20.DurationParameterName = "";
dose20.EventMode = "stop";
dose20.LagParameterName = "";
dose20.RateUnits = "";
dose20.TargetName = "Blinc_ug";
dose20.TimeUnits = "";
dose20.Name = "Blin 15 IV 4w";
dose20.Notes = "";
dose20.Tag = "";

dose21 = adddose(m1, "Blin 5-15-60 IV 3w", "schedule");
dose21.Amount = [35;105;420];
dose21.Rate = [5;15;60];
dose21.Time = [0;7;14];
dose21.Active = false;
dose21.AmountUnits = "";
dose21.DurationParameterName = "";
dose21.EventMode = "stop";
dose21.LagParameterName = "";
dose21.RateUnits = "";
dose21.TargetName = "Blinc_ug";
dose21.TimeUnits = "";
dose21.Name = "Blin 5-15-60 IV 3w";
dose21.Notes = "";
dose21.Tag = "";

dose22 = adddose(m1, "Blin 60 IV 3w", "schedule");
dose22.Amount = 1260;
dose22.Rate = 60;
dose22.Time = 0;
dose22.Active = false;
dose22.AmountUnits = "";
dose22.DurationParameterName = "";
dose22.EventMode = "stop";
dose22.LagParameterName = "";
dose22.RateUnits = "";
dose22.TargetName = "Blinc_ug";
dose22.TimeUnits = "";
dose22.Name = "Blin 60 IV 3w";
dose22.Notes = "";
dose22.Tag = "";

dose23 = adddose(m1, "8q3w inj", "schedule");
dose23.Amount = [1;1;1;1;1;1;1;1];
dose23.Rate = [0;0;0;0;0;0;0;0];
dose23.Time = [0;21;42;63;84;105;126;147];
dose23.Active = false;
dose23.AmountUnits = "";
dose23.DurationParameterName = "";
dose23.EventMode = "stop";
dose23.LagParameterName = "";
dose23.RateUnits = "";
dose23.TargetName = "injection_effect";
dose23.TimeUnits = "";
dose23.Name = "8q3w inj";
dose23.Notes = "";
dose23.Tag = "";

dose24 = adddose(m1, "8q3w inj drug", "schedule");
dose24.Amount = [1;1;1;1;1;1;1;1];
dose24.Rate = [0;0;0;0;0;0;0;0];
dose24.Time = [0;21;42;63;84;105;126;147];
dose24.Active = false;
dose24.AmountUnits = "";
dose24.DurationParameterName = "";
dose24.EventMode = "stop";
dose24.LagParameterName = "";
dose24.RateUnits = "";
dose24.TargetName = "drug_effect";
dose24.TimeUnits = "";
dose24.Name = "8q3w inj drug";
dose24.Notes = "";
dose24.Tag = "";

dose25 = adddose(m1, "TDB 0.05 mg 8q3w", "schedule");
dose25.Amount = [0.714;0.714;0.714;0.714;0.714;0.714;0.714;0.714];
dose25.Rate = [4.284;4.284;4.284;4.284;4.284;4.284;4.284;4.284];
dose25.Time = [0;21;42;63;84;105;126;147];
dose25.Active = false;
dose25.AmountUnits = "";
dose25.DurationParameterName = "";
dose25.EventMode = "stop";
dose25.LagParameterName = "";
dose25.RateUnits = "";
dose25.TargetName = "TDBc_ugperkg";
dose25.TimeUnits = "";
dose25.Name = "TDB 0.05 mg 8q3w";
dose25.Notes = "";
dose25.Tag = "";

dose26 = adddose(m1, "TDB 0.2 mg 8q3w", "schedule");
dose26.Amount = [2.857;2.857;2.857;2.857;2.857;2.857;2.857;2.857];
dose26.Rate = [17.142;17.142;17.142;17.142;17.142;17.142;17.142;17.142];
dose26.Time = [0;21;42;63;84;105;126;147];
dose26.Active = false;
dose26.AmountUnits = "";
dose26.DurationParameterName = "";
dose26.EventMode = "stop";
dose26.LagParameterName = "";
dose26.RateUnits = "";
dose26.TargetName = "TDBc_ugperkg";
dose26.TimeUnits = "";
dose26.Name = "TDB 0.2 mg 8q3w";
dose26.Notes = "";
dose26.Tag = "";

dose27 = adddose(m1, "TDB 0.8 mg 8q3w", "schedule");
dose27.Amount = [11.43;11.43;11.43;11.43;11.43;11.43;11.43;11.43];
dose27.Rate = [68.58;68.58;68.58;68.58;68.58;68.58;68.58;68.58];
dose27.Time = [0;21;42;63;84;105;126;147];
dose27.Active = false;
dose27.AmountUnits = "";
dose27.DurationParameterName = "";
dose27.EventMode = "stop";
dose27.LagParameterName = "";
dose27.RateUnits = "";
dose27.TargetName = "TDBc_ugperkg";
dose27.TimeUnits = "";
dose27.Name = "TDB 0.8 mg 8q3w";
dose27.Notes = "";
dose27.Tag = "";

dose28 = adddose(m1, "TDB 1.6 mg 8q3w", "schedule");
dose28.Amount = [22.86;22.86;22.86;22.86;22.86;22.86;22.86;22.86];
dose28.Rate = [137.16;137.16;137.16;137.16;137.16;137.16;137.16;137.16];
dose28.Time = [0;21;42;63;84;105;126;147];
dose28.Active = false;
dose28.AmountUnits = "";
dose28.DurationParameterName = "";
dose28.EventMode = "stop";
dose28.LagParameterName = "";
dose28.RateUnits = "";
dose28.TargetName = "TDBc_ugperkg";
dose28.TimeUnits = "";
dose28.Name = "TDB 1.6 mg 8q3w";
dose28.Notes = "";
dose28.Tag = "";

dose29 = adddose(m1, "TDB 3.2 mg 8q3w", "schedule");
dose29.Amount = [45.71;45.71;45.71;45.71;45.71;45.71;45.71;45.71];
dose29.Rate = [274.26;274.26;274.26;274.26;274.26;274.26;274.26;274.26];
dose29.Time = [0;21;42;63;84;105;126;147];
dose29.Active = false;
dose29.AmountUnits = "";
dose29.DurationParameterName = "";
dose29.EventMode = "stop";
dose29.LagParameterName = "";
dose29.RateUnits = "";
dose29.TargetName = "TDBc_ugperkg";
dose29.TimeUnits = "";
dose29.Name = "TDB 3.2 mg 8q3w";
dose29.Notes = "";
dose29.Tag = "";

dose30 = adddose(m1, "TDB 6.4 mg 8q3w", "schedule");
dose30.Amount = [91.43000000000001;91.43000000000001;91.43000000000001;91.43000000000001;91.43000000000001;91.43000000000001;91.43000000000001;91.43000000000001];
dose30.Rate = [548.58;548.58;548.58;548.58;548.58;548.58;548.58;548.58];
dose30.Time = [0;21;42;63;84;105;126;147];
dose30.Active = false;
dose30.AmountUnits = "";
dose30.DurationParameterName = "";
dose30.EventMode = "stop";
dose30.LagParameterName = "";
dose30.RateUnits = "";
dose30.TargetName = "TDBc_ugperkg";
dose30.TimeUnits = "";
dose30.Name = "TDB 6.4 mg 8q3w";
dose30.Notes = "";
dose30.Tag = "";

dose31 = adddose(m1, "TDB 12.8 mg 8q3w", "schedule");
dose31.Amount = [182.86;182.86;182.86;182.86;182.86;182.86;182.86;182.86];
dose31.Rate = [1097.16;1097.16;1097.16;1097.16;1097.16;1097.16;1097.16;1097.16];
dose31.Time = [0;21;42;63;84;105;126;147];
dose31.Active = false;
dose31.AmountUnits = "";
dose31.DurationParameterName = "";
dose31.EventMode = "stop";
dose31.LagParameterName = "";
dose31.RateUnits = "";
dose31.TargetName = "TDBc_ugperkg";
dose31.TimeUnits = "";
dose31.Name = "TDB 12.8 mg 8q3w";
dose31.Notes = "";
dose31.Tag = "";

dose32 = adddose(m1, "TDB 20 mg 8q3w", "schedule");
dose32.Amount = [285.71;285.71;285.71;285.71;285.71;285.71;285.71;285.71];
dose32.Rate = [1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26];
dose32.Time = [0;21;42;63;84;105;126;147];
dose32.Active = false;
dose32.AmountUnits = "";
dose32.DurationParameterName = "";
dose32.EventMode = "stop";
dose32.LagParameterName = "";
dose32.RateUnits = "";
dose32.TargetName = "TDBc_ugperkg";
dose32.TimeUnits = "";
dose32.Name = "TDB 20 mg 8q3w";
dose32.Notes = "";
dose32.Tag = "";

dose33 = adddose(m1, "TDB 30 mg 8q3w", "schedule");
dose33.Amount = [428.57;428.57;428.57;428.57;428.57;428.57;428.57;428.57];
dose33.Rate = [2571.42;2571.42;2571.42;2571.42;2571.42;2571.42;2571.42;2571.42];
dose33.Time = [0;21;42;63;84;105;126;147];
dose33.Active = false;
dose33.AmountUnits = "";
dose33.DurationParameterName = "";
dose33.EventMode = "stop";
dose33.LagParameterName = "";
dose33.RateUnits = "";
dose33.TargetName = "TDBc_ugperkg";
dose33.TimeUnits = "";
dose33.Name = "TDB 30 mg 8q3w";
dose33.Notes = "";
dose33.Tag = "";

dose34 = adddose(m1, "3qw + 7q3w inj", "schedule");
dose34.Amount = [1;1;1;1;1;1;1;1;1;1];
dose34.Rate = [0;0;0;0;0;0;0;0;0;0];
dose34.Time = [0;7;14;21;42;63;84;105;126;147];
dose34.Active = false;
dose34.AmountUnits = "";
dose34.DurationParameterName = "";
dose34.EventMode = "stop";
dose34.LagParameterName = "";
dose34.RateUnits = "";
dose34.TargetName = "injection_effect";
dose34.TimeUnits = "";
dose34.Name = "3qw + 7q3w inj";
dose34.Notes = "";
dose34.Tag = "";

dose35 = adddose(m1, "3qw + 7q3w inj drug", "schedule");
dose35.Amount = [1;1;1;1;1;1;1;1;1;1];
dose35.Rate = [0;0;0;0;0;0;0;0;0;0];
dose35.Time = [0;7;14;21;42;63;84;105;126;147];
dose35.Active = false;
dose35.AmountUnits = "";
dose35.DurationParameterName = "";
dose35.EventMode = "stop";
dose35.LagParameterName = "";
dose35.RateUnits = "";
dose35.TargetName = "drug_effect";
dose35.TimeUnits = "";
dose35.Name = "3qw + 7q3w inj drug";
dose35.Notes = "";
dose35.Tag = "";

dose36 = adddose(m1, "TDB 20 mg 3qw + 7q3w", "schedule");
dose36.Amount = [22.86;142.86;142.86;285.71;285.71;285.71;285.71;285.71;285.71;285.71];
dose36.Rate = [137.16;857.16;857.16;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26];
dose36.Time = [0;7;14;21;42;63;84;105;126;147];
dose36.Active = false;
dose36.AmountUnits = "";
dose36.DurationParameterName = "";
dose36.EventMode = "stop";
dose36.LagParameterName = "";
dose36.RateUnits = "";
dose36.TargetName = "TDBc_ugperkg";
dose36.TimeUnits = "";
dose36.Name = "TDB 20 mg 3qw + 7q3w";
dose36.Notes = "";
dose36.Tag = "";

dose37 = adddose(m1, "2qw + 7q3w inj", "schedule");
dose37.Amount = [1;1;1;1;1;1;1;1;1];
dose37.Rate = [0;0;0;0;0;0;0;0;0];
dose37.Time = [0;7;21;42;63;84;105;126;147];
dose37.Active = false;
dose37.AmountUnits = "";
dose37.DurationParameterName = "";
dose37.EventMode = "stop";
dose37.LagParameterName = "";
dose37.RateUnits = "";
dose37.TargetName = "injection_effect";
dose37.TimeUnits = "";
dose37.Name = "2qw + 7q3w inj";
dose37.Notes = "";
dose37.Tag = "";

dose38 = adddose(m1, "2qw + 7q3w inj drug", "schedule");
dose38.Amount = [1;1;1;1;1;1;1;1;1];
dose38.Rate = [0;0;0;0;0;0;0;0;0];
dose38.Time = [0;7;21;42;63;84;105;126;147];
dose38.Active = false;
dose38.AmountUnits = "";
dose38.DurationParameterName = "";
dose38.EventMode = "stop";
dose38.LagParameterName = "";
dose38.RateUnits = "";
dose38.TargetName = "drug_effect";
dose38.TimeUnits = "";
dose38.Name = "2qw + 7q3w inj drug";
dose38.Notes = "";
dose38.Tag = "";

dose39 = adddose(m1, "TDB 20 mg 2qw + 7q3w", "schedule");
dose39.Amount = [22.86;285.71;285.71;285.71;285.71;285.71;285.71;285.71;285.71];
dose39.Rate = [137.16;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26];
dose39.Time = [0;7;21;42;63;84;105;126;147];
dose39.Active = false;
dose39.AmountUnits = "";
dose39.DurationParameterName = "";
dose39.EventMode = "stop";
dose39.LagParameterName = "";
dose39.RateUnits = "";
dose39.TargetName = "TDBc_ugperkg";
dose39.TimeUnits = "";
dose39.Name = "TDB 20 mg 2qw + 7q3w";
dose39.Notes = "";
dose39.Tag = "";

dose40 = adddose(m1, "TDB 20 mg 3qw + 7q3w - flat frac", "schedule");
dose40.Amount = [95.23999999999999;95.23999999999999;95.23999999999999;285.71;285.71;285.71;285.71;285.71;285.71;285.71];
dose40.Rate = [571.4400000000001;571.4400000000001;571.4400000000001;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26;1714.26];
dose40.Time = [0;7;14;21;42;63;84;105;126;147];
dose40.Active = false;
dose40.AmountUnits = "";
dose40.DurationParameterName = "";
dose40.EventMode = "stop";
dose40.LagParameterName = "";
dose40.RateUnits = "";
dose40.TargetName = "TDBc_ugperkg";
dose40.TimeUnits = "";
dose40.Name = "TDB 20 mg 3qw + 7q3w - flat frac";
dose40.Notes = "";
dose40.Tag = "";

dose41 = adddose(m1, "Blin 5-15-60 IV 8w", "schedule");
dose41.Amount = [35;105;2520];
dose41.Rate = [5;15;60];
dose41.Time = [0;7;14];
dose41.Active = false;
dose41.AmountUnits = "";
dose41.DurationParameterName = "";
dose41.EventMode = "stop";
dose41.LagParameterName = "";
dose41.RateUnits = "";
dose41.TargetName = "Blinc_ug";
dose41.TimeUnits = "";
dose41.Name = "Blin 5-15-60 IV 8w";
dose41.Notes = "";
dose41.Tag = "";

end
