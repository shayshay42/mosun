repo_root = '/Users/shayanhajhashemi/genentech_tce_vpop_translation';
cd(fullfile(repo_root,'assets','Supp Matlab Code'));
s = sbioloadproject('TDBr26_6_paper.sbproj');
c = struct2cell(s); model = c{1};
for ii=1:length(model.variant); model.variant(ii).Active=false; end
vid=[5,9,14,20,24,25,27,28];
for ii=1:length(vid); model.variant(vid(ii)).Active=true; end
cs = getconfigset(model,'active');
cs.StopTime = 2;
d = sbiodose('d1','schedule');
d.TargetName='TDBc_ugperkg';
d.Time=[0 1];
d.Amount=[10 10];
d.Rate=[0 0];
simData = sbiosimulate(model, cs, [], d);
[T,X]=selectbyname(simData, {'TDBc_ugperkg','IL6combo'});
disp([length(T), size(X,1), size(X,2)]);
disp([T(1), T(end)]);
