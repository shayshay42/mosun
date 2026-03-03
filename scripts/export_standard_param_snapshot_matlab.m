repo_root = '/Users/shayanhajhashemi/genentech_tce_vpop_translation';
assets_dir = fullfile(repo_root,'assets');
matlab_dir = fullfile(assets_dir,'Supp Matlab Code');
cd(matlab_dir);
s = sbioloadproject('TDBr26_6_paper.sbproj'); c = struct2cell(s); model = c{1};
for ii=1:length(model.variant), model.variant(ii).Active=false; end
vid=[5,9,14,20,24,25,27,28];
for ii=1:length(vid), model.variant(vid(ii)).Active=true; end
pt = readtable(fullfile(assets_dir,'params_41540_2020_145_MOESM2_ESM.xlsx'),'Sheet','Sheet1');
for ii=1:height(pt)
  pname=string(pt.NAME(ii)); pval=pt.DLBCL(ii);
  if ~isnan(pval)
    pobj=sbioselect(model,'Type','parameter','Name',char(pname));
    if ~isempty(pobj), pobj(1).Value=pval; end
  end
end
% patient overrides
ovr = {'Bpbo_perml',250000;'Bpbref_perml',250000;'Trpbo_perml',500000;'Trpbref_perml',500000;'KBptumor',6500;'KTrptumor',125;'kBtumorprolif',0.025;'PKflag',1;'fvalidation',0;'VPid',1;'end_time',84};
for ii=1:size(ovr,1)
  pobj=sbioselect(model,'Type','parameter','Name',ovr{ii,1});
  if ~isempty(pobj), pobj(1).Value=ovr{ii,2}; end
end
n=length(model.Parameters);
name = strings(n,1); value = zeros(n,1);
for i=1:n
  name(i)=string(model.Parameters(i).Name);
  value(i)=model.Parameters(i).Value;
end
T=table(name,value);
out=fullfile(repo_root,'generated','phase1_standard_matlab','param_snapshot_matlab.csv');
writetable(T,out);
fprintf('Wrote %s\n',out);
