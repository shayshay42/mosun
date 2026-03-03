repo_root='/Users/shayanhajhashemi/genentech_tce_vpop_translation';
cd(fullfile(repo_root,'assets','Supp Matlab Code'));
s=sbioloadproject('TDBr26_6_paper.sbproj'); c=struct2cell(s); model=c{1};
for ii=1:length(model.variant), model.variant(ii).Active=false; end
vid=[5,9,14,20,24,25,27,28];
for ii=1:length(vid), model.variant(vid(ii)).Active=true; end
for ii=vid
  fprintf('idx=%d name=%s active=%d ncontent=%d\n',ii,model.variant(ii).Name,model.variant(ii).Active,size(model.variant(ii).Content,1));
end
check={'tissue2on','tissue3on','fapop_v24','fBexit','Km_tdb','Vm_tdb'};
for i=1:length(check)
  p=sbioselect(model,'Type','parameter','Name',check{i});
  fprintf('%s value=%g\n',check{i},p.Value);
end
