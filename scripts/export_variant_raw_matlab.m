repo_root='/Users/shayanhajhashemi/genentech_tce_vpop_translation';
cd(fullfile(repo_root,'assets','Supp Matlab Code'));
s=sbioloadproject('TDBr26_6_paper.sbproj'); c=struct2cell(s); model=c{1};
rows={};k=1;
for i=1:length(model.variant)
  C=model.variant(i).Content;
  for j=1:size(C,1)
    e=C{j};
    rows(k,:)={i,string(model.variant(i).Name),string(e{1}),string(e{2}),string(e{3}),string(e{4})};
    k=k+1;
  end
end
T=cell2table(rows,'VariableNames',{'idx','variant_name','f1','f2','f3','f4'});
out=fullfile(repo_root,'generated','model_tables','variants_raw_from_matlab.tsv');
writetable(T,out,'FileType','text','Delimiter','\t');
fprintf('Wrote %s\n',out);
