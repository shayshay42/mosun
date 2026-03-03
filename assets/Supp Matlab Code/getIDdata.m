function [Tdata,Xdata]=getIDdata(data,ID,datavec)

pos = (data.ID == ID);
Tdata = data.Time(pos);
Xdata = [];
for i=1:length(datavec);
    Xdata=[Xdata data.(sprintf('%s', datavec{i}))(pos)];
end

end

