function [ Out ] = xlsreadstruct( filename )

[A,B]=xlsread(filename);

header_names=B(1,:);
fs=length(header_names);
ss=length(A(1,:));

for i=1:(fs-ss)
    Out.(sprintf('%s', header_names{i})) = B(2:end,i);
end

for i=(fs-ss+1):(fs)
    Out.(sprintf('%s', header_names{i})) = A(1:end,i-fs+ss);
end

end

