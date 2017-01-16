addpath ../../raw_data
data = importdata('../../raw_data/ULIRNAseq.count.on.gene.csv');
data.collabels = data.textdata(1,2:end);
data.rowlabels = data.textdata(2:end,1);

s = strmatch('__',data.rowlabels);
data.data(s,:) = [];
data.rowlabels(s) = [];

for ii = 1:length(data.collabels);
    data.collabels{ii}(1) = [];
    data.collabels{ii}(end) = [];
end
data.sample = data.collabels;
im2 = textscan(fopen('rnaseq_sample_names.txt'),'%s%s%s%s','delimiter','\t');
fx = lookupNew(data.collabels,im2{1});
data.collabels = im2{2}(fx);
s = strmatch('passed',im2{4}(fx),'exact');
data.data = data.data(:,s);
data.collabels = data.collabels(s);
data.sample = data.sample(s);

% remove lowly expressed genes
[r,c] = strtok(data.collabels,'#');
u = unique(r);
for ii = 1:length(u);
    s = strmatch(u{ii},r,'exact');
    mi(:,ii) = min(data.data(:,s),[],2);
end

s = sum(mi>4,2);
f = find(s<2);
data.data(f,:) = [];
data.rowlabels(f) = [];

% keep samples with >500000 reads

f = find(sum(data.data)>500000);
data2.data = data.data(:,f);
data2.collabels = data.collabels(f);
data2.rowlabels = data.rowlabels;


Q = quantilenorm(log2(2+data2.data));




[r,c] = strtok(data2.collabels,'#');



u = unique(r);
Rn.data = zeros(length(data2.rowlabels),length(u));
for ii = 1:length(u);
    s = strmatch(u{ii},r,'exact');
    Rn.data(:,ii) = mean(Q(:,s),2);
end
Rn.collabels = u;
Rn.rowlabels = data2.rowlabels;





