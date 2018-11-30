function [ hpv_status ] = read_TCGA_HPV( file,p2i,indClust)
%READ_TCGA_HPV Summary of this function goes here
%   Detailed explanation goes here
[pat,hpv] = textread(file,'%s%d');
filter = isKey(p2i,pat);

pat = pat(filter);
hpv = hpv(filter);

pid = cell2mat(values(p2i,pat));
npat = length(p2i);
hpv_status = zeros(1,npat)-1;
hpv_status(pid) = hpv;
if exist('indClust','var')
nclst = max(indClust);

for c=1:nclst
    px = zeros(1,npat);
    px(indClust==c) = 1;
    for hpv_key = [0,1]
        gx = zeros(1,npat);
        gx((hpv_status==hpv_key)) = 1;
         ct = crosstab(px,gx);
        [~,pv,~] = fishertest(ct,'Tail','right');
        ratio_pos = length(intersect(find(gx==1),find(px==1)))/length(find(px==1));
        ratio_neg = length(intersect(find(gx==0),find(px==1)))/length(find(px==1));
%         fprintf('clust:%d hpv:%d %f %f pv %f\n',c,hpv_key,ratio_pos,ratio_neg,pv);
    end
    
end
end
end

