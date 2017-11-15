function [indClust,p2p] = aggregate_cluster( cluster_matrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[npat,nsample] = size(cluster_matrix);
if nsample == 1
    indClust =  cluster_matrix;
    return
end
p2p = zeros(npat,npat);
for i=1:npat
%     if mod(i,100)==0
%         fprintf('finished %f\n',i/npat);
%     end
    for j=1:npat
        diff_vec = cluster_matrix(i,:) - cluster_matrix(j,:);
        p2p(i,j) = p2p(i,j) + length(find(diff_vec==0));
%         for m=1:nsample
%             if cluster_matrix(i,m) == cluster_matrix(j,m)
%                 p2p(i,j) = p2p(i,j) + 1;
%             end
%         end
    end
end
nclst = max(cluster_matrix(:));
dist = squareform(pdist(p2p));
dist = full(dist);
[indClust,~,sumD] = kmeans(full(dist),nclst,'replicates',100);
% [indClust,~,sumD] = kmeans(full(dist),nclst,'distance','cosine','replicates',100);
%       Z = linkage(dist,'average');
%                                 indClust = cluster(Z,'maxclust',nclst);
end

