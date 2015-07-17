function [n1,n2,partition_dist] = partition_dist(L1,L2)
% PARTITION_DIST computes the partition distance between the two
% partitions.
%
% The partition distance is the number of pixels that must be deleted from
% L1 and L2 such that L1 exactly matching L2. If the number of non-zero
% pixels is the same in L1 and L2 (a typical partition on the same data),
% then the partition distance is: n - matching_cost, where n is the number
% of elements in the dataset.
%
% In our case, however, the number of non-zero pixels can be different in
% L1 and L2. Instead, we return: n1 + n2 - 2*matching_cost, where n1 and n2
% are the number of non-zero pixels in L1 and L2, respectively.
%
% From: Partition-Distance: A Problem and Class of Perfect Graphs Arising
%                           in Clustering
%       Dan Gusfield
%       Information Processing Letters, May 2002
%
% Example:
%
%   L1 = 1 1 0 0 2      L2 = 1 1 0 0 2      L3 = 1 1 1 1 1
%        0 0 0 2 2           1 0 0 2 2           1 0 0 1 1
%        3 3 0 0 0           0 0 0 0 0           0 0 0 0 0
%        3 3 0 4 0           3 3 0 0 4           2 2 0 0 3
%        3 3 0 4 4           3 3 0 4 4           2 2 0 3 3
%
%  >> partition_dist(L1,L2) = 5
%  >> partition_dist(L1,L3) = 11
%
% This distance is symmetric; the lower the better.
% Interpretation: partition_dist / (n1+n2) is the percentage of total nodes
% that have to be deleted.
%

M = compute_m(L1,L2,4)


% we want min-matching, so transform the problem.
max_val = max(max(M));
Mt = abs(M-max_val);


%% solve the matching problem in Mt.
[rowsol,cost,v,u,costMat] = lapjv(Mt);

% rowsol contains the indices of the matched partners. go through and
% comput the cost wrt M.
cost = 0;
for i=1:size(rowsol,2)
    if rowsol(i) == 0
        continue
    else
        % rowsol is wrt the larger of the two assignments halves.
        if size(Mt,1) > size(Mt,2)
            cost = cost + M(rowsol(i),i);
        else
            cost = cost + M(i,rowsol(i));
        end
    end
end

% compute the partition distance.
n1 = size(find(L1>0),1);
n2 = size(find(L2>0),1);
partition_dist = n1 + n2 - 2*cost;

end

function M = compute_m(L1,L2,connsize)
%COMPUTE_M computes the matching matrix, M.
%
% M(i,j) is the intersection of region i from L1 and region j from L2.
%

[L1,num_regions1] = bwlabel(L1,4);
[L2,num_regions2] = bwlabel(L2,connsize);

Region1 = regionprops(L1,'PixelIdxList');
Region2 = regionprops(L2,'PixelIdxList');

M = zeros(num_regions1,num_regions2);

fprintf('\t\t# of regions: L1=%i, L2=%i\n',num_regions1, num_regions2)


%% compute the M matrix.
for i=1:num_regions1
    for j=1:num_regions2
        M(i,j) = size(intersect_sorted(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
        %M(i,j) = size(intersect(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
    end
end

end
