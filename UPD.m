%========================================================================%
%    UNDER-SEGMENTATION PARTITION DISTANCE                               %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : June 15, 2011                                               %
%                                                                        %
%========================================================================%


function dist = UPD(L1,L2)
% UPD computes the under-segmentation partition distance between two partitions.
%
% The under-segmentation error is the total amount of "bleeding" caused by
% superpixels that overlap each ground truth segment, normalized by the 
% segment's area.
%
% Example:
%
%  >> L1 = [1 1 0 0 2; 0 0 0 2 2; 3 3 0 0 0; 3 3 0 4 0; 3 3 0 4 4];
%  >> L2 = [1 1 0 0 2; 1 0 0 2 2; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];
%  >> L3 = [1 1 1 1 1; 1 0 0 1 1; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];
%
%   L1 = 1 1 0 0 2      L2 = 1 1 0 0 2      L3 = 1 1 1 1 1
%        0 0 0 2 2           1 0 0 2 2           1 0 0 1 1
%        3 3 0 0 0           0 0 0 0 0           0 0 0 0 0
%        3 3 0 4 0           3 3 0 0 4           2 2 0 0 3
%        3 3 0 4 4           3 3 0 4 4           2 2 0 3 3
%
%  >> UPD(L1,L2) = ??? [fill in later]
%  >> UPD(L2,L1) = ???
%  >> UPD(L1,L3) = ???
%  >> UPD(L3,L1) = ???
%  >> UPD(L2,L3) = ???
%  >> UPD(L3,L2) = ???
%
% In our case, L1 is the ground truth and L2 is the algorithm's partition.
%


dist = Turbo_UPD(L1,L2);
%dist = SLIC_UPD(L1,L2);

end



function dist_turbo = Turbo_UPD(L1,L2)
%TURBO_UPD computes under segmentation error as defined by the TurboPixels paper.
%
% U = 1/M [ \sum_{i=1}^M (([ \sum_{sj | sj \cup gi > 0} |sj| ] - |gi|) / |gi|) ],
%       where gi = the ith groundtruth segment (there are M total)
%             sj = the jth superpixel
%
% From: TurboPixels: Fast Superpixels Using Geometric Flows.
%       PAMI, 2009.

[L1,num_regions1] = bwlabel(L1,4);
[L2,num_regions2] = bwlabel(L2,4);

fprintf('\t\t# of regions: L1=%i, L2=%i\n',num_regions1, num_regions2)

Region1 = regionprops(L1,'PixelIdxList','Area');
Region2 = regionprops(L2,'PixelIdxList','Area');

dist_turbo = 0;
for i=1:num_regions1
    dist_gi = 0;
    for j=1:num_regions2
        I = size(intersect_sorted(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
        if I > 0
            dist_gi = dist_gi + Region2(j).Area;%size(Region2(j).PixelIdxList,1);                        
        end
    end
    dist_gi = (dist_gi - Region1(i).Area) / Region1(i).Area;
    %dist_gi = (dist_gi - size(Region1(i).PixelIdxList,1)) / size(Region1(i).PixelIdxList,1);
    dist_gi
    dist_turbo = dist_turbo + dist_gi;    
end

dist_turbo = dist_turbo / num_regions1;


end


function dist_slic = SLIC_UPD(L1,L2)
%SLIC_UPD computes under segmentation error as defined by the SLIC paper.
%
% U = 1/N [ \sum_{i=1}^M ( \sum_{sj | sj \cup gi > B} |sj| ) - N ],
%       where N = the size of the image in pixels.
%             gi = the ith groundtruth segment (there are M total)
%             sj = the jth superpixel
%             B = .05 * |sj|, to account for error in the ground truth.
%
% From: SLIC Superpixels.
%       EPFL Technical  Report, 2010.


N = size(L1,1) * size(L1,2); % size of the image in pixels.

[L1,num_regions1] = bwlabel(L1,4);
[L2,num_regions2] = bwlabel(L2);

Region1 = regionprops(L1,'PixelIdxList');
Region2 = regionprops(L2,'PixelIdxList');

fprintf('\t\t# of regions: L1=%i, L2=%i\n',num_regions1, num_regions2)

dist_slic = 0;
for i=1:num_regions1
    for j=1:num_regions2
        I = size(intersect_sorted(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
        B = .05 * size(Region2(j).PixelIdxList,1); % 5% of the size of j to account for errors.
        if I > B
            dist_slic = dist_slic + size(Region2(j).PixelIdxList,1);                        
        end
    end
end

dist_slic = dist_slic - N;
dist_slic = dist_slic / N;

end


function c = intersect_sorted(a,b)
%INTERSECT_SORTED     Set intersection between sorted sets.
% INTERSECT_SORTED(A,B) when A and B are vectors returns the values common
% to both A and B.  A and B must be sorted and unique, and the result will be
% sorted and unique.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

c = a(ismember_sorted(a,b));

end



function [tf,loc] = ismember_sorted(a,s)
%ISMEMBER_SORTED   True for member of sorted set.
% ISMEMBER_SORTED(A,S) for the vector A returns an array of the same size as A
% containing 1 where the elements of A are in the set S and 0 otherwise.
% A and S must be sorted and cannot contain NaN.
%
% [TF,LOC] = ISMEMBER_SORTED(A,S) also returns an index array LOC where
% LOC(i) is the index in S which matches A(i) (highest if there are ties)
% or 0 if there is no such index.
%
% See also ISMEMBER, MATCH_SORTED, INTERSECT_SORTED, SETDIFF_SORTED, UNION_SORTED.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

% The internal function ismembc comes from ismember.m
% It requires non-sparse arrays.
a = full(a);
s = full(s);
if nargout < 2
  tf = ismembc(a,s);
else
  loc = ismembc2(a,s);
  tf = (loc > 0);
end

end