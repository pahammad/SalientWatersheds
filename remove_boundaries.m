%========================================================================%
%    REMOVE LEFTOVER BOUNDARIES                                          %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%

function L = remove_boundaries(L)
%REMOVE_BOUNDARIES  Removes leftover boundary pixels.
%
% REMOVE_BOUNDARIES(L) removes boundary (0) pixels that lie sandwiched
% between two same nonzero regions.

% Algorithm: Because the perimeter of every region is surrounded in all 8
% sides by a zero, we only need to check if a boundary pixel is surrounded
% on the north/south or west/east side by the same region. In other words,
% for the following configuration, we do not want to remove any boundary
% pixels because none are dominated by the same region.
%
%   1 0 0 0
%   0 0 0 0 
%   0 0 2 0
%
% Every true configuration has to look something like this:
%
%   0 0 0
%   0 1 0
%   0 0 0
%   2 2 0
%
% Example:
%
% img =
%
%     1     1     1     0     3
%     0     0     0     0     0
%     4     0     0     2     0
%     4     4     0     2     2
%     0     0     0     0     0
%     4     4     4     4     4
% 
% >> img = remove_boundaries(img)
%
% img =
%
%     1     1     1     0     3
%     0     0     0     0     0
%     4     0     0     2     0
%     4     4     0     2     2
%     4     4     0     0     0
%     4     4     4     4     4
%


[num_rows, num_cols] = size(L);

% repeat this procedure until nothing changes.
%   -- runs about 15 times.
changed = 1;
while (changed)
    changed = 0;
    for x=1:num_rows,
        for y=1:num_cols,

            % ignore non-boundary pixels.
            if L(x,y) ~= 0
                continue
            end

            % north and south.
            if (x-1 > 0) && (x+1 <= num_rows)
                %   not border     same region north & south  no salient
                if (L(x-1,y) ~= 0) && (L(x-1,y) == L(x+1,y)) %&& (S(x,y) == 0)
                    L(x,y) = L(x-1,y);    
                    changed = 1;
                end
            end
            
            % west and east.
            if (y-1 > 0) && (y+1 <= num_cols)
                %   not border     same region west and east
                if (L(x,y-1) ~= 0) && (L(x,y-1) == L(x,y+1)) %&& (S(x,y) == 0)
                    L(x,y) = L(x,y-1);     
                    changed = 1;
                end                        
            end

        end
    end
end
    

end