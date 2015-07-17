%========================================================================%
%    REMOVES SMALL REGIONS                                               %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%

function [L,reg_sizes] = remove_small_regions(L,xx)
%REMOVE_SMALL_REGIONS   Removes regions smaller than xx from L.
%
% REMOVE_SMALL_REGIONS(L,xx) removes regions that consist of less than xx 
% pixels. The small regions are merged with a bordering region (arbitrarily
% picked in case there are multiple bordering regions).
%
% Algorithm: 
% 1. Use regionprops to find the size of each region.
% 2. Identify the regions which have size < xx.
% 3. For each such region, look at its two-hop neighborhood to find a
% neighboring region (b).
% 4. Relabel all the small region's pixels to b.
% 5. Relabel boundary pixels to b.
%
% Returns: new label matrix.
%

% TODO: faster way to do this that doesn't require the two-step process?
% TODO: use template from run_merging and call the appropriate functions to
% merge small regions with their most similar neighbor.

[num_rows, num_cols] = size(L);

reg_sizes = regionprops(L,'Area');
for reg=1:size(reg_sizes,1)
    
    if (reg_sizes(reg).Area == 0) % region doesn't exist.
        continue
    elseif (reg_sizes(reg).Area > xx) % region is large enough.
        continue
    else % region too small, let's do it.
        
        [reg_rows, reg_cols] = find(L==reg);
        
        % find the bordering region to replace the small region with.
        bndry_region = -1;

        for i=1:size(reg_rows)
            x = reg_rows(i);
            y = reg_cols(i);
            
            % two to the west.
            if y-2 > 0 % not on the image border.       
                %  pass thru boundary && hit new region
                if (L(x,y-1) == 0) && (L(x,y-2) ~= 0) && (L(x,y-2) ~= reg)                   
                    bndry_region = L(x,y-2);
                    break;
                end
            end

            % two to the east.
            if y+2 <= num_cols
                if (L(x,y+1) == 0) && (L(x,y+2) ~= 0) && (L(x,y+2) ~= reg)
                    bndry_region = L(x,y+2);
                    break;
                end
            end

            % two north.
            if x-2 > 0
                if (L(x-1,y) == 0) && (L(x-2,y) ~= 0) && (L(x-2,y) ~= reg)
                    bndry_region = L(x-2,y);
                    break;
                end
            end

            % two south.
            if x+2 <= num_rows
                if (L(x+1,y) == 0) && (L(x+2,y) ~= 0) && (L(x+2,y) ~= reg)
                    bndry_region = L(x+2,y);
                    break;
                end
            end 
            
        end
        
        
        if bndry_region == -1
            continue;
            % leave out for now because a region might only be diagonally
            % adjacent to another region, in which case its adjacency
            % won't be detected.
            %error('Boundary region (%i) does not exist for %i',bndry_region,reg)
        end       
        

        % now, replace.
        for i=1:size(reg_rows)
            x = reg_rows(i);
            y = reg_cols(i);
                       
            % two to the west.
            if y-2 > 0 % not on the image border.       
                %  pass thru boundary && hit new region
                if (L(x,y-1) == 0) && (L(x,y-2) == bndry_region)
                    
                    % change the 0 to reg -- i could change it to
                    % bndry_region, but i'll need to convert all reg pixels
                    % to bndry_region later anyways.
                    L(x,y-1) = reg;
                                       
                    % update the count for the new region.
                    reg_sizes(bndry_region).Area = reg_sizes(bndry_region).Area + 1;
                end
            end

            % two to the east.
            if y+2 <= num_cols
                if (L(x,y+1) == 0) && (L(x,y+2) == bndry_region)
                    
                    L(x,y+1) = reg;
                    reg_sizes(bndry_region).Area = reg_sizes(bndry_region).Area + 1;                    
                end
            end

            % two north.
            if x-2 > 0
                if (L(x-1,y) == 0) && (L(x-2,y) == bndry_region)
                    
                    L(x-1,y) = reg;                    
                    reg_sizes(bndry_region).Area = reg_sizes(bndry_region).Area + 1;                    
                end
            end

            % two south.
            if x+2 <= num_rows
                if (L(x+1,y) == 0) && (L(x+2,y) == bndry_region)
                    
                    L(x+1,y) = reg;
                    reg_sizes(bndry_region).Area = reg_sizes(bndry_region).Area + 1;                    
                end
            end             
        end
               
        % change all small region to the larger region.
        L(L==reg) = bndry_region;
        
        % recompute region areas.
        reg_sizes(bndry_region).Area = reg_sizes(bndry_region).Area + reg_sizes(reg).Area;
        reg_sizes(reg).Area = 0;

        
    end
    
end

end