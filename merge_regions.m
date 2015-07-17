%========================================================================%
%    Performing logic for merging two regions                            %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%


function [L,Boundary,Regions] = merge_regions(L,Boundary,Regions,old_label,new_label)
%MERGE_REGIONS merges regions in the label matrix.
%
% MERGE_REGIONS(L,Boundary,Regios,old_label,new_label) Merges 'old_label'
% and 'new_label' in the label matrix 'L' by renaming all 'old_label'
% pixels to 'new_label' and removing the boundary between them. Also updates
% the list of 'new_label' pixels in 'Regions' and the updated 'Boundary'
% pixels.
%
% Note: 'new_label' exists and should be the min of the two labels being
% merged.
%

% TODO: 
%
%   1 1 1 1
%   1 0 0 1
%   1 0 0 1
%   1 1 1 1
%
% This can happen at the boundary of 4 regions (although only one is shown
% here). Right now, these blocks are being removed by another
% post-processing script, but it would be nice if region merging also took
% this into account explicitly. Doing it via 8hood has issues.
%


[num_rows num_cols] = size(L);

% get the boundary pixels for both regions.
merged_boundaries = [Boundary(old_label).Pixels; Boundary(new_label).Pixels];

% store new_label boundaries; initially over-estimate (duplicates + removed boundaries).
Boundary(new_label).Pixels = merged_boundaries;

% boundaries/rows of merged_boundaries to delete: 1 if delete; 0 to keep.
rows_to_delete = zeros(size(merged_boundaries,1),1);

% change old_label to new_label.
L(Regions(old_label).PixelIdxList) = new_label;


%% recursively go through boundaries of both regions while things change.
changed = 1;
while (changed)
    changed = 0;
    
    for i=1:size(merged_boundaries,1)
        
        % get boundary pixel.
        x = merged_boundaries(i,1);
        y = merged_boundaries(i,2);
        
        % already processed.
        if L(x,y) ~= 0 % or == new_label or 'rows_to_delete(i) == 1'
            %if (L(x,y) ~= new_label)
            %    error('Boundary pixel changed to something (%i) other than new_label (%i)',L(x,y),new_label)
            %end
            
            rows_to_delete(i) = 1; % ensure it's set to 1 to remove duplicates.
            continue;
        end

        % check if the boundary is surrounded by new_label.

        % north and south.
        if (x-1 > 0) && (x+1 <= num_rows) && (L(x-1,y) == new_label) && (L(x+1,y) == new_label)
            L(x,y) = new_label;
            rows_to_delete(i) = 1;
            changed = 1;
        end

        % west and east.
        if (y-1 > 0) && (y+1 <= num_cols) && (L(x,y-1) == new_label) && (L(x,y+1) == new_label)
            L(x,y) = new_label;
            rows_to_delete(i) = 1;
            changed = 1;
        end             
        
    end
    
end

% update the list of pixels for new_label.
old_boundary_pixels = Boundary(new_label).Pixels(rows_to_delete==1,:);
old_boundary_pixels = sub2ind(size(L), old_boundary_pixels(:,1), old_boundary_pixels(:,2)); % convert 2d to 1d units.
%if size(old_boundary_pixels,1) < 1, 
%    error('No boundary pixels were changed when merging regions %i and %i',old_label,new_label)
%end
                                    % old region pixels              removed boundaries      new_label_pixels     
Regions(new_label).PixelIdxList = [Regions(old_label).PixelIdxList; old_boundary_pixels; Regions(new_label).PixelIdxList];

% delete the duplicate/merged boundary pixels.
Boundary(new_label).Pixels(rows_to_delete==1,:) = [];


end
