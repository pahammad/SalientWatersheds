%========================================================================%
%    Measures for computing the distance between two regions             %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%


function dist = reg_distance(reg1,reg2,EMD_dist,Reg2Intensity,Reg2Texture,Boundary_size,Regions)
%REG_DISTANCE The distance between the region properties of two regions.
%
% REG_DISTANCE(reg1,reg2,EMD_dist,Reg2Intensity,Reg2Texture) computes the
% distance between 'reg1' and 'reg2' in feature-space ('Reg2Intensity' and
% 'Reg2Texture') using the EMD with cost matrix 'EMD_dist'. 
%
% Let:
%    d = EMD(Int(reg1),Int(reg2)) +  \alpha * \sum_{i=1}^8 EMD(Text(reg1,i),Text(reg2,i)) / 8),
%       where alpha is the weight of the texture component.
%
% Then, the function returns: exp(-d)
%
% If Reg2Intensity or Reg2Texture equal '[]', they will not be used as part
% of the distance calculation.
%

alpha = 0.1; % weight of the texture component.

dist = 0;
if size(Reg2Intensity,1) ~= 0
    
    F1 = Reg2Intensity(reg1,:);
    F2 = Reg2Intensity(reg2,:);    
    %F1 = F1 ./ sum(F1);
    %F2 = F2 ./ sum(F2);
   
    dist = emd_hat_gd_metric_mex(F1',F2',EMD_dist);

    % Old distance measures.
    %dist = dist * exp( -emd_hat_gd_metric_mex(F1',F2',EMD_dist) / (2*1^2));
    %dist = dist * exp(-relative_entropy(F1,F2,0)/20);
    %dist = dist * exp(-entropy(F1,F2));
    %dist = dist * exp( -sum( ((F1-F2).^2) )); WHY NEGATIVE?
    %dist = dist * exp( (-sqrt(sum((F1-F2).^2))) / 20); WHY NEGATIVE?
end

if size(Reg2Texture,1) ~= 0
    
    F1 = Reg2Texture(reg1,:,:);
    F2 = Reg2Texture(reg2,:,:);    
    t = 0;
    
    for i=1:8
        F1i = F1(1,:,i);% ./ sum(F1(1,:,i));
        F2i = F2(1,:,i);% ./ sum(F2(1,:,i));
        
        t = t + emd_hat_gd_metric_mex(F1i',F2i',EMD_dist);
                
        % Old distance measures.
        %t = t + sum(abs(F1i-F2i))/2; % EMD-equivalent.
        %t = t + relative_entropy(F1i,F2i,0);
        %t = t + entropy(F1i,F2i);
        %t = t + -sum( ((F1i-F2i).^2) ); WHY NEGATIVE?
        %t = t + -sqrt( sum( ((F1i-F2i).^2) ) ); WHY NEGATIVE?
    end
    dist = dist + (alpha*t/8);
end

% Ignore for now...
if size(Boundary_size,1) ~= 0
    dist = dist + 1/Boundary_size(reg1,reg2);
    %dist = dist + (size(Regions(reg1).PixelIdxList,1)+size(Regions(reg2).PixelIdxList,1) / Boundary_size(reg1,reg2)) / (100);
end 

dist = exp(-dist);

end




function H = entropy(F1,F2)
%ENTROPY The entropy of the merged region.
%

F12 = F1 + F2 + 1e-10; % make nonzero.
F12_prob = F12 ./ sum(F12);
H = -sum(F12_prob.*log(F12_prob));

if H < 0 || H > log(size(F1,2))
    error('Entropy is invalid: n=%i,H=%.3f',size(F1,2),H);
end

end


function RE = relative_entropy(F1,F2,merged)
%RELATIVE_ENTROPY The relative entropy of each to each other.
%
% If 'merged' is 1 then the RE is from F1 to F12 and F2 to F12. 
% Otherwise, the RE is from F1 to F2 and F2 to F1.
 
% TODO: add error check for bounds.

F1 = F1 + 1e-10;
F2 = F2 + 1e-10;
    
F1_prob = F1 ./ sum(F1);
F2_prob = F2 ./ sum(F2);

if merged
    F12 = F1 + F2 + 1e-10;
    F12_prob = F12 ./ sum(F12);
    
    RE1 = sum(F1_prob.*log(F1_prob./F12_prob));
    RE2 = sum(F2_prob.*log(F2_prob./F12_prob));
    
    RE = RE1 + RE2;
else
    RE1 = sum(F1_prob.*log(F1_prob./F2_prob));
    RE2 = sum(F2_prob.*log(F2_prob./F1_prob));
    RE = RE1 + RE2;
end

end