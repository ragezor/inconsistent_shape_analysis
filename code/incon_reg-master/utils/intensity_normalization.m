function [normalized_intensity, v_min, v_max] = intensity_normalization(intensity, varargin)
if isempty(varargin)
    v_min = min(intensity(:)); v_max = max(intensity(:));
else
    v_min = varargin{1}; v_max = varargin{2};
end

normalized_intensity = (intensity - v_min) ./ (v_max - v_min);


end