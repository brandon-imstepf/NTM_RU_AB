function data_new = DataToCCF(data_old, studyname, matdir_)
% DataToCCF  Re-index experimental tau data from study ROIs onto CCF atlas indices.
%
% Mouse tauopathy data is collected in experiments that measure tau in a
% specific set of brain regions (study-specific ROIs, often ~50-100 regions).
% The NTM model, however, operates on the Allen Common Coordinate Framework
% (CCF) atlas, which has 426 brain regions.  This function bridges the two:
% it maps study ROI values onto their corresponding CCF region indices.
%
% INPUTS
%   data_old  – experimental data in study-ROI space: [n_study_rois x nt]
%               If empty, the data is loaded directly from mousedata_struct
%               (useful when you want the raw data without pre-processing).
%   studyname – string key into mousedata_struct selecting the study
%               (e.g. 'Hurtado', 'IbaP301S').
%   matdir_   – path to the folder containing Mouse_Tauopathy_Data_HigherQ.mat
%
% OUTPUT
%   data_new  – data re-indexed to CCF space: [426 x nt]
%               Entries for CCF regions not covered by the study are NaN.
%               Multiple CCF regions can map to the same study ROI
%               (when a coarse study ROI corresponds to several CCF sub-regions).

% -------------------------------------------------------------------------
% Load mouse tauopathy data.
% mousedata_struct is a struct with one field per study (e.g. .Hurtado).
% Each study has fields:
%   .data    – measured tau values [n_study_rois x nt]
%   .regions – cell array mapping study ROI indices to CCF indices
%   .seed    – which study ROI(s) were initially injected (or NaN if unknown)
% -------------------------------------------------------------------------
load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'], 'mousedata_struct');

% If no data was passed in, fall back to the raw data stored in the struct.
if isempty(data_old)
    data_old = mousedata_struct.(studyname).data;
end

% -------------------------------------------------------------------------
% Build the CCF-space output array.
% rois is a cell array of length n_study_rois; rois{i} is a vector of CCF
% region indices that correspond to study ROI i.
% -------------------------------------------------------------------------
rois = mousedata_struct.(studyname).regions(:, 2);

% Pre-fill with NaN: CCF regions with no data from this study stay NaN.
data_new = NaN(426, size(data_old, 2));   % 426 = total regions in the CCF atlas

% -------------------------------------------------------------------------
% For each study ROI i, copy its data value into ALL corresponding CCF regions.
% A single study ROI can span several CCF sub-regions (many-to-one mapping).
% -------------------------------------------------------------------------
for i = 1:length(rois)
    roi_inds_i = rois{i};           % vector of CCF indices for study ROI i
    for j = 1:length(roi_inds_i)
        data_new(roi_inds_i(j), :) = data_old(i, :);   % copy row
    end
end

end
