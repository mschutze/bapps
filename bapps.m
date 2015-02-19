% Batch PET Preprocessing in SPM8 (bapps)
% Manuel Schütze - Feb 2015
%
% Description:
%   This batch scripts runs different SPM8 preprocessing routines on
%   pre-aligned nifti PET images in a folder.
% How to use:
%   First trasform your dicom images into nifti images (eg using MRIcron).
%   Rename images to patient IDs. Using SPM8's display function, align 
%   images so that the AC is at the origin and the PC is about the same 
%   axis. Put all images in the same folder. Navigate to that folder in
%   Matlab. Modify options in this script and run it by typing bapps.
%   Output folder with different steps of the preprocessing are created
%   automatically. Make sure to use check reg to check the results.

%% Config
%location of SPM8 folder
spm8 = '/Users/manuel/Documents/MATLAB/spm8/';
%template used for normalization
template = '/Users/manuel/Documents/MATLAB/spm8/templates/PET_FDG.nii';
