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
%   axis. Put all images in the same folder. Navigate to the project folder
%   in Matlab. Modify options in this script and run it by typing bapps.
%   Output folder with different steps of the preprocessing are created
%   automatically according to defined names. Make sure to use check reg 
%   to check the results.

function [] = bapps()

    clc;
    fprintf('Batch PET Preprocessing in SPM8 (bapps)\n\n');

    %% Config
    % location of SPM8 folder
    spm8 = '/Users/manuel/Documents/MATLAB/spm8/';
    % project folder (set by running this script from inside)
    root = pwd; 
    % run the following routines (in this order)
    normalize = true;
    segment = true;
    make_gm_mask = true;
    smooth = true;
    correct = true;
    % folder names (will be created if don't exist)
    f_aligned = 'c.Nifti_alinhado';
    f_normalized = 'd.Normalized';
    f_segmented = 'e.Segmented'; 
    f_smoothed = 'f.Smoothed';
    f_corrected = 'g.Corrected';
    % normalization: template
    template = '/Users/manuel/Documents/MATLAB/spm8/templates/PET_FDG.nii';
    % gm mask: threshold (any voxel with probability higher than
    % set here will be considered grey matter and included in the mask)
    gmth = '0.9';
    % gm mask: name
    gm_name = 'gm_mask_90.nii';
    % smoothing: FWHM filter
    fwhm = '12';
    % correct global signal: mask (mask used to define brain voxels for
    % global signal correction - all voxels in image are divided by the mean
    % value inside the mask - eg whole brain mask, cerebellum mask, etc)
    corr_mask = '/Users/manuel/Documents/MATLAB/mask.nii';
    % correct global signal: alterntive method (set true to first calculate
    % the mean of the whole image, then use mean/8 as a theshold to calculate
    % a second mean on the basis that any voxel with a value below is 
    % likely to be outside the brain)
    corr_alternative = false;
    % correct global signal: prefix
    corr_pref = 'nGB';
    % end config
    spm_jobman('initcfg');
    spm('defaults', 'PET');

    %% Normalization
    if(normalize)
        fprintf('-> Starting normalization routine\n\n');

        % check if input folder exists
        if(exist(f_aligned,'dir') == 0) 
            error(['The folder ',f_aligned,' does not exist']);
        end

        % find nifti files
        cd(f_aligned)
        files = dir('*.nii');
        nf = length(files);
        if(nf == 0)
            error(['No Nifti files found in ',f_aligned','\n']);
        else
            fprintf(['Found ',int2str(nf),' Nifti files in ',f_aligned,'\n\n']);
        end

        % run normalization for each image
        for i=1:nf
            % define variables for individual subjects
            curFile = files(i).name;
            curSubj = strcat(root,'/',f_aligned,'/',curFile);
            % run matlabbatch job
            fprintf(['Processing subject ',curFile,' (',int2str(i),' of ',int2str(nf),')']);
            clear matlabbatch
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {curSubj};
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {curSubj};
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {template};
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
            matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-90 -126 -72
                                                                         90 90 108];
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
            spm_jobman('run', matlabbatch);
        end

        %move files to f_normalized and cleanup
        movefile(strcat('w*.nii'),strcat('../',f_normalized));
        delete('*_sn.mat');
        cd(root);
    end
    
    %% Segmentation
    if(segment)
        fprintf('-> Starting segmentation routine\n\n');

        % check if input folder exists
        if(exist(f_normalized,'dir') == 0) 
            error(['The folder ',f_normalized,' does not exist']);
        end

        % find nifti files
        cd(f_normalized)
        files = dir('*.nii');
        nf = length(files);
        if(nf == 0)
            error(['No Nifti files found in ',f_normalized','\n']);
        else
            fprintf(['Found ',int2str(nf),' Nifti files in ',f_normalized,'\n\n']);
        end

        % run segmentation for each image
        for i=1:nf
            % define variables for individual subjects
            curFile = files(i).name;
            curSubj = strcat(root,'/',f_normalized,'/',curFile);
            % run matlabbatch job
            fprintf(['Processing subject ',curFile,' (',int2str(i),' of ',int2str(nf),')']);
            clear matlabbatch
            matlabbatch{1}.spm.spatial.preproc.data = {curSubj};
            matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
            matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 0];
            matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
            matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
            matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
            matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                                                           [spm8 'tpm/grey.nii']
                                                           [spm8 'tpm/white.nii']
                                                           [spm8 'tpm/csf.nii']
                                                           };
            matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
                                                             2
                                                             2
                                                             4];
            matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
            matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
            matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
            matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
            matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
            spm_jobman('run', matlabbatch);
        end

        % move files to f_segmented and cleanup
        movefile(strcat('c1w*.nii'),strcat('../',f_segmented));
        delete('mw*.nii');
        delete('*_sn.mat');
        cd(root);
    end
    
    %% Make GM mask
    if(make_gm_mask)
        fprintf('-> Creating GM mask\n\n');

        % check if input folder exists
        if(exist(f_segmented,'dir') == 0) 
            error(['The folder ',f_segmented,' does not exist']);
        end

        % find nifti files
        cd(f_segmented)
        c1 = dir('c1*.nii');
        nf = length(c1);
        maps = cell(1,nf);
        if(nf == 0)
            error(['No c1* Nifti files found in ',f_segmented','\n']);
        end

        % generate soft mean expression to be used by ImCalc
        global exp;
        exp = '(';
        for k=1:nf
            if k==1
                exp = strcat(exp,'i',int2str(k));
            else
                exp = strcat(exp,'+i',int2str(k));
            end
        end
        exp = strcat(exp,')./(eps');
        for l=1:nf
            exp = strcat(exp,'+(i',int2str(l),'~=0)');
        end
        exp = strcat(exp,')');

        expression = exp; %something like 'i1+iN./(eps+(i1~=0)+(iN~=0))'
        gmth_exp = strcat('i1>',gmth);

        % run ImCalc in matlabbatch job
        for i=1:nf
            curFile = c1(i).name; %get c1 (GM) images
            maps{1,i} = curFile;
        end
        temp_mean = strcat(root,'/',f_segmented,'/','gm_mean.nii');
        clear matlabbatch
        % Soft mean (SPM > Image Calculator)
        matlabbatch{1}.spm.util.imcalc.input = maps;
        matlabbatch{1}.spm.util.imcalc.output = temp_mean;
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = expression;
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 2; %saving as 8bit UINT
        % Threshold (SPM > Image Calculator)
        matlabbatch{2}.spm.util.imcalc.input = {temp_mean};
        matlabbatch{2}.spm.util.imcalc.output = gm_name;
        matlabbatch{2}.spm.util.imcalc.outdir = {''};
        matlabbatch{2}.spm.util.imcalc.expression = gmth_exp;
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 2; %saving as 8bit UINT
        spm_jobman('run', matlabbatch);
        % move files to f_segmented and cleanup
        cd(root);
    end
    
    %% Smoothing
    if(smooth)
        fprintf('-> Starting smoothing routine\n\n');

        % check if input folder exists
        if(exist(f_normalized,'dir') == 0) 
            error(['The folder ',f_normalized,' does not exist']);
        end

        % find nifti files
        cd(f_normalized)
        files = dir('*.nii');
        nf = length(files);
        if(nf == 0)
            error(['No Nifti files found in ',f_normalized','\n']);
        else
            fprintf(['Found ',int2str(nf),' Nifti files in ',f_normalized,'\n\n']);
        end

        % run smoothing for each image
        for i=1:nf
            % define variables for individual subjects
            curFile = files(i).name;
            curSubj = strcat(root,'/',f_normalized,'/',curFile);
            % run matlabbatch job
            fprintf(['Processing subject ',curFile,' (',int2str(i),' of ',int2str(nf),')']);
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = {curSubj};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [str2num(fwhm) str2num(fwhm) str2num(fwhm)];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = strcat('s',fwhm);
            spm_jobman('run', matlabbatch);
        end

        % move files to f_segmented and cleanup
        movefile(strcat('s',fwhm,'w*.nii'),strcat('../',f_smoothed));
        cd(root);
    end
    
    %% Global signal correction
    if(correct)
        fprintf('-> Starting global signal correction\n\n');

        % check if input folder exists
        if(exist(f_smoothed,'dir') == 0) 
            error(['The folder ',f_smoothed,' does not exist']);
        end

        % find nifti files
        cd(f_smoothed)
        files = dir('s*.nii');
        nf = length(files);

        for i=1:nf
            fprintf(['Processing subject ',curFile,' (',int2str(i),' of ',int2str(nf),') ']);
            % read in image volume for current subject
            img_header = spm_vol(files(i).name);
            img = spm_read_vols(img_header);
            
            % read in mask and transform into logical array (1=True / 0=False)
            % or use alternative method of thresholding voxel values
            if(corr_alternative)
                firstpassmean  = mean(img);
                mask = find(img > (firstpassmean/8));
            else
                mask = ismember(spm_read_vols(spm_vol(corr_mask)),1);
            end

            %Divide each voxel by the mean count inside mask
            mean_val = mean(img(mask));
            img_norm = img/mean_val;

            %Save normalized image
            img_header.fname = [corr_pref files(i).name];
            spm_write_vol(img_header,img_norm);
            fprintf('Done!\n');
        end
        movefile(strcat(corr_pref,'*.nii'),strcat('../',f_corrected));
        cd(root);
    end
    fprintf('\nFinished bapps');
end