function mpragelike = tbx_scfg_hmri_MPRAGElike
% Configuration file for the MPRAGElike module,
% when including the tool into the hMRI toolbox
%__________________________________________________________________________
% Copyright (C) 2026 Cyclotron Research Centre

% Written by C. Phillips,
% Cyclotron Research Centre, University of Liege, Belgium

%--------------------------------------------------------------------------
% IMAGES
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% imgt1w T1w Image
%--------------------------------------------------------------------------
imgt1w         = cfg_files;
imgt1w.tag     = 'imgt1w';
imgt1w.name    = 'T1w Image';
imgt1w.help    = {'This is the image that is placed in the numerator.'};
imgt1w.filter  = 'image';
imgt1w.ufilter = '.*';
imgt1w.num     = [1 1];
imgt1w.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% imgmtw MTw Image
%--------------------------------------------------------------------------
imgmtw         = cfg_files;
imgmtw.tag     = 'imgmtw';
imgmtw.name    = 'MTw Image';
imgmtw.val     = {{''}};
imgmtw.help    = {'This is (one of) the image(s) that is placed in the denumerator.'};
imgmtw.filter  = 'image';
imgmtw.ufilter = '.*';
imgmtw.num     = [1 1];
imgmtw.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% imgpdw PDw Images
%--------------------------------------------------------------------------
imgpdw         = cfg_files;
imgpdw.tag     = 'imgpdw';
imgpdw.name    = 'PDw Image';
imgpdw.val     = {{''}};
imgpdw.help    = {'This is (one of) the image(s) that is placed in the denumerator.'};
imgpdw.filter  = 'image';
imgpdw.ufilter = '.*';
imgpdw.num     = [1 1];
imgpdw.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% OPTIONS
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% lambda Lambda for regularisation
%--------------------------------------------------------------------------
lambda         = cfg_entry;
lambda.tag     = 'lambda';
lambda.name    = 'Lambda';
lambda.help    = {
    'The regularisation parameter ''lambda''.'
    ''
    'If set to NaN then an automatic procedure is used to estimate the value.'
    'It can also be a real value (or a vector) to test a (or several) value(s) of lambda.'
    }';
lambda.strtype = 'e';
lambda.num     = [1 Inf];
lambda.val     = {NaN}; % Estimate lambda 'on the fly' from the images

%--------------------------------------------------------------------------
% indivimg Build individidual maps if both PDw & MTw are provided, or not
%--------------------------------------------------------------------------
indivimg         = cfg_menu;
indivimg.tag     = 'indivimg';
indivimg.name    = 'Build individual maps from pair of MTw & PDw';
indivimg.help    = {
    'If both MTw and PDw images are provided,an MPRAGElike image is calculate with the average of these.'
    'Still, one could be interested in building an MRPAGElike image from each of these individually.'
    'Thus this options allows to get 3 MPRAGElike images from MTw and PDw images.'
    }';
indivimg.labels  = {
                    'Just mean of MTw and PDw'
                    'Mean and individual from MTw and PDw'
}';
indivimg.values  = {
                    false
                    true
}';
indivimg.val    = {false} ; % No individual images by default

%--------------------------------------------------------------------------
% thresh Threshold for intensity of resulting MPRAGElike image
%--------------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Intensity thresholds';
thresh.help    = {
    'Threshold for intensity of resulting MPRAGElike image'
    'Acceptable [min max] values for the reconstructed image'
    }';
thresh.strtype = 'r';
thresh.num     = [1 2];
thresh.val     = {[0 500]};

%--------------------------------------------------------------------------
% coreg Coregister the MTw/PDw images onto the T1w image, or not
%--------------------------------------------------------------------------
coreg         = cfg_menu;
coreg.tag     = 'coreg';
coreg.name    = 'Coregister the MTw/PDw onto the T1w';
coreg.help    = {
    'Coregister the MTw/PDw images onto the T1w image.'
    'If the T1w, MTw, PDw are not in the same space, then one shoudl coregister them before calculating the MPRAGElike image.'
        }';
coreg.labels  = {
                    'No coregistration'
                    'Coregister MTw/PDw onto T1w'
}';
coreg.values  = {
                    false
                    true
}';
coreg.val     = {false} ; % No coregistration by default

%--------------------------------------------------------------------------
% bidsform BIDS formating of the data
%--------------------------------------------------------------------------
bidsform         = cfg_menu;
bidsform.tag     = 'bidsform';
bidsform.name    = 'BIDS formating';
bidsform.help    = {
    'Follow BIDS convention, or not.'
    'If data are BIDS complient, then one should preserve the BIDS organisation.'
        }';
bidsform.labels  = {
                    'Not BIDS compliant data'
                    'Follow BIDS organisation'
}';
bidsform.values  = {
                    false
                    true
}';
bidsform.val     = {false} ; % Not BIDS by default

%--------------------------------------------------------------------------
% options Estimation Options
%--------------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Creation Options';
options.val     = {lambda indivimg thresh coreg bidsform};
options.help    = {'Various creation options for the MPRAGE-like image.'};

%--------------------------------------------------------------------------
% estimate Coreg: Estimate
%--------------------------------------------------------------------------
mpragelike         = cfg_exbranch;
mpragelike.tag     = 'mpragelike';
mpragelike.name    = 'MPRAGE-like';
mpragelike.val     = {imgt1w imgmtw imgpdw options};
mpragelike.help    = {
    'MPRAGE-like image creation.'
    ''
    'The method used here is based on work by Fortin et al. 2025 (https://doi.org/10.1002/mrm.30453).'
    ''
    'The idea is to take the regularized ratio between a T1w image and a PDw or MTw (or mean of these 2) image(s). This results in a well-contrasted image looking like an actual MPRAGE acquisition.'
    ''
    'A few options are possible, esp. regarding the definition/estimation of the regularisation parameter.'
    }';
mpragelike.prog = @hmri_run_mpragelike;
mpragelike.vout = @vout_estimate;

end
%% OUTPUT and RUN functions

%==========================================================================
function dep = vout_estimate(job) %#ok<*INUSD>
dep(1)            = cfg_dep;
dep(1).sname      = 'MPRAGElike Image(s)';
dep(1).src_output = substruct('.','images');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

end

%==========================================================================
function out = hmri_run_mpragelike(job)

% Collect input
fn_in = char(job.imgt1w{1} , job.imgmtw{1} , job.imgpdw{1});

params = struct( ...
    'lambda', job.options.lambda, ...
    'indiv', job.options.indivimg, ...
    'thresh', job.options.thresh, ...
    'coreg', job.options.coreg, ...
    'BIDSform', job.options.bidsform);

% Call the processing function
fn_out = hmri_MPRAGElike(fn_in,params);

% Collect output
out.images = {fn_out};

end