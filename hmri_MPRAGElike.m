function fn_out = hmri_MPRAGElike(fn_in,params)

% Function to calculate an MPRAGE-like image from a set of (2 or 3) images
% acquired with the MPM protocol.
% This procedure is directly derived from M-A Fortin's work and is only a
% Matlab reimplemenation of his work. Check the reference here under.
%
% FORMAT
% fn_out = hmri_MPRAGElike(fn_in,params)
%
% INPUT
% fn_in : char array of filenames of 2 or 3 input images,
%           1st one should be T1w (numerator),
%           2nd (+3rd if provided) should be MTw and/or PDw (denominator)
% params : structure with some parameters
%   .lambda : regularisation parameter(s) [100, def]
%             If several values are passed, i.e. in a vector, then 1 image
%             is created per value. These are labdeled 'l100' and 'l200'
%             for lambada 100 and 200 for example.
%   .indiv  : a binary flag, to decide whether individual images are
%             created for each 2nd and 3rd input filename when 3 images are
%             passed in fn_in. These images will be labelled 'i1' and 'i2'.
%             [false, def.]
%   .thresh : threshold for [min max]range in MPRAGE-like image as in paper
%             but if left empty, NO thresholding applied. [[0 500], def.]
%   .BIDSform : a binary flag, to indicate if BIDS format is followed.
%               [false, def.]
%
% OUTPUT
% fn_out : char array of filenames of generated image(s)
%          Depending on the input parameters there could be up to 3 images
%          per lambda value passed.
%
% NOTE
% A) using spm_imcalc vs loading images
% Loading the images seems very appealing but then this assumes the images
% are already coregistered, as we would work with the same voxel grid. Thus
% to allow the automatic coregistration of images (i.e. updating the
% voxel-to-realworld mappin) one should be using spm_imcalc as it use the
% realworld space of the 1st image!
% Loading the MPRAGE-like images makes it easier for the thresholding/
% fixing of the image though
% B) Thresholding of MPRAGE-like image
% In the MPRAGE-like image, capping values of extremely high values seems a
% good idea but not so sure for the negative ones...
% A better idea could be to 1/ take their absolute value, then 2/ divide by
% lambda. This way we keep some positive but low-level signal
% C) FOllowing BIDS format
% This remains to be implemented...
% - Filename should be suffixed with 'MPRAGElike' instead of 'MPM'
% - if/when iamges are coregistered, then one should do the job in a
%   separate derivatives folder
%
%
% REFERENCCE
% Fortin M.-A. et al., 2025: https://doi.org/10.1002/mrm.30453]
% Original repository https://github.com/mafortin/mprage-like
%_______________________________________________________________________
% Copyright (C) 2025 Cyclotron Research Centre

% Written by C. Phillips,
% Cyclotron Research Centre, University of Liege, Belgium

% Set defaults & check input
params_def = struct(...
    'lambda', 100, ...
    'indiv', false, ...
    'thresh', [0 500], ...
    'BIDSform', false);
if nargin<2, params = params_def; end

Nimg_in = size(fn_in,1);
if Nimg_in > 3
    error('Too many input images');
elseif Nimg_in < 2
    error('Too few input images');
end
if Nimg_in==2 && params.indiv
    params.indiv = false ; % No individual images if only 2 input images
end

% Prepare output filenames and check nr of files to be created
Nlambda = numel(params.lambda);
if Nlambda>1
    % Plan saving the different lambdas in file name
    Nd_lambda = ceil(log10(max(params.lambda(:))+1));
    % Number of digits to write the lambdas
end

pth_in = spm_file(fn_in(1,:),'fpath');
pth_out = pth_in;
fn_basename = spm_file(fn_in(1,:),'basename');
fn_tmp = fullfile(pth_out, fn_basename);

% Map input volumes
V_in = spm_vol(fn_in);
% Read in all the data
% -> more memory use but more flexible treatement than with spm_imcalc
% val_in = spm_read_vols(V_in); sz = size(val_in);
% vval_in = reshape(val_in,[prod(sz(1:3)) sz(4)])'; % row-vectorized images

% Prepare output volume(s)
V_out = V_in(1);
V_out.dt(1) = 16; % use floats!
V_out.descrip = 'MPRAGE-like image';

% Define spm_imcalc flag
ic_flags = struct( ...
    'dmtx', true, ... % load data matrix
    'dtype', 16, ...  % use floats!
    'interp', 4, ...  % 4th degree B-splines as for normalization
    'descrip', 'MPRAGE-like image');

% Get the job done
for ii=1:Nlambda
    lambda = params.lambda(ii);
    if Nlambda==1 % just one lambda
        fn_out_ii = [fn_tmp,'_MPRAGElike.nii'];
    else
        %         fn_out_ii = [sprintf('%s_l%4d',fn_tmp,round(params.lambda(ii))),'_MPRAGElike.nii'];
        fn_out_ii = sprintf( ... 
            sprintf('%%s_MPRAGElike_l%%0%dd.nii',Nd_lambda), ... % Adding the right number of 0's
            fn_tmp,lambda);
        
    end
    V_out.fname = fn_out_ii;
    %     vval_MPRlike = (vval_in(1,:)-lambda)./(mean(vval_in(2:end,:),1)+lambda);
    spm_imcalc(V_in, V_out, ...
        '(X(1,:)-lambda)./(mean(X(2:end,:),1)+lambda)', ic_flags, lambda);
    
    % Check thresholding
    if ~isempty(params.thresh) % apply thresholding
        % load MPRAGElike image into memory
        val_MPRlike = spm_read_vols(V_out); 
        % Deal with "thresholding"
        %         vval_MPRlike(vval_MPRlike<params.thresh(1)) = 0;
        val_MPRlike(val_MPRlike(:)<params.thresh(1)) = ...
            abs(val_MPRlike(val_MPRlike(:)<params.thresh(1)))/lambda;
        val_MPRlike(val_MPRlike(:)>params.thresh(2)) = 0;
        % Write out fixed volume
        spm_write_vol(V_out,val_MPRlike);
    end
    
end

fn_out = fn_out_ii;

end