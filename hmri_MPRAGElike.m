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
%   .coreg  : a binary flag, to decide whether or not the input images
%             should be coregistered to the 1st one. [false, def.]
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
% C) Coregistration
% In order to keep the original data untouched, when the input images have
% to be coregistered to the 1st one, then the 2nd (and 3rd image if
% provided) are fist copied in a temporary folder. These coregistered
% temporary images are then deleted at the end of the process.
% D) FOllowing BIDS format
% This remains to be implemented...
% - Filename should be suffixed with 'MPRAGElike' instead of 'MPM'
% - if/when images are coregistered, then one should do the job in a
%   separate derivatives folder
%
% TO-DO's
% 1) deal with BIDS organized data
% 2) add (optional) coregsitration step in processing
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
    'coreg', false, ...
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
    N_MPRcreate = 1; % creating just 1 MPRAGElike image
else
    N_MPRcreate = 3; % creating 3 MPRAGElike images: 1 combined + 2 individuals
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
fn_basename = fullfile(pth_out,spm_file(fn_in(1,:),'basename'));
% fn_tmp = fullfile(pth_out, fn_basename);

% Deal with coregistration, if requested
% To preserve the original images, the 2nd and 3rd (if provided) image(s) 
% are simply copied in a temporary folder, coregistered to the 1st one,
% then after processing the temporary volumes are deleted !
if params.coreg
    fn_in_orig = fn_in; %#ok<*NASGU>
    fn_in_c = cellstr(fn_in);
    % temporary folder is created within root folder, this might be an 
    % issue on some systems... -> deal with it later if/when necessary
   pth_tmp = fullfile(pth_in,['tmp_',datestr(now,'yyyymmddTHHMMSS')]);
   if ~exist(pth_tmp,'dir'), mkdir(pth_tmp), end
   % copy original files as 'img2.nii' and 'img3.nii', 
   % then coregister by changing the mapping only
   flags_coreg = spm_get_defaults('coreg.estimate');
   flags_coreg.params   = [0 0 0  0 0 0]; % Starting estimates
   flags_coreg.graphics = false; % no display
   for ii=2:Nimg_in
       fn_tmp = fullfile(pth_tmp,sprintf('img%d.nii',ii));
       copyfile(deblank(fn_in(ii,:)), fn_tmp);
       x = spm_coreg(fn_in(1,:), fn_tmp, flags_coreg);
       % coregistration matrix & vx-to-mm matrix for moving image
       M  = spm_matrix(x); MM = spm_get_space(fn_tmp);
       % update vx-to-mm matrix for moving image
       spm_get_space(fn_tmp, M\MM)
       % keep name of temporary moving image
       fn_in_c{ii} = fn_tmp;
   end
   fn_in = char(fn_in_c);
end

% Map input volumes
V_in = spm_vol(fn_in);

% Prepare output volume structure
V_out = V_in(1);
V_out.dt(1) = 16; % use floats!
V_out.descrip = 'MPRAGE-like image';

% Define spm_imcalc flag
ic_flags = struct( ...
    'dmtx', true, ... % load data matrix
    'dtype', 16, ...  % use floats!
    'interp', 4);  % 4th degree B-splines as for normalization


% Get the job done
fn_out_c = cell(Nlambda,N_MPRcreate);
for ii=1:Nlambda
    % Check lambda value
    lambda = params.lambda(ii);
    if isnan(lambda) % Automatic definition of lambda from images
        lmabda = estimate_lambda(fn_in);
    end
    if Nlambda==1 % just one lambda
        fn_out_ii = [fn_basename,'_MPRAGElike.nii'];
    else % adding lambda alue as suffix if multiple values
        fn_out_ii = sprintf( ...
            sprintf('%%s_MPRAGElike-l%%0%dd.nii',Nd_lambda), ... % Adding the right number of 0's
            fn_basename,lambda);        
    end
    % Save lambda in image header (description field)
    ic_flags.descrip = sprintf('MPRAGE-like image, lambda %d',round(lambda)) ;
    for jj=1:N_MPRcreate % Looping on nr of images to create per lambda
        if jj>1 % add suffix for individual images
            fn_out_ii_jj = spm_file(fn_out_ii,'suffix', ...
                sprintf('-i%d',jj-1)); % add '_i1' or '_i2'
        else
            fn_out_ii_jj = fn_out_ii;
        end
        V_out.fname = fn_out_ii_jj;
        switch jj
            case 1
                spm_imcalc(V_in, V_out, ...
                    '(X(1,:)-lambda)./(mean(X(2:end,:))+lambda)', ...
                    ic_flags, lambda);
            case 2
                spm_imcalc(V_in, V_out, ...
                    '(X(1,:)-lambda)./(X(2,:)+lambda)', ...
                    ic_flags, lambda);
            case 3
                spm_imcalc(V_in, V_out, ...
                    '(X(1,:)-lambda)./(X(3,:)+lambda)', ...
                    ic_flags, lambda);                
        end
        
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
        % Collect output files
        fn_out_c{ii,jj} = fn_out_ii_jj;
    end
end

% Cleanup temporary folder & images for gorecstration, if used
if params.coreg
    rmdir(pth_tmp,'s');
end

% Collect output
fn_out = char(fn_out_c);

end

%% SUBFUNCTION
% Estimate the value of lambda from the input image intensities

function lambda = estimate_lambda(fn_in)

% Number of images 
Nfn_in = size(fn_in,1);

% Get "global" estimate, using SPM's function
v_global = spm_global(spm_vol(fn_in));

% Get the voxel values from all images
val_in = spm_read_vols(spm_vol(fn_in));
sz_img = size(val_in);
vval_in = reshape(val_in,[prod(sz(1:3)) sz(4)]);

% Create a "brain mask" of voxels > global across images
mask_glob = vval_in(:,1)>thr_global(1);
for ii=2:Nfn_in
    mask_glob = mask_glob | vval_in(:,ii)>thr_global(ii);
end

% Take lambda as the mean across images of the median values all 
% within-mask voxel values of each image.
median_vval = zeros(1,Nfn_in);
for ii=1:Nfn_in
    median_vval = median(vval_in(mask_glob,ii));
end
lambda = mean(median_vval);

end
