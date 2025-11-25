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
%            [false, def.]
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
    Nd_lambda = ceil(log10(max(params.lambda(:))));
    % Number of digits to write the lambdas
end

pth_in = spm_file(fn_in(1,:),'fpath');
pth_out = pth_in;
fn_basename = spm_file(fn_in(1,:),'basename');
fn_tmp = fullfile(pth_out, fn_basename);

% get the job done
ic_flags = struct('dmtx',true);

% Prepare output volume(s)
V_out = V_in(1);
V_out.dt(1) = 16; % use floats!
V_out.descrip = 'MPRAGE-like image';
% get the job done
for ii=1:Nlambda
   if Nlambda==1 % just one lambda
       fn_out_ii = [fn_tmp,'_MPRAGEl.nii'];
   else
       fn_out_ii = [sprintf('%s_l%4d',fn_tmp,round(params.lambda(ii))),'_MPRAGEl.nii'];
%        fn_out_ii =
%        sprintf(sprintf('%%s_l%%%dd',Nd_lambda),fn_tmp,params.lambda(ii));
%        Trying to add the right number of 0's
   end
   V_Allin = spm_vol(fn_in);
%    ic_flags.dtype = V_MRl.dt(1); % keeping the same data type
   ic_flags.dtype = 16; % use floats
   lambda = params.lambda(ii);
   spm_imcalc(fn_in,fn_out_ii,'(X(1,:)-lambda)./(mean(X(2:end,:),1)+lambda)',ic_flags,lambda)
    % Check thresholding
    if ~isempty(params.thresh) % apply thresholding
        vval_MPRlike(vval_MPRlike<params.thresh(1)) = 0;
        vval_MPRlike(vval_MPRlike>params.thresh(2)) = 0;
    end
    spm_write_vol(V_out,reshape(vval_MPRlike',sz(1:3)));
    
end

fn_out = fn_out_ii;

end