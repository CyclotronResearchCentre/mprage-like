% Small test scripts
% 
% Note: We assume that SPM is on Matlab path
%_______________________________________________________________________
% Copyright (C) 2025 Cyclotron Research Centre

% Written by C. Phillips, 
% Cyclotron Research Centre, University of Liege, Belgium

% Pick up some "old" test data
pth_data_top = 'D:\ccc_DATA\Test_MPRAGElike';

% pth_data_MS7T = fullfile( pth_data_top , 'MS7T');
% fn_T1w = spm_select('FPList',pth_data_MS7T,'^sub.*-T1w_.*_MPM\.nii$');
% fn_PDw = spm_select('FPList',pth_data_MS7T,'^sub.*-PDw_.*_MPM\.nii$');
% fn_MTw = spm_select('FPList',pth_data_MS7T,'^sub.*-MTw_.*_MPM\.nii$');

% IRONSLEEP data
pth_data_IRONS = fullfile( pth_data_top , 'IRONSLEEP');
fn_T1w = spm_select('FPList',pth_data_IRONS,'^sub.*-T1w_.*_MPM\.nii$');
fn_PDw = spm_select('FPList',pth_data_IRONS,'^sub.*-PDw_.*_MPM\.nii$');
fn_MTw = spm_select('FPList',pth_data_IRONS,'^sub.*-MTw_.*_MPM\.nii$');

% SCAIFIELD data
% CRC
pth_data_SF = fullfile( pth_data_top,'SCAIFIELD','CRC','Pilot-04_S1','NIfTI');
fn_T1w = spm_select('FPList',pth_data_SF,'^s.*-0012-.*-01-.*\.nii$');
fn_PDw = spm_select('FPList',pth_data_SF,'^s.*-0014-.*-01-.*\.nii$');
fn_MTw = spm_select('FPList',pth_data_SF,'^s.*-0016-.*-01-.*\.nii$');
% NTN
pth_data_SF = fullfile( pth_data_top,'SCAIFIELD','NTN','24_09_02_BVSS_C_10_RESCAN_MPM_02','NIfTI');
fn_T1w = spm_select('FPList',pth_data_SF,'^s.*-0005-.*-01-.*\.nii$');
fn_PDw = spm_select('FPList',pth_data_SF,'^s.*-0007-.*-01-.*\.nii$');
fn_MTw = spm_select('FPList',pth_data_SF,'^s.*-0009-.*-01-.*\.nii$');
% DZNE
pth_data_SF = fullfile( pth_data_top,'SCAIFIELD','DZNE','24_08_22-TQMRI_42644_1c','NIfTI');
fn_T1w = spm_select('FPList',pth_data_SF,'^s.*-0008-.*-01-.*\.nii$');
fn_PDw = spm_select('FPList',pth_data_SF,'^s.*-0006-.*-01-.*\.nii$');
fn_MTw = spm_select('FPList',pth_data_SF,'^s.*-0019-.*-01-.*\.nii$');

% fn_in =char(fn_T1w,fn_MTw);
fn_in =char(fn_T1w,fn_MTw,fn_PDw);

fn_out = hmri_MPRAGElike(fn_in)

params = struct(...
    'lambda', 100, ...
    'indiv', false, ...
    'thresh', [0 500], ...
    'coreg', false, ...
    'BIDSform', false);

params.lambda = [50 100 200]; % multiple lambda's
params.indiv = true;
params.coreg = true;

params.lambda = [NaN 100];
params.coreg = true;
fn_out = hmri_MPRAGElike(fn_in,params)

%% 
% Check "global signal" of all images
% __________________________________________________________________________
%  
%   spm_global returns the mean counts integrated over all the slices from
%   the volume.
%  
%   The mean is estimated after discounting voxels outside the object using
%   a criteria of greater than > (global mean)/8.
%  __________________________________________________________________________

v_global = spm_global(spm_vol(fn_in));
% v_global = % MS7T
%   218.5170 123.9108 149.8330 
% v_global = % IRONSLEEP
%   390.5394 479.7623 476.1641 
% v_global = % SCAIFIELD
%   716.4834 706.7789 756.9217 % CRC
%   929 1025 922 % NTN
%   909.6434 812.4403  839.2374 % DZNE

% Check range, mean/mode, histogram
val_in = spm_read_vols(spm_vol(fn_in));
sz = size(val_in)
vval_in = reshape(val_in,[prod(sz(1:3)) sz(4)]);
mean_all = mean(vval_in)
% mean_all =  % MS7T
%   194.5161  123.1467  146.8916
% mean_all =  % IRONSLEEP
%   192.3433  222.7650  229.1573
% mean_all =  % SCAIFIELD
%   342.9494  353.9259  356.7661 % CRC
%   509.4290  582.0244  505.3355 % NTN
%   477.6616  423.5554  429.9221 % DZNE

thr_global = mean_all*.8
% thr_global =  % MS7T
%   155.6129   98.5174  117.5133
% -> quite OK but some extra bits for T1w...
% thr_global = % IRONSLEEP
%   153.8746  178.2120  183.3258
% thr_global = % SCAIFIELD
%   274.3595  283.1407  285.4129 % CRC
%   407.5432  465.6195  404.2684 % NTN
%   382.1293  338.8443  343.9377 % DZNE

% Crete a mask based on the union of the 3 thresholded images
fn_mask = spm_file(fn_in(1,:),'suffix','_mask');
mask_glob = vval_in(:,1)>thr_global(1);
for ii=2:numel(thr_global)
    mask_glob = mask_glob | vval_in(:,ii)>thr_global(ii);
end
Vmsk = spm_vol(fn_in(1,:));
Vmsk.fname = fn_mask;
Vmsk.dt(1) = 2;
spm_write_vol(Vmsk,reshape(mask_glob,sz(1:3)))

spm_check_registration(char(fn_in,fn_mask))

% Look at histograms
figure, hist(vval_in,100)
figure,
% all values
subplot(3,2,1), histogram(vval_in(:,1))
subplot(3,2,3), histogram(vval_in(:,2))
subplot(3,2,5), histogram(vval_in(:,3))
% only in mask
subplot(3,2,2), histogram(vval_in(mask_glob,1))
subplot(3,2,4), histogram(vval_in(mask_glob,2))
subplot(3,2,6), histogram(vval_in(mask_glob,3))



% Check median values
% All voxels
median_all = median(vval_in)
% median_all = % MS7T
%     45    36    38 % -> pretty useless !
% median_all = % IRONSLEEP 
%     22    20    24 % -> pretty useless !
% median_all =  % SCAIFIELD
%     34    41    34 % CRC
%     90   133    87 % NTN
%     67    57    54 % DZNE

mediam_msk = [median(vval_in(mask_glob,1)) median(vval_in(mask_glob,2)) median(vval_in(mask_glob,3))]
% mediam_msk = % MS7T
%    552   334   436 % Seems very interesting!
% mediam_msk = % IRONSLEEP 
%    470   553   581
% mediam_msk = % SCAIFIELD
%          963         956        1054 % CRC
%         1055        1275        1095 % NTN
%         1130         970        1042 % DZNE

mean(mediam_msk)
% 534.6667 % IRONSLEEP 
% SCAIFIELD -> take it divided by 10 ???
% 991 % CRC
% 1142 % NTN
% 1047 % DZNE

% Check range as percentile
prctile(vval_in,[5 95])
%     21    20    20 % MS7T
%    813   525   650
%            2           2           2 % IRONSLEEP 
%          776         982        1001
%            5           5           5 % SCAIFIELD
%         1555        1525        1621 % CRC
%            5           5           5 % NTN
%         2016        2071        1983
%            5           5           5 % DZNE
%         1901        1799        1847
        
prctile([vval_in(mask_glob,1) vval_in(mask_glob,2) vval_in(mask_glob,3)],[5 95])
%          162          64          66 % MS7T
%         1033         685         848
%          154          80          86 % IRONSLEEP 
%          971        1256        1275
%          151         286         150 % SCAIFIELD
%         1940        1851        1988 % CRC
%            5           5           5 %NTN
%         2016        2071        1983
%          341         190         173 % DZNE
%         2394        2360        2334

%% Lambda validation ?
% Check Structural Similarity Index Metric (SSIM)

