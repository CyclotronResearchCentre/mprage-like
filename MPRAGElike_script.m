% Small test scripts
% 
% Note: We assume that SPM is on Matlab path
%_______________________________________________________________________
% Copyright (C) 2025 Cyclotron Research Centre

% Written by C. Phillips, 
% Cyclotron Research Centre, University of Liege, Belgium

% Pick up some "old" test data
pth_data = 'D:\ccc_DATA\Test_MPRAGElike';
fn_T1w = spm_select('FPList',pth_data,'^sub.*-T1w_.*\.nii$');
fn_PDw = spm_select('FPList',pth_data,'^sub.*-PDw_.*\.nii$');
fn_MTw = spm_select('FPList',pth_data,'^sub.*-MTw_.*\.nii$');

fn_in =char(fn_T1w,fn_MTw);

