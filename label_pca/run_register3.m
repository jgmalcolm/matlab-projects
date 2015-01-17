function run_register3(pose)
  
  % use first map
%   maps = load('images/3D', 'maps');
%   maps = maps.maps;
%   m = maps{1};
%   clear maps;
%   maps{1} = m;
%   
%   % transform to get second map
%   label_cnt = numel(unique(maps{1}));
%   lmap = label_map(maps{1}, label_cnt);
%   lmap = reshape(lmap', [size(maps{1}) label_cnt-1]);
%   maps{2} = f_inv(transform3(lmap, pose));
  
   %for loading the NRRD reader
  %addpath('/projects/schiz/guest/kquintus/scripts/matlab/');
  addpath('/projects/schiz/software/Matlab'); 
  addpath /projects/schiz/pi/marc/matlab/nrrdUtilities/;
  
   %load the first map
   %vol1 = nrrdZipLoad('/projects/schiz/strct/fe/mclean/ROI/total_vol/ROI_final_jeong/case00747/case00747_fROI_resample.nhdr');
   %vol2 = nrrdZipLoad('/projects/schiz/strct/fe/mclean/ROI/total_vol/ROI_final_jeong/case00749/case00749_fROI.nhdr');
   files = dir('../../../labelPCA/labelRegistration/relabeled_jeong/*.nhdr');
   
%    vol1 = nrrdZipLoad('../../../labelPCA/labelRegistration/relabeled_jeong/case51_fROI.nhdr');
%    vol2 = nrrdZipLoad('../../../labelPCA/labelRegistration/relabeled_jeong/case3_fROI.nhdr');
%    vol3 = nrrdZipLoad('../../../labelPCA/labelRegistration/relabeled_jeong/case0581_fROI.nhdr');
% 
%    vol1(vol1~=0)=1;
%    vol2(vol2~=0)=1;
%    vol3(vol3~=0)=1;
%    maps{1} = vol1;
%    maps{2} = vol2;
%    maps{3} = vol3;
%    clear vol1;
%    clear vol2;
%    clear vol3;

   for i=1:length(files)
       fname = sprintf('../../../labelPCA/labelRegistration/relabeled_jeong/%s',files(i).name);
       maps{i} = nrrdZipLoad(fname)+1;
       %vol = maps{i};
       %vol(vol~=0)=1;
       %maps{i} = vol;
   end
   clear vol;
   
  % register
  [maps_ pose_] = register3(maps);

  for i =1:length(files)
      fname = sprintf('../../../labelPCA/labelRegistration/jeong_registered/%s',files(i).name);
      nrrdSave(fname,maps_{i});
  end
  fname = sprintf('../../../labelPCA/labelRegistration/jeong_registered/mean_map.nhdr');
  nrrdSave(fname,maps_{i+1});
  
%   nrrdSave('../../../labelPCA/labelRegistration/jeong_registered/case51_fROI.nhdr',maps_{1});
%   nrrdSave('../../../labelPCA/labelRegistration/jeong_registered/case3_fROI.nhdr',maps_{2});
%   nrrdSave('../../../labelPCA/labelRegistration/jeong_registered/case0581_fROI.nhdr',maps_{3});
%   save '../../../labelPCA/labelRegistration/jeong_registered/51_3_pose.mat' 'pose_';
end

function map = f_inv(lmap)
  sz = size(lmap);
  if numel(sz)==3, sz(4) = 1; end
  lmap = reshape(lmap, [prod(sz([1:3])) sz(4)]);
  lmap = lmap'; % TODO: rearrange these lines
  label_cnt = size(lmap,1) + 1;
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:3)));
end
