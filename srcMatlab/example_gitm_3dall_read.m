% close all
% clear all
% clc


%% =======================================

%% =======================================
function [utime_all,lat_all,lon_all,alt_all,Tn_all,rho_all,VnUp_all]=GITM_3DALL_read(fpath)

% fpath='/Users/xianjing/Desktop/no_backup/Gitm/Runs/Variability/by-05+0.0/data/';

% % datasub includes all the model output file  size (39x76x76x54)
% % the name of the the parameter is in vars

fnames=dir([fpath '3DALL*.bin']);

% data_all=[];
utime_all=[];
on2_all=[];
n2_all=[];
o_all=[];

for flen_index=1:length(fnames)-1
%  for flen_index=20:30
    
%      flen_index1=flen_index
    
fnames_each=[fpath fnames(flen_index).name ];
fileID = fopen(fnames_each);

% version=3.13
A = fread(fileID,[1 1],'long');
version = fread(fileID,[1 1],'double');
A = fread(fileID,[1 1],'long');

% nLons=76, nLats=76, nAlts=54
A = fread(fileID,[1 1],'long');
nLons = fread(fileID,[1 1],'long');
nLats = fread(fileID,[1 1],'long');
nAlts = fread(fileID,[1 1],'long');
A = fread(fileID,[1 1],'long');

% nVars=39
A = fread(fileID,[1 1],'long');
nVars = fread(fileID,[1 1],'long');
A = fread(fileID,[1 1],'long');

% 
vars = [];
for iVar=1:nVars

    A = fread(fileID,[1 1],'long');
    Var1 = char(fread(fileID,[40],'char'));
    A = fread(fileID,[1 1],'long');
% Var1 is a matrix, make Var1 into string
    svar=[];
    for i=1:length(Var1)
        svar=[svar Var1(i)];
    end
 % put the name of 39 variables inside a matrix  
    vars = [vars;svar];
end

A = fread(fileID,[1 1],'long');
iTime = fread(fileID,[7],'long');
iTime=iTime';
A = fread(fileID,[1 1],'long');

% iTime
doy=datenum(iTime(1), iTime(2), iTime(3))-datenum(iTime(1), 0, 0);
utime=doy+iTime(4)/24+iTime(5)/60/24+iTime(6)/60/60/24;
utime_all=[utime_all utime];

for iVar=1:nVars
    A = fread(fileID,[1 1],'long');
    for iAlt = 1:nAlts
        tmp = fread(fileID,[nLons nLats],'double');
        datasub(iVar,:,:,iAlt) = tmp;
    end
    A = fread(fileID,[1 1],'long');
end
fclose(fileID);

% datasub includes all the model output file  size (39x76x76x54)
% the name of the the parameter is in vars

lon = squeeze(datasub(1,:,:,:));
lon_deg=lon/pi*180;
lat = squeeze(datasub(2,:,:,:));
lat_deg=lat/pi*180;
alt = squeeze(datasub(3,:,:,:));
Tn = squeeze(datasub(16,:,:,:));
rho = squeeze(datasub(4,:,:,:));
VnUp = squeeze(datasub(19,:,:,:));




lat_all(flen_index,:,:,:)=lat_deg(:,:,:); 
lon_all(flen_index,:,:,:)=lon_deg(:,:,:); 
alt_all(flen_index,:,:,:)=alt(:,:,:); 
Tn_all(flen_index,:,:,:)=Tn(:,:,:); 
rho_all(flen_index,:,:,:)=rho(:,:,:); 
VnUp_all(flen_index,:,:,:)=VnUp(:,:,:); 


end



