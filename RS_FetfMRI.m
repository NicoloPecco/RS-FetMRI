%% Article: ....


%% Code Description

% RS-FetfMRI semiautomatic processing (RS-fMRIsp) is an open-source package (source code available at http://.....)
% totally integrated in SPM (Statistical Parametric Maps). 

% It is supported by MacOS, Linux and Windiws operating systems if MATLAB is available. RS-fMRIsp is divided into 6 modules where the first
% three, from M1 to M3, works ‘Within Session’ (WS) while the last four, from M4 to cM7, work ‘Between Session’ (BS).

%% Requirements:

%      1) Matlab2013 or above;
%      2) SPM12 or above;
%      3) RS-FetfMRI package: 'art', 'Template' folder and
%                             'Create_Template.m' Matlab function

%% Initialization:

%      1) Set path for the ‘art’ folder in Matlab;
%      2) Inside the ‘art’ folder there are three ‘.cfg’ text files. In
%         these files, three paths need to be modified (‘image_dir’,
%         ‘motion_dir’ and ‘mask_file’): 

%               - File_1Srub.cfg:  Update each path maintaining this
%                 folder at the end:

%                       1. image_dir: /M2_WS_03_Scrub/9999/ 
%                       2. motion_dir: /M2_WS_03_Scrub/9999/ 
%                       3. mask_file: /M2_WS_03_Scrub/9999/meanm9999_001.nii

%               - File_2Scrb_orig.cfg AND File_2Scrb.cfg:  Update each path but maintain this folder at the end:

%                       1. image_dir: /M4_BS_02_Scrub/ 
%                       2. motion_dir: /M4_BS_02_Scrub/Motion_file/  
%                       3. mask_file: /M4_BS_02_Scrub/9999/r2mr1m9999_0001.nii   


%% Follow the MANUAL for more detailed code explanation.


%% INITIALIZATION: Path creation and Setting Paths

 close all
 clc
 clear all

% General Path set up:
global path_general;

path=cd();
path_to_cffg_file=[path,'/art/'];
path_general=[path,'/'];

addpath(path_to_cffg_file);
cd(path_general)
folder2create={'M1_PP_01_OrigVol','M1_PP_02_4Dto3D','M1_PP_03_Reorient','M1_PP_04_Rename','M2_WS_01_Mask','M2_WS_02_Realign','M2_WS_03_Scrub','M2_WS_04_Rename','M3_WS_01_SegRefVols','M3_WS_02_MaskRefVols','M3_WS_03_RealignMaskRefVols','M3_WS_04_Rename','M4_BS_01_Realign_Reslice','M4_BS_02_Scrub','M4_BS_03_Rename','M5_BS_01_SegMeanRefVol','M5_BS_02_MaskMeanRefVol','M5_BS_03_NormMaskMeanRefVol','M6_BS_01_MaskAllVols','M6_BS_02_NormMaskAllVols'};
for num_folder=1:size(folder2create,2)
    if ~exist(folder2create{1,num_folder}, 'dir')
           mkdir(folder2create{1,num_folder});
    else
        disp([folder2create{1,num_folder},' already exists.'])
    end
end

% Paths name creation:
%Module1
path_to_M1_PP_01_OrigVol=[path_general,'M1_PP_01_OrigVol'];
path_to_M1_PP_02_4Dto3D=[path_general,'M1_PP_02_4Dto3D'];
path_to_M1_PP_03_Reorient=[path_general,'M1_PP_03_Reorient'];
path_to_M1_PP_04_Rename=[path_general,'M1_PP_04_Rename'];
%Module2
path_to_M2_WS_01_Mask=[path_general,'M2_WS_01_Mask'];
path_to_M2_WS_02_Realign=[path_general,'M2_WS_02_Realign'];
path_to_M2_WS_03_Scrub=[path_general,'M2_WS_03_Scrub'];
path_to_M2_WS_04_Rename=[path_general,'M2_WS_04_Rename'];
%Module3
path_to_M3_WS_01_SegRefVols=[path_general,'M3_WS_01_SegRefVols'];
path_to_M3_WS_02_MaskRefVols=[path_general,'M3_WS_02_MaskRefVols'];
path_to_M3_WS_03_RealignMaskRefVols=[path_general,'M3_WS_03_RealignMaskRefVols'];
path_to_M3_WS_04_Rename=[path_general,'M3_WS_04_Rename'];
%Module4
path_to_M4_BS_01_Realign_Reslice=[path_general,'M4_BS_01_Realign_Reslice'];
path_to_M4_BS_02_Scrub=[path_general,'M4_BS_02_Scrub'];
path_to_M4_BS_03_Rename=[path_general,'M4_BS_03_Rename'];
%Module5
path_to_M5_BS_01_SegMeanRefVol=[path_general,'M5_BS_01_SegMeanRefVol'];
path_to_M5_BS_02_MaskMeanRefVol=[path_general,'M5_BS_02_MaskMeanRefVol'];
path_to_M5_BS_03_NormMaskMeanRefVol=[path_general,'M5_BS_03_NormMaskMeanRefVol'];
%Module6
path_to_M6_BS_01_MaskAllVols=[path_general,'M6_BS_01_MaskAllVols'];
path_to_M6_BS_02_NormMaskAllVols=[path_general,'M6_BS_02_NormMaskAllVols'];
%Template
path_original_tmp=[path_general,'Template/Template_for_session'];

cd(path_general)
save('All_paths','path_general','path_original_tmp','path_to_cffg_file','path_to_M1_PP_01_OrigVol','path_to_M1_PP_02_4Dto3D','path_to_M1_PP_03_Reorient','path_to_M1_PP_04_Rename','path_to_M2_WS_01_Mask','path_to_M2_WS_02_Realign','path_to_M2_WS_03_Scrub','path_to_M2_WS_04_Rename','path_to_M3_WS_01_SegRefVols','path_to_M3_WS_02_MaskRefVols','path_to_M3_WS_03_RealignMaskRefVols','path_to_M3_WS_04_Rename','path_to_M4_BS_01_Realign_Reslice','path_to_M4_BS_02_Scrub','path_to_M4_BS_03_Rename','path_to_M5_BS_01_SegMeanRefVol','path_to_M5_BS_02_MaskMeanRefVol','path_to_M5_BS_03_NormMaskMeanRefVol','path_to_M6_BS_01_MaskAllVols','path_to_M6_BS_02_NormMaskAllVols');
spm_input('Please copy your original volumes folder into the folder -M1_PP_01_OrigVol- and press yes to start the processing:','-1','bd','Yes');

%% MODULO 1:
disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||       //||  ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||          // ||  ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||         ||  ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||             ||  ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||       __||__')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance imaging fMRI -------------------')
disp('-------------------------- Within Session Processing  ---------------------------')
disp('---------------------From original session to reorientation----------------------')
disp('---------------------------------------------------------------------------------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                       .\( / / . ) _)_)_______/                                  ')
disp('                        .( ( (  (__/__/  .                                       ')
disp('                          \._)_). )                                              ')
disp('                           ._)_)./                                               ')
disp('                           ._)_)./                                               ')
disp('                                                                                 ')


%% Creation of folders for each session. Renaming nifti files with default session name inside M1_PP_01_OrigVol.
cd(path_general)
diary Module1_logfile

disp('First of all, Good Day! Thanks for using this program! ')
disp('Creating folder and Renaming 4D file with default name.')

cd(path_to_M1_PP_01_OrigVol);
list=dir();
% Get vector to recognize folders only:
list(ismember({list.name},{'.','..'}))=[];
Session_name = cell(length(list),1);
for i=1:length(list)
       Session_name{i,1}=num2str(i,'%02d');
       Session_name{i,1}=['10',Session_name{i,1}];
       disp(['Default name for session ',num2str(i,'%01d'),' : ',Session_name{i,1}]);
end

for i=1:size(Session_name,1)
    cd([path_to_M1_PP_01_OrigVol,'/',list(i).name])
    name_orig=dir('*.nii');
    name_orig(ismember({name_orig.name},{'.','..'}))=[];
    if size(name_orig,1)==1
        copyfile([path_to_M1_PP_01_OrigVol,'/',list(i).name,'/',name_orig.name],[path_to_M1_PP_02_4Dto3D,'/',Session_name{i,1},'.nii'])
    else
     spm_input(['Error: More than one .nii detected in folder (',list(i).name,') . Inside each folder there must be only one Nifti file. Check the name_orig variable and rerun the code. '],'-1','bd','Yes');   
     disp('More than one .nii detected in folder. Inside each folder there must be only one Nifti file. See the name_orig variable.');
     return
    end
     disp(['Naming for Session ',Session_name{i,1},' completed!!'])
end

%% Converting session from 4D to 3D; From M1_PP_01_OrigVol to M1_PP_02_4Dto3D
disp('Converting session from 4D to 3D using SPM:')

for i=1:size(Session_name,1)
    cd(path_to_M1_PP_02_4Dto3D)
    mkdir(Session_name{i,1})
end

% 4D to 3D using SPM:
clear matlabbatch
for i=1:size(Session_name,1)
    input=[path_to_M1_PP_02_4Dto3D,'/',Session_name{i,1},'.nii'];
    output=[path_to_M1_PP_02_4Dto3D,'/',Session_name{i,1}];
    matlabbatch{1}.spm.util.split.vol = {input};
    matlabbatch{1}.spm.util.split.outdir = {output};
    spm_jobman('run',matlabbatch);
    disp(['Session ',Session_name{i,1},' 4D to 3D conversion completed!!'])
end

% Volume for session:
for i=1:size(Session_name,1)
    clear list_tmp
    list_tmp=length(dir([path_to_M1_PP_02_4Dto3D,'/',Session_name{i,1},'/',Session_name{i,1},'*.nii']));
    vol4session_tmp(i)=list_tmp;
    vol4session=min(vol4session_tmp);
end

%% From M1_PP_02_4Dto3D to M1_PP_03_Reorient

disp('Copying files from M1_PP_02_4Dto3D to M1_PP_03_Reorient')
%disp('This step is one of the most important..')

for i=1:size(Session_name,1)
    from=[path_to_M1_PP_02_4Dto3D,'/',Session_name{i,1}];
    to=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1}];
    copyfile(from,to)
end

cd(path_to_M1_PP_03_Reorient)
file2show = cell(6,1);
REF_vol=cell(1,size(Session_name,1));
for i=1:size(Session_name,1)
    clear input file2show
    for j=1:6
        %file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(j*5+(j-1)*5,'%05d'),'.nii'];
        file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(((floor(vol4session/6)*j)),'%05d'),'.nii'];
        answ_tmp(j)=((floor(vol4session/6)*j));
    end
    spm_check_registration(file2show{1,:},file2show{2,:},file2show{3,:},file2show{4,:},file2show{5,:},file2show{6,:});
    answ={[num2str(answ_tmp(1)),'|',num2str(answ_tmp(2)),'|',num2str(answ_tmp(3)),'|',num2str(answ_tmp(4)),'|',num2str(answ_tmp(5)),'|',num2str(answ_tmp(6))]};
    promt=['Choose the Reference for session ',Session_name{i,1}, ' : '];
    spm_input(['This step is one of the most important!  Six images will be displayed per session: ',Session_name{i,1},'. 1) Choose a reference image; 2) Reorient the chosen reference volume *be sure to set the origin at the level of Anterior Commissure* (See Manual); 3) Apply the transformation parameters to all volumes of the current session.'],'-1','bd','Yes');
    REF_vol1=spm_input(promt,'-1','m',answ,[answ_tmp(1),answ_tmp(2),answ_tmp(3),answ_tmp(4),answ_tmp(5),answ_tmp(6)]);
    REF_vol(i)={num2str(REF_vol1,'%05d')};
    spm_image('Display',[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',REF_vol{1,i},'.nii']);

         if i== size(Session_name,1)
                    try
                    spm_input('Reorienting Ref. image:','-1','b','Continue|--',[1,0],1);
                    catch
                    end
         else
                    try
                    spm_input('Reorient and press continue::','-1','b','Continue|--',[1,0],1);
                    catch
                    end
         end
end

for i=1:size(Session_name,1)
    clear input file2show answ_tmp
    for j=1:6 
        %file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(j*5+(j-1)*5,'%05d'),'.nii'];
        file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(((floor(vol4session/6)*j)),'%05d'),'.nii'];
    end
    spm_check_registration(file2show{1,:},file2show{2,:},file2show{3,:},file2show{4,:},file2show{5,:},file2show{6,:});
    if i==1
    spm_input(['These are the reoriented images. Here you have the chance to see and, if needed, to return to the reference image orientation section for the session: ',Session_name{i,1}],'-1','bd','Yes');
    end
    Promt='Check reoriented images:';
    m = spm_input(Promt,'-1','Continue|Reorient',['Y','N']);
    if m=='Y'
        continue
        else
        clear input file2show
        for j=1:6
            %file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(j*5+(j-1)*5,'%05d'),'.nii'];
            file2show{j,1}=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(((floor(vol4session/6)*j)),'%05d'),'.nii'];
            answ_tmp(j)=((floor(vol4session/6)*j));
        end
        answ={[num2str(answ_tmp(1)),'|',num2str(answ_tmp(2)),'|',num2str(answ_tmp(3)),'|',num2str(answ_tmp(4)),'|',num2str(answ_tmp(5)),'|',num2str(answ_tmp(6))]};
        spm_check_registration(file2show{1,:},file2show{2,:},file2show{3,:},file2show{4,:},file2show{5,:},file2show{6,:});
        promt=['Choose the Reference for session ',Session_name{i,1}, ' : '];
        REF_vol1=spm_input(promt,'-1','m',answ,[answ_tmp(1),answ_tmp(2),answ_tmp(3),answ_tmp(4),answ_tmp(5),answ_tmp(6)]);
        REF_vol(i)={num2str(REF_vol1,'%05d')};
        spm_image('Display',[path_to_M1_PP_03_Reorient,'/',Session_name{i,1},'/',Session_name{i,1},'_',REF_vol{1,i},'.nii']);

         if i== size(Session_name,1)
                    try
                    spm_input('Reorienting Ref. image:','-1','b','Continue|--',[1,0],1);
                    catch
                    end
         else
                    try
                    spm_input('Reorient and press continue::','-1','b','Continue|--',[1,0],1);
                    catch
                    end
         end
    end
end

disp('Copying files from M1_PP_03_Reorient to M1_PP_04_Rename')

for i=1:size(Session_name,1)
    disp(['Copying files from M1_PP_03_Reorient to M1_PP_04_Rename for session: ',Session_name{i,1}])
    from=[path_to_M1_PP_03_Reorient,'/',Session_name{i,1}];
    to=[path_to_M1_PP_04_Rename,'/',Session_name{i,1}];
    copyfile(from,to)
end
disp('Copying files from M1_PP_03_Reorient to M1_PP_04_Rename Completed')

%% Placing reference Volume of each session at first place 

for i=1:size(Session_name,1)
    f = [path_to_M1_PP_04_Rename,'/',Session_name{i,1}];
    cd(f)
    ref_path = [path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/',Session_name{i,1},'_',REF_vol{1,i},'.nii'];
    new_name=[path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/',Session_name{i,1},'_00000.nii'];
    movefile(ref_path,new_name); % chiamo la ref 00015.nii come 'Ref.nii'
    files = dir('*.nii');
    ii=1;
          for id = 1:length(files)
             num=[Session_name{i,1},'_',num2str(id,'%03d'),'.nii'];
             movefile([path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/',files(ii).name],[path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/',num]);
             ii=ii+1;
          end
        clear files
        disp(['Reference Volume: ',REF_vol{1,i}, ' was placed at the beginning of the Session: ', Session_name{i,1},'. Completed!!'])
end

% Saving Variables of module 1:
cd(path_general)
save('Variables_M1','REF_vol','Session_name','vol4session')
diary off
disp('End of Module 1!')

spm_input('Orientation of all of the Fetal Functional images completed! Proceed with 1st-pass Masking, Realignment and 1st-pass Scrubbing procedure in module 2.','-1','bd','Ok!');

%% Module 2

disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||      //-\\  ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||         //  //  ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||        //   ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||           //    ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||      //___  ')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance Imaging fMRI -------------------')
disp('-------------------------- Within Session Processing ----------------------------')
disp('---From 1-pass Masking with session-specific reoriented fetal-specific mask -----')
disp('-------------to Realignment and 1st pass Scrubbing procedure with ART------------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                       .\( / / . ) _)_)_______/                                  ')
disp('                        .( ( (  (__/__/  .                                       ')
disp('                          \._)_). )                                              ')
disp('                           ._)_)./                                               ')
disp('                           ._)_)./                                               ')
disp('                                                                                 ')

%% Copying files inside fMRI masking
disp('Copying files from M1_PP_04_Rename to M2_WS_01_Mask.')

for i=1:size(Session_name,1)
    from=[path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/*.nii'];
    to=[path_to_M2_WS_01_Mask,'/',Session_name{i,1}];
    copyfile(from,to)
    cd(to)
    mkdir([Session_name{i,1},'_template'])
    mkdir([Session_name{i,1},'template_anp'])
end

%% Creating and centering all template for each session:
cd(path_general)
diary Module2_logfile

spm_input(['Creating session-specific template mask in the user image space for each session (voxel resolution, image dimensions and orientation). Template creation will take approximately ',num2str((((size(Session_name,1)*110))*2)/60),' minutes (for 1.5T).'],'-1','bd','Yes');

cd(path_general)
Create_Template(path_original_tmp,path_to_M1_PP_04_Rename,path_to_M2_WS_01_Mask,Session_name);

% Translation of template
num_template=17;
for i=1:size(Session_name,1)
    disp(['Template mask translation for Session: ',Session_name{i,1}])
    %Calculating origin for reference Volume:
        VV=spm_vol([path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_001.nii']);
        Sess_origin = VV.mat\[0 0 0 1]';

    for j=1:num_template

    %Calculating origin for each template:
        clear input output option cmd matlabbatch cmdout_tmp

        V=spm_vol([path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/scGW_',num2str(j+20),'.nii']);
        temp_origin = V.mat\[0 0 0 1]';

     % Calculating origins differences:
        x_diff=(Sess_origin(1,1) - temp_origin(1,1));
        y_diff=(Sess_origin(2,1) - temp_origin(2,1));
        z_diff=(Sess_origin(3,1) - temp_origin(3,1));

     %Translating original templates:
        matlabbatch{1}.spm.util.reorient.srcfiles = {[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/scGW_',num2str(j+20),'.nii']};%{}
        matlabbatch{1}.spm.util.reorient.transform.transM = [1 0 0 x_diff
                                                             0 1 0 y_diff
                                                             0 0 1 z_diff
                                                             0 0 0 1];
        matlabbatch{1}.spm.util.reorient.prefix = Session_name{i,1};
        spm_jobman('run',matlabbatch);

      %Translating smoothed templates:
        for jjj=1:3
        matlabbatch{1}.spm.util.reorient.srcfiles = {[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'template_anp/banp_',num2str(jjj),'_',num2str(j+20),'.nii']};%{}
        matlabbatch{1}.spm.util.reorient.transform.transM = [1 0 0 x_diff
                                                             0 1 0 y_diff
                                                             0 0 1 z_diff
                                                             0 0 0 1];
        matlabbatch{1}.spm.util.reorient.prefix = [Session_name{i,1}];
        spm_jobman('run',matlabbatch);
        end
    end
    cd([path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template'])
    delete s*
    cd([path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'template_anp'])
    delete ban*

    disp(['Template masks translation for Session ',Session_name{i,1},' Completed!'])

end

%% Insert Gestational Week and choose one Template:
template2show = cell(1);

clear input
k=1;
    for j=21:37
        template2show{k,1}=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_template/',Session_name{1,1},'scGW_',num2str(j),'.nii'];
        k=k+1;
    end
ref2show=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_001.nii'];
spm_check_registration(ref2show,template2show{:,:});

spm_input(['The chosen session-specific template mask should include the entire fetal brain and a limited quantity of the surrounding anatomical structures (black edges). Please use the zoom function for a more accurate choice.'],'-1','bd','Yes');
GW2beused=spm_input('Choose the session-specific template mask:','-1','m','21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36|37',[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37],30);

% First choice
clear input template2show ref2show temp2show
template2show = cell(1);
    for j=1:3
        template2show{j,1}=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'template_anp/',Session_name{1,1},'banp_',num2str(j),'_',num2str(GW2beused),'.nii'];
    end
ref2show=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_001.nii'];
temp2show=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_template/',Session_name{1,1},'scGW_',num2str(GW2beused),'.nii'];
spm_check_registration(ref2show,temp2show,template2show{:,:});

% Final choice

spm_input('For the session-specific template mask you have chosen, select the one that covers the whole brain in each direction (es. L/R, A/P; I/S).','-1','bd','Ok!');
template2bused=spm_input('Insert a number between 1 and 4: ','-1','m','1|2|3|4',[1,2,3,4]);

%% Applying the session-specific template mask ( 1-pass Masking):

disp('Masking images with the chosen session-specific template mask (1-pass Masking): ');

clear matlabbatch list_to_mask
for i=1:size(Session_name,1)
    cd([path_to_M2_WS_01_Mask,'/',Session_name{i,1}]);
    list_to_mask=dir('*.nii');
    for id=1:length(list_to_mask)
        clear matlabbatch list_to_mask
        matlabbatch{1}.spm.util.imcalc.input{1,1} = char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_',num2str(id,'%03d'),'.nii,1'));
        if template2bused==1
             matlabbatch{1}.spm.util.imcalc.input{2,1} = char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/',Session_name{i,1},'scGW_',num2str(GW2beused),'.nii,1'));
        else
             matlabbatch{1}.spm.util.imcalc.input{2,1} = char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'template_anp/',Session_name{i,1},'banp_',num2str(template2bused-1),'_',num2str(GW2beused),'.nii,1'));
        end
        matlabbatch{1}.spm.util.imcalc.output = char(strcat('m',Session_name{i,1},'_', num2str(id,'%03d'), '.nii'));
        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1})};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
    end
    tmp=dir([Session_name{i,1},'*nii']);
    delete(tmp.name);
    clear tmp
end

for i=1:size(Session_name,1)
clear input  ref2show file2show temp2beshow
       for j=1:6
            %file2show{j,1}=[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/m',Session_name{i,1},'_',num2str(j*5+(j-1)*5,'%03d'),'.nii'];
            file2show{j,1}=[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/m',Session_name{i,1},'_',num2str(((floor(vol4session/6)*j)),'%03d'),'.nii'];
       end
        ref2show=[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/m',Session_name{i,1},'_001.nii'];
        nomask2show=[path_to_M1_PP_04_Rename,'/',Session_name{i,1},'/',Session_name{i,1},'_001.nii'];
        if template2bused==1
             temp2beshow=char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/',Session_name{i,1},'scGW_',num2str(GW2beused),'.nii,1'));
        else
             temp2beshow= char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'template_anp/',Session_name{i,1},'banp_',num2str(template2bused-1),'_',num2str(GW2beused),'.nii,1'));
        end
        spm_check_registration(ref2show,nomask2show,temp2beshow,file2show{:,:});
        if i<size(Session_name,1)
            spm_input('1st-pass mask review:','-1','b','Continue|--',[1,0]);
        else
            spm_input('1st-pass mask review:','-1','b','Proceed|--',[1,0]);
        end
end

%% From 05_0_fMRI_Masking to 05_fMRI_Realign_to_First_Within_Session
disp('Copying files from M2_WS_01_Mask to M2_WS_02_Realign.')

for i=1:size(Session_name,1)
    from=[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/m*'];
    to=[path_to_M2_WS_02_Realign,'/',Session_name{i,1}];
    copyfile(from,to)
end

%% Realignment WS to first volume with SPM:
clear answ_tmp answ
% Visualization 
whg2show=cell(GW2beused-21,1);
k=1;
clear answ
for j=21:GW2beused
     whg2show{k,1}=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_template/',Session_name{1,1},'scGW_',num2str(j),'.nii'];
     answ_tmp(k)=j;
     k=k+1;
end
ref2show=[path_to_M2_WS_02_Realign,'/',Session_name{1,1},'/m',Session_name{1,1},'_001.nii'];
spm_check_registration(ref2show,whg2show{:,:});
answw='';
for hh=1:length(answ_tmp)
answ=num2str(answ_tmp(hh));
    if hh<length(answ_tmp)
        answw=[answw,answ,'|'];
    else
        answw=[answw,answ];
    end
end

spm_input('Now you should select the Tissue-Weighting mask that will be used to calculate Within Session motion parameters in the Realignment algorithm. Therefore it should completely fit within the brain. ','-1','bd','Ok!');
wgh_winner=spm_input('Choose the Tissue-Weighting mask: ','-1','m',answw,answ_tmp);
wgh_winner=num2str(wgh_winner);
wgh_winner=cellstr(wgh_winner);

clear matlabbatch list list_tmpp path_tmp list_tmpp lista_volumi path_folder_tmp

        for ii=1:size(Session_name,1)
            clear matlabbatch list 
            path_folder_tmp=[path_to_M2_WS_02_Realign,'/',Session_name{ii,1},'/']; %Prima cartella
            cd(path_folder_tmp);
            lista_volumi=dir('*.nii');u=0;% lista volumi nella cartella
                for j=1:length(lista_volumi)
                    TTFF = length(regexpi(lista_volumi(j).name,'nii'));
                    if TTFF
                      u=u+1;
                      list_tmpp(u,:)=lista_volumi(j).name;
                    end
                end
        list=cellstr(list_tmpp);
            for iii=1:length(list)
            path_folder_ttmp(iii,:)=[path_folder_tmp,list(iii)];% Creating path for each volume
            path_tmp(iii,:)=[path_folder_ttmp{iii,:}];
            end
            path=cellstr(path_tmp); %All session volumes for realignment

% Realign: Estimate
        clear matlabbatch

        matlabbatch{1}.spm.spatial.realign.estwrite.data = {path}';
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {[path_to_M2_WS_01_Mask,'/',Session_name{ii,1},'/',Session_name{ii,1},'_template/',Session_name{ii,1},'scGW_',wgh_winner{1,1},'.nii,1']};
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r1';
        spm_jobman('run',matlabbatch);

        clear list path_tmp list_tmpp lista_volumi path_folder_tmp
        disp(['Realignment of ',Session_name{ii,1},' volumes to reference volume completed!!'])
        end

%Renaming:

        for i=1:size(Session_name,1)
            cd([path_to_M2_WS_02_Realign,'/',Session_name{i,1}]);
            list=dir(['m',Session_name{i,1},'*']);
            for ii=1:size(list,1)
                current_name=list(ii).name;
                movefile([path_to_M2_WS_02_Realign,'/',Session_name{i,1},'/',current_name],[path_to_M2_WS_02_Realign,'/',Session_name{i,1},'/r1',current_name]);
            end
        end

spm_input('Within each Session folder found in the M2_WS_02_Realign folder you can find the: 1) output of each realigned volume with -r1- as the prefix; 2) mean volume calculated from the average of all session volumes; 3) a txt file (rp*.txt) and a ps file (*.ps) reporting translation and rotation motion parameters and charts of each volume with respect to the reference image. ','-1','bd','Ok!');


%% From M2_WS_02_Realign to M2_WS_03_Scrub

disp('Copy from M2_WS_02_Realign to M2_WS_03_Scrub!!')

for i=1:size(Session_name,1)
    from=[path_to_M2_WS_02_Realign,'/',Session_name{i,1}];
    to=[path_to_M2_WS_03_Scrub,'/',Session_name{i,1}];
    copyfile(from,to)
    cd(to)
    delete('*ps')
end

disp('Copying from M2_WS_02_Realign to M2_WS_03_Scrub Completed!!')

%% ART

disp('1st-pass scrubbing procedure with Artifact Detection Tool is running:')

%Configuring text file:
for i=1:size(Session_name,1)
cd(path_to_cffg_file)
fid = fopen('File_1Srub.cfg','rt');
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;
% replace string S1 with string S2
% (replacing '9999' string with correct session name)
folder_out=path_to_M2_WS_03_Scrub;
if i==1
    Y = strrep(X, '9999', Session_name{i,1});
else
    Y = strrep(X, Session_name{i-1,1}, Session_name{i,1});
end
fid2 = fopen('File_1Srub.cfg','wt'); %sovrascrivo il file ogni volta;
fwrite(fid2,Y) ;
fclose (fid2) ;
disp(['ART for ', Session_name{i,1},' is Running!!'])
cd(path_to_cffg_file)
art('sess_file','File_1Srub.cfg');
disp(['ART for ', Session_name{i,1},' Completed!!'])
end
% Restoring text file. File File_da_usare originale:
%metti 9999 al posto della sessione
cd(path_to_cffg_file)
fid = fopen('File_1Srub.cfg','rt');
X = fread(fid);
fclose(fid);
X = char(X.');
% replace string S1 with string S2
Y = strrep(X,Session_name{end,1}, '9999');
fid2 = fopen('File_1Srub.cfg','wt');%sovrascrivo il file ogni volta;
fwrite(fid2,Y);
fclose (fid2);
disp('ART for all sessions Completed!!')
spm_input('Check the Artifact Detection Tool outlier graphs. This inspection could help you in the next section where you will be ask to remove/keep outlier images. ','-1','bd','Ok!');
%spm_input('Press Continue to proceed:','-1','b','Continue|--',[1,0]);
close all

%% From M2_WS_03_Scrub to M2_WS_04_Rename

disp('Copying files from M2_WS_03_Scrub to M2_WS_04_Rename and moving ART detected ouliers into a separate subfolder:')

for i=1:size(Session_name,1)
cd([path_to_M2_WS_03_Scrub,'/',Session_name{i,1}])
mkdir('scrubbed_volumes');
load('art_regression_outliers.mat')
index = and(any(R == 1, 2), any(R == 1, 2));
outliers=find(index); %outliers variable contains scrubbed volumes for ART
%Keeping first volume if scrubbed
if sum(outliers==1)==1
    outliers(1)=[];
    disp('Reference Volume was detected as an outliers, but it was kept.');
end
list_totale=dir('r1*');
disp(['Deleting volumes number ',num2str(outliers') ,' from Session ',Session_name{i,1},'.'])
lista_to_delete=list_totale(outliers,:);
    if size(lista_to_delete,1)==0
        disp(['No detected outliers for session ',Session_name{i,1},'.'])
    else
        for jj=1:size(lista_to_delete,1)
        from=([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/',lista_to_delete(jj).name]);
        to=([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes']);
        movefile(from,to);
        end
    end
    clear list_totale lista_to_delete outliers index R
end

disp('1-pass scrubbed procedure Completed!')
close all


%% Scrubbed volumes visualization and choice to keep the volumes:

spm_input('SPM display window will show scrubbed volume(s) for each session. Instructions for removing detected outliers: ''Yes'' to scrub all of the volumes, ''No'' to keep all of the volumes or ''OnebyOne'' button to select the volume(s) you want to keep for the current session. ','-1','bd','Ok!');
clear ref2show
for i=1:size(Session_name,1)
    ref2show{1,1}=[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii,1'];
    clear img2show list
    cd([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes'])
    list=dir('*nii');
    if size(list,1)>=23
        part_disp=floor(size(list,1)/9);
        disp(['Scrubbed volumes for session ',Session_name{i,1},' will be splitted in ',part_disp, ' Display.'])
        resto = mod(size(list,1),9);
            if resto==0
              %part_disp=part_disp;
            else
              part_disp=part_disp+1;
            end
        for d=1:part_disp

                if (resto~=0 && d==part_disp)
                    clear list_current
                    for nnn=1:resto
                        list_current{nnn,1}=list((9*(d-1)+1)+(nnn-1)).name;
                    end
                else
                    clear list_current
                    for nnn=1:9
                        list_current{nnn,1}=list((9*(d-1)+1)+(nnn-1)).name;
                    end
                end

                clear img2show
                for j=1:size(list_current,1)
                 img2show{j,1}=[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list_current{j,1},',1'];%mask
                end
                spm_check_registration(ref2show{1,1},img2show{:,:});
                m=spm_input(['Scrub for Session ',Session_name{i,1}],'-1','b','Yes|No|OnebyOne',['Y','N','o']); m=char(m);
                %m=input(['Do you want to scrubb all volumes for session ',Session_name{i,1},'? Press Y to scrubb all, Press N to keep all or any button to select volumes you want to keep: '],'s');
                if m=='Y'

                    continue

                elseif m=='N'

                    for jj=1:size(list_current,1)
                        movefile([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list_current{jj,1}],[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/', list_current{jj,1}]);
                    end

                else

                        for jjj=1:size(list_current,1)
                        %prmt=['Do you want to keep ',list_current{jjj,1} ,' volume? (Press 1 to keep, 0 to scrub): '];
                        spm_check_registration(ref2show{1,1},img2show{jjj,:});
                        takebackscrub(jjj)=spm_input(['Keep ',list(jj).name ,'?'],'-1','y/n',[1,0]);
                            if takebackscrub(jjj)
                                 movefile([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/',list_current{jjj,1}],[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/',list_current{jjj,1}]);
                            end
                        end

                end

        end

    elseif size(list,1)==0
        continue
    else
        clear img2show
        for j=1:size(list,1)
             img2show{j,1}=[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list(j).name,',1'];%mask
        end
        spm_check_registration(ref2show{1,1},img2show{:,:});
        m=spm_input(['Scrub for Session ',Session_name{i,1}],'-1','b','Yes|No|OnebyOne',['Y','N','o']); m=char(m);
        if m=='Y'
            continue

        elseif m=='N'

                    for jj=1:size(list,1)
                        movefile([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list(jj).name],[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/', list(jj).name]);
                    end
        else
            for jj=1:size(list,1)
            prmt=['Do you want to keep ',list(jj).name ,' volume? (Press 1 to keep, 0 to scrub): '];
            spm_check_registration(ref2show{1,1},img2show{jj,:});
            takebackscrub(jj)=spm_input(['Keep ',list(jj).name ,'?'],'-1','y/n',[1,0]);
                if takebackscrub(jj)
                movefile([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/',list(jj).name],[path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/',list(jj).name]);
                end
            end
        end
    end
end

%Eliminating session if number of volumes is less than 1/3%

flag=zeros(size(Session_name,1),1);
spm_input('Sessions with less than 2/3 of the initial number of volumes will be deleted due to excessive Within Session movement.','-1','bd','Ok!');

disp('Sessions with less than 2/3 of the initial number of volumes will be deleted due to excessive Within Session movement.');

   for i=1:size(Session_name,1)
       check_volumes=dir([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/*nii']);
       if length(check_volumes)>=floor(vol4session*2/3)
          flag(i)=0;
          disp(['Session ',Session_name{i,1},' was deleted.'])
       else
          flag(i)=1;
          disp(['Session ',Session_name{i,1},' was kept.'])
       end
   end

Session_name=Session_name(logical(flag));   
spm_input(['Remaining session(s) are: ', num2str(length(Session_name)),' out of ',num2str(length(flag)),'.'],'-1','bd','Ok!');

%% Rename files inside M2_WS_04_Rename:

disp('Renaming each session after 1-pass scrubbing procedure:')

cd(path_to_M2_WS_04_Rename)
for i=1:size(Session_name,1)
    mkdir(Session_name{i,1})
end

for i=1:size(Session_name,1)
    cd([path_to_M2_WS_03_Scrub,'/',Session_name{i,1}]);
    files = dir('r1*');
    ii=1;
         for id = 1:length(files)
             num=['r1m',Session_name{i,1},'_',num2str(id,'%03d'),'.nii'];
             copyfile([path_to_M2_WS_03_Scrub,'/',Session_name{i,1},'/',files(ii).name],[path_to_M2_WS_04_Rename,'/',Session_name{i,1},'/',num]);
             ii=ii+1;
         end
        clear files
        disp(['Renaming volumes of Session ',Session_name{i,1},', DONE!'])
end

close all
disp('Renaming Completed!')
disp('You have COMLETED MODULE 2!!')

cd(path_general)
save('Variables_M2.mat','template2bused','GW2beused','wgh_winner','Session_name','vol4session')
diary off

cd(path_to_M1_PP_01_OrigVol);
list=dir();
% Get vector to recognize folders only:
list(ismember({list.name},{'.','..'}))=[];
if length(list)>1
   
spm_input('You have more than one session therefore you will proceed with the standard procedure.','-1','bd','Ok!');
   
%% MODULO 3:
disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||       //°// ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||             //  ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||         \\  ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||              )) ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||       __//  ')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance Imaging fMRI -------------------')
disp('-------------------------- Within Session Processing  ---------------------------')
disp('----------From Session-specific Functional References Segmentation --------------')
disp('--------to Mean inner-brain Masked Functional References Calculation-------------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                       .\( / / . ) _)_)_______/                                  ')
disp('                        .( ( (  (__/__/  .                                       ')
disp('                          \._)_). )                                              ')
disp('                           ._)_)./                                               ')
disp('                           ._)_)./                                               ')
disp('                                                                                 ')

%% Creation of Mask_ref_Vol folder and creation od session folders inside 07_2_Segment_Ref_Volumes:

cd(path_to_M3_WS_01_SegRefVols);
for i=1:size(Session_name,1)
    mkdir(Session_name{i,1})
end

disp('Insertion of each reference volume inside the respective folder:')
for i=1:size(Session_name,1)
    from=[path_to_M2_WS_04_Rename,'/',Session_name{i,1},'/','r1m',Session_name{i,1},'_001.nii'];
    to=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1}];
    copyfile(from,to)
    disp(['Inserting reference of ',Session_name{i,1},' inside ',Session_name{i,1} ,' folder. Done!'])
end

%% Session-Specific Functional References Segmentation for inner-brain mask calculation:

    % Option to remask images:
    
cd(path_general)
diary Module3_logfile
for i=1:size(Session_name,1)
        clear temp2show file2show
        disp(['These are the 1-pass session-specific masked images for session: ',Session_name{i,1}]);

           % Reference of session i
        file2show{1,1}=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii'];
        % Template:
        k=1;
        for j=21:GW2beused
             temp2show{k,1}=[path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/',Session_name{i,1},'scGW_',num2str(j),'.nii'];
             answ_tmp(k)=j;
             k=k+1;
        end
        spm_check_registration(file2show{1,1},temp2show{:,:});
        answw='';
        for hh=1:length(answ_tmp)
        answ=num2str(answ_tmp(hh));
            if hh<length(answ_tmp)
                answw=[answw,answ,'|'];
            else
                answw=[answw,answ];
            end
        end
        clear input
        m=spm_input('Mask the Ref. image again?','-1','y/n',['N','Y']);m=char(m);
            if m=='Y'
                continue
            else
                if i==1
            spm_input(['The chosen session-specific template mask should include the entire fetal brain and a minimal amount of the skull (black edges).'],'-1','bd','Yes');
                end
            GW_tmp=spm_input('Choose the session-specific template mask','-1','m',answw,answ_tmp,30);

                % Specific choices
                clear input template2show ref2show temp2show
                    for j=1:3
                        template2show{j,1}=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'template_anp/',Session_name{1,1},'banp_',num2str(j),'_',num2str(GW_tmp),'.nii'];
                    end
                file2show=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii'];
                temp2show=[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_template/',Session_name{1,1},'scGW_',num2str(GW_tmp),'.nii'];
                spm_check_registration(file2show,temp2show,template2show{:,:});

                spm_input('For the mask you have chosen, select the one that covers the whole brain in each direction (es. L/R, A/P; I/S).','-1','bd','Ok!');
                template2bused=spm_input('Insert a number between 1 and 4: ','-1','m','1|2|3|4',[1,2,3,4]);
%Applying mask:

                    cd([path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1}]);
                        clear matlabbatch list_to_mask

                                matlabbatch{1}.spm.util.imcalc.input{1,1} = char(strcat(path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii,1'));
                                if template2bused==1
                                     matlabbatch{1}.spm.util.imcalc.input{2,1} = char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'_template/',Session_name{i,1},'scGW_',num2str(GW_tmp),'.nii,1'));
                                else
                                     matlabbatch{1}.spm.util.imcalc.input{2,1} = char(strcat(path_to_M2_WS_01_Mask,'/',Session_name{i,1},'/',Session_name{i,1},'template_anp/',Session_name{i,1},'banp_',num2str(template2bused-1),'_',num2str(GW_tmp),'.nii,1'));
                                end
                        matlabbatch{1}.spm.util.imcalc.output = char(strcat('r1m',Session_name{i,1},'_001.nii'));
                        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1})};
                        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
                        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
                        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                        spm_jobman('run',matlabbatch);
                        GW_tmp=GW2beused;
            end
end

    % Segmentation and Gestational Week selection

spm_input(['In this section Session-Specific Functional References Segmentation of each reference image is processed. It is suggested to use the GW that was used in module 2: ',num2str(GW2beused)],'-1','bd','Ok!');
Decision_GW = spm_input(['Do you want to proceed with ',num2str(GW2beused),' GW for classes C1:C7 and inner and outer brain?'],'-1','bd','yes|no',[1,0]);

if Decision_GW==1
    GW_winner=GW2beused;
    GWin_out_winner=GW2beused;
else

    GW_winner=spm_input('Select Gestation Week for classes C1:C7:','-1','m','21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36|37',[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37],30);
    GWin_out_winner=spm_input('Select Gestation Week for Inner/Outer brain:','-1','m','21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36|37',[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37],30);
end

%Creation of path for masked references images:

disp('Creating Session-Specific Functional References masks.');
clear input path3D matlabbatch
input= cell(length(Session_name),1);
for i=1:size(Session_name,1)
input(i,:)=cellstr([path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii,1']);% input images
end
matlabbatch{1}.spm.spatial.preproc.channel.vols = {
                                                   ''
                                                   };
matlabbatch{1}.spm.spatial.preproc.channel.biasreg =  0.01;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc1_34GW_CorticalPlate_Cerebellum_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc2_34GW_WM_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc3_34GW_CSF_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc4_34GW_DGM_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc5_34GW_Hippocampi_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc6_34GW_Amygdala_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc7_34GW_Brainstem_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(7).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(7).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(8).tpm = {strcat(path_general,'Template/Template_Priors_Seg/Inbrains_and_Outbrains/34GW/src8_InnerBrain_34GW.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(8).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(8).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(8).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(9).tpm = {strcat(path_general,'Template/Template_Priors_Seg/Inbrains_and_Outbrains/34GW/src9_OuterBrain_34GW.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(9).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(9).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(9).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'none';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{1}.spm.spatial.preproc.channel.vols = input;
for i=1:9
    TMp=matlabbatch{1}.spm.spatial.preproc.tissue;
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm=strrep(TMp(i).tpm,'34',num2str(GW_winner));
    if (i==8 || i==9) % In and Outer brain
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm=strrep(TMp(i).tpm,'34',num2str(GWin_out_winner));
    end
end
disp(['Segmentation is running with GW for C1:C7: ',num2str(GW_winner)]);
disp(['and for inner and outer brain with GW: ',num2str(GWin_out_winner)]);
spm_jobman('run',matlabbatch);
disp('Segmentation COMPLETED!')

%% Creating inner-brain mask for all sessions:

disp('Creating inner-brain Session-Specific Functional References masks for the reference image of each Session!');
clear matlabbatch input
matlabbatch{1}.spm.util.imcalc.input = '';
matlabbatch{1}.spm.util.imcalc.output = 'c9_mask';
matlabbatch{1}.spm.util.imcalc.outdir = '';
matlabbatch{1}.spm.util.imcalc.expression = '1-i1';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: c9_mask', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

for i=1:size(Session_name,1)
input=cellstr([path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/c9r1m',Session_name{i,1},'_001.nii,1']);% input images
output=cellstr([path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1}]);% input images
matlabbatch{1}.spm.util.imcalc.input = input;
matlabbatch{1}.spm.util.imcalc.outdir=output;
disp(['Creation of inner-brain Session-Specific Functional References masks for Session number: ',Session_name{i,1},' COMPLETED!']);
spm_jobman('run',matlabbatch);
end

%% Visual sc9 inspection

spm_input(['Now the inner-brain Session-Specific Functional Reference mask will be displayed for each session. Segmentation was processed with : ',num2str(GW_winner) ,' Gestational Week.'],'-1','bd','Ok!');

% Comment this section to skip the visualization;
clear ref2show mask2show
for i=1:size(Session_name,1)
disp(['Display session: ',Session_name{i,1}, ' and relative session-specific inner-brain mask.'])
ref2show{1,1}=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/r1m',Session_name{i,1},'_001.nii,1'];% Nomask images
mask2show{1,1}=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/sc9_mask.nii'];%mask
spm_check_registration(ref2show{1,1},mask2show{1,1});
clear ref2show mask2show
if i== size(Session_name,1)
spm_input('Press continue to proceed:','-1','b','Continue|--',[1,0]);
else
spm_input('Press continue to see next:','-1','b','Continue|--',[1,0]);
end
end

% Deleting temporarly file

for i=1:size(Session_name,1)
   cd([path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1}])
   file2bedeleted=dir();
   file2bedeleted(ismember({file2bedeleted.name},{'.','..','sc9_mask.nii','c9_mask.nii'}))=[];
   delete(file2bedeleted.name)
end


%% Copying files

cd(path_to_M3_WS_02_MaskRefVols)
for i=1:size(Session_name,1)
mkdir(Session_name{i,1})
end
% Copio file rinominati scrubbati
for i=1:size(Session_name,1)
    from=[path_to_M2_WS_04_Rename,'/',Session_name{i,1}];
    to=[path_to_M3_WS_02_MaskRefVols,'/',Session_name{i,1}];
    copyfile(from,to)
    disp(['Copying volumes for session ',Session_name{i,1},' COMPLETED!']);
end
% Copying sc9 masks
for i=1:size(Session_name,1)
    from=[path_to_M3_WS_01_SegRefVols,'/',Session_name{i,1},'/sc9_mask.nii'];
    to=[path_to_M3_WS_02_MaskRefVols,'/',Session_name{i,1}];
    copyfile(from,to)
    disp(['Copying session-specific inner-brain mask for session ',Session_name{i,1},' COMPLETED!']);
end

%% Applying mask to all volumes:
clear matlabbatch
for iii=1:size(Session_name,1)
    disp(['Applying session-specific inner-brain mask to all volumes of ',Session_name{iii,1},' COMPLETED!']);
    cd([path_to_M3_WS_02_MaskRefVols,'/',Session_name{iii,1}]);
    list_to_mask=dir('r1m*');
    niiPrefix='r1m';
    for i=1:length(list_to_mask)

        matlabbatch{1}.spm.util.imcalc.input = {
                                                strcat(path_to_M3_WS_02_MaskRefVols,'/',Session_name{iii,1},'/',char(niiPrefix),Session_name{iii,1},'_',num2str(i,'%03d'),'.nii,1')
                                                strcat(path_to_M3_WS_02_MaskRefVols,'/',Session_name{iii,1},'/sc9_mask.nii,1')
                                                };
        matlabbatch{1}.spm.util.imcalc.output = char(strcat('m',niiPrefix,Session_name{iii,1},'_', num2str(i,'%03d'), '.nii'));
        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(path_to_M3_WS_02_MaskRefVols,'/',Session_name{iii,1})};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
    end
    clear matlabbatch list_to_mask
end

for i=1:size(Session_name,1)
     cd([path_to_M3_WS_02_MaskRefVols,'/',Session_name{i,1}]);
     delete('r1*')
     disp(['Deleting r1* input from: ' Session_name{i,1}]);
end

%% Realign masked ref.
% Insertion of reference volume inside M3_WS_03_RealignMaskRefVols :
for i=1:size(Session_name,1)
    from=[path_to_M3_WS_02_MaskRefVols,'/',Session_name{i,1},'/','mr1m',Session_name{i,1},'_001.nii'];
    to=path_to_M3_WS_03_RealignMaskRefVols;
    copyfile(from,to)
    disp(['Insertion of Functional Reference volume of ',Session_name{i,1},' inside M3_WS_03_RealignMaskRefVols.']);
end
        clear matlabbatch
        matlabbatch{1, 1}.spm.spatial.realign.estwrite.data='';
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {[path_to_M2_WS_01_Mask,'/',Session_name{1,1},'/',Session_name{1,1},'_template/',Session_name{1,1},'scGW_',wgh_winner{1,1},'.nii,1']};
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % do not reslice images only mean
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1; % naerest N.
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
        for  i=1:size(Session_name,1)
         matlabbatch{1, 1}.spm.spatial.realign.estwrite.data{i,1}={
                                                         strcat(path_to_M3_WS_03_RealignMaskRefVols,'/mr1m',Session_name{i,1},'_001.nii,1')
                                                         };
        end
        spm_jobman('run',matlabbatch);
        disp('Realignment Done and mean reference created!.');

movefile([path_to_M3_WS_03_RealignMaskRefVols,'/','meanmr1m',Session_name{1,1},'_001.nii'],[path_to_M3_WS_03_RealignMaskRefVols,'/Mean_fMRI_Ref_Vol.nii']);

%% Visual Mean Reference inspection

clear ref2show
ref2show = cell(length(Session_name),1);
for i=1:size(Session_name,1)
ref2show{i,:}=[path_to_M3_WS_03_RealignMaskRefVols,'/mr1m',Session_name{i,1},'_001.nii,1'];% Nomask images
end
mean2show{1,1}=[path_to_M3_WS_03_RealignMaskRefVols,'/Mean_fMRI_Ref_Vol.nii'];%mask
spm_check_registration(ref2show{:,:},mean2show{:,:});

%% 08_fMRI_Rename: 

for i=1:size(Session_name,1)
    %disp(['Copying masked volumes of session ', Session_name{i,1}, ' from M3_WS_02_MaskRefVols to M3_WS_04_Rename folder.' ]);
    from=[path_to_M3_WS_02_MaskRefVols,'/',Session_name{i,1}];
    to=[path_to_M3_WS_04_Rename,'/',Session_name{i,1}];
    copyfile(from,to)
    cd(to)
    disp('DONE');
    delete sc*
    if i==1
    copyfile([path_to_M3_WS_03_RealignMaskRefVols,'/','Mean_fMRI_Ref_Vol.nii'],to) % sara da chiamare 000
    movefile([path_to_M3_WS_04_Rename,'/',Session_name{1,1},'/Mean_fMRI_Ref_Vol.nii'],[path_to_M3_WS_04_Rename,'/',Session_name{1,1},'/mr1m',Session_name{1,1},'_000.nii'])
    disp(['Positioning Mean Reference at first place of session: ' Session_name{i,1}]);
    end
end

% Rename:
disp('Renaming all volumes ');
count=0;
for i=1:size(Session_name,1)
    f = [path_to_M3_WS_04_Rename,'/',Session_name{i,1}];
    cd(f)
    files = dir('*.nii');
    ii=1;
         for id = (1+count):(length(files)+count)
             num=['mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
             movefile([path_to_M3_WS_04_Rename,'/',Session_name{i,1},'/',files(ii).name],[path_to_M3_WS_04_Rename,'/',Session_name{i,1},'/',num]);
             ii=ii+1;
         end
        count=count+length(files);
        clear files
        disp(['Renaming volumes of Session ',Session_name{i,1},', DONE!'])
end
disp('Renaming Completed!!');

spm_input('Mean inner-brain Masked Funcrtional Reference volume was calculated and named Mean_fMRI_Ref_Vol.nii and it was placed as the first volume of the first session as required for the Between Session processing. Press ok to proceed with Realignment of inner-brain masked Functional images in module 4: ','-1','bd','Ok!');

% GWin_out_winner GW_winner
cd(path_general)
save('Variables_M3.mat','GW_winner','GWin_out_winner','Session_name','vol4session','GW2beused')
disp('Copying DONE!')
disp('Module 3 Completed!')
diary off

%% MODULO 4:
disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||       //    ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||          //     ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||     //      ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||        //__||_  ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||        ||   ')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance Imaging fMRI -------------------')
disp('-------------------------- Between Session Processing  --------------------------')
disp('-----Realignment of all inner-brain Session-Specific Masked Functional scans ----')
disp('---------to Mean BS Functional Reference and 2-pass Scrubbing procedure----------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                      .\( / / . ) _)_)_________/                                 ')
disp('                        .( ( (  (__/__/.../                                      ')
disp('                           . \._)_). )                                           ')
disp('                             ._)_)./                                             ')
disp('                             ._)_)./                                             ')
disp('                                                                                 ')


%% Copying file from previous module

disp('Copying file from M3_WS_04_Rename to M4_BS_01_Realign_Reslice:')
for i=1:size(Session_name,1)
    from=[path_to_M3_WS_04_Rename,'/',Session_name{i,1}];
    to=[path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1}];
    copyfile(from,to)
end

cd(path_general)
diary Module4_logfile

resl=spm_input('If you proceed to single subject statistical analyses, you should reslice Realigned inner-brain masked Functional images. For Group analysis reslicing is not recommanded. Press ''Yes'' to reslice or ''No'' to not reslice:' ,'-1','bd','Yes|No',[1,0]);

if resl==1
    resl=resl+1;
end

clear matlabbatch
        matlabbatch{1}.spm.spatial.realign.estwrite.data = ''; 
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [resl 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r2';

for i=1:size(Session_name,1)
       cd([path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1}]);
       clear list input
       list=dir('*nii');
        for j=1:length(list)
              input(j,:)=cellstr([path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1},'/',list(j).name,',1']);% input images
              matlabbatch{1}.spm.spatial.realign.estwrite.data{1, i}=input;
        end
end
spm_jobman('run',matlabbatch);

%Renaming
for i=1:size(Session_name,1)
    cd(path_to_M4_BS_02_Scrub)
    mkdir(Session_name{i,1})
end


if resl==0
count=0;
for i=1:size(Session_name,1)
    cd([path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1}])
    files = dir('mr1m*');
    ii=1;
         for id = (1+count):(length(files)+count)
             num=['r2mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
             copyfile([path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1},'/',files(ii).name],[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',num]);
             ii=ii+1;
         end
         count=count+length(files);
         clear files
end

else
    
    for i=1:size(Session_name,1)
    from=[path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1},'/r2*'];
    to=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1}];
    copyfile(from,to)
    end
    

end

% Visualization of Realign of
disp('Visual inspection of random BS Realigned inner-brain masked Functional images.')
clear ref2show mean2show img2show
mean2show{1,1}=[path_to_M4_BS_01_Realign_Reslice,'/',Session_name{1,1},'/mr1m',Session_name{1,1},'_0001.nii,1'];
c=0;
for i=1:size(Session_name,1)
    clear randomnum 
    cd([path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1}]);
    all_files=dir('mr1*');
        if size(all_files,1)~=0
            randomnum=round(1+rand(1,2)'*(size(all_files,1)-1));
        end
    for j=1:2
         img2show{j+2*c,1}=[path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1},'/',all_files(randomnum(j)).name,',1'];%mask
    end
    c=c+1;
end
spm_check_registration(mean2show{1,1},img2show{:,:});


%% SCRUB

spm_input('Realignment of all inner-brain Session-Specific Masked Functional Scans completed. Press ok to proceed with the 2nd scrubbing procedure with ART, DVARS and FD. ','-1','bd','Ok!');

cd(path_to_M4_BS_02_Scrub)
mkdir('Motion_file')
mkdir('art_single_sessions')

for i=1:size(Session_name,1)
    from=[path_to_M4_BS_01_Realign_Reslice,'/',Session_name{i,1},'/','*txt'];
    to=[path_to_M4_BS_02_Scrub,'/Motion_file/'];
    copyfile(from,to)
end

% ART 2
%Configuring test file:
disp('Updating the Configuration file for ART.')
cd(path_to_cffg_file)
fid = fopen('File_2Scrb_orig.cfg','rt');
X = fread(fid) ;
fclose(fid) ;
X = char(X.');
Y = strrep(X, '9999', Session_name{1,1});

% Change X file:
Y(11)=num2str(length(Session_name));
for i=1:size(Session_name,1)
Imagestobeinserted{i,:}=(['session ', num2str(i), ' image ',Session_name{i,1},'/', 'r2mr1m',Session_name{i,1},'_????.nii']);
end

cd([path_to_M4_BS_02_Scrub,'/Motion_file/'])
list=dir('*txt');
        for nn=1:size(list,1)
        list_uff{nn,1}=strcat(list(nn).name);
        Utils(nn,:)=regexp(list_uff{nn,1},'_','split');
        end

util=Utils(:,3);
[ordinare,indicix]=sort(util);
list_finale=list(indicix,:);

for i=1:size(Session_name,1)
Motiontobeinserted{i,:}=(['session ', num2str(i), ' motion ', list_finale(i).name]);
end

textlength=length(Y); % inizio parte interessante del file cfg
counter=textlength+1;

for i=1:size(Session_name,1)
   l1=length(char(Imagestobeinserted(i,:)));
   Y(counter:counter+l1-1)=char(Imagestobeinserted(i,:));
   counter=counter+l1;
   Y(counter:counter+4)='     ';
   counter=counter+4;
end

for i=1:size(Session_name,1)
   l1=length(char(Motiontobeinserted(i,:)));
   Y(counter:counter+l1-1)=char(Motiontobeinserted(i,:));
   counter=counter+l1;
   Y(counter:counter+4)='     ';
   counter=counter+4;
   if i==length(Session_name)
   fine='end';
   Y(counter:counter+length(fine)-1)=char(fine);
   end
end

cd(path_to_cffg_file)
fid2 = fopen('File_2Scrub.cfg','wt');
fwrite(fid2,Y);
fclose (fid2);
i=1;
folder_out=path_to_M4_BS_02_Scrub;
cd(path_to_cffg_file)
disp('Updating of Configuration file for ART Completed (File_2Scrub.cfg)')
disp('ART in progress!')
cd(path_to_cffg_file)
art('sess_file','File_2Scrub.cfg');

%clear variable:
clear Y counter Motiontobeinserted Utils util list Imagestobeinserted X art

%close all

%Caricare in cartella 'motion_outliers' i volumi 4D di ogni sessione (copio da 00_Original..).
cd(path_to_M4_BS_02_Scrub)
mkdir('motion_outliers');

%SMP 3D to 4D
disp('Creating 4D files!')
for i=1:size(Session_name,1)
    clear matlabbatch list input
    pathtosession=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1}];
    cd(pathtosession)
    list=dir('*.nii');
        for j=1:length(list)
              input(j,:)=cellstr([pathtosession,'/',list(j).name,',1']);% input images
        end
        matlabbatch{1}.spm.util.cat.vols=input;
    matlabbatch{1}.spm.util.cat.name = [Session_name{i,1},'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.RT = NaN;
    spm_jobman('run',matlabbatch);
    movefile([pathtosession,'/',Session_name{i,1},'.nii'],[path_to_M4_BS_02_Scrub,'/motion_outliers']);
end

disp('Creation of 4D files for DVARS and FD Completed!')

%% Dvars FD for each session:

%Creo una cartella per ogni sessione e dentro metto le cartelle DVARS e FD:
for i=1:size(Session_name,1)
    cd([path_to_M4_BS_02_Scrub,'/motion_outliers'])
    mkdir(Session_name{i,1})
    cd([path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1}])
    mkdir('DVARS')
    mkdir('FD')
end

%% DVARS && FD

lliist_rp=dir([path_to_M4_BS_02_Scrub,'/Motion_file/*.txt']);

for i=1:size(Session_name,1)
    
clear func_fn func_spm Y Tmean mask_tmean brainmed_tmp DY Dvars Dvars_outliers brain_value_mediana
clear R_MP T_MP Tdiff Rdiff FD_outliers FD Rdiff_fin Tdiff_fin

cd([path_to_M4_BS_02_Scrub,'/motion_outliers'])

disp(['Running DVARS and FD for session : ',Session_name{i,1}]);

% Loading 4D of relative session:
func_fn = [path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1},'.nii'];
func_spm = spm_vol(func_fn);
V2 = spm_read_vols(func_spm);
X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); T0 = size(V2,4);
I0 = prod([X0,Y0,Z0]);
Y  = reshape(V2,[I0,T0]); clear V2 V1;
Y = double(Y); % In ogni colonna ho un immagine e nelle righe i voxel

% Temporal mean:
Tmean=mean(Y,2);
Y1=Y(Y~=0);
thr2 = prctile(Y1,2);
thr98 = prctile(Y1,98);
robthr=thr2 + 0.1 * (thr98 - thr2 );
% calculating mask (vector)
mask_tmean=Tmean>robthr; 
maskmean=sum(mask_tmean==1)/length(mask_tmean);% calculating fraction of 1-intensity voxels (vector)
for iii=1:T0
    brainmed_tmp(:,iii)=Y(:,iii).*mask_tmean;
end
brain_value_mediana=median(brainmed_tmp(brainmed_tmp~=0)); 

% Dvars calculation
DY=(diff(Y,1,2));
for iii=1:T0-1
    DY(:,iii)=DY(:,iii).*mask_tmean;
end
DY=mean(DY.^2);
DY=DY/maskmean;
DY=sqrt(DY);
Dvars=(DY/brain_value_mediana)*1000;
thr25 = prctile(Dvars,25);
thr75 = prctile(Dvars,75);
threshv = thr75 + 1.5 * (thr75-thr25);
Dvars_value=[0,Dvars]; % calcolerei la soglia senza lo 0 Dvars_value=Dvars;
% Dvars Outliers detection
Dvars_outliers=(Dvars_value>threshv);

% Savings output

dlmwrite([path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1},'/DVARS/output_vals.txt'],Dvars_value')
dlmwrite([path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1},'/DVARS/output.txt'],Dvars_outliers')

%% Frame-Wise Displacement

MP=load([path_to_M4_BS_02_Scrub,'/Motion_file/',lliist_rp(i).name]);
drflag = [0 0 0 1 1 1]; % 0 trans 1 rotation
r_Idx   = find(drflag);
t_idx   = find(~drflag);
R_MP=MP(:,r_Idx);
T_MP=MP(:,t_idx);
%Differenziale e aggiungo l'ultima riga che rimane uguale(tolgo 0)
Tdiff=diff(T_MP,1,1);
Tdiff=[Tdiff;-T_MP(end,:)]; %ok
%Differenziale e aggiungo l'ultima riga che rimane uguale(tolgo 0)
Rdiff=diff(R_MP,1,1);
Rdiff=[Rdiff;-R_MP(end,:)]; %ok
%process rotazioni
Rdiff=abs(Rdiff)*50;
Rdiff_fin=mean(Rdiff,2)*3;
%process traslazioni
Tdiff_fin=(mean(abs(Tdiff),2))*3;
%merge
FD=Rdiff_fin+Tdiff_fin;
% FD Outliers detection
thr25 = prctile(FD,25);
thr75 = prctile(FD,75);
threshv = thr75 + 1.5 * (thr75-thr25);
FD_outliers=[0;(FD>threshv)];
FD_outliers(end)=[];
FD=[0;FD];
FD(end)=[];

% Savings output

dlmwrite([path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1},'/FD/output_vals.txt'],FD)
dlmwrite([path_to_M4_BS_02_Scrub,'/motion_outliers/',Session_name{i,1},'/FD/output.txt'],FD_outliers)

disp(['DVARS and FD for session : ',Session_name{i,1}, ' COMPLETED!']);
end

disp('Concatenating Outliers of DVARS, FD and ART!');

% acquisire elenco delle sottocartelle in motion_outliers
cd(path_to_M4_BS_02_Scrub)
d = dir('motion_outliers');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% Creo vettore in cui memorizzo il numero di volumi in ogni sottocartella
numFilesInFolders = zeros(size(nameFolds, 1), 1);
% e lo riempio basandomi sulle info contenute nei file "motion_outliers/*/FD/output_vals.txt"
for i = 1:length(nameFolds)
  currentFile = fullfile('motion_outliers', nameFolds(i), 'FD', 'output_vals.txt');
  A = importdata(char(currentFile));
  numFilesInFolders(i) = size(A, 1);
end

% creare vettore colonna di output. Lungo 60 * numero_sessioni, pieno di zeri.
outputVect = [];

for i = 1:length(nameFolds)
   tempZeros = zeros(numFilesInFolders(i), 1);
   dummyZeros = zeros(numFilesInFolders(i), 1);
   % aprire DVARS, cercare file output.txt
   %currentFile = strcat('motion_outliers/', nameFolds(i), '/DVARS/output.txt');
   currentFile = fullfile('motion_outliers', nameFolds(i), 'DVARS', 'output.txt');
   tempCol = [];
   if exist(char(currentFile), 'file') == 2
       %  disp(currentFile)
       %  se non esiste, pace
       %  se esiste, leggere il contenuto con importdata('output.txt')
       A = importdata(char(currentFile));
       A = [A dummyZeros];
       A = A';
       %  merge degli 1 su ununico vettore colonna
       tempCol = sum(A);
       tempCol = tempCol';
       tempZeros = tempZeros + tempCol;
       
   end

   % ripetere con FD/output.txt
   %currentFile = strcat('motion_outliers/', nameFolds(i), '/FD/output.txt');
   currentFile = fullfile('motion_outliers', nameFolds(i), 'FD', 'output.txt');
   if exist(char(currentFile), 'file') == 2
       disp(currentFile)
       % se non esiste, pace
       % se esiste, leggere il contenuto con importdata('output.txt')
       A = importdata(char(currentFile));
       A = [A dummyZeros];
       A = A';
       %  merge degli 1 su ununico vettore colonna
       tempCol = sum(A);
       tempCol = tempCol';
       tempZeros = tempZeros + tempCol;
   end
   tempZeros = tempZeros > 0;

   outputVect = [outputVect; tempZeros];
end

cd(path_to_M4_BS_02_Scrub)
load(fullfile('art_single_sessions', 'art_outliers.mat'));
art = out_idx;
for i = 1:length(art)
   outputVect(art(i)) = 1;
end

% outputIdx contains outputVect values that differs form 0:
outputIdx = [];
for i = 1:length(outputVect)
   if outputVect(i) == 1
      outputIdx = [outputIdx, i];
   end
end

R = zeros(sum(numFilesInFolders), length(outputIdx));
for i = 1:length(outputIdx)
    R( outputIdx(i), i ) = 1;
end

% disp('Outliers index: ')
outliersIdx = outputIdx';

% saving R regress_dvars_fd_art file
save('regress_dvars_fd_art', 'R');
% saving  outliers index in indici_outliers
save('indici_outliers', 'outliersIdx');
% saving  same index in a csv.
csvwrite('indici_outliers.csv', outliersIdx);

% Moveoutliers:
% Creating scrubbed_volumes folders
for i=1:size(Session_name,1)
    cd([path_to_M4_BS_02_Scrub,'/',Session_name{i,1}])
    mkdir('scrubbed_volumes')
end

disp('Placing outlier volumes in the relative session ''scrubbed_volumes'' folder.');

%muovo gli outliers
cd(path_to_M4_BS_02_Scrub)
load('indici_outliers.mat') % outliersIdx
for i=1:length(outliersIdx)
    current_out = outliersIdx(i);
    daconfrontare=num2str(current_out,'%04d');
    for ii=1:size(Session_name,1)
        current_Session = Session_name{ii,1};
        cd([path_to_M4_BS_02_Scrub,'/',Session_name{ii,1}])
        list=dir('*nii');
            for j=1:length(list)
                session_tmp=regexp(list(j).name,'_','split');
                volume_nii=char(session_tmp(2));
                volume_nii(end-3:end)='';
                 % confronto i due nomi ( outlier e file corrente) Se sono
                 % uguali sposto il file in scrubbed volume
                if (sum(daconfrontare == volume_nii)==4 && sum(daconfrontare == '0001')~=4)
                    from=([path_to_M4_BS_02_Scrub,'/',Session_name{ii,1},'/',list(j).name]);
                    to=([path_to_M4_BS_02_Scrub,'/',Session_name{ii,1},'/scrubbed_volumes/',list(j).name]);
                    movefile(from,to)
                    disp(['Removing outlier ',num2str(current_out),' from session ',Session_name{ii,1},'.']);
                end
            end
    end
end

%% Scrubbed volumes visualization and choice to keep the volumes:

spm_input('SPM display window will show scrubbed volume(s) for each session. Instructions for removing detected outliers: ''Yes'' to scrub all of the volumes, ''No'' to keep all of the volumes or ''OnebyOne'' button to select volumes you want to keep for the current session. ','-1','bd','Ok!');

for i=1:size(Session_name,1)
    clear img2show list lista_temporanea ref2show
    cd([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes']);
    lista_temporanea=dir([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/*nii']);
    list=dir('*nii');
    if size(list,1)>=23
        ref2show{1,1}=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',lista_temporanea(1).name,',1'];
        disp(['Scrubbed volumes for session ',Session_name{i,1},' will be splitted in ',part_disp, ' Display.'])
        part_disp=floor(size(list,1)/9);
        resto = mod(size(list,1),9);
            if resto==0
              %part_disp=part_disp;
            else
              part_disp=part_disp+1;
            end
        for d=1:part_disp

                if (resto~=0 && d==part_disp)
                    clear list_current
                    for nnn=1:resto
                        list_current{nnn,1}=list((9*(d-1)+1)+(nnn-1)).name;
                    end
                else
                    clear list_current
                    for nnn=1:9
                        list_current{nnn,1}=list((9*(d-1)+1)+(nnn-1)).name;
                    end
                end

                clear img2show
                for j=1:size(list_current,1)
                 img2show{j,1}=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list_current{j,1},',1'];%mask
                end
                spm_check_registration(ref2show{1,1},img2show{:,:});
                clear input
                m=spm_input(['Scrub for Session ',Session_name{i,1}],'-1','b','Yes|No|OnebyOne',['Y','N','o']); m=char(m);
                if m=='Y'

                    continue

                elseif m=='N'

                    for jj=1:size(list_current,1)
                        movefile([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list_current{jj,1}],[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/', list_current{jj,1}]);
                    end

                else

                        for jjj=1:size(list_current,1)
                        prmt=['Do you want to keep ',list_current{jjj,1} ,' volume? (Press 1 to keep, 0 to scrub): '];
                        spm_check_registration(ref2show{1,1},img2show{jjj,:});
                         takebackscrub(jjj)=spm_input(['Keep ',list(jj).name ,'?'],'-1','y/n',[1,0]);
                            if takebackscrub(jjj)
                                 movefile([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/',list_current{jjj,1}],[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',list_current{jjj,1}]);
                            end
                        end

                end

        end

    elseif size(list,1)==0
        continue
    else
        ref2show{1,1}=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',lista_temporanea(1).name,',1'];
        for j=1:size(list,1)
             img2show{j,1}=[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list(j).name,',1'];%mask
        end
         spm_check_registration(ref2show{1,1},img2show{:,:});
         clear input
        m=spm_input(['Scrub for Session ',Session_name{i,1}],'-1','b','Yes|No|OnebyOne',['Y','N','o']); m=char(m);        
        if m=='Y'
            continue

        elseif m=='N'

                    for jj=1:size(list,1)
                        movefile([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/', list(jj).name],[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/', list(jj).name]);
                    end
        else
            for jj=1:size(list,1)
            %prmt=['Do you want to keep ',list(jj).name ,' volume? (Press 1 to keep, 0 to scrub): '];
            spm_check_registration(ref2show{1,1},img2show{jj,:});
             takebackscrub(jj)=spm_input(['Keep ',list(jj).name ,'?'],'-1','y/n',[1,0]);
                if takebackscrub(jj)
                movefile([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/',list(jj).name],[path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',list(jj).name]);
                end
            end
        end
    end
end

%Eliminating session if number of volumes is less than 33%
clear flag check_volumes
flag=zeros(size(Session_name,1),1);
disp('Sessions with less than 2/3 of the initial number of volumes will be deleted due to excessive movement Between Session.')
   for i=1:size(Session_name,1)
       check_volumes=dir([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/scrubbed_volumes/*nii']);
       if length(check_volumes)>=floor(vol4session*2/3)
          flag(i)=0;
          disp(['Session ',Session_name{i,1},' was deleted.'])
       else
          flag(i)=1;
          disp(['Session ',Session_name{i,1},' was kept.'])
       end
   end

Session_name=Session_name(logical(flag));
spm_input(['Remaining session(s) are: ', num2str(length(Session_name)),' out of ',num2str(length(flag)),'.'],'-1','bd','Ok!');
   

%% Copy and renaming file from M4_BS_02_Scrub to M4_BS_03_Rename
cd(path_to_M4_BS_03_Rename)

for i=1:size(Session_name,1)
    mkdir(Session_name{i,1});
end

disp('Renaming all volumes: ');
count=0;
for i=1:size(Session_name,1)
    f = [path_to_M4_BS_02_Scrub,'/',Session_name{i,1}];
    cd(f)
    files = dir('*.nii');
    ii=1;
         for id = (1+count):(length(files)+count)
             num=['r2mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
             copyfile([path_to_M4_BS_02_Scrub,'/',Session_name{i,1},'/',files(ii).name],[path_to_M4_BS_03_Rename,'/',Session_name{i,1},'/',num]);
             ii=ii+1;
         end
        count=count+length(files);
        clear files
        disp(['Renaming volumes of Session ',Session_name{i,1},', DONE!'])
end
disp('Renaming Completed!!');

disp('Copying file in 17_1_Segment folder!');
disp('End of Module 4!!');

cd(path_general)
save('Variables_M4','Session_name')

diary off

spm_input('Between Session processing completed. Press ok to proceed with Segmentation of the Mean BS Functional Reference volume for the inner-brain BS mask calculation in module 5: ','-1','bd','Ok!');

%% MODULO 5:
disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||       ____  ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||          //     ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||     //___   ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||            //   ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||    ___//    ')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance imaging fMRI -------------------')
disp('--------------------------- Between Session Processing---------------------------')
disp('-----------From Segmentation of the Mean BS Functional Reference-----------------')
disp('--------to application of inner-brain BS mask and Deformation parameters---------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                       .\( / / . ) _)_)_______/                                  ')
disp('                        .( ( (  (__/__/  .                                       ')
disp('                          \._)_). )                                              ')
disp('                           ._)_)./                                               ')
disp('                           ._)_)./                                               ')
disp('                                                                                 ')

%% DECISION ABOUT GROUP ANALYSIS:
% In case of One session the processing start again here afer module 2.

else
spm_input('You have only one session therefore M3 and M4 modules will be skipped and their folders will be deleted.','-1','bd','Ok!');    
%copio il primo file della prima sessione dentro M5_BS_01_SegMeanRefVol
from=[path_to_M2_WS_04_Rename,'/',Session_name{1,1},'/','r1m',Session_name{1,1},'_001.nii'];
to=[path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii'];
copyfile(from,to)
% Deleting useless folder:
fol2delete={'M3_WS_01_SegRefVols','M3_WS_02_MaskRefVols','M3_WS_03_RealignMaskRefVols','M3_WS_04_Rename','M4_BS_01_Realign_Reslice','M4_BS_02_Scrub','M4_BS_03_Rename'};
for  i=1:length(fol2delete)
    cd(path_general)
    rmdir(fol2delete{1,i});
end
end

cd(path_general)
diary Module5_logfile
GA=spm_input('Choose if you want to proceed with:' ,'-1','bd','Group Analysis|Single subject',[1,0]);

if GA==1

 folder2beremoved={'M5_BS_01_SegMeanRefVol','M5_BS_02_MaskMeanRefVol','M5_BS_03_NormMaskMeanRefVol','M6_BS_01_MaskAllVols','M6_BS_02_NormMaskAllVols'}; 
 folder2becreated={'M5_GR_BS_01_SegMeanRefVol','M5_GR_BS_02_MaskMeanRefVol','M5_GR_BS_03_NormMaskMeanRefVol','M6_GR_BS_01_MaskAllVols','M6_GR_BS_02_NormMaskAllVols'};

    for i=1:5
            movefile(['/home/dati/Nico_auto_15T/Dati_fmri/',folder2beremoved{1,i}],['/home/dati/Nico_auto_15T/Dati_fmri/',folder2becreated{1,i}])
    end
        
    % Changing path

            %Module5
            path_to_M5_BS_01_SegMeanRefVol=[path_general,'M5_GR_BS_01_SegMeanRefVol'];
            path_to_M5_BS_02_MaskMeanRefVol=[path_general,'M5_GR_BS_02_MaskMeanRefVol'];
            path_to_M5_BS_03_NormMaskMeanRefVol=[path_general,'M5_GR_BS_03_NormMaskMeanRefVol'];
            %Module6
            path_to_M6_BS_01_MaskAllVols=[path_general,'M6_GR_BS_01_MaskAllVols'];
            path_to_M6_BS_02_NormMaskAllVols=[path_general,'M6_GR_BS_02_NormMaskAllVols'];
            
            cd(path_to_M1_PP_01_OrigVol);
            clear list
            list=dir();
            list(ismember({list.name},{'.','..'}))=[];
            if length(list)>1
            from=[path_to_M4_BS_03_Rename,'/',Session_name{1,1},'/','r2mr1m',Session_name{1,1},'_0001.nii'];
            to=[path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii'];
            copyfile(from,to)
            end
            
spm_input('You chose the Group analysis!','-1','bd','Ok!');

else
            cd(path_to_M1_PP_01_OrigVol);
            clear list
            list=dir();
            list(ismember({list.name},{'.','..'}))=[];
            if length(list)>1    from=[path_to_M4_BS_03_Rename,'/',Session_name{1,1},'/','r2mr1m',Session_name{1,1},'_0001.nii'];
                to=[path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii'];
                copyfile(from,to)
            end
    
spm_input('You chose the single subject analysis!','-1','bd','Ok!');

end

%% Segment Mean Ref and ask for GW e GW inner outer:

spm_input('Next you will be asked for the Gestational Week in order to proceed with the Segmentation of the Mean BS Functional Reference volume for the inner-brain BS mask and Deformation parameters calculation. ','-1','bd','Ok!');

m='Y';
while(1)

if m=='Y'

GW_winner=spm_input('Select Gestation Week for classes C1:C7:','-1','m','21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36|37',[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37],30);
GWin_out_winner=spm_input('Select Gestation Week for Inner/Outer brain:','-1','m','21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36|37',[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37],30);
GW_winner=num2str(GW_winner);
GW_winner=cellstr(GW_winner);
GWin_out_winner=num2str(GWin_out_winner);
GWin_out_winner=cellstr(GWin_out_winner);

clear input path3D matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc1_34GW_CorticalPlate_Cerebellum_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc2_34GW_WM_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc3_34GW_CSF_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc4_34GW_DGM_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc5_34GW_Hippocampi_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc6_34GW_Amygdala_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).tpm = {strcat(path_general,'Template/Template_Priors_Seg/34GW/04_Smoothing_Classi/sc7_34GW_Brainstem_Gholipour_2017.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(7).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(7).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(8).tpm = {strcat(path_general,'Template/Template_Priors_Seg/Inbrains_and_Outbrains/34GW/src8_InnerBrain_34GW.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(8).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(8).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(8).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(9).tpm = {strcat(path_general,'Template/Template_Priors_Seg/Inbrains_and_Outbrains/34GW/src9_OuterBrain_34GW.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(9).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(9).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(9).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'none';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 2;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
input=([path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii,1']);% input images
matlabbatch{1}.spm.spatial.preproc.channel.vols{1,1} = input;

% Adding Gestation Week to get the right Template
for i=1:9
    TMp=matlabbatch{1}.spm.spatial.preproc.tissue;
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm=strrep(TMp(i).tpm,'34',GW_winner{1,1});
    if (i==8 || i==9) % In and Outer brain
    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm=strrep(TMp(i).tpm,'34',GWin_out_winner{1,1});
    end
end
disp(['Segmenting reference volume session ',Session_name{1,1},' .']);
spm_jobman('run',matlabbatch);
disp(['Segmentation reference volume of session ',Session_name{1,1},' Completed!!']);

%% Inverto la maschera c9 in sc9

disp('Inverting c9 mask: ');

% copio il file c9:
from=[path_to_M5_BS_01_SegMeanRefVol,'/c9*'];
to=path_to_M5_BS_02_MaskMeanRefVol;
copyfile(from,to);

% copio il file :
from=[path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii'];
to=path_to_M5_BS_03_NormMaskMeanRefVol;
copyfile(from,to)

clear matlabbatch input
input=cellstr([path_to_M5_BS_02_MaskMeanRefVol,'/c9r2mr1m',Session_name{1,1},'_0001.nii,1']);% input images
output=cellstr(path_to_M5_BS_02_MaskMeanRefVol);% input images
matlabbatch{1}.spm.util.imcalc.input = {'Text to change'};
matlabbatch{1}.spm.util.imcalc.output = 'c9_mask';
matlabbatch{1}.spm.util.imcalc.outdir = {'Text to change'};
matlabbatch{1}.spm.util.imcalc.expression = '1-i1';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: c9_mask', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{1}.spm.util.imcalc.input = input;
matlabbatch{1}.spm.util.imcalc.outdir=output;

spm_jobman('run',matlabbatch);

disp('Inversion of c9 mask Completed!!');
cd(path_to_M5_BS_02_MaskMeanRefVol)
delete ('c9r2*')

clear list
cd(path_to_M5_BS_02_MaskMeanRefVol)
list=dir('*nii');
spm_check_registration(list.name,[path_to_M5_BS_01_SegMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii']);
disp('Check the inner-brain mask for the Mean Functional Reference Image.')
spm_input('Mean Funct. Ref. Mask review:','-1','b','Continue|--',[1,0]);

%% Applico la masschera sc9 alla Mean:

% Copying sc9 mask file:
from=[path_to_M5_BS_02_MaskMeanRefVol,'/sc9*'];
to=path_to_M5_BS_03_NormMaskMeanRefVol;
copyfile(from,to);

% Copying y*.nii file:
from=[path_to_M5_BS_01_SegMeanRefVol,'/y*'];
to=path_to_M5_BS_03_NormMaskMeanRefVol;
copyfile(from,to);

clear matlabbatch input output

matlabbatch{1}.spm.util.imcalc.input = {'Text to change'};
matlabbatch{1}.spm.util.imcalc.output = 'strip_mean_ref_funct';
matlabbatch{1}.spm.util.imcalc.outdir = {'Text to change'};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.spatial.normalise.write.subj.def = {'Text to change'};
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Image Calculator: ImCalc Computed Image: strip_mean_ref_funct', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
                                                          NaN NaN NaN];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [NaN NaN NaN];
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

input1=([path_to_M5_BS_03_NormMaskMeanRefVol,'/r2mr1m',Session_name{1,1},'_0001.nii,1']);
input2=([path_to_M5_BS_03_NormMaskMeanRefVol,'/sc9_mask.nii,1']);
output=path_to_M6_BS_01_MaskAllVols;
output2=([path_to_M5_BS_03_NormMaskMeanRefVol,'/y_r2mr1m',Session_name{1,1},'_0001.nii']);
matlabbatch{1}.spm.util.imcalc.input{1,1} = input1;
matlabbatch{1}.spm.util.imcalc.input{2,1} = input2;
matlabbatch{1}.spm.util.imcalc.outdir=output;
matlabbatch{1, 2}.spm.spatial.normalise.write.subj.def{1, 1}=output2;
spm_jobman('run',matlabbatch);

%Checkreg Wstrip Template pulito della GW_winner scelta
clear Wstrip2show temp2show
disp('Check the Normalized Mean Functional Reference Image with subject GW template space.')
Wstrip2show=[path_to_M5_BS_02_MaskMeanRefVol,'/wstrip_mean_ref_funct.nii'];
temp2show=[path_general,'Template/Template_orig/STA',GW_winner{1,1},'.nii'];
spm_check_registration(Wstrip2show,temp2show);
%chk=input('Press 1 to proceed with next module or press 0 to Segment again: ');
chk=spm_input('Perform Segmentation again?','-1','b','Yes|No',[0,1]);
                    if (chk==1)
                        m ='N';
                        break
                    else
                        m ='Y';
                    end
else
break
end
end

disp('Copying Completed!');
disp('End of Module 5!');

% GWin_out_winner GW_winner
cd(path_general)
save('Variables_M5.mat','GW_winner','GWin_out_winner','Session_name')
diary off

spm_input('Application of inner-brain BS mask and deformation parameters to the Mean Functional Reference volume Completed. Press ok to proceed with application of deformation parameters to all of the functional scans in module 6: ','-1','bd','Ok!');

%% MODULO 6:
disp('---------------------------------------------------------------------------------')
disp('---------------------------------------------------------------------------------')
disp('|||\\\\///|||    //°°°\\     |||\\    |||     |||  |||       |||||||       ____  ')
disp('||| \\\///|||   //     \\    ||| \\   |||     |||  |||       |||          //     ')
disp('|||  \\// |||  ((       ))   |||  ||  |||     |||  |||       |||||||     //___   ')
disp('|||   \/  |||   \\     //    ||| //   |||     |||  |||       |||        //  //   ')
disp('|||       |||    \\___//     |||//    |||||||||||  |||||||   |||||||   //__//    ')
disp('---------------------------------------------------------------------------------')
disp('------------------ Functional Magnetic Resonance imaging fMRI -------------------')
disp('-------------------------- Between Session Processing  --------------------------')
disp('---From application of deformation parameters to all of the functional scans-----')
disp('-----------to smoothing of all inner-brain BS masked functional scans------------')
disp('------------------------------END BETWEEN SESSION--------------------------------')
disp('                           . . . . . . . . .                                     ')
disp('                         ./ ( (   ) . ).(   )).\.                                ')
disp('                       .)  ) ) (  . (  . ) (  ..\.                               ')
disp('                      .() ( (   ) . ) . (  ))  ( |                               ')
disp('                      .(( ) ) ) (  . ( ..)  ()__/                                ')
disp('                       .\( / / . ) _)_)_______/                                  ')
disp('                        .( ( (  (__/__/  .                                       ')
disp('                          \._)_). )                                              ')
disp('                           ._)_)./                                               ')
disp('                           ._)_)./                                               ')
disp('                                                                                 ')
%%

disp('Copying r2* files from 12_2_fMRI_Rename to 20_Strip_Norm_Funct')
cd(path_to_M1_PP_01_OrigVol);
list=dir();
% Get vector to recognize folders only:
list(ismember({list.name},{'.','..'}))=[];

if length(list)>1
for i=1:size(Session_name,1)
    from=[path_to_M4_BS_03_Rename,'/',Session_name{i,1}];
    to=[path_to_M6_BS_01_MaskAllVols,'/',Session_name{i,1}];
    copyfile(from,to);
end
disp('Copying Completed!');
else
mkdir([path_to_M6_BS_01_MaskAllVols,'/',Session_name{1,1}])
list_one_sess=dir([path_to_M2_WS_04_Rename,'/',Session_name{1,1},'/*nii']);
  for i=1:length(list_one_sess)
    from=[path_to_M2_WS_04_Rename,'/',Session_name{1,1},'/',list_one_sess(i).name];
    to=[path_to_M6_BS_01_MaskAllVols,'/',Session_name{1,1},'/r2mr1m',Session_name{1,1},'_',num2str(i,'%04d'),'.nii'];
    copyfile(from,to);
  end
disp('Copying Completed!');
end

disp('Copying sc9_mask and strip_meaa_func from M5_BS_03_NormMaskMeanRefVol to M6_BS_01_MaskAllVols')
from=[path_to_M5_BS_03_NormMaskMeanRefVol,'/s*'];
to=path_to_M6_BS_01_MaskAllVols;
copyfile(from,to);
disp('Copying Completed!');

disp('Copying wstrip_mean_ref_funct.nii from M5_BS_02_MaskMeanRefVol to M6_BS_02_NormMaskAllVols')
from=[path_to_M5_BS_02_MaskMeanRefVol,'/w*'];
to=path_to_M6_BS_02_NormMaskAllVols;
copyfile(from,to);
disp('Copying Completed!');

disp('Copying y_file from M5_BS_03_NormMaskMeanRefVol to M6_BS_02_NormMaskAllVols')
from=[path_to_M5_BS_03_NormMaskMeanRefVol,'/y*'];
to=path_to_M6_BS_02_NormMaskAllVols;
copyfile(from,to);


%% Applying mask to all volume:

cd(path_general)
diary Module6_logfile
cd(path_to_M6_BS_02_NormMaskAllVols);
for i=1:size(Session_name,1)
mkdir(Session_name{i,1});
end

count=0;
clear matlabbatch
for iii=1:size(Session_name,1)

    cd([path_to_M6_BS_01_MaskAllVols,'/',Session_name{iii,1}]);
    list_to_mask=dir('r2*');
    niiPrefix='r2mr1m';
    for id = (1+count):(length(list_to_mask)+count)

        matlabbatch{1}.spm.util.imcalc.input = {
                                                strcat(path_to_M6_BS_01_MaskAllVols,'/',Session_name{iii,1},'/',char(niiPrefix),Session_name{iii,1},'_',num2str(id,'%04d'),'.nii,1');
                                                strcat(path_to_M6_BS_01_MaskAllVols,'/sc9_mask.nii,1')
                                                };
        matlabbatch{1}.spm.util.imcalc.output = char(strcat('m',niiPrefix,Session_name{iii,1},'_', num2str(id,'%04d'), '.nii'));
        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{iii,1})};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
    end
    count=count+length(list_to_mask);
    clear matlabbatch list_to_mask
end

%Visualization of masked images:
disp('Visual inspection of the Mean Functional Reference image and random masked BS images.')
clear ref2show masked2show img2show
ref2show{1,1}=[path_to_M5_BS_02_MaskMeanRefVol,'/strip_mean_ref_funct.nii,1'];
c=0;
for i=1:size(Session_name,1)
    clear randomnum
    cd([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1}]);
    all_files=dir('mr2*');
        if size(all_files,1)~=0
            randomnum=round(1+rand(1,2)'*(size(all_files,1)-1));
        end
    for j=1:2
         img2show{j+2*c,1}=[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',all_files(randomnum(j)).name,',1'];%mask
    end
    c=c+1;
end
spm_check_registration(ref2show{1,1},img2show{:,:});
spm_input('Check masked volumes:','-1','b','Continue|--',[1,0]);

%Applying batch:
clear matlabbatch input output
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
                                                          NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

matlabbatch{1, 1}.spm.spatial.normalise.write.subj.resample='';

% y file..
matlabbatch{1, 1}.spm.spatial.normalise.write.subj.def{1, 1}= strcat(path_to_M6_BS_02_NormMaskAllVols,'/y_r2mr1m',Session_name{1,1},'_0001.nii');

%All volumes..
count=0;
for iii=1:size(Session_name,1)

    cd([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{iii,1}]);
    files=dir('mr2*');

    for id = (1+count):(length(files)+count)

        matlabbatch{1, 1}.spm.spatial.normalise.write.subj.resample{id, 1}=strcat(path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{iii,1},'/mr2mr1m',Session_name{iii,1},'_',num2str(id,'%04d'),'.nii,1');

    end
    count=count+length(files);
    clear  files
end

disp('Writing all smoothed normalized masked images...')
spm_jobman('run',matlabbatch);


%% Renaming final dataset.

cd([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{1,1}])
file2bedel= dir([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{1,1},'/*_0001.nii']);

if length(file2bedel)==3

disp(['Delating swmr2mr1m',Session_name{1,1},'_0001.nii wmr2mr1m',Session_name{1,1},'_0001.nii mr2mr1m',Session_name{1,1},'_0001.nii']);
    
    delete(file2bedel.name)
    disp('Renaming all sw* volumes : ');
    count=0;
    clear id num files
    for i=1:size(Session_name,1)
        f = [path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1}];
        cd(f)
        files = dir('sw*');
        ii=1;
             for id = (1+count):(length(files)+count)
                 num=['swmr2mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
                 movefile([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',files(ii).name],[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',num]);
                 ii=ii+1;
             end
            count=count+length(files);
            clear files
            disp(['Renaming sw* volumes of Session ',Session_name{i,1},', DONE!'])
    end
    disp('Renaming Completed!!');

        disp('Renaming all wm* volumes : ');
    count=0;
    clear id num files
    for i=1:size(Session_name,1)
        f = [path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1}];
        cd(f)
        files = dir('wm*');
        ii=1;
             for id = (1+count):(length(files)+count)
                 num=['wmr2mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
                 movefile([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',files(ii).name],[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',num]);
                 ii=ii+1;
             end
            count=count+length(files);
            clear files
            disp(['Renaming wm* volumes of Session ',Session_name{i,1},', DONE!'])
    end
    disp('Renaming Completed!!');
    
        disp('Renaming all m* volumes : ');
    count=0;
    clear id num files
    for i=1:size(Session_name,1)
        f = [path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1}];
        cd(f)
        files = dir('m*');
        ii=1;
             for id = (1+count):(length(files)+count)
                 num=['mr2mr1m',Session_name{i,1},'_',num2str(id,'%04d'),'.nii'];
                 movefile([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',files(ii).name],[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',num]);
                 ii=ii+1;
             end
            count=count+length(files);
            clear files
            disp(['Renaming m* volumes of Session ',Session_name{i,1},', DONE!'])
    end
    else 
    disp('Mean_ref_vol was not present in session 1.')
end
    disp('Renaming Completed!!');

%Visualization of masked images:
disp('Visual inspection of normalized masked images:')
clear ref2show masked2show img2show
Temp2show{1,1}=[path_general,'Template/Template_orig/STA',GW_winner{1,1},'.nii,1'];
wstripmean2show{1,1}=[path_to_M6_BS_02_NormMaskAllVols,'/wstrip_mean_ref_funct.nii'];
c=0;
    for i=1:size(Session_name,1)
        clear randomnum
        cd([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1}]);
        all_files=dir('wmr2*');
            if size(all_files,1)~=0
                randomnum=round(1+rand(1,2)'*(size(all_files,1)-1));
            end
        for j=1:2
             img2show{j+2*c,1}=[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',all_files(randomnum(j)).name,',1'];%mask
        end
        c=c+1;
    end
    spm_check_registration(Temp2show{1,1},wstripmean2show{1,1},img2show{:,:});
    spm_input('Check normalized volumes:','-1','b','Continue|--',[1,0]);
    
%Creating list fot visualization:
disp('Creation of the 4D timeseries including all smoothed normalized images retained for stats!!')
disp('4D visualization on SPM!')
FileNames_tot='';
for i=1:size(Session_name,1)
    clear list2show_sw FileNames
    list2show_sw=dir([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/sw*']);
    for k=1:length(list2show_sw)
       FileNames{k,:}=[path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{i,1},'/',list2show_sw(k).name];
    end
    FileNames_tot=[FileNames_tot;FileNames];
end
spm_ov_browser('ui',char(FileNames_tot))
spm_input('Check final 4D movie:','-1','b','Continue|--',[1,0]);

% 4D timeseries creation
%SMP 3D to 4D

    disp('Creating 4D files!')
    clear matlabbatch
    cd(path_to_M6_BS_02_NormMaskAllVols)
    matlabbatch{1}.spm.util.cat.vols=FileNames_tot;
    matlabbatch{1}.spm.util.cat.name = '4D.nii';
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.RT = NaN;
    spm_jobman('run',matlabbatch);
    delete('*.mat')
    movefile([path_to_M6_BS_02_NormMaskAllVols,'/',Session_name{1,1},'/4D.nii'],path_to_M6_BS_02_NormMaskAllVols)
    
disp('4D timeseries creation completed!!')

disp('Well done!!! You have completed all modules!')
disp('You are ready to analyze your images!!!')
disp('Now it is time for a beer, cheers!')
diary off

spm_input('Well done!!! You have completed all modules! You are ready to analyze your images!! ','-1','bd','Ok!');
