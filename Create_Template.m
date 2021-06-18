function [] = Create_Template(path_original_tmp,path_to_M1_PP_04_Rename,path_to_M2_WS_01_Mask,Session_name)

% This function allows to create fMRI template for each Gestational Week in the User specific space (image dimension and voxel resolution).
% INPUT:
    %1) path_original_tmp: path of original template; Example ('/home/dati/Nico_auto_15T/Dati_fmri/Template_4_bet/Original_bin_Template' )
    %2) Session_name.
% Reference image is taken from folder 04_fMRI_Rename. It will be the
% first volume of the first session.

cd(path_original_tmp);
for i=1:size(Session_name,1)
    mkdir(Session_name{i,1})
    from=[path_original_tmp,'/*.nii'];
    to=[path_original_tmp,'/',Session_name{i,1}];
    copyfile(from,to)
end

for hh=1:size(Session_name,1)
clear input output option cmd matlabbatch llist
cd([path_original_tmp,'/',Session_name{hh,1}]);
matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr([path_to_M1_PP_04_Rename,'/',Session_name{hh,1},'/',Session_name{hh,1},'_001.nii']);
    for i=1:17
            matlabbatch{1}.spm.spatial.coreg.write.source = cellstr([path_original_tmp,'/',Session_name{hh,1},'/GW_',num2str(i+20),'.nii']);
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'c';
            spm_jobman('run',matlabbatch);
    end
end

    %smooth isotropico
for hh=1:size(Session_name,1)
     clear input output option cmd matlabbatch
     cd([path_original_tmp,'/',Session_name{hh,1}]);
    for i=1:17
        clear input output option cmd matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = cellstr([path_original_tmp,'/',Session_name{hh,1},'/cGW_',num2str(i+20),'.nii,1']);
        matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2]; % Prima era 2 2 2 (modificato in data: 12 Febb)
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 't';
        spm_jobman('run',matlabbatch);
    end
cd([path_original_tmp,'/',Session_name{hh,1}]);
%delete('*cGW')
end

for hh=1:size(Session_name,1)
     clear llist input output option cmd matlabbatch
    %Binarizzo dopo lo smooth
    cd([path_original_tmp,'/',Session_name{hh,1}]);
    for i=1:17
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[path_original_tmp,'/',Session_name{hh,1},'/tcGW_',num2str(i+20),'.nii']};
        matlabbatch{1}.spm.util.imcalc.output = ['scGW_',num2str(i+20)];
        matlabbatch{1}.spm.util.imcalc.outdir = {[path_to_M2_WS_01_Mask,'/',Session_name{hh,1},'/',Session_name{hh,1},'_template']};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1>',num2str(0.2)];%Diminuire thr per ingrandire le maschere;
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
    end
    delete('*tcGW')
end

%% 2 part
    %smooth anisotropico
    ker=[4 4 4;6 6 6;8 8 8];
for hh=1:size(Session_name,1)
     clear input output option cmd matlabbatch
     cd([path_original_tmp,'/',Session_name{hh,1}]);
    for i=1:17
        for jj=1:3
        clear input output option cmd matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = cellstr([path_original_tmp,'/',Session_name{hh,1},'/cGW_',num2str(i+20),'.nii,1']);
        matlabbatch{1}.spm.spatial.smooth.fwhm = ker(jj,:); % Prima era 2 2 2 (modificato in data: 12 Febb)
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = ['anp_',num2str(jj),'_'];
        spm_jobman('run',matlabbatch);
        end
    end
%cd([path_original_tmp,'/',Session_name{hh,1}]);
end

for hh=1:size(Session_name,1)
     clear llist input output option cmd matlabbatch
    %Binarizzo dopo lo smooth
    cd([path_original_tmp,'/',Session_name{hh,1}]);
    for i=1:17
        for jj=1:3
        clear matlabbatch
        matlabbatch{1}.spm.util.imcalc.input = {[path_original_tmp,'/',Session_name{hh,1},'/anp_',num2str(jj),'_cGW_',num2str(i+20),'.nii']};
        matlabbatch{1}.spm.util.imcalc.output = ['banp_',num2str(jj),'_',num2str(i+20)];
        matlabbatch{1}.spm.util.imcalc.outdir = {[path_to_M2_WS_01_Mask,'/',Session_name{hh,1},'/',Session_name{hh,1},'template_anp']};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1>',num2str(0.3)];
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch);
        end
    end
end

for hh=1:size(Session_name,1)
cd([path_original_tmp,'/',Session_name{hh,1}]);
delete('*')
end
cd(path_original_tmp)
for hh=1:size(Session_name,1)
rmdir(Session_name{hh,1})
end
end















