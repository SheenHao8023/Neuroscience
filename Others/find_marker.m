clc;clear;

pathn = 'D:\matlab workspace\fnirs_XXX\fnirsdata1_nirs2mat'
cd(pathn);
filename_all=dir(fullfile(pathn,'*.mat'))
for sub=1:45
   load(filename_all(sub).name)
   for m=1:4
   mark1_all{sub,m}=find(nirs_data.vector_onset==m)
   mark1_size{sub,m}=size(mark1_all{sub,m})
   end
end
 save mark1_all
 
 %{[92;1231;2077;3360;4159;4887;5610],[399;1465;2289;3559;4433;5142;5890],[571;1588;2700;3832;4567;5261;6059],[1018;1816;3136;3948;4647;5382;6137];[57;882;1652;2432;3225;4004;4833],[319;1150;1998;2768;3555;4320;5137],[441;1263;2140;2930;3686;4459;5290],[655;1775;2233;3036;3787;4590;5385]}