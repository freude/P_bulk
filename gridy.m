function b=gridy()

c=textread(strcat(pwd,'/dis_scr/mesh_sample.mesh'),'%s','delimiter','\n');
I1=find(~cellfun(@isempty,strfind(c,'Vertices')))+2;
I2=find(~cellfun(@isempty,strfind(c,'Tetrahedra')))-2;
b=cell2mat(cellfun(@str2num,c(I1:I2), 'UniformOutput', false));

b=b(:,1:3);





