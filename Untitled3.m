srcFiles = dir('C:\Users\BERLIN\Desktop\Bone fracture detection\dataset\non-fracture\*.jpg');  % the folder in which ur images exists
for i = 1 : length(srcFiles)
filename = strcat('C:\Users\BERLIN\Desktop\Bone fracture detection\dataset\non-fracture\',srcFiles(i).name);
im = imread(filename);
k=imresize(im,[300,300]);
newfilename=strcat('C:\Users\BERLIN\Desktop\Bone fracture detection\dataset\non-fracture\Non-fracture\',srcFiles(i).name);
imwrite(k,newfilename,'jpg');
end