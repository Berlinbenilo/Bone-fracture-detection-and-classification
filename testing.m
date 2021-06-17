load detector1.mat
[filename,filepath] = uigetfile('file selector');
fullpath = strcat(filepath,filename);
img = imread(fullpath);
[bboxes,scores] = detect(acfDetector,img);

flag=0;
for i = 1:length(scores)
    if scores(i,1) >= 85.2 && scores(i,1)<= 91
        annotation = sprintf('Fracture: %.1f',scores(i));
        img = insertObjectAnnotation(img,'rectangle',bboxes(i,:),annotation,'LineWidth',4);
%         flag = flag+1;
%         bbox1(flag,1:4) = bboxes(i,1:4);
%         score1(flag,1) = scores(i,1);
%         boundingBoxArea(flag,1) = bbox1(flag,3) * bbox1(flag,4);
    end
end

% for bb = 1:length(boundingBoxArea)
%     if boundingBoxArea(bb,1) == max(boundingBoxArea)
%         annotation = sprintf('Fracture: %.1f',score1(bb));
%         img = insertObjectAnnotation(img,'rectangle',bbox1(bb,:),annotation,'LineWidth',4);
%     end
% end

figure
imshow(img)