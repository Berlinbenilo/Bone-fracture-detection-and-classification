fracture = gTruth(:,1:2);
fracture.imageFilename = fullfile(fracture.imageFilename);
acfDetector = trainACFObjectDetector(fracture,'NegativeSamplesFactor',2);
save detector1.mat