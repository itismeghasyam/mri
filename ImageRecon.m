directory=pwd; 
cd s6190;%enter file directory for dicom images
for i=1:172
    str=strcat('*.MRDC.',num2str(i));
    listing(i) = dir(str);
end

dImage=zeros(256,256,length(listing));
for i=1:length(listing)
    dImage(:,:,i)=dicomread(listing(i).name);
end
cd (directory);

figure(1);

imshow3Dfull(dImage);