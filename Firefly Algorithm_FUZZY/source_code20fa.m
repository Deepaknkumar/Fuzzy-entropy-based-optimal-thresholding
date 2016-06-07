clear all;
% clc;
% close all;
I = imread('peppers_color.jpg');        
figure,imshow(I),title('Original Image')

Number_of_levels = 5;       %%% Input Number of Levels Required(Currently valid for two only)

sgrays = double(I);
size_of_image = size(sgrays(:,:,1));
ssize_of_image = size(sgrays(:,:,1));

threshlevels = [];
mum = [];
mum = [];
mub = [];

n = 15;
for k=1:3
    E(k,:) = imhist(I(:,:,k));
%     E(k,:) = E(k,:)/sum(E(k,:));
    [Entp(k,:),values(k,:),convplot(k,:)] = fireflysc20(n,Number_of_levels,E(k,:));       %%% threshold levels and minimum Entropy obtained
        
        for i=1:Number_of_levels
           threshlevels(k,i) = threshExtractersc20(values(k,((i-1)*3)+1:3*i));
        end
        
    pE = E(k,:);
    figure,plot(pE),title('Probability Distribution Curve')
    xlabel('Pixel Values')
    ylabel('Respective values')
    hold on
    stem(threshlevels(k,:),(ones(1,length(threshlevels(k,:)))),'fill','r')
    hold off
end

for i=1:3
    figure,plot(convplot(i,:)),title('Convergence Plot')
    xlabel('Iterations')
    ylabel('Fitness Function Value')
end

fprintf('\nThe threshold values for %d levels:  \n',Number_of_levels)

threshlevels
%%------------------Segmenting out Image --------------------------%%
segimage = [];
for iv = 1:3
for i=1:ssize_of_image(1)
    for j=1:ssize_of_image(2)
       if  sgrays(i,j,iv) <= threshlevels(iv,1)
           segimage(i,j,iv) = 0;
       end
    end
end

for t=2:Number_of_levels
    for i=1:ssize_of_image(1)
        for j=1:ssize_of_image(2)
            if  sgrays(i,j,iv) > threshlevels(iv,t-1) && sgrays(i,j,iv) <= threshlevels(iv,t)
              segimage(i,j,iv) = threshlevels(iv,t-1);
            end
        end
    end
end
        
for i=1:ssize_of_image(1)
    for j=1:ssize_of_image(2)
       if  sgrays(i,j,iv) > threshlevels(iv,Number_of_levels) && sgrays(i,j,iv)< 256
           segimage(i,j,iv) = threshlevels(iv,Number_of_levels);
       end
    end
end
end

main_image = uint8(segimage);
figure,imshow(main_image),title('Segmented Image')

for i=1:3
D = ((segimage(:,:,i)) -(sgrays(:,:,i))).^2;
mse(i) = sum(D(:))/numel(main_image(:,:,i));
psnr(i) = 10*log10(255*255/mse(i));

K=[0.6 0.6];
 [mssim(i), ssim_map] = ssim_index(segimage(:,:,i),sgrays(:,:,i), K);  
 [FSIM(i), FSIMc] = FeatureSIM(segimage(:,:,i),sgrays(:,:,i));
end

fprintf('\nMean Square Error = %f\n',sum(mse)/3);
fprintf('Peak Signal to Noise Ratio = %f\n',sum(psnr)/3);
fprintf('MSSIM = %f\n',sum(mssim)/3);
fprintf('FSIM = %f\n',sum(FSIM)/3);
fprintf('Objective Function Value = %f %f %f\n',Entp(1),Entp(2),Entp(3));