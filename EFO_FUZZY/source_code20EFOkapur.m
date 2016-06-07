clear all;
% clc;
% close all;
I = imread('mandril_color.jpg');        
figure,imshow(I),title('Original Image')

Number_of_levels = 8;       %%% discrete values

sgrays = double(I);
size_of_image = size(sgrays(:,:,1));
ssize_of_image = size(sgrays(:,:,1));

threshlevels = [];

N_var=Number_of_levels;
N_emp=15;
Max_gen=700;
minval= 0;
maxval= 255;
R_rate=0.3;
Ps_rate=0.2;
P_field=0.1;
N_field=0.45;

for k=1:3
    E(k,:) = imhist(I(:,:,k));
%     E(k,:) = E(k,:)/sum(E(k,:));
    [Entp(k,:),threshlevels(k,:),convplot(k,:)] = EFOsc20kapur(N_var,N_emp,Max_gen,minval,maxval,R_rate,Ps_rate,P_field,N_field,E(k,:));       %%% threshold levels and minimum Entropy obtained
        
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