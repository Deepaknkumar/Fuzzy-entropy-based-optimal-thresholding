clear all;
% clc;
% close all;
I = imread('mandril_color.jpg');        
figure,imshow(I),title('Original Image')

Number_of_levels = 4;       %%% discrete

sgrays = double(I);
size_of_image = size(sgrays(:,:,1));
ssize_of_image = size(sgrays(:,:,1));

caddterm = ssize_of_image(1)*ssize_of_image(2)*8;

for q=1:3
for i=1:256
    bp = sgrays(:,:,q)>i-1;      %%% Mark every Value if it is greater than pixel value i-1 
    negmat = bp-1;       %%% Obtain -1 matrix for every value equals zero
    b = negmat+bp;       %%% Put +1 if Value is greater than i-1 else -1
    
    b(1,:) = 0;
    b(end,:) = 0;
    b(:,1) = 0;
    b(:,end) = 0;
    
    pienergy = 0;          %%% To store Energy of Every pixel value
    for j=2:size_of_image(1)-1
        for k=2:size_of_image(2)-1
            pienergy = pienergy + (b(j,k).^2) - sum(sum(b(j,k).*b(j-1:j+1,k-1:k+1)));
        end
    end
    E(q,i) = pienergy + caddterm;
    clc;
    fprintf(' Matrix %d Energy Density Function %d %% Completed......\n',q,uint8((i*100)/256))
end
end

threshlevels = [];
 n=15;          %number of Particles
 nv = Number_of_levels*3;        %number of variables
 lim = [zeros(nv,1) 255.*ones(nv,1)]; %lower and upper bound of variables
 vcf = 2;       %velocity clamping factor
 cc = 3;        %cognitive constant
 sc = 2;        %social constant
 miniw = .4;    %Min Inertia weight
 maxiw = .9;    %Max Inertia weight
 numiter = 700;
 
for k=1:3
    E(k,:) = imhist(I(:,:,k));
%     E(k,:) = E(k,:)/sum(E(k,:));
    [Entp(k,:),values(k,:),convplot(k,:)] = Particle_Swarm_Optimizationasc20(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,numiter,E(k,:));       %%% threshold levels and minimum Entropy obtained
    
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
fprintf('Objective Function Value = %f\n',Entp(1));