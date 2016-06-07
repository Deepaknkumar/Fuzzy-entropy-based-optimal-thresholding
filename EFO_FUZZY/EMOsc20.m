%Diego Oliva, Erik Cuevas, Gonzalo Pajares, Daniel Zaldivar, Valentín Osuna. 
%A Multilevel Thresholding algorithm using electromagnetism optimization
%Universidad Complutense de Madrid / Universidad de Guadalajara

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The algorithm was published as:
%Diego Oliva, Erik Cuevas, Gonzalo Pajares, Daniel Zaldivar, Valentín Osuna. 
%A Multilevel Thresholding algorithm using electromagnetism optimization, 
%Journal of Neurocomputing, 139, (2014), 357-381.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Intructions:
% I -> Original Image, could be RGB or Gray Scale
% level -> Number of threshold to find
% This version works with KAPUR as fitness function.



close all
clear all
% Se carga la imagen RGB o escala de grises
I = imread('test1.jpg');
% Se asigna el nivel de segmentacion 
    level = 2;


% Se obtienen los histogramas si la imagen es RGB uno por cada canal si es
% en escala de grises solamente un historgrama.
if size(I,3) == 1 %grayscale image
    [n_countR, x_valueR] = imhist(I(:,:,1));
elseif size(I,3) == 3 %RGB image
    %histograma para cada canal RGB
    [n_countR, x_valueR] = imhist(I(:,:,1));
    [n_countG, x_valueG] = imhist(I(:,:,2));
    [n_countB, x_valueB] = imhist(I(:,:,3));
end
Nt = size(I,1) * size(I,2); %Cantidad total de pixeles en la imagen RENG X COL
%Lmax niveles de color a segmentar 0 - 256
Lmax = 256;   %256 different maximum levels are considered in an image (i.e., 0 to 255)

% Distribucion de probabilidades de cada nivel de intensidad del histograma 0 - 256 
for i = 1:Lmax
    if size(I,3) == 1 
        %grayscale image
        probR(i) = n_countR(i) / Nt;
    elseif size(I,3) == 3 
        %RGB image    
        probR(i) = n_countR(i) / Nt;
        probG(i) = n_countG(i) / Nt;
        probB(i) = n_countB(i) / Nt;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametros del problema de segmentacion
%N_PAR = level  - 1; %number of thresholds (number of levels-1) (dimensiones)
N_PAR = level;
dim = N_PAR;  


for ii = 1:1  % for para pruebas estadisticas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros de la poblacion
%maximo de iteraciones
MAXITER = 200;
%m cantidad de puntos, n dimensiones en las cuales se trabaja
m = 50; %cantidad de miembros de la poblacion
n = dim; %dimensiones del problema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del espacio de busqueda
%Crea el espacio de busqueda
%u -> vector de limites superiores de cada dimension
%l -> vector de limites inferiores de cada dimension
%xR, xG, xB -> poblaciones inicializadas en cero

if size(I,3) == 1
    % Imagen en escala de grises
    u = ones(1,dim) * Lmax;
    l = ones(1,dim);
    xR = zeros(m,n);

elseif size(I,3) == 3
    % Imagen RGB
    u = ones(1,dim) * Lmax;
    l = ones(1,dim);

    xR = zeros(m,n);

 
    xG = zeros(m,n);
    

    xB = zeros(m,n);
end

    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %EMO Original               %%Version 4%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Generar y evaluar la poblacion que sera evolucionada durante la optimizacion
    %se generan m valores aleatorios y se evaluan en la funcion de fitness Ec. 4       
    if size(I,3) == 1
        % Imagen en escala de grises
        % Inicializacion aleatoria
        xR = incializa(m,n,u,l,xR);
        for si=1:length(xR)
           xR(si,:)=sort(xR(si,:)); 
        end
        % Evaluar poblacion en la funcion de fitness
        %[fitR, fitBestR] = fitnessIMG(I, m, Lmax, level, xR, probR);
        fitR = Kapur(m,n,level,xR,probR);
        
        % Elige el mejor elemento de la poblacion en base al fitness
        [aR, bR] = max(fitR); %Maximiza
        %gBestR = xR(bR, :);
        %gbestvalueR = fitR(bR);
        xBestR = xR(bR, :);
        fxBestR = fitR(bR);
        
    elseif size(I,3) == 3
        % Imagen RGB
        % Inicializacion aleatoria para cada canal R, G, B
        xR = incializa(m,n,u,l,xR);
        xG = incializa(m,n,u,l,xG);
        xB = incializa(m,n,u,l,xB);
        
        % Evalua la poblacion de cada canal en la funcion de fitness
        %[fitR, fitBestR, fitG, fitBestG, fitB, fitBestB] = fitnessIMG(I, m, Lmax, level, xR, probR, xG, probG, xB, probB);
        fitR = Kapur(m,n,level,xR,probR);
        fitG = Kapur(m,n,level,xG,probG);
        fitB = Kapur(m,n,level,xB,probB);
        
        % Se elige el mejor elemento de cada poblacion en base al fitness 
        [aR,bR] = max(fitR); %maximiza
        %gBestR=xR(bR,:);
        %gbestvalueR = fitR(bR);
        xBestR = xR(bR, :);
        fxBestR = fitR(bR);
        
        [aG,bG] = max(fitG); % maximiza
        %gBestG=xG(bG,:);
        %gbestvalueG = fitG(bG);
        xBestG = xG(bG, :);
        fxBestG = fitG(bG);
        
        [aB,bB] = max(fitB); % maximiza
        %gBestB = xR(bB,:);
        %gbestvalueB = fitB(bB);
        xBestB = xB(bR, :);
        fxBestB = fitR(bR);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %delta -> valor de vecindad de busqueda
    %LSITER -> valor  maximo de iteraciones para la busqueda
    delta = 0.025;
    LSITER = 4;
    cc = 0;
    it = 0;
    
    while it < MAXITER
    %while cc < MAXITER * 0.01
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if size(I,3) == 1
            % Imagen escala de grises
            [bR,xR,fitR,xBestR,fxBestR] = local(I,Lmax,level,LSITER,delta,m,n,l,u,xR,fitR, probR);
            FxR = calcFO(m,n,xR,fitR,xBestR,fxBestR);
            xR = moveO(FxR,xR,m,n,bR,l,u);
            for si=1:length(xR)
                xR(si,:)=sort(xR(si,:)); 
            end
            
            
            fitR = Kapur(m,n,level,xR,probR);
            BestvalR_fit(it + 1) = fxBestR;
            
        elseif size(I,3) == 3
            % Imagen RGB
            %[bR,xR,fitR,xBestR,fxBestR,bG,xG,fitG,xBestG,fxBestG,bB,xB,fitB,xBestB,fxBestB]=localO(I,Lmax,level,LSITER,delta,m,n,l,u,xR,fitR, probR, xG,fitG, probG, xB, fitB, probB);
            
            [bR,xR,fitR,xBestR,fxBestR] = local(I,Lmax,level,LSITER,delta,m,n,l,u,xR,fitR,probR);
            FxR = calcFO(m,n,xR,fitR,xBestR,fxBestR);
            xR = moveO(FxR,xR,m,n,bR,l,u);
            
            [bG,xG,fitG,xBestG,fxBestG] = local(I,Lmax,level,LSITER,delta,m,n,l,u,xG,fitG,probG);
            FxG = calcFO(m,n,xG,fitG,xBestG,fxBestG);
            xG = moveO(FxG,xG,m,n,bR,l,u);
            
            [bB,xB,fitB,xBestB,fxBestB] = local(I,Lmax,level,LSITER,delta,m,n,l,u,xB,fitB,probB);
            FxB = calcFO(m,n,xB,fitB,xBestB,fxBestB);
            xB = moveO(FxB,xB,m,n,bB,l,u);
            
            %[fitR, fitBestR, fitG, fitBestG, fitB, fitBestB] = fitnessIMG(I, m, Lmax, level, xR, probR, xG, probG, xB, probB);
            fitR = Kapur(m,n,level,xR,probR);
            fitG = Kapur(m,n,level,xG,probG);
            fitB = Kapur(m,n,level,xB,probB);
                       
            
            
            BestvalR(it + 1,:) = xBestR;
            BestvalG(it + 1,:) = xBestG;
            BestvalB(it + 1,:) = xBestB;
        end     
        
        it = it + 1;
        
    end
    
    
    if size(I,3) == 1
            BestvalR(ii) = fxBestR;
    elseif size(I,3) == 3    
            BestvalR(ii,:) = xBestR;
            BestvalG(ii,:) = xBestG;
            BestvalB(ii,:) = xBestB;
    end
    
    
    
end
    

        
    
    
if size(I,3) == 1 %grayscale image
    gBestR = sort(xBestR);
    Iout = imageGRAY(I,gBestR);
    Iout2 = mat2gray(Iout);
    intensity = gBestR(1:n-1) %best threshold
    STDR =  std(BestvalR)      %standar deviation of fitness
    MEANR = mean(BestvalR)     %mean of fitness
    PSNRV = PSNR(I, Iout)      %PSNR bestween the segmented and the original image
    fxBestR                    %Best fitness
    %Plot results in image
    figure
    imshow(Iout)
    figure
    imshow(I)
    
    %plot the thresholds over the histogram
    figure 
    plot(probR)
    hold on
    vmax = max(probR);
    for i = 1:n-1
        line([intensity(i), intensity(i)],[0 vmax],[1 1],'Color','r','Marker','.','LineStyle','-')

        hold on
    end
    hold off
    figure
    plot(BestvalR_fit)

    
elseif size(I,3) == 3 %RGB image
    gBestR = sort(xBestR);
    gBestG = sort(xBestG);
    gBestB = sort(xBestB);
    Iout = imageRGB(I,gBestR,gBestG,gBestB);
    intensity = [gBestR(1:n-1);gBestG(1:n-1);gBestB(1:n-1)]
    STDR =  std(BestvalR)
    STDG =  std(BestvalG)
    STDB =  std(BestvalB)
    figure
    imshow(Iout)
end


function [fmin,best,fminval] = EMOsc20(nub,E)

 