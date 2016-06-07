function z = fitnessfuncsc20otsu(u,E,muT)
            threshlevel = length(u);
            sz = 0:1:255;
            %u = round(u);
            %%%---------------- For First threshold level--------------%%%
            omega = sum(E(1:u(1)));         %%% Calculation of omega for first level
            mu = sum(sz(1:u(1)).*(E(1:u(1))))/omega;
            F(1) = omega*((mu-muT).^2);
                for ij = 2:threshlevel          %%% Calculation of Entropy of each level
                   omega = sum(E(u(ij-1)+1:u(ij)));         %%% Calculation of omega for first level
                   mu = sum(sz(u(ij-1)+1:u(ij)).*(E(u(ij-1)+1:u(ij))))/omega;
                   F(ij) = omega*((mu-muT).^2);
                end
            omega = sum(E(u(ij)+1:256));         %%% Calculation of omega for first level
            mu = sum(sz(u(ij)+1:256).*(E(u(ij)+1:256)))/omega;
            F(ij+1) = omega*((mu-muT).^2);
            z = -1*sum(F);         %%% Sum of entropies of all levels