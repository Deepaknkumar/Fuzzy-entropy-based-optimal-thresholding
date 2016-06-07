function mu = mucreator(values)

    a1 = values(1);
    b1 = values(2);
    c1 = values(3);
    a2 = values(4);
    b2 = values(5);
    c2 = values(6);
    
    mu = zeros(1,256);
    for k=1:256
        if k <= a1
            mu(k) = 0;
        elseif k > a1 && k <= b1
            mu(k) = ((k-a1)^2)/((c1-a1)*(b1-a1));
        elseif k > b1 && k <= c1
            mu(k) = 1 - (((k-c1)^2)/((c1-a1)*(c1-b1)));
        elseif k > c1 && k <= a2
            mu(k) = 1;
        elseif k > a2 && k <= b2
            mu(k) = 1 - ((k-a2)^2)/((c2-a2)*(b2-a2));
        elseif k > b2 && k <= c2
            mu(k) = (((k-c2)^2)/((c2-a2)*(c2-b2)));
        elseif k > c2
            mu(k) = 0;
        end
    end