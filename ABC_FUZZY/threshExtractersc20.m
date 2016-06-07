function threshvalue = threshExtractersc20(uValues)
    a1 = uValues(1);
    b1 = uValues(2);
    c1 = uValues(3);
    
    if b1 >= (a1+c1)/2
        threshvalue = a1 + sqrt((c1-a1)*(b1-a1)/2);
    elseif b1 < (a1+c1)/2
        threshvalue = c1 - sqrt((c1-a1)*(c1-b1)/2);
    end