function O=histogrameq(I)
    O=zeros(size(I));

    freq=zeros([1 256]);
    
    for i=0:255
       load=numel(I(I==i));
       for j=i+1:256
            freq(j)=freq(j)+load;
       end
    end
    
    total=numel(I);
    freq=freq./total;

    freq=freq.*255;
    freq=round(freq);
    
    for i=0:255
        O(I==i)=freq(i+1);
    end
end