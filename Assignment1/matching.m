function O=matching(I,ref)
    O=zeros(size(I));
    M=zeros(256,1);
    cdf1=cdf(I);
    cdf2=cdf(ref);
    for i=1:256
        [alpha,idx]=min(abs(cdf1(i)-cdf2));
        M(i)=idx-1;
    end
    for i=0:255
        O(I==i)=M(i+1);
    end
    O=uint8(O);
end