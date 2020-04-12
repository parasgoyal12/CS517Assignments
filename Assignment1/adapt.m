function O=adapt(I,windowsize)
O=zeros(size(I));

n=floor(windowsize/2);

I=padarray(I,[n n],'symmetric','both');

for i=1+n:size(I,1)-n
    for j=1+n:size(I,2)-n
       
            fromrow=i-n;
            torow=i+n;
            fromcol=j-n;
            tocol=j+n;
            
            neighbour=I(fromrow:torow,fromcol:tocol);
            cdf1=cdf(neighbour);
            M=round(cdf1.*255);
            
            O(i-n,j-n)=M(I(i,j)+1);
        
    end
end
O=uint8(O);
end