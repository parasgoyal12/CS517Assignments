%Name-Paras Goyal
%References->Stack Overflow(stackoverflow.com)
%2.http://fourier.eng.hmc.edu/e161/lectures/contrast_transform/node3.html
%3. MatLab Docs for various functions

function val=A1_Paras_2018CSB1111_2019_CS517(qID,fname_inp1,fname_inp2,frame_out,prmts,toshow)
    
    val=[];
    if(qID==1), val=nearest(fname_inp1,prmts(1),prmts(2),toshow,frame_out);end
    if(qID==2), val=bilinear(fname_inp1,prmts(1),prmts(2),toshow,frame_out);end
    if(qID==3), val=rotation(fname_inp1,prmts,toshow,frame_out);end
    if(qID==4), val=bitplane(fname_inp1,prmts,toshow,frame_out);end
    if(qID==5), val=regis(fname_inp1,fname_inp2,prmts,toshow,frame_out);end
    if(qID==6), val=histogrameq(fname_inp1,toshow,frame_out);end
    if(qID==7), val=matching(fname_inp1,fname_inp2,toshow,frame_out);end
    if(qID==8), val=adapt(fname_inp1,toshow,frame_out);end
end
function xx=nearest(img,h,w,toshow,frame_out)
    img=imread(img);
    O=zeros(h,w);
    [hold,wold]=size(img);
    rh = hold/h;
    rw = wold/w;
    for i=1:w
        x=round((i-0.5)*rw + 0.5);
        for j=1:h
            y=round((j-0.5)*rh+0.5);
            if x==0
                x=1;
            elseif x== wold+1
                x=hold;
            end
            if y==0
                y=1;
            elseif y==hold+1
                y=hold;
            end
            O(j,i)=img(y,x);
        end
    end
    O=uint8(O);
    imwrite(O,sprintf('%s.tif',frame_out));
    B=imresize(img,[h w],'nearest');
    xx=rmse(O,B);
    if(toshow==1)
        figure;subplot(2,2,1);
    imshow(img);
    title("Input Image");
    subplot(2,2,2);
    imshow(B);
    title(sprintf("InbuiltNearestN.Mean=%0.2f .RMSE=%0.2f ",mean(B(:)),xx));
    subplot(2,2,3);
    imshow(O);
    title(sprintf("Custom NearestN. Mean=%0.2f RMSE=%0.2f",mean(O(:)),xx));
    subplot(2,2,4);
    B=imresize(img,[h w],'bicubic');
    imshow(B);
    title(sprintf("Bicubic. Mean=%0.2f RMSE=%0.2f",mean(B(:)),rmse(B,O)));
    end
end
function xx=bilinear(A,h,w,toshow,frame_out)
    img=imread(A);
    A=imread(A);
    [hold,wold]=size(A);
    O=zeros(h,w);
    rh=hold/h;
    rw=wold/w;
    A=double(A);
    for i=1:h
        x=(i-0.5)*rh+0.5;
        for j=1:w
            y=(j-0.5)*rw+0.5; 
            x1=floor(x);
            x2=ceil(x);
            y1=floor(y);
            y2=ceil(y);
            if x1==0
                x1=1;
            end
            if x2==hold+1
               x2=hold;
            end
            if y1==0
                y1=1;
            end
            if y2==wold+1
                y2=wold;
            end            
            alpha1=abs(x-x1);
            alpha2=abs(y-y1);            
            if x1==x2
                alpha1=0.5;
            end
            if y1==y2
                alpha2=0.5;
            end                
            h1=alpha2*A(x1,y2)+(1-alpha2)*A(x1,y1);
            h2=alpha2*A(x2,y2)+(1-alpha2)*A(x2,y1);
            O(i,j)=(1-alpha1)*h1+(alpha1)*h2;
        end
    end
    
    O=uint8(O);
    imwrite(O,sprintf('%s.tif',frame_out));
    B=imresize(img,[h w],'bilinear');
    xx=rmse(O,B);
    if(toshow==1)
        figure;subplot(2,2,1);
    imshow(img);
    title("Input Image");
    subplot(2,2,2);
    imshow(B);
    title(sprintf("InbuiltBiliner.Mean=%0.2f .RMSE=%0.2f ",mean(B(:)),xx));
    subplot(2,2,3);
    imshow(O);
    title(sprintf("CustomBilinear. Mean=%0.2f RMSE=%0.2f",mean(O(:)),xx));
    subplot(2,2,4);
    B=imresize(img,[h w],'bicubic');
    imshow(B);
    title(sprintf("Bicubic. Mean=%0.2f RMSE=%0.2f",mean(B(:)),rmse(B,O)));
    end    
end
function xr=rotation(I,angle,toshow,frame_out)
    I=imread(I);
    p=0.5;
    [m,n]=size(I);
     a=pi*angle/180;
    %a=angle;
    om=ceil(m*abs(cos(a))+n*abs(sin(a)));
    on=ceil(m*abs(sin(a))+n*abs(cos(a)));
    if(angle==90)
        om=n;on=m;
    end
    if(angle==270)
        om=n;on=m;
    end
    if(angle==180)
        om=m;on=n;
    end
    O=uint8(zeros([om,on]));
    xo=ceil(m/2)-p;
    yo=ceil(n/2)-p;
    midx=ceil(om/2)-p;
    midy=ceil(on/2)-p;
    for i=1:om
        for j=1:on
    %         x=round((i-midx)*cos(a)+(j-midy)*sin(a))+xo;
    %         y=round(-(i-midx)*sin(a)+(j-midy)*cos(a))+yo;
            x=(i-midx-p)*cos(a)+(j-midy-p)*sin(a)+xo+p;
            y=-(i-midx-p)*sin(a)+(j-midy-p)*cos(a)+yo+p;
            x1=floor(x);
            x2=ceil(x);
            y1=floor(y);
            y2=ceil(y);
            alpha1=abs(x-x1);
            alpha2=abs(y-y1);
            if(x1==0),x1=1;x2=1;end
            if(x2==(m+1)),x2=m;x1=m;end
            if(y1==0),y1=1;y2=1;end
            if(y2==(n+1)),y2=n;y1=n;end
            if(x1==x2),alpha1=0.5;end
            if(y2==y1),alpha2=0.5;end
    %             if(x>0&&x<=m&&y>0&&y<=n)
            if(x1>0 && y1>0 && x2<=m && y2<=n)
                t1=alpha2*I(x1,y2)+(1-alpha2)*I(x1,y1);
                t2=alpha2*I(x2,y2)+(1-alpha2)*I(x2,y1);
                O(i,j)=alpha1*(t2)+(1-alpha1)*t1;
    %                O(i,j)=I(x,y);
            end
        end
    end
    O=uint8(O);
    B=imrotate(I,angle,'bilinear');
    error=rmse(B,O);
    xr=[om,on,error];
    imwrite(O,sprintf('%s.tif',frame_out));
    if(toshow==1)
        figure;
        subplot(1,3,1);
        imshow(I);
        title("Input Image");
        subplot(1,3,2);
        imshow(O);
        title(sprintf("CustomRotateBilinear. Mean=%0.2f RMSE=%0.2f size=[%d %d]",mean(O(:)),error,om,on));
        subplot(1,3,3);
        imshow(B);
        title(sprintf("ImRotateBilinear. Mean=%0.2f Size=[%d %d]",mean(B(:)),size(B,1),size(B,2)));
    end
end
function rms=bitplane(I,vales,toshow,frame_out)
    I=imread(I);
    I=double(I);
    O=zeros([size(I,1),size(I,2),numel(vales)]);
    
    for i=1:numel(vales)
       O(:,:,i)=bitand(I,vales(i)); 
    end
    
    O=uint8(O);
    rms=zeros(1,numel(vales));
    
    I=uint8(I);
    if(toshow)
    figure;
    subplot(2,2,1);
    imshow(I);
    title("Input Image");
    end
    for i=1:numel(vales)
        rms(i)=rmse(I,O(:,:,i));
        imwrite(O(:,:,i),sprintf("%s_%d.tif",frame_out,i));
        if(toshow)
            subplot(ceil(numel(vales)/2),ceil(numel(vales)/2),i+1);
            B=O(:,:,i);
            imshow(B);     
            title(sprintf("Img corrs. =%d Mean=%.2f RMSE=%0.2f",vales(i),mean(B(:)),rms(i)));
        end
    end
end
function xr=regis(I1,I2,M,toshow,frame_out)
    outputm=688;
    outputn=688;
    I=imread(I2);
    I=double(I);
%    Uncomment the code below to enable automatic x and y calculation
%    Cprime=genMatrix(R+0.5,P+0.5);
%    [m,n]=size(I);
%    xprime=[m-0.5 n-0.5 (m-0.5)*(n-0.5) 1]*Cprime(1:4)+0.5;
%    yprime=[m-0.5 n-0.5 (m-0.5)*(n-0.5) 1]*Cprime(5:8)+0.5;
%    outputm=round(xprime);
%    outputn=round(yprime);
    P=M(:,3:4)-0.5;
    R=M(:,1:2)-0.5;
    mat=zeros(8,8);
    mat(1,1:2)=R(1,:);
    mat(1,4)=1;
    mat(1,3)=R(1,1)*R(1,2);
    mat(2,5:8)=mat(1,1:4);
    mat(3,1:2)=R(2,:);
    mat(3,4)=1;
    mat(3,3)=R(2,1)*R(2,2);
    mat(4,5:8)=mat(3,1:4);
    mat(5,1:2)=R(3,:);
    mat(5,4)=1;
    mat(5,3)=R(3,1)*R(3,2);
    mat(6,5:8)=mat(5,1:4);
    mat(7,1:2)=R(4,:);
    mat(7,4)=1;
    mat(7,3)=R(4,1)*R(4,2);
    mat(8,5:8)=mat(7,1:4);
    A=[P(1,1) P(1,2) P(2,1) P(2,2) P(3,1) P(3,2) P(4,1) P(4,2)];
    A=A';
    C=mat\A;
    O=zeros(outputm,outputn);
    c1=C(1:4,:);
    c2=C(5:8,:);
    [m,n]=size(I);
    for i=1:outputm
        for j=1:outputn
            x=[i-0.5 j-0.5 (i-0.5)*(j-0.5) 1]*c1+0.5;
            y=[i-0.5 j-0.5 (i-0.5)*(j-0.5) 1]*c2+0.5;
            x1=floor(x);
            x2=ceil(x);
            y1=floor(y);
            y2=ceil(y);            
            if(x1<=0)
                x1=1;
            end
            if(x2<=0)
                x2=1;
            end
            if(x1>m)
                x1=m;
            end
            if(x2>m)
                x2=m;
            end            
            alpha1=abs(x-x1);
            alpha2=abs(y-y1);
            if(y1<=0)
                y1=1;
            end
            if(y2<=0)
                y2=1;
            end
            if(y1>n)
                y1=n;
            end
            if(y2>n)
                 y2=n;
            end            
             if(x1==x2)
                 alpha1=0.5;
             end
             if(y1==y2)
                 alpha2=0.5;
             end
             t1=alpha2*I(x1,y2)+(1-alpha2)*I(x1,y1);
             t2=alpha2*I(x2,y2)+(1-alpha2)*I(x2,y1);            
             O(i,j)=alpha1*t2+(1-alpha1)*t1;           
        end
    end
    O=uint8(O);
    imwrite(O,sprintf("%s.tif",frame_out));
    I1=imread(I1);
    rms=rmse(I1,O);
    xr=[outputm,outputn,rms];
    if(toshow)
        figure;
        subplot(2,2,1);
        imshow(I1);
        imshow(I1);
        title("Input Image 1");
        subplot(2,2,2);
        imshow(I);
        title("Input Image 2");
        subplot(2,2,3);
        imshow(O);
        title(sprintf("Generated RMS=%.2f Mean=%.2f Size=[%d %d]",rms,mean(O(:)),outputm,outputn));
        subplot(2,2,4);
        O=double(O)./255;
        I1=double(I1)./255;
        err=sqrt((O-I1).^2);
        imshow(err);
        title("Error Image");
    end
end
function rms=histogrameq(I,toshow,frame_out)
    I=imread(I);
    I=double(I);
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
    O=uint8(O);
    imwrite(O,sprintf('%s.tif',frame_out));
    B=histeq(I);
    rms=rmse(B,O);
    if(toshow)
        figure;
        subplot(2,3,1);
        imshow(I);
        title(sprintf("Input Image Mean=%.2f",mean(I(:))));
        subplot(2,3,2);
        imhist(I);
        title(sprintf("Histogram of Input"));
        subplot(2,3,3);
        imshow(O);
        title(sprintf("Histogram equalized Image Mean=%.2f Rms=%.2f",mean(O(:)),rms));
        subplot(2,3,4);
        imhist(O);
        title("Histogram of My Output");
        subplot(2,3,5);
        imshow(B);
        title(sprintf("System equalized Image Mean=%.2f Rms=%.2f",mean(O(:)),rms));
        subplot(2,3,6);
        imhist(B);
        title("Histogram of System Output");    
    end
end
function rms=matching(I,ref,toshow,frame_out)
    I=imread(I);
    I=double(I);
    ref=imread(ref);
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
    imwrite(O,sprintf("%s.tif",frame_out));
    B=imhistmatch(I,ref);
    rms=rmse(B,O);
    if(toshow)
        figure;
        subplot(2,4,1);
        imshow(I);
        title(sprintf("Input Image Mean=%.2f",mean(I(:))));
        subplot(2,4,2);
        imhist(I);
        title(sprintf("Histogram of Input"));
        subplot(2,4,3);
        imshow(ref);
        title(sprintf("Reference Image Mean=%.2f",mean(ref(:))));
        subplot(2,4,4);
        imhist(ref);
        title("Histogram of Reference");
        subplot(2,4,5);
        imshow(O);
        title(sprintf("Histogram Matched Image Mean=%.2f Rms=%.2f",mean(O(:)),rms));
        subplot(2,4,6);
        imhist(O);
        title("Histogram of My Output");
        subplot(2,4,7);
        imshow(B);
        title(sprintf("System Matched Image Mean=%.2f Rms=%.2f",mean(B(:)),rms));
        subplot(2,4,8);
        imhist(B);
        title("Histogram of System Matched");    
    end
end
function rms=adapt(I,toshow,frame_out)
    I=double(imread(I));
    I1=I;
    windowsize=8;
    O=zeros(size(I));
    n=floor(windowsize/2);
    I=padarray(I,[n n],'symmetric','both');
    for i=1+n:size(I,1)-n
        disp(i);
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
    
    B=adapthisteq(uint8(I1),'NumTiles',[windowsize windowsize],'ClipLimit',0);
    rms=rmse(B,O);
    imwrite(O,sprintf("%s.tif",frame_out));
    if(toshow)
            figure;
            subplot(2,3,1);
            imshow(uint8(I1));
            title(sprintf("Input Image Mean=%.2f",mean(I1(:))));
            subplot(2,3,2);
            imhist(uint8(I1));
            title(sprintf("Histogram of Input"));
            subplot(2,3,3);
            imshow(O);
            title(sprintf("Histogram equalized Image Mean=%.2f Rms=%.2f",mean(O(:)),rms));
            subplot(2,3,4);
            imhist(O);
            title("Histogram of My Output");
            subplot(2,3,5);
            imshow(B);
            title(sprintf("System equalized Image Mean=%.2f Rms=%.2f",mean(B(:)),rms));
            subplot(2,3,6);
            imhist(B);
            title("Histogram of System Output");    
    end
end

function xx=rmse(A,B)
    t1=min(size(A,1),size(B,1));
    t2=min(size(A,2),size(B,2));
    A=double(A);
    B=double(B);
    err=double(A(1:t1,1:t2)-B(1:t1,1:t2));
    xx=sqrt(sum(sum(double(err.*err)))/(numel(err)));
end
function freq=cdf(I)
    freq=zeros([1 256]);    
    for i=0:255
       load=numel(I(I==i));
       for j=i+1:256
            freq(j)=freq(j)+load;
       end
    end    
    total=numel(I);
    freq=freq/total; 
end
function M=genMatrix(P,R)
    P=P-0.5;
    R=R-0.5;
    mat=zeros(8,8);
    mat(1,1:2)=R(1,:);
    mat(1,4)=1;
    mat(1,3)=R(1,1)*R(1,2);
    mat(2,5:8)=mat(1,1:4);
    mat(3,1:2)=R(2,:);
    mat(3,4)=1;
    mat(3,3)=R(2,1)*R(2,2);
    mat(4,5:8)=mat(3,1:4);
    mat(5,1:2)=R(3,:);
    mat(5,4)=1;
    mat(5,3)=R(3,1)*R(3,2);
    mat(6,5:8)=mat(5,1:4);
    mat(7,1:2)=R(4,:);
    mat(7,4)=1;
    mat(7,3)=R(4,1)*R(4,2);
    mat(8,5:8)=mat(7,1:4);
    A=[P(1,1) P(1,2) P(2,1) P(2,2) P(3,1) P(3,2) P(4,1) P(4,2)];
    A=A';
    M=mat\A;
end