function [O,O2]=jpeg(fname_inp,toshow,toplot)
    QX=[ 16 11 10 16 24 40 51 61;
         12 12 14 19 26 58 60 55;
         14 13 16 24 40 57 69 56;
         14 17 22 29 51 87 80 62; 
         18 22 37 56 68 109 103 77;
         24 35 55 64 81 104 113 92;
         49 64 78 87 103 121 120 101;
         72 92 95 98 112 100 103 99];
    quality=50;
    block_size=4;
    if quality > 50
        QX = round(QX.*(ones(8)*((100-quality)/50)));
        
    elseif quality < 50
        QX = round(QX.*(ones(8)*(50/quality)));
        
    end
    QX=QX(1:block_size,1:block_size);
    n=size(QX,1);
    input_image=rgb2gray(imread(fname_inp));
    dct_quantized=encode_jpeg(fname_inp,QX,n);
    O=cell(n,1);
    rms=zeros(n,1);
    for i=1:n
        O{i}=decode_jpeg(dct_quantized,QX,n,i);
        rms(i)=compute_imgs_rmse(input_image,O{i});
    end
    O2=cell(n,1);
    fig=figure;
    for i=1:n
        imshow(O{i});
        title(sprintf("RMS=%0.2f",rms(i)));
        drawnow
        frame=getframe(fig);
        O2{i}=frame2im(frame);
    end
    close;
    if toshow==1
        a=ceil(sqrt(n));
        figure('Name','Images');
        for i=1:n
            subplot(a,a,i);
            imshow(O{i});
            title(sprintf('Feature Length=%d,RMSE=%0.2f',i,rms(i)));
        end
    end
    if toplot==1
        figure('Name','Plot of rms vs feature length');
        plot(1:n,rms,'r-');
        xlabel('Feature Lenght');
        ylabel('RMSE');
    end
end
function O=encode_jpeg(fname_inp,QX,n)
    I=imread(fname_inp);
    I=double(rgb2gray(I));
    I=padarray(I,[ceil(size(I,1)/n)*n-size(I,1),ceil(size(I,2)/n)*n-size(I,2)],'post');
    [row,col]=size(I);
    I=I-128;
    QX=double(QX);
    dct_quantized=zeros(row,col);
    for i1=1:n:row
     for i2=1:n:col
         zBLOCK=I(i1:i1+n-1,i2:i2+n-1);
         win1=dct2(zBLOCK);
         win2=round(win1./QX);
         dct_quantized(i1:i1+n-1,i2:i2+n-1)=win2;
     end
    end
    O=dct_quantized;
end
function O=decode_jpeg(dct_quantized,QX,n,features)
    feature=zeros(n);
    [row,col]=size(dct_quantized);
    feature(1:features,1:features)=ones(features);
    dct_restored=zeros(row,col);
    for i1=1:n:row
     for i2=1:n:col
         win2=dct_quantized(i1:i1+n-1,i2:i2+n-1);
         win2=win2.*feature;
         win3=win2.*QX;
         dct_restored(i1:i1+n-1,i2:i2+n-1)=idct2(win3)+128;
     end
    end
    O=uint8(dct_restored);
end
function xx = compute_imgs_rmse(A,B)
    t1=min(size(A,1), size(B,1));
    t2=min(size(A,2), size(B,2));
    err=double(A(1:t1,1:t2)-B(1:t1,1:t2));
    xx=sqrt( sum(sum(double(err.*err))) / (numel(err)) );
end

 