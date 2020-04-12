A=rgb2gray(imread('peppers.png'));

angle=30;
[m,n]=size(imrotate(A,30));
O=zeros([m n]);
tx=ceil(size(A,1)/2)-0.5;
ty=ceil(size(A,2)/2)-0.5;
%matrix=[1 0 ceil(m/2)-0.5;0 1 ceil(n/2)-0.5;0 0 1]*[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1]*[1 0 -tx;0 1 -ty;0 0 1];
matrix=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1];
alpha=1;


for lambda=0:0.05:1
    %matrix=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1]*[1 0 -tx;0 1 -ty;0 0 1];
    matrix=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1];
    matrix=matrix.*lambda;
    matrix(1,1)=matrix(1,1)-lambda+1;
    matrix(2,2)=matrix(2,2)-lambda+1;
    matrix(3,3)=1;
    alpha1=matrix*[1 1 1];
    alpha2=matrix*[]
    matrix=[1 0 ceil(m/2)-0.5;0 1 ceil(n/2)-0.5;0 0 1]*matrix;
    for i=1:m
        for j=1:n
            B=matrix\[i-0.5;j-0.5;1];
            x=B(1);
            y=B(2);
            if(x>0 && x<=size(A,1) && y>0 && y<=size(A,2))
                x=ceil(x);
                y=ceil(y);
                O(i,j)=A(x,y);
            end
        end
    end
    
    O=uint8(O);
    imshow(O);
    im{alpha}=O;
    alpha=alpha+1;
end

ff='test.gif';
for idx=1:numel(0:0.05:1)
    A=im{idx};

    if idx==1
        imwrite(A,gray,ff,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,gray,ff,'gif','WriteMode','append','DelayTime',0.1);
    end
end