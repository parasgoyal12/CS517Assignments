img=imread('mozart.jpg');
angle=30;
img_ref=imref2d(size(img));
sz=size(imrotate(img,angle));
 var=0:0.025:1;
 mat1=[-2,1,0;-1,2,0;300,344,1];
 t=affine2d(mat1);
 sz=size(imwarp(img,t));
i=1;
fig=figure;
background=zeros(2000);
im=cell(numel(var));
for lambda=var
    %mat=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1];
    
    mat=mat1.*lambda;
    mat(1,1)=mat(1,1)-lambda+1;
    mat(2,2)=mat(2,2)-lambda+1;
    mat(3,3)=1;

    tform=affine2d(mat);

    [outp,ref]=imwarp(img,img_ref,tform,'cubic');
%     %outp=padarray(outp,[ceil(sz(1)*0.5-size(outp,1)*0.5),ceil(sz(2)*0.5-size(outp,2)*0.5)],'both');
%     %outp=imresize(outp,[sz(1) sz(2)]);
%     a=size(outp,1);
%     b=size(outp,2);
%     a=sz(1)*0.5-a*0.5;
%     b=sz(2)*0.5-b*0.5;
% 
%     outp=padarray(outp,[ceil(a),ceil(b)],'pre');
%     outp=padarray(outp,[floor(a),floor(b)],'post');
    imshowpair(outp,ref,background,imref2d(size(background)));
    drawnow
    frame=getframe(fig);
    im{i}=frame2im(frame);
    i=i+1;
end

% ff='test.gif';
% delay=0.1;
% for idx=1:numel(var)
%      [A,map]=rgb2ind(im{idx},256);
% 
%     if idx==1
%         imwrite(A,map,ff,'gif','LoopCount',Inf,'DelayTime',delay);
%     else
%         imwrite(A,map,ff,'gif','WriteMode','append','DelayTime',delay);
%     end
% end
% for idx=numel(var):-1:1
%     [A,map]=rgb2ind(im{idx},256);
%     imwrite(A,map,ff,'gif','WriteMode','append','DelayTime',delay);
% end
% close;