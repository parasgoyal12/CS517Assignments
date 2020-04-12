function A2a_Paras_2018CSB1111_2020_CS517(image_name,dummy,angle,dummy2,afnTM,output,toshow)    
    img=(imread(image_name));
    frame_rate=10;
    frames=40;
    lambda=1/frames;
    var=0:lambda:1;
    if(dummy==1),mat1=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1];
        else,mat1=afnTM;
    end
%     sz=size(imrotate(img,angle));
    
%     mat1=[0,-1,0;1,0,0;0,0,1];
    t=affine2d(mat1);
    szt=size(imwarp(img,t));
    sz(1)=max([szt(1),size(img,1)]);
    sz(2)=max([szt(2),size(img,2)]);
    i=1;
    fig=figure;
    im=cell([numel(var),1]);
    for lambda=var
        %mat=[cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0;0 0 1];
        
        mat=mat1.*lambda;
        mat(1,1)=mat(1,1)-lambda+1;
        mat(2,2)=mat(2,2)-lambda+1;
        mat(3,3)=1;
        if(det(mat)==0),continue;end
        tform=affine2d(mat);

        outp=imwarp(img,tform,'cubic');
        %outp=padarray(outp,[ceil(sz(1)*0.5-size(outp,1)*0.5),ceil(sz(2)*0.5-size(outp,2)*0.5)],'both');
        %outp=imresize(outp,[sz(1) sz(2)]);
        a=size(outp,1);
        b=size(outp,2);
        a=abs(sz(1)*0.5-a*0.5);
        b=abs(sz(2)*0.5-b*0.5);
        if(a<0),a=0;end
        if(b<0),b=0;end
        outp=padarray(outp,[ceil(a),ceil(b)],'pre');
        outp=padarray(outp,[floor(a),floor(b)],'post');
        imshow(outp);
        drawnow
        frame=getframe(fig);
        im{i}=frame2im(frame);
        i=i+1;
    end
    
    ff=sprintf('%s.gif',output);
    delay=1/frame_rate;
    for idx=1:numel(var)
        
        if(size(im{idx},3)==1)
             [A,map]=gray2ind(im{idx},256);
        else
            [A,map]=rgb2ind(im{idx},256);
         end
         if(isempty(A)),continue;end 
        if idx==1
            imwrite(A,map,ff,'gif','LoopCount',Inf,'DelayTime',delay);
        else
            imwrite(A,map,ff,'gif','WriteMode','append','DelayTime',delay);
        end
    end
    for idx=numel(var):-1:1
%         [A,map]=rgb2ind(im{idx},256);
         if(size(im{idx},3)==1)
             [A,map]=gray2ind(im{idx},256);
         else
            [A,map]=rgb2ind(im{idx},256);
         end
         if(isempty(A)),continue;end
        imwrite(A,map,ff,'gif','WriteMode','append','DelayTime',delay);
    end
    close;

    if(toshow~=0)
        figure;
        len=ceil(sqrt(frames));
        for idx=1:numel(var)
            subplot(len,len,idx);
            imshow(im{idx});
        end
    end
    
end