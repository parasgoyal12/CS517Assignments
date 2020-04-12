function generate_gif(I,output,frame_rate)
    ff=sprintf('%s.gif',output);
    delay=1/frame_rate;
    im=I;
    for idx=1:max(size(I,1),size(I,2))
        
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
    close;
end    