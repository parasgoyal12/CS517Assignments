function A2b_Paras_2018CSB1111_2020_CS517(fname_inp1, fname_inp2, tpts, fname_out, toshow )
    im1=imread(fname_inp1);
    im2=imread(fname_inp2);
    im1_pts=tpts(:,1:2);
    im2_pts=tpts(:,3:4);
    frames=21;
    frame_rate=10;
    wrap_frac=linspace(0,1,frames);
    im=morph(im1,im2,im1_pts,im2_pts,wrap_frac);
        ff=sprintf('%s.gif',fname_out);
    delay=1/frame_rate;
    for idx=1:frames
        
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
    if (toshow~=0)
        figure;
        len=ceil(sqrt(frames));
        for idx=1:frames
            subplot(len,len,idx);
            imshow(im{idx});
        end
    end
end
function O=morph(im1,im2,im1_pts,im2_pts,ratio)
    m1= size(im1,1);
    n1=size(im1,2);
    m2= size(im2,1);
    n2=size(im2,2);
    
    O = cell(length(ratio), 1);
    im1_pts = [double(im1_pts); 0, 0; size(im1, 2) + 1, 0; 0, size(im1, 1) + 1; size(im1, 2) + 1, size(im1, 1) + 1];
    im2_pts = [double(im2_pts); 0, 0; size(im2, 2) + 1, 0; 0, size(im2, 1) + 1; size(im2, 2) + 1, size(im2, 1) + 1];
    [~,b,~]=unique(im1_pts,"rows");
    [~,c,~]=unique(im2_pts,"rows");
    if size(c,1)<size(b,1)
          im1_pts=im1_pts(c,:);
          im2_pts=im2_pts(c,:);
    else
        im1_pts=im1_pts(b,:);
        im2_pts=im2_pts(b,:);
    end

    for t=1:length(ratio)

        alpha=round((1-ratio(t)).*m1+ratio(t).*m2);
        beta=round((1-ratio(t)).*n1+ratio(t).*n2);
        morphed_pts = im1_pts * (1 - ratio(t)) + im2_pts * ratio(t);
        TR=delaunayTriangulation(morphed_pts(:,1),morphed_pts(:,2));
        tri=TR.ConnectivityList;

        [XX,YY]=meshgrid(1:beta,1:alpha);
        XI = [XX(:), YY(:)];
        T=pointLocation(TR,XI);
        cood_bary = zeros(size(XI, 1), 3);
        cood1 = zeros(size(XI, 1), 2);
        cood2 = zeros(size(XI, 1), 2);
        
        cood1_mat = cell(size(tri, 1), 1);
        cood2_mat = cell(size(tri, 1), 1);
        for i = 1: size(tri, 1)
            cood1_mat{i} = [[im1_pts(tri(i, 1), :), 1]' [im1_pts(tri(i, 2), :), 1]' [im1_pts(tri(i, 3), :), 1]'];
            cood2_mat{i} = [[im2_pts(tri(i, 1), :), 1]' [im2_pts(tri(i, 2), :), 1]' [im2_pts(tri(i, 3), :), 1]'];
        end

        
        for i = 1: size(XI, 1)
            cood_bary(i,:)=cartesianToBarycentric(TR,T(i),XI(i,:))';
            cood_homo1 = (cood1_mat{T(i)} * cood_bary(i, :)');
            cood1(i, :) = cood_homo1(1:2)'/cood_homo1(3);
            cood_homo2 = cood2_mat{T(i)} * cood_bary(i, :)';
            cood2(i, :) = (cood_homo2(1:2))'/cood_homo2(3);
        end

        O1=interp(reshape(cood1(:, 1), alpha, beta), reshape(cood1(:, 2), alpha, beta), im1);
        O2 = interp( reshape(cood2(:, 1), alpha, beta), reshape(cood2(:, 2), alpha, beta), im2 );
        O{t} = uint8(O1 * (1 - ratio(t)) + O2 * ratio(t));

        [W1, H1] = meshgrid(linspace(1, beta, size(im2, 2)), linspace(1, alpha, size(im2, 1)));
        O{t} = interp(W1, H1, O{t});
    end
end

function O  = interp( x, y, im )
    a = zeros(size(x, 1), size(x, 2), 3);
    b = zeros(size(x, 1), size(x, 2), 3);
    c = zeros(size(x, 1), size(x, 2), 3);
    d = zeros(size(x, 1), size(x, 2), 3);
    [m,n,~]=size(im);
    frx = floor(x);
    cx = ceil(x);
    frx(frx < 1) = 1;
    frx(frx > n) = n;
    cx(cx < 1) = 1;
    cx(cx > n) = n;    
    x1 = x - frx;
    x1(x1 <= 0) = 1;
    x2 = cx - x;
    x2(x2 < 0) = 1;
    
    fry = floor(y);
    cy = ceil(y);
    fry(fry < 1) = 1;
    fry(fry > m) = m;
    cy(cy < 1) = 1;
    cy(cy >m) = m;    
    y1 = y - fry;
    y1(y1 <= 0) = 1;
    y2 = cy - y;
    y2(y2 < 0) = 1;
    
    for i = 1: size(x, 1)
        for j = 1: size(x, 2)
            a(i, j, :) = im(fry(i, j), frx(i, j), :);
            b(i, j, :) = im(fry(i, j), cx(i, j), :);
            c(i, j, :) = im(cy(i, j), frx(i, j), :);
            d(i, j, :) = im(cy(i, j), cx(i, j), :);
        end
    end


    O= (b.*x1 + a.*x2).*y2 + (d.*x1 + c.* x2).*y1;
    O= uint8(O);

end
