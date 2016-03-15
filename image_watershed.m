function [IM2d_blob]=image_watershed(i)
% [ image watershed transform ]
subsegmethod=2;

load result5.mat;
BW_obj=1;
thed=0.9;
thed2=10;
B=blob{i};
IM2d=B.Image;BW=IM2d;
IM3d=zeros(size(IM2d));
for k=1:B.blobArea
    IM3d(B.PixelList(k,2)-B.BoundingBox(1,2)+0.5,B.PixelList(k,1)-B.BoundingBox(1,1)+0.5)=B.PixelValues(k,1);
end

% BW: binary image
% BW_obj=1(object=1), =0(object=0)
% D: distance transform of the complement 
% L: labeling of distance transform of the complement 
    figure;imagesc(IM2d);axis equal;colorbar;
    figure;imagesc(IM3d);axis equal;colorbar;
if subsegmethod==1; % 1: watershed method
    if BW_obj==1
        D=bwdist(~BW);
        D=-D;
        D(~BW)=-Inf; 
        L=watershed(D);
        L=watershed_fix(L,IM3d);
        rgb=label2rgb(L);
        figure;imshow(rgb,'InitialMagnification','fit');
    elseif BW_obj==0
        D=bwdist(BW);
        D=-D;D(~BW)=-Inf;  
        L=watershed(D);rgb=label2rgb(L);
        figure;imshow(rgb,'InitialMagnification','fit');
    end
    
elseif subsegmethod==2; % 2: adhesion cut
    [BW2,BWedge]=edge_add(IM2d);
    k2=conv_point(BWedge,thed,IM2d);
    IM2d_blob=link_neck(IM2d,k2,BWedge,thed2);
    figure;imagesc(IM2d_blob);axis equal;colorbar;
    
    % draw neck points
    IM2d_add=uint8(IM2d);
    for ki=1:size(k2,2)
        IM2d_add(BWedge(k2(1,ki),2),BWedge(k2(1,ki),1))=2;
    end
    figure;imagesc(IM2d_add);axis equal;colorbar;
end
    
end

function Lfix=watershed_fix(L,IM3d)
r=1;
threh=20;

[n,m]=size(L);Lfix=L;
for i=1:n
    for j=1:m
        if L(i,j)==0           
            Lfix(i,j)=water_fix(IM3d,L,i,j,r,threh,n,m);
        end        
    end
end

end
function output=water_fix(IM3d,L,i,j,r,threh,n,m)
    region_t2=IM3d(max(i-r,1):min(i+r,n),max(j-r,1):min(j+r,m));
    region_t1=L(max(i-r,1):min(i+r,n),max(j-r,1):min(j+r,m));
    region_t=[region_t1(:)';region_t2(:)'];region_t=(sortrows(region_t'))';
    tabL=tabulate(region_t(1,:));
    tabl_out=[];ii=0;
    for k=1:size(tabL,1)
        if tabL(k,2)~=0
            tabl_out=[tabl_out;[tabL(k,1) mean(mean(region_t(2,ii+1:ii+tabL(k,2))))]];
            ii=ii+tabL(k,2);
        end
    end
    if tabl_out(1,1)==0
        tabl_out=tabl_out(2:end,:);
    end
    if std(tabl_out(:,2))<threh
        [a,b]=min(abs(tabl_out(:,2)-IM3d(i,j)));
        output=tabl_out(b,1);
    else
        output=0;
    end
end
function [BW2,BWedge]=edge_add(BW1)
BW2=edge(BW1);
for i=1:size(BW1,1)
    for j=1:size(BW1,2)     
        if i==1||i==size(BW1,1)||j==1||j==size(BW1,2)
            if BW1(i,j)==1
                BW2(i,j)=1;
            end
        end     
    end
end
[xy(:,2),xy(:,1)]=find(BW2);n=size(xy,1);
xy=[xy zeros(n,1)];
BWedge=zeros(n,2);BWedge(1,:)=xy(1,1:2);xy(1,3)=1;
for k=2:n
    xy=xy(xy(:,3)~=1,:);
    A=distance(xy(:,1:2),BWedge(k-1,:));
    [a,b]=min(A);xy(b,3)=1;
    BWedge(k,:)=xy(b,1:2);
end
BWedge=BWedge(BWedge(:,1)~=0,:);
end
function k2=conv_point2(BWedge,thed)
% Obtain neck points k2 from BWedge xy coordinates with thed distance threshold

    % convex hull of BWedge: k
    k=convhull(BWedge);k2=[];
    if sum(diff(k)./abs(diff(k)))<0 % check order of k equal to BWedge
        k=k(end:-1:1);
    end
    % find neck points on cell boundary: k2
    for i=2:size(k,1) % for each convex region
        xh=linspace(BWedge(k(i-1),1),BWedge(k(i),1),50); % linspace convex hull xy into 50 pieces: (xh,yh) 
        yh=linspace(BWedge(k(i-1),2),BWedge(k(i),2),50);
        xw=BWedge(k(i-1):k(i),1);yw=BWedge(k(i-1):k(i),2);w_list=k(i-1):k(i); % wall xy: (xw,yw)
        conv_d=zeros(1,size(xw,1));
        for j=1:size(xw,1)
            conv_d(1,j)=min(distance([xh;yh]',[xw(j,1) yw(j,1)])); % find distance to hull for each wall point            
        end
        [a,b]=max(conv_d); % max distance a to hull in each convex region, neck point: b
        if a>thed
            k2=[k2 w_list(b)]; % k2: neck points
        end
    end
end
function k2=conv_point(BWedge,thed,IM2d)
% Obtain neck points k2 from BWedge xy coordinates with thed distance threshold
    n=2;
    [im_a,im_b]=size(IM2d);
    ded_vec=angle_line(BWedge,n);
    k=find(ded_vec<thed);k2=[];
    for iik=1:size(k,1)
        if isempty(find(1:n==BWedge(k(iik),1)))==0||isempty(find(1:n==BWedge(k(iik),2)))==0||isempty(find(im_b-n+1:im_b==BWedge(k(iik),1)))==0||isempty(find(im_a-n+1:im_a==BWedge(k(iik),2)))==0
            k(iik)=0;
        end
    end
    k=k(k~=0);
    if sum(diff(k)./abs(diff(k)))<0 % check order of k equal to BWedge
        k=k(end:-1:1);
    end
    % find neck points on cell boundary: k2
    k0=k;
    %
    for j=1:size(k,1)
       if isempty(find(k-k(j)<5&k-k(j)>0))==0||k(j)~=0
           b=find(k-k(j)<5&k-k(j)>0);
           if isempty(b)==0
               b_list=[[k(j) k(b)'];[ded_vec(k(j)) ded_vec(k(b))']];
               [~,a]=min(b_list(2,:));
               k(k==b_list(1,a))=0;
           end
       end
    end
       %
k2=k0';
    function ded_vec=angle_line(BWedge,n)
    % calculate direction vector of xy in window of boundary
    a=size(BWedge,1);vec=zeros(a,2); %vec: mean vector series on wboundary with window size n
    ded_vec=zeros(a,1); % vecd: dot product between i, i+1 mean vector
    for ii=1:a
        if ii-1<n
            vec(ii,1:2)=mean(diff([BWedge(a-n+ii:a,:);BWedge(1:ii+n,:)]));
        elseif ii>a-n
            vec(ii,1:2)=mean(diff([BWedge(ii-n:a,:);BWedge(1:n-a+ii,:)]));
        else
            vec(ii,1:2)=mean(diff(BWedge(ii-n:ii+n,:)));
        end
        % calculate dplot
        if ii>1
            ded_vec(ii,1)=dot(vec(ii-1,1:2)/sqrt(vec(ii-1,1)^2+vec(ii-1,2)^2),vec(ii,1:2)/sqrt(vec(ii,1)^2+vec(ii,2)^2));
        end
    end
    ded_vec(1,1)=dot(vec(1,1:2)/sqrt(vec(1,1)^2+vec(1,1)^2),vec(a,1:2)/sqrt(vec(a,1)^2+vec(a,1)^2));   
    end
end
function IM2d_blob=link_neck(IM2d,k2,BWedge,thed2)
% Seperate attached cells into blobs using neck points k2

if isempty(k2)==1
    IM2d_blob=IM2d; % if k2 is empty
elseif size(k2,2)==1
    IM2d_blob=IM2d; % if only one k2
else
    nn=size(k2,2); %number of neck points:nn
    k2points=[k2;zeros(2,nn)]; %if #k2 >1, start connect neck points
    for i=1:nn
        k2points(2:3,i)=[BWedge(k2(1,i),1);BWedge(k2(1,i),2)]; % neck point xy coordinates
    end
    % distance matrix of neck points: k2dist=[starti;endi;dist ..]
    k2dist0=[];
    for i=1:nn-1
        k2dist0=[k2dist0 [i*ones(1,nn-i);i+1:nn]];        
    end     
    k2dist=[k2dist0;pdist(k2points(2:3,:)')];
    % linkage from distance matrix k2dist
    %link_pair=zeros(3,fix(nn/2));
    link_pair=[];lpend=size(k2dist,2);
    for k=1:lpend
        [a,b]=min(k2dist(3,:));
        if a>thed2
            link_pair=[link_pair k2dist(:,b)];kde=size(link_pair,2);
            for kk=1:2 %delete matrix contains link_pair elements
                k2dist=k2dist(:,k2dist(kk,:)~=link_pair(1,kde));k2dist=k2dist(:,k2dist(kk,:)~=link_pair(2,kde));
            end
        else
            if b==1
                k2dist=k2dist(:,2:end);
            elseif b==size(k2dist,2)
                k2dist=k2dist(:,1:end-1);
            else
                k2dist=[k2dist(:,1:b-1) k2dist(:,b+1:end)];
            end
        end
    end
    IM2d_blob=separate_bob(IM2d,k2points,link_pair);
end
    function IM2d_blob=separate_bob(IM2d,k2points,link_pair)
        for j=1:size(link_pair,2)
            x=[k2points(2,link_pair(1,j)) k2points(2,link_pair(2,j))];
            y=[k2points(3,link_pair(1,j)) k2points(3,link_pair(2,j))];
            if abs(x(2)-x(1))>=abs(y(2)-y(1))
                if x(end)-x(1)>=0
                    X=x(1):x(end);
                else
                    X=x(1):-1:x(end);
                end
                Y=interp1(x,y,X);Y=round(Y);
            else
                if y(end)-y(1)>=0
                    Y=y(1):y(end);
                else
                    Y=y(1):-1:y(end);
                end
                X=interp1(y,x,Y);X=round(X);
            end           
            for jj=1:size(X,2)
                IM2d(Y(1,jj),X(1,jj))=0;
            end
        end
        IM2d_blob=IM2d;
    end
end