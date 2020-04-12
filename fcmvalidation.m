clc;
clf;
clear;

SR=imread('l3d44o0401apr14/final.tif');
I1=SR(:,:,1);
I1s=I1;   
I2=SR(:,:,2);
I3=SR(:,:,3);
I4=SR(:,:,4);


% I1 = imresize(I1,[256 256]);
% I2 = imresize(I2,[256 256]);
% I3 = imresize(I3,[256 256]);
% I4 = imresize(I4,[256 256]);

climage=zeros(size(I1,1),size(I1,2),3);



fid1 = fopen('features.txt','wt');
for ii = 1:(size(I1,1)-1)
    for jj = 1:(size(I1,2))
        fprintf(fid1,'%d,%d,%d,%d\n',I1(ii,jj),I2(ii,jj),I3(ii,jj),I4(ii,jj));
    end    
    %fprintf(fid1,'%d',I4(ii,size(I1,2))); 
end



for jj2 = 1:(size(I1,2)-1)
    fprintf(fid1,'%d,%d,%d,%d\n',I1(size(I1,1),jj2),I2(size(I1,1),jj2),I3(size(I1,1),jj2),I4(size(I1,1),jj2));
end

fprintf(fid1,'%d,%d,%d,%d',I1(size(I1,1),size(I1,2)),I2(size(I1,1),size(I1,2)),I3(size(I1,1),size(I1,2)),I4(size(I1,1),size(I1,2)));
%fprintf(fid1,'\n');

%I1(1)
img8a = uint8(I1);
%imshow(img8a);

img8b = uint8(I2);
%figure,imshow(img8b);

img8c = uint8(I3);
%figure,imshow(img8c);

img8d = uint8(I4);
%figure,imshow(img8d);

data=load('features.txt');

%  data = load('fcmdata.dat');  % load some sample data
n_clusters = 4;         
options = [2.0 25 0.001 1];
% number of clusters
[center,U,obj_fcn] = fcm(data, n_clusters,options);

 maxU = max(U);
 
cluster1 = find(U(1,:) == maxU);
cluster2 = find(U(2,:) == maxU);
cluster3 = find(U(3,:) == maxU);
cluster4 = find(U(4,:) == maxU);
% cllabels=cluster(U,data);
idx=maxU;

cllabels=zeros(1,65536);

for cli1 = 1:(numel(cluster1))
    cllabels(cluster1(cli1))=1;
end

for cli2 = 1:(numel(cluster2))
    cllabels(cluster2(cli2))=2;
end

for cli3 = 1:(numel(cluster3))
    cllabels(cluster3(cli3))=3;
end


for cli4 = 1:(numel(cluster4))
    cllabels(cluster4(cli4))=4;
end
% % % % % % % % % % % % % % 
% cluster1 = (cllabels == 1); % |1| for cluster 1 membership
% cluster2 = (cllabels == 2); % |2| for cluster 2 membership
% cluster3 = (cllabels == 3);
% cluster4 = (cllabels == 4);
% figure
% gscatter(data(:,1),data(:,2),cllabels,'rb','+o')
% legend('Cluster 1','Cluster 2','Location','best')

%cllabels=cllabels.';
% 

NDVI = double(I3 - I2) ./ double(I3+I2);
NDWI = double(I1-I3) ./ double(I1+I3);
NDBI = double(I4-I3) ./ double(I4+I3);

cl1wvalue = double(0.0);
cl2wvalue = double(0.0);
cl3wvalue = double(0.0);
cl4wvalue = double(0.0);

cl1wcount=0;
cl2wcount=0;
cl3wcount=0;
cl4wcount=0;

for iim2 = 1:(size(I1,1))
    for jjm2 = 1:(size(I1,2))
        if cllabels((iim2-1)*size(I1,2)+jjm2)==1 && ~(isnan(NDWI(iim2,jjm2)))
            cl1wvalue = double(cl1wvalue) + double(NDWI(iim2,jjm2));
            cl1wcount=cl1wcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==2 && ~(isnan(NDWI(iim2,jjm2)))
            cl2wvalue = double(cl2wvalue) + double(NDWI(iim2,jjm2));
            cl2wcount=cl2wcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==3 && ~(isnan(NDWI(iim2,jjm2)))
            cl3wvalue = double(cl3wvalue) + double(NDWI(iim2,jjm2));
            cl3wcount=cl3wcount+1;   
        
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==4 && ~(isnan(NDWI(iim2,jjm2)))
            cl4wvalue = double(cl4wvalue) + double(NDWI(iim2,jjm2));
            cl4wcount=cl4wcount+1;   
        end    
    end
end

cl1wvalue = cl1wvalue / cl1wcount;
cl2wvalue = cl2wvalue / cl2wcount;
cl3wvalue = cl3wvalue / cl3wcount;
cl4wvalue = cl4wvalue / cl4wcount;

cl1wvalue = cl1wvalue*100000;
cl2wvalue = cl2wvalue*100000;
cl3wvalue = cl3wvalue*100000;
cl4wvalue = cl4wvalue*100000;

cl1vvalue = double(0.0);
cl2vvalue = double(0.0);
cl3vvalue = double(0.0);
cl4vvalue = double(0.0);

cl1vcount=0;
cl2vcount=0;
cl3vcount=0;
cl4vcount=0;

for iim2 = 1:(size(I1,1))
    for jjm2 = 1:(size(I1,2))
        if cllabels((iim2-1)*size(I1,2)+jjm2)==1 && ~(isnan(NDVI(iim2,jjm2)))
            cl1vvalue = double(cl1vvalue) + double(NDVI(iim2,jjm2));
            cl1vcount=cl1vcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==2 && ~(isnan(NDVI(iim2,jjm2)))
            cl2vvalue = double(cl2vvalue) + double(NDVI(iim2,jjm2));
            cl2vcount=cl2vcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==3 && ~(isnan(NDVI(iim2,jjm2)))
            cl3vvalue = double(cl3vvalue) + double(NDVI(iim2,jjm2));
            cl3vcount=cl3vcount+1;   
        
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==4 && ~(isnan(NDVI(iim2,jjm2)))
            cl4vvalue = double(cl4vvalue) + double(NDVI(iim2,jjm2));
            cl4vcount=cl4vcount+1;   
        end    
    end
end

cl1vvalue = cl1vvalue / cl1vcount;
cl2vvalue = cl2vvalue / cl2vcount;
cl3vvalue = cl3vvalue / cl3vcount;
cl4vvalue = cl4vvalue / cl4vcount;

cl1vvalue = cl1vvalue*100000;
cl2vvalue = cl2vvalue*100000;
cl3vvalue = cl3vvalue*100000;
cl4vvalue = cl4vvalue*100000;


cl1bvalue = double(0.0);
cl2bvalue = double(0.0);
cl3bvalue = double(0.0);
cl4bvalue = double(0.0);

cl1bcount=0;
cl2bcount=0;
cl3bcount=0;
cl4bcount=0;

for iim2 = 1:(size(I1,1))
    for jjm2 = 1:(size(I1,2))
        if cllabels((iim2-1)*size(I1,2)+jjm2)==1 && ~(isnan(NDBI(iim2,jjm2)))
            cl1bvalue = double(cl1bvalue) + double(NDBI(iim2,jjm2));
            cl1bcount=cl1bcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==2 && ~(isnan(NDBI(iim2,jjm2)))
            cl2bvalue = double(cl2bvalue) + double(NDBI(iim2,jjm2));
            cl2bcount=cl2bcount+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==3 && ~(isnan(NDBI(iim2,jjm2)))
            cl3bvalue = double(cl3bvalue) + double(NDBI(iim2,jjm2));
            cl3bcount=cl3bcount+1;   
        
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==4 && ~(isnan(NDBI(iim2,jjm2)))
            cl4bvalue = double(cl4bvalue) + double(NDBI(iim2,jjm2));
            cl4bcount=cl4bcount+1;   
        end    
    end
end

cl1bvalue = cl1bvalue / cl1bcount;
cl2bvalue = cl2bvalue / cl2bcount;
cl3bvalue = cl3bvalue / cl3bcount;
cl4bvalue = cl4bvalue / cl4bcount;

cl1bvalue = cl1bvalue*100000;
cl2bvalue = cl2bvalue*100000;
cl3bvalue = cl3bvalue*100000;
cl4bvalue = cl4bvalue*100000;




cluster1c=1111;
cluster2c=1111;
cluster3c=1111;
cluster4c=1111;





if cl1vvalue >= cl2vvalue && cl1vvalue >= cl3vvalue && cl1vvalue >= cl4vvalue 
    cluster1c=0010;
elseif cl2vvalue >= cl1vvalue && cl2vvalue >= cl3vvalue && cl2vvalue >= cl4vvalue    
    cluster2c=0010;
elseif cl3vvalue >= cl1vvalue && cl3vvalue >= cl2vvalue  && cl3vvalue >= cl4vvalue   
    cluster3c=0010;
elseif cl4vvalue >= cl1vvalue && cl4vvalue >= cl2vvalue  && cl4vvalue >= cl3vvalue   
    cluster4c=0010;
end


if cl1bvalue >= cl2bvalue && cl1bvalue >= cl3bvalue && cl1bvalue >= cl4bvalue 
    cluster1c=0100;
elseif cl2bvalue >= cl1bvalue && cl2bvalue >= cl3bvalue && cl2bvalue >= cl4bvalue    
    cluster2c=0100;
elseif cl3bvalue >= cl1bvalue && cl3bvalue >= cl2bvalue  && cl3bvalue >= cl4bvalue   
    cluster3c=0100;
elseif cl4bvalue >= cl1bvalue && cl4bvalue >= cl2bvalue  && cl4bvalue >= cl3bvalue   
    cluster4c=0100;
end

if cl1wvalue >= cl2wvalue && cl1wvalue >= cl3wvalue && cl1wvalue >= cl4wvalue 
    cluster1c=0001;
elseif cl2wvalue >= cl1wvalue && cl2wvalue >= cl3wvalue && cl2wvalue >= cl4wvalue    
    cluster2c=0001;
elseif cl3wvalue >= cl1wvalue && cl3wvalue >= cl2wvalue  && cl3wvalue >= cl4wvalue   
    cluster3c=0001;
elseif cl4wvalue >= cl1wvalue && cl4wvalue >= cl2wvalue  && cl4wvalue >= cl3wvalue   
    cluster4c=0001;
end

if cluster1c == 1111
    cluster1c=1000;
elseif cluster2c == 1111
    cluster2c=1000;
elseif cluster3c == 1111
    cluster3c=1000;
elseif cluster4c == 1111
    cluster4c=1000;
end

%climage=zeros(size(I1,1),size(I1,2),3);

for iim = 1:(size(I1,1))
    for jjm = 1:(size(I1,2))
        if cllabels((iim-1)*size(I1,2)+jjm)==1
            if cluster1c == 0100 || cluster1c == 1111
                climage(iim,jjm,1)=1.0;
            elseif cluster1c == 0010
                climage(iim,jjm,2)=1.0;
            elseif cluster1c == 0001
                climage(iim,jjm,3)=1.0;
            elseif cluster1c == 1000
                climage(iim,jjm,1)=1.0;
                climage(iim,jjm,2)=1.0;
              
            
            end
            

        elseif cllabels((iim-1)*size(I1,2)+jjm)==2
            if cluster2c == 0100 || cluster2c==1111
                climage(iim,jjm,1)=1.0;
            elseif cluster2c == 0010
                climage(iim,jjm,2)=1.0;
            elseif cluster2c == 0001
                climage(iim,jjm,3)=1.0;
            elseif cluster2c == 1000
                climage(iim,jjm,1)=1.0;
                climage(iim,jjm,2)=1.0;
            
            end
            
        elseif cllabels((iim-1)*size(I1,2)+jjm)==3
            if cluster3c == 0100 || cluster3c==1111
                climage(iim,jjm,1)=1.0;
            elseif cluster3c == 0010
                climage(iim,jjm,2)=1.0;
            elseif cluster3c == 0001
                climage(iim,jjm,3)=1.0;
            elseif cluster3c == 1000
                climage(iim,jjm,1)=1.0;
                climage(iim,jjm,2)=1.0;
            end
            
        elseif cllabels((iim-1)*size(I1,2)+jjm)==4
            if cluster4c == 0100 || cluster4c==1111
                climage(iim,jjm,1)=1.0;
            elseif cluster4c == 0010
                climage(iim,jjm,2)=1.0;
            elseif cluster4c == 0001
                climage(iim,jjm,3)=1.0;
            elseif cluster4c == 1000
                climage(iim,jjm,1)=1.0;
                climage(iim,jjm,2)=1.0;
            end
            
        end
    end
end


% for iim = 1:(size(I1,1))
%     for jjm = 1:(size(I1,2))
%         if cllabels((iim-1)*size(I1,2)+jjm)==1
%             climage(iim,jjm,1)=1.0;
%         elseif cllabels((iim-1)*size(I1,2)+jjm)==4
%             climage(iim,jjm,2)=1.0;
%         elseif cllabels((iim-1)*size(I1,2)+jjm)==2
%             climage(iim,jjm,3)=1.0;
%         elseif cllabels((iim-1)*size(I1,2)+jjm)==3
%            climage(iim,jjm,1)=1.0;
%            climage(iim,jjm,2)=1.0;
%         end
%     end
% end


imshow(climage);
%imwrite(climage,'Predicted256.jpg')

% index1 = find(U(1,:) == maxU);
% index2 = find(U(2,:) == maxU);
% plot(data(index1,1),data(index1,2),'ob')
% hold on
% plot(data(index2,1),data(index2,2),'or')
% plot(center(1,1),center(1,2),'xb','MarkerSize',15,'LineWidth',3)
% plot(center(2,1),center(2,2),'xr','MarkerSize',15,'LineWidth',3)
% hold off
% opt = genfisOptions('FCMClustering');
% opt.NumClusters = Nc;
% opt.Exponent = options(1);
% opt.MaxNumIteration = options(2);
% opt.MinImprovement = options(3);
% opt.Verbose = options(4);
% inputData = data(:,1:M);
% outputData = data(:,index1+1:end);
% figure();
% imshow(outputData);
% imshow(data);



I1 = double(I1);
I2 = double(I2);
I3 = double(I3);
I4 = double(I4);  





cl1totala=0;
cl1totalb=0;
cl1totalc=0;
cl1totald=0;
cl1n=0;

cl2totala=0;
cl2totalb=0;
cl2totalc=0;
cl2totald=0;
cl2n=0;

cl3totala=0;
cl3totalb=0;
cl3totalc=0;
cl3totald=0;
cl3n=0;

cl4totala=0;
cl4totalb=0;
cl4totalc=0;
cl4totald=0;
cl4n=0;



for iim2 = 1:(size(I1,1))
    for jjm2 = 1:(size(I1,2))
        if cllabels((iim2-1)*size(I1,2)+jjm2)==1
            cl1totala=cl1totala+I1(iim2,jjm2);
            cl1totalb=cl1totalb+I2(iim2,jjm2);
            cl1totalc=cl1totalc+I3(iim2,jjm2);
            cl1totald=cl1totald+I4(iim2,jjm2);
            cl1n=cl1n+1;
            
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==2
            cl2totala=cl2totala+I1(iim2,jjm2);
            cl2totalb=cl2totalb+I2(iim2,jjm2);
            cl2totalc=cl2totalc+I3(iim2,jjm2);
            cl2totald=cl2totald+I4(iim2,jjm2);
            cl2n=cl2n+1;
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==3
            cl3totala=cl3totala+I1(iim2,jjm2);
            cl3totalb=cl3totalb+I2(iim2,jjm2);
            cl3totalc=cl3totalc+I3(iim2,jjm2);
            cl3totald=cl3totald+I4(iim2,jjm2);
            cl3n=cl3n+1;
            
        elseif cllabels((iim2-1)*size(I1,2)+jjm2)==4
            cl4totala=cl4totala+I1(iim2,jjm2);
            cl4totalb=cl4totalb+I2(iim2,jjm2);
            cl4totalc=cl4totalc+I3(iim2,jjm2);
            cl4totald=cl4totald+I4(iim2,jjm2);
            cl4n=cl4n+1;
            
        end    
    end
end

cl1cena=cl1totala/cl1n;
cl1cenb=cl1totalb/cl1n;
cl1cenc=cl1totalc/cl1n;
cl1cend=cl1totald/cl1n;


cl2cena=cl2totala/cl2n;
cl2cenb=cl2totalb/cl2n;
cl2cenc=cl2totalc/cl2n;
cl2cend=cl2totald/cl2n;


cl3cena=cl3totala/cl3n;
cl3cenb=cl3totalb/cl3n;
cl3cenc=cl3totalc/cl3n;
cl3cend=cl3totald/cl3n;


cl4cena=cl4totala/cl4n;
cl4cenb=cl4totalb/cl4n;
cl4cenc=cl4totalc/cl4n;
cl4cend=cl4totald/cl4n;


cl1dist=0;
cl2dist=0;
cl3dist=0;
cl4dist=0;

for iim3 = 1:(size(I1,1))
    for jjm3 = 1:(size(I1,2))
        if cllabels((iim3-1)*size(I1,2)+jjm3)==1
            cl1dist = cl1dist + sqrt((I1(iim3,jjm3)-cl1cena)^2 + (I2(iim3,jjm3)-cl1cenb)^2 +(I3(iim3,jjm3)-cl1cenc)^2 +(I4(iim3,jjm3)-cl1cend)^2);
            
        elseif cllabels((iim3-1)*size(I1,2)+jjm3)==2
            cl2dist = cl2dist + sqrt((I1(iim3,jjm3)-cl2cena)^2 + (I2(iim3,jjm3)-cl2cenb)^2 +(I3(iim3,jjm3)-cl2cenc)^2 +(I4(iim3,jjm3)-cl2cend)^2);
        elseif cllabels((iim3-1)*size(I1,2)+jjm3)==3
            cl3dist = cl3dist + sqrt((I1(iim3,jjm3)-cl3cena)^2 + (I2(iim3,jjm3)-cl3cenb)^2 +(I3(iim3,jjm3)-cl3cenc)^2 +(I4(iim3,jjm3)-cl3cend)^2);
        elseif cllabels((iim3-1)*size(I1,2)+jjm3)==4
            cl4dist = cl4dist + sqrt((I1(iim3,jjm3)-cl4cena)^2 + (I2(iim3,jjm3)-cl4cenb)^2 +(I3(iim3,jjm3)-cl4cenc)^2 +(I4(iim3,jjm3)-cl4cend)^2);
            
        end    
    end
end

cl1dist = cl1dist/cl1n;
cl2dist = cl2dist/cl2n;
cl3dist = cl3dist/cl3n;
cl4dist = cl4dist/cl4n;

tdistance = cl1dist+cl2dist+cl3dist+cl4dist;