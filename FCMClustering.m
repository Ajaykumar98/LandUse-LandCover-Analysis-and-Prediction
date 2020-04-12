clc;
clf;
clear;


SR=imread('Region2finalPrediction.tif');
I1=SR(:,:,1);
I1s=I1;   
I2=SR(:,:,2);
I3=SR(:,:,3);
I4=SR(:,:,4);

%Extracting the four bands from the stacked image
 I1 = imresize(I1,[256 256]);
 I2 = imresize(I2,[256 256]);
 I3 = imresize(I3,[256 256]);
 I4 = imresize(I4,[256 256]);

climage=zeros(size(I1,1),size(I1,2),3);



fid1 = fopen('features.txt','wt');
for ii = 1:(size(I1,1)-1)
    for jj = 1:(size(I1,2))
        fprintf(fid1,'%d,%d,%d,%d\n',I1(ii,jj),I2(ii,jj),I3(ii,jj),I4(ii,jj));
    end    
    %fprintf(notfid1,'%d',I4(ii,size(I1,2))); 
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

n_clusters = 4;         
options = [2.0 25 0.001 1];

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

Nume1=double(int16(I3)-int16(I2));
Nume2=double(int16(I1)-int16(I3));
Nume3=double(int16(I4)-int16(I3));



Denom1=double(int16(I3)+int16(I2));
Denom2=double(int16(I1)+int16(I3));
Denom3=double(int16(I4)+int16(I3));


NDVI = Nume1 ./ Denom1;
NDWI = Nume2 ./ Denom2;
NDBI = Nume3 ./ Denom3;

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

clbarray=zeros(1,4);
clvarray=zeros(1,4);
clwarray=zeros(1,4);

clbarray(1)=cl1bvalue;
clbarray(2)=cl2bvalue;
clbarray(3)=cl3bvalue;
clbarray(4)=cl4bvalue;

clvarray(1)=cl1vvalue;
clvarray(2)=cl2vvalue;
clvarray(3)=cl3vvalue;
clvarray(4)=cl4vvalue;

clwarray(1)=cl1wvalue;
clwarray(2)=cl2wvalue;
clwarray(3)=cl3wvalue;
clwarray(4)=cl4wvalue;

classign= zeros(1,4);
classignvalues=zeros(1,4);


clvflag=0;
clwflag=0;
clbflag=0;

clflagarray=[clvflag clwflag clbflag];

while check(clflagarray)

    if clflagarray(3) == 0
    
        %For Builtups
        tindex=retmaxi(clbarray);
        %clbarray(tindex)=-inf;
        %clbflag=1;

        if classign(tindex)==0
            classign(tindex)=0100;
            classignvalues(tindex)=clbarray(tindex);
            clbarray(tindex)=-inf;
            clflagarray(3)=1;
        end
    end
    
    if clflagarray(2)==0
        
         %For Water
        tindex=retmaxi(clwarray);
        %clwarray(tindex)=-inf;
        %clwflag=1;

        if classign(tindex)==0
            classign(tindex)=0001;
            classignvalues(tindex)=clwarray(tindex);
            clwarray(tindex)=-inf;
            clflagarray(2)=1;

        elseif classign(tindex)==0100
            if clwarray(tindex) > classignvalues(tindex)
                classign(tindex)=0001;
                classignvalues(tindex)=clwarray(tindex);
                clwarray(tindex)=-inf;
                clflagarray(2)=1;
                clflagarray(3)=0;
            else
                clflagarray(2)=0;
                clbarray(tindex)=-inf;
            end     

        elseif classign(tindex)==0010
            if clwarray(tindex) > classignvalues(tindex)
                classign(tindex)=0001;
                classignvalues(tindex)=clwarray(tindex);
                clwarray(tindex)=-inf;
                clflagarray(2)=1;
                clflagarray(1)=0;
            else
                clflagarray(2)=0;
                clbarray(tindex)=-inf;
            end
        end
    end
    
    %For Vegetation
    if clflagarray(1)==0
          tindex=retmaxi(clvarray);
        %clwarray(tindex)=-inf;
        %clwflag=1;

        if classign(tindex)==0
            classign(tindex)=0010;
            classignvalues(tindex)=clvarray(tindex);
            clvarray(tindex)=-inf;
            clflagarray(1)=1;

        elseif classign(tindex)==0100
            if clvarray(tindex) > classignvalues(tindex)
                classign(tindex)=0010;
                classignvalues(tindex)=clvarray(tindex);
                clvarray(tindex)=-inf;
                clflagarray(1)=1;
                clflagarray(3)=0;
            else
                clflagarray(1)=0;
                clvarray(tindex)=-inf;
            end

        elseif classign(tindex)==0001
            if clvarray(tindex) > classignvalues(tindex)
                classign(tindex)=0010;
                classignvalues(tindex)=clvarray(tindex);
                clvarray(tindex)=-inf;
                clflagarray(1)=1;
                clflagarray(2)=0;
            else
                clflagarray(1)=0;
                clvarray(tindex)=-inf;
            end
        end
    end
        
        
end

for iim =1:4
    if classign(iim)==0
        classign(iim)=1000;
        break;
    end
end

cluster1c=classign(1);
cluster2c=classign(2);
cluster3c=classign(3);
cluster4c=classign(4);

% 
% if cl1vvalue >= cl2vvalue && cl1vvalue >= cl3vvalue && cl1vvalue >= cl4vvalue 
%     cluster1c=0010;
% elseif cl2vvalue >= cl1vvalue && cl2vvalue >= cl3vvalue && cl2vvalue >= cl4vvalue    
%     cluster2c=0010;
% elseif cl3vvalue >= cl1vvalue && cl3vvalue >= cl2vvalue  && cl3vvalue >= cl4vvalue   
%     cluster3c=0010;
% elseif cl4vvalue >= cl1vvalue && cl4vvalue >= cl2vvalue  && cl4vvalue >= cl3vvalue   
%     cluster4c=0010;
% end
% 
% 
% if cl1bvalue >= cl2bvalue && cl1bvalue >= cl3bvalue && cl1bvalue >= cl4bvalue 
%     cluster1c=0100;
% elseif cl2bvalue >= cl1bvalue && cl2bvalue >= cl3bvalue && cl2bvalue >= cl4bvalue    
%     cluster2c=0100;
% elseif cl3bvalue >= cl1bvalue && cl3bvalue >= cl2bvalue  && cl3bvalue >= cl4bvalue   
%     cluster3c=0100;
% elseif cl4bvalue >= cl1bvalue && cl4bvalue >= cl2bvalue  && cl4bvalue >= cl3bvalue   
%     cluster4c=0100;
% end
% 
% if cl1wvalue >= cl2wvalue && cl1wvalue >= cl3wvalue && cl1wvalue >= cl4wvalue 
%     cluster1c=0001;
% elseif cl2wvalue >= cl1wvalue && cl2wvalue >= cl3wvalue && cl2wvalue >= cl4wvalue    
%     cluster2c=0001;
% elseif cl3wvalue >= cl1wvalue && cl3wvalue >= cl2wvalue  && cl3wvalue >= cl4wvalue   
%     cluster3c=0001;
% elseif cl4wvalue >= cl1wvalue && cl4wvalue >= cl2wvalue  && cl4wvalue >= cl3wvalue   
%     cluster4c=0001;
% end

% if cluster1c == 1111
%     cluster1c=1000;
% elseif cluster2c == 1111
%     cluster2c=1000;
% elseif cluster3c == 1111
%     cluster3c=1000;
% elseif cluster4c == 1111
%     cluster4c=1000;
% end

%climage=zeros(size(I1,1),size(I1,2),3);

for iim = 1:(size(I1,1))
    for jjm = 1:(size(I1,2))
        if cllabels((iim-1)*size(I1,2)+jjm)==1
            if cluster1c == 0100 
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
            if cluster2c == 0100 
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
            if cluster3c == 0100 
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
            if cluster4c == 0100 
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

%imwrite(climage,'Predicted256png.png')

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


function output=check(ipvector)

flag=0;
for ii=1:size(ipvector,2)
    if ipvector(ii)==0
        flag=1;
        break
    end
end

if flag==1
    output=1;
else 
    output=0;
end
end

function output2=retmaxi(ipvector2)

temp=-inf;
tindex=-inf;
for ii2=1:size(ipvector2,2)
    if ipvector2(ii2)>temp
        temp=ipvector2(ii2);
        tindex=ii2;
    end
end
output2=tindex;

end


