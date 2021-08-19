function [hyperbolaLocations]= ShearletPreprocess(GPRData_org)
%[rebarSignalImage, hyperbolaLocations]= cplx_auto_ShearletPreprocess(GPRData)

%GPRData = imresize(GPRData, 2);

addpath('shearlet_toolbox')

GPRData = imresize(GPRData_org(:,:,1),[512,512]);

% Load image
%x=double(imread('barbara.gif'));
x=double(GPRData);

[L L]=size(x);


% Create noisy image
sigma=20;
x_noisy=x+sigma.*randn(L,L);

% setup parameters for shearlet transform
lpfilt='maxflat';
% .dcomp(i) indicates there will be 2^dcomp(i) directions 
shear_parameters.dcomp =[ 3  3  4  4];
% .dsize(i) indicate the local directional filter will be
% dsize(i) by dsize(i)
shear_parameters.dsize =[32 32 16 16];
% 
%Tscalars determine the thresholding multipliers for
%standard deviation noise estimates. Tscalars(1) is the
%threshold scalar for the low-pass coefficients, Tscalars(2)
%is the threshold scalar for the band-pass coefficients, 
%Tscalars(3) is the threshold scalar for the high-pass
%coefficients. 

Tscalars=[0 3 4];

%There are three possible ways of implementing the 
%local nonsubsampled shearlet transform (nsst_dec1e,
%nsst_dec1, nsst_dec2). For this demo, we have created 
%a flag called shear_version to choose which one to
%test.

%shear_version=0; %nsst_dec1e
shear_version=1; %nsst_dec1
%shear_version=2; %nsst_dec2

% compute the shearlet decompositon
if shear_version==0
  [dst,shear_f]=nsst_dec1e(x_noisy,shear_parameters,lpfilt);
elseif shear_version==1
  [dst,shear_f]=nsst_dec1(x_noisy,shear_parameters,lpfilt);
elseif shear_version==2
  [dst,shear_f]=nsst_dec2(x_noisy,shear_parameters,lpfilt);
end

% Determines via Monte Carlo the standard deviation of
% the white Gaussian noise for each scale and 
% directional component when a white Gaussian noise of
% standard deviation of 1 is feed through.
if shear_version==0
   dst_scalars=nsst_scalars_e(L,L,shear_f,lpfilt);
else
   dst_scalars=nsst_scalars(L,L,shear_f,lpfilt);
end



% apply hard threshold to the shearlet coefficients
dst=nsst_HT(dst,sigma,Tscalars,dst_scalars);



display_flag=0;
if display_flag==1
   figure('name', 'preprocess', 'visible','off');
   imagesc(dst{1})
   %imshow(uint8(dst{1}))
   for i=1:length(dst)-4
       l=size(dst{i+1},3);
       JC=ceil(l/2);
       JR=ceil(l/JC);
       figure('name', 'preprocess', 'visible','on');
       %ttt=dst{i+1}(:,:,1);
       for k=1:l
           subplot(JR,JC,k)
           %ttt=ttt+dst{i+1}(:,:,k);
           imagesc(abs(dst{i+1}(:,:,k)))
           axis off
           axis image
       end 
       %{
       ttt=uint8(ttt);
       x=[0 8 16 180 200 255];
        Y=[0 70000 8000 0 3139 0];
        xi=0:255;
        yi = interp1q(x',Y',xi');
        ImgHistMatch1 = histeq(ttt, yi);
        threshold = graythresh(ImgHistMatch1);
        bwImage1 = imbinarize(ImgHistMatch1, threshold);
        
       imshow(bwImage1)
       axis off
       axis image
       %}
   end
end % display_flag 

%{
figure('name', 'preprocess', 'visible','on');
imagesc(abs(dst{2}(:,:,1)));axis off
figure('name', 'preprocess', 'visible','on');
imagesc(abs(dst{2}(:,:,4)));axis off
figure('name', 'preprocess', 'visible','on');
imagesc(abs(dst{2}(:,:,5)));axis off
figure('name', 'preprocess', 'visible','on');
imagesc(abs(dst{2}(:,:,8)));axis off
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate hyperbola locations (vertival and horizontal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1 = 10*uint8(abs(dst{2}(:,:,1)));
part4 = 10*uint8(abs(dst{2}(:,:,4)));
part5 = 10*uint8(abs(dst{2}(:,:,5)));
part8 = 10*uint8(abs(dst{2}(:,:,8)));
[horiLocation1] = cplx_auto_horizontalLocation(part1); % initial 5 horizontal positions
[horiLocation4] = cplx_auto_horizontalLocation(part4);
[horiLocation5] = cplx_auto_horizontalLocation(part5);
[horiLocation8] = cplx_auto_horizontalLocation(part8);







horiLocation_l = [horiLocation4;horiLocation5];  % left
horiLocation_l = sort(horiLocation_l);
rsrtH_l = [horiLocation_l(1)];
k=1;
for i=2:length(horiLocation_l)
    if horiLocation_l(i)-rsrtH_l(k)<60               % combine closed points left
        rsrtH_l(k) = (rsrtH_l(k)+horiLocation_l(i))/2;
        continue;
    end
    rsrtH_l = [rsrtH_l;horiLocation_l(i)];
    k=k+1;
end


horiLocation_r = [horiLocation1;horiLocation8];
horiLocation_r = sort(horiLocation_r);
rsrtH_r = [horiLocation_r(1)];
k=1;
for i=2:length(horiLocation_r)
    if horiLocation_r(i)-rsrtH_r(k)<60                % combine closed points right
        rsrtH_r(k) = (rsrtH_r(k)+horiLocation_r(i))/2;
        continue;
    end
    rsrtH_r = [rsrtH_r;horiLocation_r(i)];
    k=k+1;
end



% combine 1 4 5 6 7 8
GPRDataRM = 10*uint8(abs(dst{2}(:,:,1)+dst{2}(:,:,4)+dst{2}(:,:,5)...
    +dst{2}(:,:,6)+dst{2}(:,:,7)+dst{2}(:,:,8)));
%subplot(2,3,2),
%figure('name', 'preprocess', 'visible','on');
%imshow(GPRDataRM);hold on

%zrv = zeros(size(rsrtH_l))+50;
%plot(rsrtH_l,zrv,'rx','MarkerSize',14,'LineWidth',1.5);hold on % left

%zrv = zeros(size(rsrtH_r))+50;
%plot(rsrtH_r,zrv,'bx','MarkerSize',14,'LineWidth',1.5);hold on % right


%{
% generate the horizontal points of each hyperbola
rbLocationH = zeros(max(length(rsrtH_l),length(rsrtH_r)),1);
if rsrtH_l(1)<rsrtH_r(1)
    for i=1:length(rbLocationH)-1
        rbLocationH(i) = (rsrtH_l(i)+rsrtH_r(i))/2;
    end
    if length(rsrtH_l)==length(rsrtH_r)
        rbLocationH(i+1) = (rsrtH_l(end)+rsrtH_r(end))/2;
    else
        rsrtH_r(length(rsrtH_r)+1)=rsrtH_r(end)+20;
        rbLocationH(i+1) = rsrtH_l(end)+20;
    end    
elseif rsrtH_l(1)>rsrtH_r(1)
    rbLocationH(1)=rsrtH_r(1)-20;
    for i=2:length(rbLocationH)
        rbLocationH(i) = (rsrtH_l(i-1)+rsrtH_r(i))/2;
    end    
end

if rsrtH_l(end)>rsrtH_r(end) && length(rsrtH_r)<5
    rbLocationH(end+1)= rsrtH_l(end)+20;
end
%}

% generate the horizontal points of each hyperbola
entrLR = [rsrtH_l;rsrtH_r];
entrLR = sort(entrLR);
newEntrLR = [];
k=1;
if sum(ismember(rsrtH_l,entrLR(1)))==1  % 1st point is left
    newEntrLR(k)=entrLR(1);k=k+1;
else                                    % 1st point is right
    newEntrLR(k)=entrLR(1)-50;k=k+1;
    newEntrLR(k)=entrLR(1);k=k+1;
end
for i=2:length(entrLR)-1
    % left right
    if sum(ismember(rsrtH_l,entrLR(i)))==1 && sum(ismember(rsrtH_r,entrLR(i+1)))==1
        newEntrLR(k)=entrLR(i);k=k+1;
    end
    % right left
    if sum(ismember(rsrtH_r,entrLR(i)))==1 && sum(ismember(rsrtH_l,entrLR(i+1)))==1
        newEntrLR(k)=entrLR(i);k=k+1;
    end
    % left left
    if sum(ismember(rsrtH_l,entrLR(i)))==1 && sum(ismember(rsrtH_l,entrLR(i+1)))==1
        newEntrLR(k)=entrLR(i);k=k+1;
        newEntrLR(k)=entrLR(i)+20;k=k+1;
    end
    % right right
    if sum(ismember(rsrtH_r,entrLR(i)))==1 && sum(ismember(rsrtH_r,entrLR(i+1)))==1
        newEntrLR(k)=entrLR(i);k=k+1;
        newEntrLR(k)=entrLR(i+1)-20;k=k+1;
    end
end
if sum(ismember(rsrtH_r,entrLR(end)))==1  % last point is right
    newEntrLR(k)=entrLR(end);
else                                      % last point is left
    newEntrLR(k)=entrLR(end)+50;k=k+1;
    newEntrLR(k)=entrLR(end);
end

% insert missing points
rbLocationH=[];
k=1;
for i=1:2:length(newEntrLR)-1
    if newEntrLR(i+1)-newEntrLR(i)>110
        rbLocationH(k)=newEntrLR(i)+20;k=k+1;
        rbLocationH(k)=newEntrLR(i+1)-20;k=k+1;
    else        
        rbLocationH(k)=(newEntrLR(i)+newEntrLR(i+1))/2;k=k+1;
    end
end





% generate the vertical points of each hyperbola
[gY,~]=find(GPRDataRM>100);
rbLocationV = mean(gY)*ones(size(rbLocationH));






hyperbolaLocations = [rbLocationH' rbLocationV'];



%plot(rbLocationH,rbLocationV,'yx','MarkerSize',14,'LineWidth',1.5);
%plot(rbLocationH,[35 35 35 35 35],'yx','MarkerSize',14,'LineWidth',1.5);










