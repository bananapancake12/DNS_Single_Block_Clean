% Read Spectra (spU and spUV)
% write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
% write(10) N
% write(10) y,dthetai,dthdy
% write(10) istat
% write(10) buffSP

% ATTENTION ribh hardcoded
% ATTENTION dudy has to be calculated using the same scheme as in the code
clear all
%fname  = 'stats001x001_old.dat';
fname  = 'stats004x008.dat';

yinterp = .70;
ptstilex = 24;
ptstilez = 24;
jribn = 1;

fid    = fopen(fname,'r');
         fread(fid,1,'int32');
t      = fread(fid,1,'float64');
Re     = fread(fid,1,'float64');
alp    = fread(fid,1,'float64');
bet    = fread(fid,1,'float64');
Lx = 2*pi/alp;
Lz = 2*pi/bet;
Ly = 2;
mpgx   = fread(fid,1,'float64');
nband  = fread(fid,1,'int32');
iter   = fread(fid,1,'int32');
         fread(fid,88,'int32');
         fread(fid,2,'int32');
Ngal   = fread(fid,[4 nband+2],'int32');
         fread(fid,2,'int32');
nyu     = Ngal(4,nband+1)-Ngal(4,1)+2;
yu      = fread(fid,nyu,'float64');
dthetai= fread(fid,1,'float64');
dthdyu  = fread(fid,nyu,'float64');
         fread(fid,2,'int32');
         
nyv     = Ngal(3,nband+1)-Ngal(3,1)+2;
yv      = fread(fid,nyv,'float64');
dthetai= fread(fid,1,'float64');
dthdyv  = fread(fid,nyv,'float64');
         fread(fid,2,'int32');
istat  = fread(fid,1,'int32');

nyp = nyu-2;

% tiles
tilex = Ngal(1,2)*2/3/ptstilex;
tilez = Ngal(2,2)*2/3/ptstilez;

% Y-coordinate of bands shifted up Ngal(3,1)+1 units.
% Y-coordinate of bands shifted up Ngal(3,1)+1 units.
nyyu          = Ngal(4,1:nband+1)-Ngal(4,1)+1;
nyyu(1    )   = 0;
nyyu(nband+1) = nyyu(nband+1)+1;
nyu1          = nyyu(2);

nyyv          = Ngal(3,1:nband+1)-Ngal(3,1)+1;
nyyv(1    )   = 0;
nyyv(nband+1) = nyyv(nband+1)+1;
nyv1          = nyyv(2);

nyyp = nyyu;
nyyp(1) = nyyp(1);
nyyp(2) = nyyp(2)-2;
nyyp(3) = nyyp(3);
nyyp(end) = nyyp(end)-2;
nyp1          = nyyp(2);

ypp = yu(2:end-1);

% istat stores the number of statistics accumulated in the file, so dividing
%  by istat we obtain the average.
% nstat stores the number of statistics accumulated times the grid size at each plane
nstatu = 0*yu+istat;
for iband = 1:nband
  nstatu(nyyu(iband)+1:nyyu(iband+1)) =  nstatu(nyyu(iband)+1:nyyu(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end

nstatv = 0*yv+istat;
for iband = 1:nband
  nstatv(nyyv(iband)+1:nyyv(iband+1)) =  nstatv(nyyv(iband)+1:nyyv(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end


nstatp = 0*ypp+istat;
for iband = 1:nband
  nstatp(nyyp(iband)+1:nyyp(iband+1)) =  nstatp(nyyp(iband)+1:nyyp(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end


% Read the general statistics
         fread(fid,2,'int32');
Um     = fread(fid,nyu,'float64')./nstatu;
         fread(fid,2,'int32');
U2m  = fread(fid,nyu,'float64')./nstatu;
         fread(fid,2,'int32');
Vm   = fread(fid,nyv,'float64')./nstatv;
         fread(fid,2,'int32');
V2m  = fread(fid,nyv,'float64')./nstatv;
         fread(fid,2,'int32');
Wm   = fread(fid,nyu,'float64')./nstatu;
         fread(fid,2,'int32');
W2m  = fread(fid,nyu,'float64')./nstatu;
         fread(fid,2,'int32');
Pm   = fread(fid,nyp,'float64')./nstatp;
         fread(fid,2,'int32');
P2m  = fread(fid,nyp,'float64')./nstatp;
         fread(fid,2,'int32');
UVm  = fread(fid,nyu,'float64')./nstatu;
       fread(fid,2,'int32');
UWm  = fread(fid,nyu,'float64')./nstatu;
       fread(fid,2,'int32');
VWm  = fread(fid,nyv,'float64')./nstatv;
       fread(fid,2,'int32');
wxm  = fread(fid,nyv,'float64')./nstatv;
       fread(fid,2,'int32');
wx2m = fread(fid,nyv,'float64')./nstatv;
       fread(fid,1,'int32');


dnx  = Ngal(1,2)/tilex*2/3;
dnz  = Ngal(2,2)/tilez*2/3;
dnyu = Ngal(4,2)-Ngal(4,1)+1;
dnyv = Ngal(3,2)-Ngal(3,1)+1;
dnyp  = Ngal(4,2)-1-(Ngal(4,1)+1)+1;

nstatuC = istat*tilex*tilez;%*2;
nstatvC = istat*tilex*tilez;%*2;
nstatpC = istat*tilex*tilez;%*2;

 fread(fid,1,'int32');

 

UmC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
UmC = reshape(UmC,dnx,dnz,dnyu);

figure(97)
surf(UmC(:,:,1))
shading flat

fread(fid,2,'int32');

U2mC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
U2mC = reshape(U2mC,dnx,dnz,dnyu);

fread(fid,2,'int32');

VmC = fread(fid,dnx*dnz*dnyv,'float64')/nstatvC;
VmC = reshape(VmC,dnx,dnz,dnyv);

figure(98)
surf(VmC(:,:,1))
shading flat

fread(fid,2,'int32');

V2mC = fread(fid,dnx*dnz*dnyv,'float64')/nstatvC;
V2mC = reshape(V2mC,dnx,dnz,dnyv);

fread(fid,2,'int32');

WmC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
WmC = reshape(WmC,dnx,dnz,dnyu);

figure(99)
surf(WmC(:,:,1))
shading flat

fread(fid,2,'int32');

W2mC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
W2mC = reshape(W2mC,dnx,dnz,dnyu);

fread(fid,2,'int32');

PmC = fread(fid,dnx*dnz*dnyp,'float64')/nstatpC;
PmC = reshape(PmC,dnx,dnz,dnyp);

fread(fid,2,'int32');

P2mC = fread(fid,dnx*dnz*dnyp,'float64')/nstatpC;
P2mC = reshape(P2mC,dnx,dnz,dnyp);

fread(fid,2,'int32');

UVmC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
UVmC = reshape(UVmC,dnx,dnz,dnyu);

fread(fid,2,'int32');

UWmC = fread(fid,dnx*dnz*dnyu,'float64')/nstatuC;
UWmC = reshape(UWmC,dnx,dnz,dnyu);

fread(fid,2,'int32');

VWmC = fread(fid,dnx*dnz*dnyv,'float64')/nstatvC;
VWmC = reshape(VWmC,dnx,dnz,dnyv);

fread(fid,2,'int32');

wxmC = fread(fid,dnx*dnz*dnyv,'float64')/nstatvC;
wxmC = reshape(wxmC,dnx,dnz,dnyv);

fread(fid,2,'int32');

wx2mC = fread(fid,dnx*dnz*dnyv,'float64')/nstatvC;
wx2mC = reshape(wx2mC,dnx,dnz,dnyv);

figure(100)
surf(PmC(:,:,1))
shading flat

fread(fid,1,'int32');
       
ribh = 0;
fclose(fid);

% dudy         = .5*dthetai*[dthdy(2:end); dthdy(end)].*[2*(Um(2)-Um(1));(Um(3:nyu)-Um(1:nyu-2));2*(Um(nyu)-Um(nyu-1))];
% dudy         = .5*dthetai*[dthdy(2:end); dthdy(end)].*[(-Um(3)+4*Um(2)-3*Um(1));(Um(3:ny)-Um(1:ny-2));(3*Um(ny)-4*Um(ny-1)+Um(ny-2))];
%dudy         = .5*dthetai*[dthdyu(2:end); dthdyu(end)].*[(-Um(3)+4*Um(2)-3*Um(1));(Um(3:nyu)-Um(1:nyu-2));(3*Um(nyu)-4*Um(nyu-1)+Um(nyu-2))];

for j=1:length(yv)
dudy(j) = dthetai*[dthdyv(j)].*[(Um(j+1)-Um(j))];
dthdyv2(j) = 1/(yu(j+1)-yu(j));
end
dudy=dudy';
for j=1:length(yv)
UVmv(j) = ((yu(j+1)-yv(j))*UVm(j)+(yv(j)-yu(j))*UVm(j+1))/(yu(j+1)-yu(j));
end
%UVmv(1)=0;
%UVmv(length(yv))=0;
UVmv=UVmv';

tot_stress   = -UVmv+dudy/Re;
indlim       = find(abs(yv)<=yinterp);
ts           = tot_stress(indlim);
ts           = .5*(ts(end:-1:1)-ts(1:end));
p1           = polyfit(yv(indlim),ts,1);
utau(jribn)  = sqrt(polyval(p1,1+.75*ribh))
Retau(jribn) = Re*utau(jribn)
yp           = Retau(jribn)*(yv+1);

%Retau=180
%utau=Retau/3250

%colours;
    colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];

%%
% figure(1)
% plot(y,tot_stress);
figure(2), clf, hold on;
plot(yp,tot_stress/utau(jribn).^2,'Color',[38  139 210]/256); %Blue
plot(yp,-UVmv/utau(jribn).^2,'Color',[133 153   0]/256);       %Green
plot(yp,dudy/Re/utau(jribn).^2,'Color',[220 50   47]/256);    %Red
legend('Total Stress','Re Stress','dudy');

load Re180.mat

figure(1)
clf
hold on
plot(yu(2:end-1),Um(2:end-1)/utau,'ko')%'Color',[38  139 210]/256,'LineWidth',2)%,'Linestyle','--')
%title('Mean velocity profile')
xlabel('y')
ylabel('U+')

figure(3)
clf
hold on 
load rgm_uv
plot(SC180_y-1,SC180_uv,'r-',1-SC180_y,-SC180_uv,'Color',[220 50   47]/256)
plot(yu,UVm/utau(jribn).^2,'Color',[38  139 210]/256)
title('u''^+v''^+ (vs torroja)')

figure(4)
clf
load rgm_v
plot(SC180_y-1,SC180_vp,'r-',1-SC180_y,SC180_vp,'Color',[220 50   47]/256)
hold on 
plot(yv,sqrt(V2m-Vm.^2)/utau(jribn),'Color',[38  139 210]/256)
title('v''^+ (vs torroja)')

figure(5)
clf
load rgm_w
plot(SC180_y-1,SC180_wp,'r-',1-SC180_y,SC180_wp,'Color',[220 50   47]/256)
hold on 
plot(yu,sqrt(W2m-Wm.^2)/utau(jribn),'Color',[38  139 210]/256)
title('w''^+ (vs torroja)')

figure(6)
clf
load rgm_u
plot(SC180_y-1,SC180_up,'r-',1-SC180_y,SC180_up,'Color',[220 50   47]/256)
hold on 
plot(yu,sqrt(U2m-Um.^2)/utau(jribn),'Color',[38  139 210]/256)
title('u''^+ (vs torroja)')

figure(7)
clf
hold on 
plot(ypp,sqrt(P2m-Pm.^2)/utau(jribn)^2)
title('p''^+')


figure(28)
clf
hold on 
plot(yv,sqrt(wx2m-wxm.^2)/utau/Retau)
title('wx''')

figure(8)
clf
%semilogx(JellyMVP(2,:),JellyMVP(1,:),'Color',[133 153   0]/256,'LineWidth',1,'Linestyle','--')
title('Mean velocity profile')
hold on
semilogx(yp(1:(length(yp)+1)/2),Um(1:(length(yp)+1)/2)/utau(jribn),'ko')
xlim([1 200])
%semilogx(yp(1:(length(yp)+1)/2),(Um(1:(length(yp)+1)/2)-Um(2))/utau(jribn),'r-')
%semilogx(yp(40:(length(yp)+1)/2-15),Um(40:(length(yp)+1)/2-15)/utau(jribn),'g-')
%semilogx([1:1:180],log([1:1:180])/0.41+5.2,'r-')
%semilogx([1:1:180],log([1:1:180])*2.5+5.5,'m-')


xlabel('y+')
ylabel('U+')


figure(9)
clf
plot(SC180_y*Retau(jribn),SC180_up,'--','LineWidth',1.5,'Color',colors(5,:))
hold on
plot(yu(1:end-1)*utau*3250+Retau(jribn),sqrt(U2m(1:end-1)-Um(1:end-1).^2)/utau(jribn),'-','LineWidth',1.5,'Color',colors(1,:))

title('u''^+ y+ (nu=1/3250)')


 for j = 1:dnyu
     uCondrms(j,1) = rms(reshape(UmC(:,:,j),dnx*dnz,1)-Um(j));
     uCondrms(j,1) = rms(reshape(UmC(:,:,j)-Um(j),dnx*dnz,1));
 end
  for j = 1:dnyv
     vCondrms(j,1) = rms(reshape(VmC(:,:,j),dnx*dnz,1)-Vm(j));
     vCondrms(j,1) = rms(reshape(VmC(:,:,j)-Vm(j),dnx*dnz,1));
  end
  for j = 1:dnyu
     wCondrms(j,1) = rms(reshape(WmC(:,:,j),dnx*dnz,1)-Wm(j));
     wCondrms(j,1) = rms(reshape(WmC(:,:,j)-Wm(j),dnx*dnz,1));
 end


upp = sqrt(U2m(1:dnyu)-Um(1:dnyu).^2)-uCondrms;

plot(yu(1:dnyu)*Retau+Retau(jribn),upp(1:dnyu)/utau(jribn),'-','LineWidth',1.5,'Color',colors(3,:))
plot(yu(1:dnyu)*Retau+Retau(jribn),uCondrms/utau(jribn),'-','LineWidth',1.5,'Color',colors(6,:))


semilogx(yu(1:dnyu)*Retau+Retau(jribn),sqrt((sqrt(U2m(1:dnyu)-Um(1:dnyu).^2)/utau(jribn)).^2-(uCondrms/utau(jribn)).^2),'-','LineWidth',1.5,'Color',colors(7,:))


axis([0 100 0 6])

figure(10)
clf
hold on
plot(SC180_y*Retau(jribn),SC180_vp,'--','LineWidth',1.5,'Color',colors(5,:))
plot(yv(1:end-1)*utau*3250+Retau(jribn),sqrt(V2m(1:end-1)-Vm(1:end-1).^2)/utau(jribn),'-','LineWidth',1.5,'Color',colors(1,:))
axis([0 100 0 1])

figure(11)
clf
hold on
plot(SC180_y*Retau(jribn),SC180_wp,'--','LineWidth',1.5,'Color',colors(5,:))
plot(yu(1:end-1)*utau*3250+Retau(jribn),sqrt(W2m(1:end-1)-Wm(1:end-1).^2)/utau(jribn),'-','LineWidth',1.5,'Color',colors(1,:))


plot(ypp*utau*3250+Retau(jribn),sqrt(P2m-Pm.^2)/utau(jribn)^2,'-','LineWidth',1.5,'Color',colors(4,:))

axis([0 100 0 2])

figure(12)
clf
hold on
axis([0 100 -1 3])
semilogx(yu(1:dnyu)*Retau+Retau(jribn),sqrt((sqrt(U2m(1:dnyu)-Um(1:dnyu).^2)/utau(jribn)).^2-(uCondrms/utau(jribn)).^2),'ko')
semilogx(yu(dnyu+1:end)*Retau+Retau(jribn),sqrt((sqrt(U2m(dnyu+1:end)-Um(dnyu+1:end).^2)/utau(jribn)).^2),'ko')

semilogx(yv(1:dnyv)*Retau+Retau(jribn),sqrt((sqrt(V2m(1:dnyv)-Vm(1:dnyv).^2)/utau(jribn)).^2-(vCondrms/utau(jribn)).^2),'ks')
semilogx(yv(dnyv+1:end)*Retau+Retau(jribn),sqrt((sqrt(V2m(dnyv+1:end)-Vm(dnyv+1:end).^2)/utau(jribn)).^2),'ks')

semilogx(yu(1:dnyu)*Retau+Retau(jribn),sqrt((sqrt(W2m(1:dnyu)-Wm(1:dnyu).^2)/utau(jribn)).^2-(wCondrms/utau(jribn)).^2),'kd')
semilogx(yu(dnyu+1:end)*Retau+Retau(jribn),sqrt((sqrt(W2m(dnyu+1:end)-Wm(dnyu+1:end).^2)/utau(jribn)).^2),'kd')

xlabel('y+')
ylabel('u''+')

%plot((Martell09_R12(:,1)+1)*Retau,-sqrt(abs(Martell09_R12(:,2))),'k-')

% semilogx([1:1:180],log([1:1:180])/0.41+5.2+(Um(15)+Um(16))/utau(jribn),'k-')
% semilogx([1:1:180],log([1:1:180])*2.50+5.5+(Um(15)+Um(16))/utau(jribn),'c-')
% bfy=log(yp(40:(length(yp)+1)/2-15));
% bfu=Um(40:(length(yp)+1)/2-15)/utau(jribn);
% Pfit=polyfit(bfy,bfu,1);
% kappa=1/Pfit(1)
% B=Pfit(2)
% 
% Us=(Um(15)+Um(16))/2/utau(jribn)
% 
% fitcoeff=25;
% Pfit=polyfit((yu((size(Um)-152)/2+1:size(Um)-(size(Um)-152)/2)),(Um((size(Um)-152)/2+1:size(Um)-(size(Um)-152)/2)),fitcoeff);
% Umfit=polyval(Pfit,yv((size(Um)-152)/2:151+(size(Um)-152)/2+1));
% figure(1)
% plot(Umfit/utau,yv((size(Um)-152)/2:151+(size(Um)-152)/2+1),'r')
% 
% max(Um)
