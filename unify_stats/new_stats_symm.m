% ATTENTION ribh hardcoded
% ATTENTION dudy has to be calculated using the same scheme as in the code
clear all

fname  = 'stats001x001.dat';

yinterp = .70;
ptstilex = 128;
ptstilez = 128;
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
nstatv = 0*yv+istat;
nstatp = 0*ypp+istat;
for iband = 1:1
  nstatu(nyyu(iband)+1:nyyu(iband+1)) =  nstatu(nyyu(iband)+1:nyyu(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
  nstatv(nyyv(iband)+1:nyyv(iband+1)) =  nstatv(nyyv(iband)+1:nyyv(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
  nstatp(nyyp(iband)+1:nyyp(iband+1)) =  nstatp(nyyp(iband)+1:nyyp(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end
for iband = 2:2
  nstatu(nyyu(iband)+1:nyyu(iband+1)) =  nstatu(nyyu(iband)+1:nyyu(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1);
  nstatv(nyyv(iband)+1:nyyv(iband+1)) =  nstatv(nyyv(iband)+1:nyyv(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1);
  nstatp(nyyp(iband)+1:nyyp(iband+1)) =  nstatp(nyyp(iband)+1:nyyp(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end
for iband = 3:3
  nstatu(nyyu(iband)+1:nyyu(iband+1)) =  nstatu(nyyu(iband)+1:nyyu(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
  nstatv(nyyv(iband)+1:nyyv(iband+1)) =  nstatv(nyyv(iband)+1:nyyv(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
  nstatp(nyyp(iband)+1:nyyp(iband+1)) =  nstatp(nyyp(iband)+1:nyyp(iband+1)) ...
                                   *Ngal(1,iband+1)*Ngal(2,iband+1)*2/3*2/3;
end

nstatvN = 0*yv+istat;
for iband = 1:3
  nstatvN(nyyv(iband)+1:nyyv(iband+1)) =  nstatvN(nyyv(iband)+1:nyyv(iband+1)) ...
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
wxm  = fread(fid,nyv,'float64')./nstatvN;
       fread(fid,2,'int32');
wx2m = fread(fid,nyv,'float64')./nstatvN;
       fread(fid,1,'int32');

% Average over top and bottom walls
for j=1:nyu/2
   Um_s(j)=(Um(j)+Um(end+1-j))/2;
   U2m_s(j)=(U2m(j)+U2m(end+1-j))/2;
   Wm_s(j)=(Wm(j)+Wm(end+1-j))/2;
   W2m_s(j)=(W2m(j)+W2m(end+1-j))/2;
   UVm_s(j)=(UVm(j)-UVm(end+1-j))/2;
   UWm_s(j)=(UWm(j)+UWm(end+1-j))/2;
end
   
for j=1:(nyv-1)/2
   Vm_s(j)=(Vm(j)+Vm(end+1-j))/2;
   V2m_s(j)=(V2m(j)+V2m(end+1-j))/2;
   VWm_s(j)=(VWm(j)+VWm(end+1-j))/2;
   wxm_s(j)=(wxm(j)+wxm(end+1-j))/2;
   wx2m_s(j)=(wx2m(j)+wx2m(end+1-j))/2;
end

for j=1:nyu/2-1
   Pm_s(j)=(Pm(j)+Pm(end+1-j))/2;
   P2m_s(j)=(P2m(j)+P2m(end+1-j))/2;
end

%Conditional statistics
dnx  = Ngal(1,2)/tilex*2/3;
dnz  = Ngal(2,2)/tilez*2/3;
dnyu = Ngal(4,2)-Ngal(4,1)+1;
dnyv = Ngal(3,2)-Ngal(3,1)+1;
dnyp  = Ngal(4,2)-1-(Ngal(4,1)+1)+1;

nstatuC = istat*tilex*tilez*2;
nstatvC = istat*tilex*tilez*2;
nstatpC = istat*tilex*tilez*2;

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


%Retau=197.5
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
plot(yp,tot_stress/utau(jribn).^2,'Color',colors(1,:))
plot(yp,-UVmv/utau(jribn).^2,'Color',colors(2,:))
plot(yp,dudy/Re/utau(jribn).^2,'Color',colors(5,:))
legend('Total Stress','Re Stress','dudy');

load Re180.mat

figure(1)
clf
hold on
plot(yu(2:end-1),Um(2:end-1)/utau,'LineWidth',1.5,'Color',colors(1,:))%'Color',[38  139 210]/256,'LineWidth',2)%,'Linestyle','--')
%title('Mean velocity profile')
xlabel('y')
ylabel('U+')

figure(101)
clf
%load Seo197_P78_MVPmUS
%plot(Seo197_P78_MVPmUS(2,:),Seo197_P78_MVPmUS(1,:),'LineWidth',1.5,'Color',colors(2,:))
set(gca,'XScale','log')
hold on
Uslip = (Um(2)*(-1-yu(1))+Um(1)*(yu(2)+1))/(yu(2)-yu(1));
Usplus = Uslip/utau
plot((yu(2:(end-1)/2)+1)*180,(Um(2:(end-1)/2)-Uslip)/utau,'LineWidth',1.5,'Color',colors(1,:))

figure(3)
clf
hold on 
load rgm_uv
plot(SC180_y-1,SC180_uv,'LineWidth',1.5,'Color',colors(5,:))
plot(yu(1:end/2),UVm_s/utau(jribn).^2,'LineWidth',1.5,'Color',colors(1,:))
title('u''^+v''^+ (vs torroja (green))')

figure(4)
clf
load rgm_v
plot(SC180_y-1,SC180_vp,'LineWidth',1.5,'Color',colors(5,:))
hold on 
%load Seo197_P78_vp
%plot(Seo197_P78_vp(:,1)/197.5-1,Seo197_P78_vp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(yv(1:(end-1)/2),sqrt(V2m_s-Vm_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
title('v''^+ (vs torroja (green) and Seo (red))')

figure(5)
clf
load rgm_w
plot(SC180_y-1,SC180_wp,'LineWidth',1.5,'Color',colors(5,:))
hold on 
%load Seo197_P78_wp
%plot(Seo197_P78_wp(:,1)/197.5-1,Seo197_P78_wp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(yu(1:end/2),sqrt(W2m_s-Wm_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
title('w''^+ (vs torroja (green) and Seo (red))')

figure(6)
clf
load rgm_u
plot(SC180_y-1,SC180_up,'LineWidth',1.5,'Color',colors(5,:))
hold on
%load Seo197_P78_up
%plot(Seo197_P78_up(:,1)/197.5-1,Seo197_P78_up(:,2),'LineWidth',1.5,'Color',colors(2,:)) 
plot(yu(1:end/2),sqrt(U2m_s-Um_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
title('u''^+ (vs torroja (green) and Seo (red))')

figure(7)
clf
hold on 
%load Seo197_P78_pp
%plot(Seo197_P78_pp(:,1)/197.5-1,Seo197_P78_pp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(ypp(1:end/2),sqrt(P2m_s-Pm_s.^2)/utau(jribn)^2,'LineWidth',1.5,'Color',colors(1,:))
title('p''^+ (vs Seo (red))')

figure(28)
clf
hold on 
plot(yv(1:(end-1)/2),sqrt(wx2m_s-wxm_s.^2)/utau/Retau,'LineWidth',1.5,'Color',colors(1,:))
title('wx''')

figure(8)
clf
title('Mean velocity profile')
hold on
semilogx(yp(1:(length(yp)+1)/2),Um(1:(length(yp)+1)/2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
xlim([1 200])
set(gca,'XScale','log')

xlabel('y+')
ylabel('U+')

  for j = 1:dnyu
     %uCondrms(j,1) = rms(reshape(UmC(:,:,j),dnx*dnz,1)-Um_s(j));
     uCondrms(j,1) = rms(reshape(UmC(:,:,j)-Um_s(j),dnx*dnz,1));
  end
  for j = 1:dnyv
     %vCondrms(j,1) = rms(reshape(VmC(:,:,j),dnx*dnz,1)-Vm_s(j));
     vCondrms(j,1) = rms(reshape(VmC(:,:,j)-Vm_s(j),dnx*dnz,1));
  end
  for j = 1:dnyu
     %wCondrms(j,1) = rms(reshape(WmC(:,:,j),dnx*dnz,1)-Wm_s(j));
     wCondrms(j,1) = rms(reshape(WmC(:,:,j)-Wm_s(j),dnx*dnz,1));
  end
  for j = 1:dnyp
     %pCondrms(j,1) = rms(reshape(PmC(:,:,j),dnx*dnz,1)-Pm_s(j));
     pCondrms(j,1) = rms(reshape(PmC(:,:,j)-Pm_s(j),dnx*dnz,1));
  end
%   for j = 1:dnyu
%      %pCondrms(j,1) = rms(reshape(UVmC(:,:,j),dnx*dnz,1)-Pm_s(j));
%      UVCondrms(j,1) = rms(reshape(UVmC(:,:,j)-UVm_s(j),dnx*dnz,1));
%   end
  
figure(9)
clf
plot(SC180_y*Retau(jribn),SC180_up,'LineWidth',1.5,'Color',colors(5,:))
hold on
%plot(Seo197_P78_up(:,1),Seo197_P78_up(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(yu(1:end/2)*utau*3250+Retau(jribn),sqrt(U2m_s-Um_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))

title('u''^+ (vs torroja (green) and Seo (red))')

upp = sqrt(U2m_s(1:dnyu)-Um_s(1:dnyu).^2)-uCondrms';

plot(yu(1:dnyu)*Retau+Retau(jribn),upp(1:dnyu)/utau(jribn),'LineWidth',1.5,'Color',colors(4,:))
plot(yu(1:dnyu)*Retau+Retau(jribn),uCondrms/utau(jribn),'LineWidth',1.5,'Color',colors(6,:))

semilogx(yu(1:dnyu)*Retau+Retau(jribn),sqrt((sqrt(U2m_s(1:dnyu)-Um_s(1:dnyu).^2)/utau(jribn)).^2-(uCondrms/utau(jribn))'.^2),'-','LineWidth',1.5,'Color',colors(7,:))
axis([0 100 0 6])

figure(10)
clf
hold on
plot(SC180_y*Retau(jribn),SC180_vp,'LineWidth',1.5,'Color',colors(5,:))
%plot(Seo197_P78_vp(:,1),Seo197_P78_vp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(yv(1:(end-1)/2)*utau*3250+Retau(jribn),sqrt(V2m_s-Vm_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
title('v''^+ (vs torroja (green) and Seo (red))')

vpp = sqrt(V2m_s(1:dnyv)-Vm_s(1:dnyv).^2)-vCondrms';

plot(yv(1:dnyv)*Retau+Retau(jribn),vpp(1:dnyv)/utau(jribn),'LineWidth',1.5,'Color',colors(4,:))
plot(yv(1:dnyv)*Retau+Retau(jribn),vCondrms/utau(jribn),'LineWidth',1.5,'Color',colors(6,:))

semilogx(yv(1:dnyv)*Retau+Retau(jribn),sqrt((sqrt(V2m_s(1:dnyv)-Vm_s(1:dnyv).^2)/utau(jribn)).^2-(vCondrms/utau(jribn))'.^2),'--','LineWidth',1.5,'Color',colors(7,:))
axis([0 100 0 1.2])

figure(11)
clf
hold on
plot(SC180_y*Retau(jribn),SC180_wp,'LineWidth',1.5,'Color',colors(5,:))
%plot(Seo197_P78_wp(:,1),Seo197_P78_wp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(yu(1:end/2)*utau*3250+Retau(jribn),sqrt(W2m_s-Wm_s.^2)/utau(jribn),'LineWidth',1.5,'Color',colors(1,:))
title('w''^+ (vs torroja (green) and Seo (red))')

wpp = sqrt(W2m_s(1:dnyu)-Wm_s(1:dnyu).^2)-wCondrms';

plot(yu(1:dnyu)*Retau+Retau(jribn),wpp(1:dnyu)/utau(jribn),'LineWidth',1.5,'Color',colors(4,:))
plot(yu(1:dnyu)*Retau+Retau(jribn),wCondrms/utau(jribn),'LineWidth',1.5,'Color',colors(6,:))

semilogx(yu(1:dnyu)*Retau+Retau(jribn),sqrt((sqrt(W2m_s(1:dnyu)-Wm_s(1:dnyu).^2)/utau(jribn)).^2-(wCondrms/utau(jribn))'.^2),'--','LineWidth',1.5,'Color',colors(7,:))
axis([0 100 0 2])

figure(13)
clf
hold on
%plot(Seo197_P78_pp(:,1),Seo197_P78_pp(:,2),'LineWidth',1.5,'Color',colors(2,:))
plot(ypp(1:end/2)*utau*3250+Retau(jribn),sqrt(P2m_s-Pm_s.^2)/utau(jribn)^2,'LineWidth',1.5,'Color',colors(1,:))
title('p''^+ (vs Seo (red))')


ppp = sqrt(P2m_s(1:dnyp)-Pm_s(1:dnyp).^2)-pCondrms';

plot(yu(1:dnyp)*Retau+Retau(jribn),ppp(1:dnyp)/utau(jribn)^2,'LineWidth',1.5,'Color',colors(4,:))
plot(yu(1:dnyp)*Retau+Retau(jribn),pCondrms/utau(jribn)^2,'LineWidth',1.5,'Color',colors(6,:))

semilogx(yu(1:dnyp)*Retau+Retau(jribn),sqrt((sqrt(P2m_s(1:dnyp)-Pm_s(1:dnyp).^2)/utau(jribn)^2).^2-(pCondrms/utau(jribn)^2)'.^2),'-','LineWidth',1.5,'Color',colors(7,:))
axis([0 100 0 7])


for j = 1:dnyu
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j),dnx*dnz,1)-UVm_s(j));
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j)-UVm_s(j),dnx*dnz,1));
   UVCondrms(j,1) = rms(reshape(UVmC(:,:,j),dnx*dnz,1));
end
  
figure(14)
clf
hold on
plot(yu(1:dnyu)*Retau+Retau(jribn),UVCondrms/utau^2,'LineWidth',1.5,'Color',colors(6,:))
plot(yu(1:(end+1)/2)*Retau+Retau(jribn),UVm_s/utau^2,'LineWidth',1.5,'Color',colors(5,:))
%plot(yu(1:(end+1)/2)*Retau+Retau(jribn),(UVm(1:(end+1)/2)-UVm_s')/utau^2,'LineWidth',1.5,'Color',colors(3,:))


figure(105)
surf(UVmC(:,:,10))
shading flat


for j = 1:dnyu
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j),dnx*dnz,1)-UVm_s(j));
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j)-UVm_s(j),dnx*dnz,1));
   UWCondrms(j,1) = rms(reshape(UWmC(:,:,j)-UWm_s(j),dnx*dnz,1));
end

figure(15)
clf
hold on
plot(yu(1:dnyu)*Retau+Retau(jribn),UWCondrms/utau^2,'LineWidth',1.5,'Color',colors(6,:))
plot(yu(1:(end+1)/2)*Retau+Retau(jribn),UWm_s/utau^2,'LineWidth',1.5,'Color',colors(5,:))


figure(106)
surf(UWmC(:,:,10))
shading flat


for j = 1:dnyv
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j),dnx*dnz,1)-UVm_s(j));
   %UVCondrms(j,1) = rms(reshape(UVmC(:,:,j)-UVm_s(j),dnx*dnz,1));
   VWCondrms(j,1) = rms(reshape(VWmC(:,:,j),dnx*dnz,1));
end

figure(16)
clf
hold on
plot(yv(1:dnyv)*Retau+Retau(jribn),VWCondrms/utau^2,'LineWidth',1.5,'Color',colors(6,:))
plot(yv(1:(end)/2)*Retau+Retau(jribn),VWm_s/utau^2,'LineWidth',1.5,'Color',colors(5,:))


figure(107)
surf(VWmC(:,:,10))
shading flat

figure(97)
surf(UmC(:,:,10))
shading flat

figure(98)
surf(VmC(:,:,10))
shading flat

figure(99)
surf(WmC(:,:,10))
shading flat


load sc_stats

figure(200)
clf
plot(sc_yu(1:end/2)+1,-(UVm_s/utau(jribn).^2-sc_UVm_s/sc_utau(jribn).^2),'LineWidth',1.5,'Color',colors(4,:))

du_actual = (Um_s(end)/utau)-(sc_Um_s(end)/sc_utau)
du2 = (Uslip/utau)
du3 = trapz(yu(1:end/2)+1,-(UVm_s/utau(jribn).^2-sc_UVm_s/sc_utau(jribn).^2))*Retau
du = du2-du3




T1_r_sc = (1-(Um_s(end)/sc_Um_s(end)))*((sc_Um_s(end)/sc_utau)/(Um_s(end)/utau))^2
T2_r_sc = -(Uslip/utau)/(sc_Um_s(end)/sc_utau)
T3_r_sc = trapz(yu(1:end/2)+1,-(UVm_s/utau(jribn).^2-sc_UVm_s/sc_utau(jribn).^2))*Retau/(sc_Um_s(end)/sc_utau)

T_r_tot_sc = T1_r_sc+T2_r_sc+T3_r_sc
T_r_23_sc =  T2_r_sc+T3_r_sc
DRcent_sc = (Um_s(end)-sc_Um_s(end))/sc_Um_s(end)


T1_r = (1-(Um_s(end)/sc_Um_s(end)))*((sc_Um_s(end)/sc_utau)/(Um_s(end)/utau))^2
T2_r = -(Uslip/utau)/(Um_s(end)/utau)
T3_r = trapz(yu(1:end/2)+1,-(UVm_s/utau(jribn).^2-sc_UVm_s/sc_utau(jribn).^2))*Retau/(Um_s(end)/utau)

T_r_tot = T1_r+T2_r+T3_r
T_r_23 =  T2_r+T3_r
mDRcent = -((Um_s(end)/utau).^2-(sc_Um_s(end)/sc_utau).^2)/(Um_s(end)/utau).^2




sl47xz_P2m_s = P2m_s;
sl47xz_Pm_s = Pm_s;
sl47xz_Retau = Retau;
sl47xz_U2m_s = U2m_s;
sl47xz_Um = Um;
sl47xz_Uslip = Uslip;
sl47xz_Um_s = Um_s;
sl47xz_utau = utau;
sl47xz_UVm_s = UVm_s;
sl47xz_UWm_s = UWm_s;
sl47xz_V2m_s = V2m_s;
sl47xz_Vm_s = Vm_s;
sl47xz_VWm_s = VWm_s;
sl47xz_W2m_s = W2m_s;
sl47xz_Wm_s = Wm_s;
sl47xz_yu = yu;
sl47xz_yv = yv;

save('../../../Code_170112_Compsl47/new_postpro/unify_stats/sl47xz_stats.mat','sl47xz_P2m_s', 'sl47xz_Pm_s', 'sl47xz_Retau', 'sl47xz_U2m_s', 'sl47xz_Um', 'sl47xz_Uslip', 'sl47xz_Um_s', 'sl47xz_utau', 'sl47xz_UVm_s', 'sl47xz_UWm_s', 'sl47xz_V2m_s', 'sl47xz_Vm_s', 'sl47xz_VWm_s', 'sl47xz_W2m_s', 'sl47xz_Wm_s', 'sl47xz_yu', 'sl47xz_yv')

