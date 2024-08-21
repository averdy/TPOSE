% TPOSE salinity budget
% nonlinear free surface, r*
% read diagnostics for each 4-month state estimate 2014-2018
% save budget terms


%----------------------------------------------------------------
% PART 1: read overlapping state estimates

addpath ~/scripts_m/
addpath ~/scripts_m/m_map/

% skip maps if savexy=0
savexy=1;

% reference salinity
Sref= 34.96;
rho = 1029;

% select area
x = 109:563; % 140E to eastern boundary
y = 64:93;   % 5S to 5N
z = 1:22;    % top 300 m

% model grid
load /home/averdy/tpose/grid_3/grid 
XC = XC(x,y);
YC = YC(x,y);
RC = RC(z);
DRF = DRF(z);
Depth = Depth(x,y);
hFacC = hFacC(x,y,z);
[nx,ny,nz] = size(hFacC);
dz = permute(repmat(DRF,[1,nx,ny]),[2,3,1]).*hFacC;

% cell volume, face areas (for flux calculations)
volume = zeros(nx,ny,nz);
areaWest = zeros(nx+1,ny,nz);
areaSouth = zeros(nx,ny+1,nz);
areaWest_noh = zeros(nx+1,ny,nz);
areaSouth_noh = zeros(nx,ny+1,nz);
areaTop = zeros(nx,ny,nz+1);
for k=1:nz
	volume(:,:,k) = hFacC(:,:,k).*RAC(x,y)*DRF(k);
	areaTop(:,:,k) = RAC(x,y);
	if x(end)==1080
		areaWest(:,:,k) = DYG([x 1],y).*DRF(k).*hFacW([x 1],y,k);
		areaWest_noh(:,:,k) = DYG([x 1],y).*DRF(k);
	else
		areaWest(:,:,k) = DYG([x x(end)+1],y).*DRF(k).*hFacW([x x(end)+1],y,k);
		areaWest_noh(:,:,k) = DYG([x x(end)+1],y).*DRF(k);
	end
	areaSouth(:,:,k) = DXG(x,[y y(end)+1]).*DRF(k).*hFacS(x,[y y(end)+1],k);
	areaSouth_noh(:,:,k) = DXG(x,[y y(end)+1]).*DRF(k);
end
areaTop(:,:,nz+1) = RAC(x,y);
area = RAC(x,y);


mon={'jan','mar','may','jul','sep','nov'};
cnt1=0; cnt2=59;

for yr=2014:2018

	for m=1:6
 
		clear adv adv_x adv_y adv_z div_x div_y div_z dif dif_h dif_v kpp airsea surf tend corr res S *_xy adv_N adv_S adv_W adv_E div_N div_S div_W div_E adv_b adv_t adv_140E adv_5N adv_5S adv_300m 
 
		% path to diagnostics
		cd(['/data/SO6/TPOSE_diags/tpose3/' char(mon(m)) num2str(yr) '/last_iter_diags/']);

		diag_budget_file = 'diag_salt_budget';

		% time stepping: 1 day
		nt = length(dir('diag_salt_budget*data'));
		dt = 48;
		ts = dt:dt:(nt*dt);

		% read diagnostics
		% calculate tendencies 
		for t=1:nt

			% SALT field 
			tmp = rdmds('diag_state',ts(t),'rec',2);
			tracer = tmp(x,y,z);
			S(t) = sum(tracer(:).*volume(:))./sum(volume(:));

			if z(1)==1
				% tendency due to E-P
				% diagnostic: oceFWflx
				% kg/m^2/s net surface Fresh-Water flux into the ocean (+=down), >0 decreases salinity
				flux = rdmds('diag_surf',ts(t),'rec',4);
				% divide by thickness of first layer 
				% (hFacC accounts for partial cells)
				tmp = flux(x,y)./(rho*dz(:,:,1));
				tmp(isnan(tmp))=0;
				%tmp=-tracer(:,:,1).*tmp;
				tmp=-Sref*tmp;
				airsea(t) = sum(sum(tmp.*volume(:,:,1)));
				%airsea(t) = sum(sum(tmp.*volume(:,:,1)./rstarfac(:,:,1)));
				if savexy
					airsea_xy(:,:,t) = tmp.*volume(:,:,1);
				end
				else
					airsea(t) = 0;
			end

			% advection
			% diagnostics: ADVx_TH, ADVy_TH, ADVr_TH
			% advective flux 
			advflux = rdmds(diag_budget_file,ts(t),'rec',1:3);

			advflux_x = advflux([x x(end)+1],y,z,1);
			advflux_y = advflux(x,[y y(end)+1],z,2);
			advflux_z = advflux(x,y,[z z(end)+1],3);
			% tmp = -(adv_x + adv_y + adv_z);
			% minus sign because advective flux is on lhs, moving to rhs
			tmpx = -diff(advflux_x,1,1);
			adv_x(t) = sum(tmpx(:));
			tmpy = -diff(advflux_y,1,2);
			adv_y(t) = sum(tmpy(:));
			tmpz = diff(advflux_z,1,3);
			adv_z(t) = sum(tmpz(:));

			% divergence
			% diagnostics: UVEL, VVEL, WVEL
			vel = rdmds('diag_state_mass',ts(t),'rec',1:3);
			if x(end)==564
				U = vel([x 1],y,z,1).*areaWest_noh;
			else
				U = vel([x x(end)+1],y,z,1).*areaWest_noh;
			end
			V = vel(x,[y y(end)+1],z,2).*areaSouth_noh;
			W = vel(x,y,[z z(end)+1],3).*areaTop; 

			% advection components
			tmpx = tracer.*diff(U,1,1);
			div_x(t) = sum(tmpx(:));
			tmpy = tracer.*diff(V,1,2);
			div_y(t) = sum(tmpy(:));

			W1=W; % needed for volume budget
			if z(1)==1
				W1(:,:,1) = 0;
			end
			tmpz = tracer.*(-diff(W1,1,3));
			div_z(t) = sum(tmpz(:));

			% lateral advection
			% u.(T-Tref), v.(T-Tref)
			adv_W(t) = ( sum(sum(advflux_x(1,:,:)))  )-( sum(sum(U(1,:,:)*Sref)) );
			adv_E(t) = ( sum(sum(advflux_x(end,:,:))) )-( sum(sum(U(end,:,:)*Sref)) );
			adv_S(t) = ( sum(sum(advflux_y(:,1,:))) )-( sum(sum(V(:,1,:)*Sref)) );
			adv_N(t) = ( sum(sum(advflux_y(:,end,:))) )-( sum(sum(V(:,end,:)*Sref)) );;

			adv_140E(:,:,t) = squeeze(advflux_x(1,:,:)-Sref*(U(1,:,:))); % 140E
			adv_5S(:,:,t) = squeeze(advflux_y(:,1,:)-Sref*(V(:,1,:))); % 5S
			adv_5N(:,:,t) = squeeze(advflux_y(:,end,:)-Sref*(V(:,end,:))); % 5N

			% vertical advection
			adv_t(t) = ( sum(sum(advflux_z(:,:,1))) )-( sum(sum(W1(:,:,1)*Sref)) );
			adv_b(t) = ( sum(sum(advflux_z(:,:,end))) )-( sum(sum(W1(:,:,end)*Sref)) );
			adv_300m(:,:,t) = squeeze(advflux_z(:,:,end)-Sref*(W1(:,:,end)));

			% mixing
			% diagnostics: DFxE_TH, DFyE_TH, DFrI_TH
			% diffusive flux 
			diffflux = rdmds(diag_budget_file,ts(t),'rec',4:7);
			if x(end)==1080
				diffflux_x = diffflux([x 1],y,z,1);
			else
				diffflux_x = diffflux([x x(end)+1],y,z,1);
			end
			diffflux_y = diffflux(x,[y y(end)+1],z,2);
			diffflux_z = diffflux(x,y,[z z(end)+1],3);
			dif_x = diff(diffflux_x,1,1);
			dif_y = diff(diffflux_y,1,2);
			dif_z = -diff(diffflux_z,1,3);

			% minus sign because diffusive flux is on lhs, moving to rhs
			tmp = -(dif_x + dif_y + dif_z);
			dif(t) = sum(tmp(:));
			if savexy
				mix_xy(:,:,t) = sum(tmp,3);
			end
			tmp = -(dif_x + dif_y);
			dif_h(t) = sum(tmp(:));
			tmp = -(dif_z);
			dif_v(t) = sum(tmp(:));

			% KPP 
			tmp = rdmds(diag_budget_file,ts(t),'rec',9);
			tmp = tmp(x,y,[z z(end)+1]);
			tmp = diff(tmp,1,3);
			kpp(t) = sum(tmp(:));
			if savexy
				mix_xy(:,:,t) = mix_xy(:,:,t)+sum(tmp,3);
			end

			% total tendency, diag
			% note that this includes the multiplication by rstarfac 
			tmp = rdmds(diag_budget_file,ts(t),'rec',10);
			tmp = tmp(x,y,z)/86400;
			tend(t) = sum(tmp(:).*volume(:));
			if savexy
				tend_xy(:,:,t) = sum(tmp.*volume,3);
			end

		end

		mix = dif+kpp;
		surf = airsea;
		adv_x = adv_x+div_x;
		adv_y = adv_y+div_y;
		adv_z = adv_z+div_z;
		res = adv_x+adv_y+adv_z+mix+surf-tend;

		if m==1 | m==3 | m==5
			 airsea_xy1(:,:,cnt1+(1:nt)) = airsea_xy;
			 tend_xy1(:,:,cnt1+(1:nt)) = tend_xy;
			 adv_140E1(:,:,cnt1+(1:nt)) = adv_140E;
			 adv_5S1(:,:,cnt1+(1:nt)) = adv_5S;
			 adv_5N1(:,:,cnt1+(1:nt)) = adv_5N;
			 adv_300m1(:,:,cnt1+(1:nt)) = adv_300m;
			 adv_x1(cnt1+(1:nt)) = adv_x;
			 adv_y1(cnt1+(1:nt)) = adv_y;
			 adv_z1(cnt1+(1:nt)) = adv_z;
			 mix1(cnt1+(1:nt)) = mix;
			 surf1(cnt1+(1:nt)) = surf;
			 tend1(cnt1+(1:nt)) = tend;
			 res1(cnt1+(1:nt)) = res;
			 S1(cnt1+(1:nt)) = S;
			 adv_W1(cnt1+(1:nt)) = adv_W;
			 adv_S1(cnt1+(1:nt)) = adv_S;
			 adv_N1(cnt1+(1:nt)) = adv_N;
			 adv_t1(cnt1+(1:nt)) = adv_t;
			 adv_b1(cnt1+(1:nt)) = adv_b;
			 cnt1=cnt1+nt;
		else
			 airsea_xy2(:,:,cnt2+(1:nt)) = airsea_xy;
			 tend_xy2(:,:,cnt2+(1:nt)) = tend_xy;
			 adv_140E2(:,:,cnt2+(1:nt)) = adv_140E;
			 adv_5S2(:,:,cnt2+(1:nt)) = adv_5S;
			 adv_5N2(:,:,cnt2+(1:nt)) = adv_5N;
			 adv_300m2(:,:,cnt2+(1:nt)) = adv_300m;
			 adv_x2(cnt2+(1:nt)) = adv_x;
			 adv_y2(cnt2+(1:nt)) = adv_y;
			 adv_z2(cnt2+(1:nt)) = adv_z;
			 mix2(cnt2+(1:nt)) = mix;
			 surf2(cnt2+(1:nt)) = surf;
			 tend2(cnt2+(1:nt)) = tend;
			 res2(cnt2+(1:nt)) = res;
			 S2(cnt2+(1:nt)) = S;
			 adv_W2(cnt2+(1:nt)) = adv_W;
			 adv_S2(cnt2+(1:nt)) = adv_S;
			 adv_N2(cnt2+(1:nt)) = adv_N;
			 adv_t2(cnt2+(1:nt)) = adv_t;
			 adv_b2(cnt2+(1:nt)) = adv_b;
			 cnt2=cnt2+nt;
		end

	end % for m

end % for yr

surf_xy1=airsea_xy1;
surf_xy2=airsea_xy2;


%----------------------------------------------------------------
% PART 2: patch overlapping state estimates

% weights for patching 

% Spatial patterns
nx=455; ny=30; nz=22;

time = datenum([2014 1 1]):datenum([2019 2 28]);
nt = length(time);
nt1 = 1827; %size(tend_xy1,3);

% overlap: 
dt = 61;
dt2 = (dt-1)/2;

clear ind
begmonth = datenum(2014, 4:2:200, 1);
begmonth(begmonth>time(end)) = [];
for i=1:length(begmonth)
	ind(i) = find(time==begmonth(i));
end
ind = [0 ind];

clear weight*
weight1(1:ind(2)-dt2-1) = 1;
weight2(1:ind(2)-dt2-1) = 0;

for i=2:2:length(ind)
	weight1(ind(i-1)+1:ind(i)-dt2-1) = 1;
	weight2(ind(i-1)+1:ind(i)-dt2-1) = 0;
	weight1(ind(i)-dt2:ind(i)+dt2) = linspace(1,0,dt);
	weight2(ind(i)-dt2:ind(i)+dt2) = linspace(0,1,dt);
end

for i=3:2:length(begmonth)
	weight2(ind(i-1)+dt2:ind(i)-dt2-1) = 1;
	weight1(ind(i-1)+dt2:ind(i)-dt2-1) = 0;
	weight2(ind(i)-dt2:ind(i)+dt2) = linspace(1,0,dt);
	weight1(ind(i)-dt2:ind(i)+dt2) = linspace(0,1,dt);
end

weight2(end+1:nt) = weight2(end);
weight1(end+1:nt) = weight1(end);

weight1=permute(repmat(weight1',[1 nx ny]),[2 3 1]);
weight2=permute(repmat(weight2',[1 nx ny]),[2 3 1]);


% initialize

surf_xy1(:,:,nt1+1:nt) = 0;
surf_xy2(:,:,1:59) = 0;
adv_140E1(:,:,nt1+1:nt) = 0;
adv_140E2(:,:,1:59) = 0;
adv_5N1(:,:,nt1+1:nt) = 0;
adv_5N2(:,:,1:59) = 0;
adv_5S1(:,:,nt1+1:nt) = 0;
adv_5S2(:,:,1:59) = 0;
adv_300m1(:,:,nt1+1:nt) = 0;
adv_300m2(:,:,1:59) = 0;


% average

surf_xy = weight1.*surf_xy1 + weight2.*surf_xy2;
adv_300m = weight1.*adv_300m1 + weight2.*adv_300m2;

weight1=weight1(:,1:nz,:);
weight2=weight2(:,1:nz,:);
adv_5N = weight1.*adv_5N1 + weight2.*adv_5N2;
adv_5S = weight1.*adv_5S1 + weight2.*adv_5S2;

weight1=weight1(1:ny,:,:);
weight2=weight2(1:ny,:,:);
adv_140E = weight1.*adv_140E1 + weight2.*adv_140E2;



% scale spatial patterns: divide by area

% select area
x = 109:563;  
y = 64:93; 
z = 1:22; 

% model grid
load /home/averdy/tpose/grid_3/grid 
XC = XC(x,y);
YC = YC(x,y);
RC = RC(z);
DRF = DRF(z);
hFacC = hFacC(x,y,z);
[nx,ny,nz] = size(hFacC);
dz = permute(repmat(DRF,[1,nx,ny]),[2,3,1]).*hFacC;

% cell volume, face areas (for flux calculations)
volume = zeros(nx,ny,nz);
areaWest = zeros(nx+1,ny,nz);
areaSouth = zeros(nx,ny+1,nz);
areaTop = zeros(nx,ny,nz+1);
for k=1:nz
	volume(:,:,k) = hFacC(:,:,k).*RAC(x,y)*DRF(k);
	areaTop(:,:,k) = RAC(x,y);
	if x(end)==1080
		areaWest(:,:,k)  = DYG([x 1],y).*DRF(k).*hFacW([x 1],y,k);
	else
		areaWest(:,:,k)  = DYG([x x(end)+1],y).*DRF(k).*hFacW([x x(end)+1],y,k);
	end
	areaSouth(:,:,k) = DXG(x,[y y(end)+1]).*DRF(k).*hFacS(x,[y y(end)+1],k);
end
areaTop(:,:,nz+1) = RAC(x,y);
area = RAC(x,y);

areaWest140E = squeeze(areaWest(1,:,:));
areaSouth5N = squeeze(areaSouth(:,end,:));
areaSouth5S = squeeze(areaSouth(:,1,:));

adv_140E = adv_140E./repmat(areaWest140E,[1 1 nt]);
adv_5N = adv_5N./repmat(areaSouth5N,[1 1 nt]);
adv_5S = adv_5S./repmat(areaSouth5S,[1 1 nt]);

surf_xy = surf_xy./repmat(areaTop(:,:,1),[1 1 nt]);
adv_300m = adv_300m./repmat(areaTop(:,:,22),[1 1 nt]);



nz=22;

% mean seas cycle
tmp=surf_xy; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,ny,365,5),4);
surf_seas = tmpS;
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
surf_anom = tmp-tmpS;

tmp=adv_300m; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,ny,365,5),4);
adv_300m_seas = tmpS;
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
adv_300m_anom = tmp-tmpS;

tmp=adv_5N; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
adv_5N_seas = tmpS;
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
adv_5N_anom = tmp-tmpS;

tmp=adv_5S; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
adv_5S_seas = tmpS;
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
adv_5S_anom = tmp-tmpS;

tmp=adv_140E; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),ny,nz,365,5),4);
adv_140E_seas = tmpS;
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
adv_140E_anom = tmp-tmpS;


% remove leap day
time1=time; time1(365*2+31+29)=[]; 

% monthly means
ndays=repmat([31 28 31 30 31 30 31 31 30 31 30 31],[1 5]);

ind=1;
for n=1:length(ndays)
	surf_monthly(:,:,n) = mean(surf_xy(:,:,ind:ind+ndays(n)),3);
	adv_300m_monthly(:,:,n) = mean(adv_300m(:,:,ind:ind+ndays(n)),3);
	adv_140E_monthly(:,:,n) = mean(adv_140E(:,:,ind:ind+ndays(n)),3);
	adv_5N_monthly(:,:,n) = mean(adv_5N(:,:,ind:ind+ndays(n)),3);
	adv_5S_monthly(:,:,n) = mean(adv_5S(:,:,ind:ind+ndays(n)),3);
	surf_anom_monthly(:,:,n) = mean(surf_anom(:,:,ind:ind+ndays(n)),3);
	adv_300m_anom_monthly(:,:,n) = mean(adv_300m_anom(:,:,ind:ind+ndays(n)),3);
	adv_140E_anom_monthly(:,:,n) = mean(adv_140E_anom(:,:,ind:ind+ndays(n)),3);
	adv_5N_anom_monthly(:,:,n) = mean(adv_5N_anom(:,:,ind:ind+ndays(n)),3);
	adv_5S_anom_monthly(:,:,n) = mean(adv_5S_anom(:,:,ind:ind+ndays(n)),3);
	ind=ind+ndays(n);
end


cd /data/averdy/tpose/budgets/
save salinity_budget_faces.mat *anom *monthly *xy *seas surf_xy *300m *140E *5N *5S


% time series

time = datenum([2014 1 1]):datenum([2019 2 28]);
nt = length(time);
nt1 = length(S1);

% overlap: 
dt = 61;
dt2 = (dt-1)/2;

clear ind
begmonth = datenum(2014, 4:2:200, 1);
begmonth(begmonth>time(end)) = [];
for i=1:length(begmonth)
	ind(i) = find(time==begmonth(i));
end
ind = [0 ind];

clear weight*
weight1(1:ind(2)-dt2-1) = 1;
weight2(1:ind(2)-dt2-1) = 0;

for i=2:2:length(ind)
	weight1(ind(i-1)+1:ind(i)-dt2-1) = 1;
	weight2(ind(i-1)+1:ind(i)-dt2-1) = 0;
	weight1(ind(i)-dt2:ind(i)+dt2) = linspace(1,0,dt);
	weight2(ind(i)-dt2:ind(i)+dt2) = linspace(0,1,dt);
end

for i=3:2:length(begmonth)
	weight2(ind(i-1)+dt2:ind(i)-dt2-1) = 1;
	weight1(ind(i-1)+dt2:ind(i)-dt2-1) = 0;
	weight2(ind(i)-dt2:ind(i)+dt2) = linspace(1,0,dt);
	weight1(ind(i)-dt2:ind(i)+dt2) = linspace(0,1,dt);
end

weight2(end+1:nt) = weight2(end);
weight1(end+1:nt) = weight1(end);


% initialize

S1(nt1+1:nt) = 0;
S2(1:59) = 0;
adv_x1(nt1+1:nt) = 0;
adv_x2(1:59) = 0;
adv_y1(nt1+1:nt) = 0;
adv_y2(1:59) = 0;
adv_z1(nt1+1:nt) = 0;
adv_z2(1:59) = 0;
surf1(nt1+1:nt) = 0;
surf2(1:59) = 0;
tend1(nt1+1:nt) = 0;
tend2(1:59) = 0;
mix1(nt1+1:nt) = 0;
mix2(1:59) = 0;
res1(nt1+1:nt) = 0;
res2(1:59) = 0;

% average

S = weight1.*S1 + weight2.*S2;
adv_x = weight1.*adv_x1 + weight2.*adv_x2;
adv_y = weight1.*adv_y1 + weight2.*adv_y2;
adv_z = weight1.*adv_z1 + weight2.*adv_z2;
mix = weight1.*mix1 + weight2.*mix2;
surf = weight1.*surf1 + weight2.*surf2;
tend = weight1.*tend1 + weight2.*tend2;
res = weight1.*res1 + weight2.*res2;

adv_W1(nt1+1:nt) = 0;
adv_W2(1:59) = 0;
adv_S1(nt1+1:nt) = 0;
adv_S2(1:59) = 0;
adv_N1(nt1+1:nt) = 0;
adv_N2(1:59) = 0;
adv_b1(nt1+1:nt) = 0;
adv_b2(1:59) = 0;
adv_t1(nt1+1:nt) = 0;
adv_t2(1:59) = 0;

adv_W = weight1.*adv_W1 + weight2.*adv_W2;
adv_S = weight1.*adv_S1 + weight2.*adv_S2;
adv_N = weight1.*adv_N1 + weight2.*adv_N2;
adv_t = weight1.*adv_t1 + weight2.*adv_t2;
adv_b = weight1.*adv_b1 + weight2.*adv_b2;

res = adv_x+adv_y+adv_z+mix+surf-tend;


% divide by volume

area(Depth(x,y)==0)=0;
vol = sum(volume(:));

adv_x = adv_x./repmat(vol,[1 nt]);
adv_y = adv_y./repmat(vol,[1 nt]);
adv_z = adv_z./repmat(vol,[1 nt]);
mix = mix./repmat(vol,[1 nt]);
surf = surf./repmat(vol,[1 nt]);
tend = tend./repmat(vol,[1 nt]);
res = res./repmat(vol,[1 nt]);

adv_W = adv_W./repmat(vol,[1 nt]);
adv_S = adv_S./repmat(vol,[1 nt]);
adv_N = adv_N./repmat(vol,[1 nt]);
adv_t = adv_t./repmat(vol,[1 nt]);
adv_b = adv_b./repmat(vol,[1 nt]);

% 30-day running mean
tend_sm = NaN*tend;
adv_x_sm = NaN*adv_x;
adv_y_sm = NaN*adv_y;
adv_z_sm = NaN*adv_z;
surf_sm = NaN*surf;
mix_sm = NaN*mix;
adv_W_sm = NaN*adv_W;
adv_S_sm = NaN*adv_S;
adv_N_sm = NaN*adv_N;
adv_t_sm = NaN*adv_t;
adv_b_sm = NaN*adv_b;
for t=16:nt-15
	 tend_sm(t)=mean(tend(t-15:t+15));
	 adv_x_sm(t)=mean(adv_x(t-15:t+15));
	 adv_y_sm(t)=mean(adv_y(t-15:t+15));
	 adv_z_sm(t)=mean(adv_z(t-15:t+15));
	 surf_sm(t)=mean(surf(t-15:t+15));
	 mix_sm(t)=mean(mix(t-15:t+15));
	 adv_W_sm(t)=mean(adv_W(t-15:t+15));
	 adv_S_sm(t)=mean(adv_S(t-15:t+15));
	 adv_N_sm(t)=mean(adv_N(t-15:t+15));
	 adv_t_sm(t)=mean(adv_t(t-15:t+15));
	 adv_b_sm(t)=mean(adv_b(t-15:t+15));
end
res_sm = tend_sm-adv_x_sm-adv_y_sm-adv_z_sm-surf_sm-mix_sm;


% mean seas cycle
% remove leap day, fill to 5 years
tmp=adv_x; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_x_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_y; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_y_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_z; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_z_anom = tmp-repmat(tmpS',[1 5]);

tmp=surf; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
surf_anom = tmp-repmat(tmpS',[1 5]);

tmp=mix; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
mix_anom = tmp-repmat(tmpS',[1 5]);

tmp=tend; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
tend_anom = tmp-repmat(tmpS',[1 5]);

tmp=res; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
res_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_W; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_W_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_S; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_S_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_N; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_N_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_b; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_b_anom = tmp-repmat(tmpS',[1 5]);

tmp=adv_t; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
adv_t_anom = tmp-repmat(tmpS',[1 5]);


% remove leap day
time1=time; time1(365*2+31+29)=[]; 

nt=length(res_anom);
if length(time1)<nt
	time1(end+1:nt)=NaN;
else
	time1=time1(1:nt);
end


% 30-day running mean
tend_anom_sm = NaN*tend_anom;
adv_x_anom_sm = NaN*adv_x_anom;
adv_y_anom_sm = NaN*adv_y_anom;
adv_z_anom_sm = NaN*adv_z_anom;
surf_anom_sm = NaN*surf_anom;
mix_anom_sm = NaN*mix_anom;
adv_W_anom_sm = NaN*adv_W_anom;
adv_S_anom_sm = NaN*adv_S_anom;
adv_N_anom_sm = NaN*adv_N_anom;
adv_t_anom_sm = NaN*adv_t_anom;
adv_b_anom_sm = NaN*adv_b_anom;
for t=16:1825-15
	 tend_anom_sm(t)=mean(tend_anom(t-15:t+15));
	 adv_x_anom_sm(t)=mean(adv_x_anom(t-15:t+15));
	 adv_y_anom_sm(t)=mean(adv_y_anom(t-15:t+15));
	 adv_z_anom_sm(t)=mean(adv_z_anom(t-15:t+15));
	 surf_anom_sm(t)=mean(surf_anom(t-15:t+15));
	 mix_anom_sm(t)=mean(mix_anom(t-15:t+15));
	 adv_W_anom_sm(t)=mean(adv_W_anom(t-15:t+15));
	 adv_S_anom_sm(t)=mean(adv_S_anom(t-15:t+15));
	 adv_N_anom_sm(t)=mean(adv_N_anom(t-15:t+15));
	 adv_t_anom_sm(t)=mean(adv_t_anom(t-15:t+15));
	 adv_b_anom_sm(t)=mean(adv_b_anom(t-15:t+15));
end
res_anom_sm = tend_anom_sm-adv_x_anom_sm-adv_y_anom_sm-adv_z_anom_sm-surf_anom_sm-mix_anom_sm;

clear *1 *2 *xy* adv_140* adv_5* adv_300*
clear surf_seas surf_monthly surf_anom_monthly

cd /data/averdy/tpose/budgets
save salinity_budgetbox *sm *anom adv_* surf tend mix res


