% TPOSE volume budget
% nonlinear free surface, r*
% read diagnostics for each 4-month state estimate 2014-2018
% save budget terms


%----------------------------------------------------------------
% PART 1: read overlapping state estimates

addpath ~/scripts_m/
addpath ~/scripts_m/m_map/

% skip maps if savexy=0
savexy=1;

% select area
x = 1:563; % all lons
y = 1:167; % all lats
z = 1:22; % top 300 m

% model grid
load /home/averdy/tpose/grid_3/grid 

% for geostrophy:
rotationPeriod=86164;
FG = 4*pi/rotationPeriod.*sind(YG);
FC = 4*pi/rotationPeriod.*sind(YC);
% equator
FG(:,79)=Inf; 
dycfg = repmat(DYC.*FG,[1 1 length(RC)]);
dxcfc = repmat(DXC.*FC,[1 1 length(RC)]);
[nx,ny,nz] = size(hFacC);
hFacC_forP = hFacC;

% crop grid
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

% constants
rho = 1035;
Cp = 3994;


mon={'jan','mar','may','jul','sep','nov'};
cnt1=0; cnt2=59;

for yr=2014:2018

	for m=1:length(mon)
 
		clear *_xy *_TP5 w_300m u_140E v_5N v_5S
 
		% path to diagnostics
		cd(['/data/SO6/TPOSE_diags/tpose3/' char(mon(m)) num2str(yr) '/last_iter_diags/']);

		% time stepping: 1 day
		nt = length(dir('diag_state_mass*data'));
		dt = 48;
		ts = dt:dt:(nt*dt);

		% read diagnostics
		% calculate tendencies 
		for t=1:nt

			% tendency due to air-sea flux 
			% diagnostic: oceFWflx
			% forcV(:,:,:,cnt)=mygrid.mskC.*mk3D(oceFWflx,mygrid.mskC)./(dzMat*rhoconst);
			oceFWflx = rdmds('diag_surf',ts(t),'rec',4);
			mask=hFacC; mask(mask>0)=1; % to identify bottom cells
			if z(1)==1
				forcV(:,:) = mask(:,:,1).*oceFWflx(x,y)./(rho*dz(:,:,1));
			else
				forcV(:,:) = 0*oceFWflx(x,y) ;
			end 
			forcV(isnan(forcV))=0;
			forcV(isinf(forcV))=0;
			surf_TP5(t) = sum(sum(forcV(109:end,64:93).*volume(109:end,64:93,1)));
			if savexy
				surf_xy(:,:,t) = sum(forcV.*volume,3);
			end

			% divergence
			% diagnostics: UVEL, VVEL, WVEL
			vel = rdmds('diag_state',ts(t),'rec',1:3);
			if savexy
				% velocities, not transport
				w_300m(:,:,t) = squeeze(vel(:,:,23,3));
				u_140E(:,:,t) = squeeze(vel(109,:,:,1)); % 140E
				v_5S(:,:,t) = squeeze(vel(:,64,:,2)); % 5S
				v_5N(:,:,t) = squeeze(vel(:,94,:,2)); % 5N
			end

			% for transport:
			% nonlin free surface, rstar:
			% UVELMASS because hFacC is now time variable
			% already has hFacC
			velmass = rdmds('diag_state_mass',ts(t),'rec',1:3);
			if x(end)==564
				U = velmass([x 1],y,z,1).*areaWest_noh;
			else
				U = velmass([x x(end)+1],y,z,1).*areaWest_noh;
			end
			V = velmass(x,[y y(end)+1],z,2).*areaSouth_noh;
			if z(end)==51
				W = velmass(x,y,z,3).*areaTop(:,:,1:end-1); 
				W(:,:,end+1) = 0*W(:,:,end); 
			else
				W = velmass(x,y,[z z(end)+1],3).*areaTop; 
			end

			% for vol:
			if z(1)==1
				W(:,:,1)=0;
			end

			Vsouth_TP5(t) = sum(sum(sum(V(109:end,64,:))));
			Vnorth_TP5(t) = sum(sum(sum(V(109:end,94,:))));
			Uwest_TP5(t) = sum(sum(sum(U(109,64:93,:))));
			W300_TP5(t) = sum(sum(sum(W(109:end,64:93,end))));

			tmpx = diff(U,1,1);
			tmpy = diff(V,1,2);
			tmpz = (-diff(W,1,3));
			tmph=tmpx+tmpy;
			div_h_TP5(t) = sum(sum(sum(tmph(109:end,64:93,:))));
			div_z_TP5(t) = sum(sum(sum(tmpz(109:end,64:93,:))));
			if savexy
				div_h_xy(:,:,t) = sum(tmpx+tmpy,3);
				div_v_xy(:,:,t) = sum(tmpz,3);
			end


			% geostrophic component
			% phihyd in m2/s2 (P/rho)
			% because r*: now using PHIHYDcR
			P = rdmds('diag_state',ts(t),'rec',6);
			% ? if eta not included: P = P + g*eta;
			P(hFacC_forP==0)=NaN;
			Ug = NaN*P;
			Ug(:,2:end,:) = -P(:,2:end,:)+P(:,1:end-1,:); % on v point 
			Ug = Ug./dycfg; % in m^2/s^2 / m * s
			Vg = NaN*P;
			Vg(2:end,:,:) = (P(2:end,:,:)-P(1:end-1,:,:)); % on u point 
			Vg = Vg./dxcfc;

			% Now average to correct points
			Ug1=Ug; Ug(isnan(Ug))=0; 
			Vg1=Vg; Vg(isnan(Vg))=0; 
			frac = (~isnan(Ug1(1:end-1,1:end-1,:)))+(~isnan(Ug1(2:end,1:end-1,:)))+(~isnan(Ug1(1:end-1,2:end,:)))+(~isnan(Ug1(2:end,2:end,:)));
			Ug(2:end,1:end-1,:) = (Ug(1:end-1,1:end-1,:) + Ug(1:end-1,2:end,:) + Ug(2:end,1:end-1,:) + Ug(2:end,2:end,:))./frac;
			frac = (~isnan(Vg1(1:end-1,1:end-1,:)))+(~isnan(Vg1(2:end,1:end-1,:)))+(~isnan(Vg1(1:end-1,2:end,:)))+(~isnan(Vg1(2:end,2:end,:)));
			Vg(1:end-1,2:end,:) = (Vg(1:end-1,1:end-1,:) + Vg(1:end-1,2:end,:) + Vg(2:end,1:end-1,:) + Vg(2:end,2:end,:))./frac;


			Ug(isnan(Ug))=0;
			Vg(isnan(Vg))=0;

			if savexy
			% velocities not transport
				ug_140E(:,:,t) = squeeze(Ug(109,:,:)); % 140E
				vg_5S(:,:,t) = squeeze(Vg(:,64,:)); % 5S
				vg_5N(:,:,t) = squeeze(Vg(:,94,:)); % 5N
			end

			Ug=Ug([x x(end)+1],y,z);
			Vg=Vg(x,[y y(end)+1],z);

			Ug = Ug.*areaWest;
			Vg = Vg.*areaSouth;

			Ug(U==0)=0;
			Vg(V==0)=0;

			hFacC1=hFacC; hFacC1(end+1,:,:)=hFacC1(end,:,:);
			hFacC2=hFacC; hFacC2(:,end+1,:)=hFacC2(:,end,:);

			Ug(hFacC1==0)=0; % make sure velocities are zero where there is land!
			Vg(hFacC2==0)=0;

			Ua=U-Ug;
			Va=V-Vg;

			Vgsouth_TP5(t) = sum(sum(sum(Vg(109:end,64,:))));
			Vgnorth_TP5(t) = sum(sum(sum(Vg(109:end,94,:))));
			Ugwest_TP5(t) = sum(sum(sum(Ug(109,64:93,:))));
			Vasouth_TP5(t) = sum(sum(sum(Va(109:end,64,:))));
			Vanorth_TP5(t) = sum(sum(sum(Va(109:end,94,:))));
			Uawest_TP5(t) = sum(sum(sum(Ua(109,64:93,:))));

			tmpx = diff(Ug,1,1);
			tmpy = diff(Vg,1,2);
			tmph=tmpx+tmpy;
			divg_h_TP5(t) = sum(sum(sum(tmph(109:end,64:93,:))));
			%divg_h_TP5(t) = sum(sum(sum(tmph)));
			if savexy
				divg_h_xy(:,:,t) = sum(tmpx+tmpy,3);
			end

			tmpx = diff(Ua,1,1);
			tmpy = diff(Va,1,2);
			tmph=tmpx+tmpy;
			diva_h_TP5(t) = sum(sum(sum(tmph(109:end,64:93,:))));
			if savexy
				diva_h_xy(:,:,t) = sum(tmpx+tmpy,3);
			end

			% total tendency, diag
			% tendV(:,:,:,cnt)=(1./mk3D(mygrid.Depth,mygrid.mskC)).*mk3D((ETAN_SNAP2−ETAN_SNAP1)/(86400),mygrid.mskC);
			clear tmp
			if t==1
				tmp = rdmds('diag_snap2d',[0 ts(t)],'rec',1);
				tmp=tmp(x,y,:);
			elseif t==nt
				tmp(:,:,1) = 0*rdmds('diag_snap2d',ts(t-1),'rec',1);
				tmp(:,:,2) = 0*tmp(:,:,1);
				tmp=tmp(x,y,:);
			else
				tmp = rdmds('diag_snap2d',ts(t-1:t),'rec',1);
				tmp=tmp(x,y,:);
			end
			tmp = diff(tmp,1,3)/86400./Depth(x,y);
			tmp(isnan(tmp))=0;
			tmp(isinf(tmp))=0;
			tend_TP5(t) = sum(sum(sum(tmp(109:end,64:93,:).*volume(109:end,64:93,:))));
			%tend_TP5(t) = sum(sum(sum(tmp.*volume)));
			if savexy
				tend_xy(:,:,t) = sum(tmp.*volume,3);
			end

		end % for t

		if savexy
			if m==1 | m==3 | m==5
				surf_xy_m1(:,:,cnt1+(1:nt)) = surf_xy;
				div_h_xy_m1(:,:,cnt1+(1:nt)) = div_h_xy;
				div_v_xy_m1(:,:,cnt1+(1:nt)) = div_v_xy;
				tend_xy_m1(:,:,cnt1+(1:nt)) = tend_xy;
				W300m_xy_m1(:,:,cnt1+(1:nt)) = w_300m;
				U_140E_xy_m1(:,:,cnt1+(1:nt)) = u_140E;
				V_5S_xy_m1(:,:,cnt1+(1:nt)) = v_5S;
				V_5N_xy_m1(:,:,cnt1+(1:nt)) = v_5N;
				Ug_140E_xy_m1(:,:,cnt1+(1:nt)) = Ug_140E_xy;
				Vg_5S_xy_m1(:,:,cnt1+(1:nt)) = Vg_5S_xy;
				Vg_5N_xy_m1(:,:,cnt1+(1:nt)) = Vg_5N_xy;
				Ua_140E_xy_m1(:,:,cnt1+(1:nt)) = Ua_140E_xy;
				Va_5S_xy_m1(:,:,cnt1+(1:nt)) = Va_5S_xy;
				Va_5N_xy_m1(:,:,cnt1+(1:nt)) = Va_5N_xy;
			else
				surf_xy_m2(:,:,cnt2+(1:nt)) = surf_xy;
				div_h_xy_m2(:,:,cnt2+(1:nt)) = div_h_xy;
				div_v_xy_m2(:,:,cnt2+(1:nt)) = div_v_xy;
				tend_xy_m2(:,:,cnt2+(1:nt)) = tend_xy;
				W300m_xy_m2(:,:,cnt2+(1:nt)) = w_300m;
				U_140E_xy_m2(:,:,cnt2+(1:nt)) = u_140E;
				V_5S_xy_m2(:,:,cnt2+(1:nt)) = v_5S;
				V_5N_xy_m2(:,:,cnt2+(1:nt)) = v_5N;
				Ug_140E_xy_m2(:,:,cnt2+(1:nt)) = Ug_140E_xy;
				Vg_5S_xy_m2(:,:,cnt2+(1:nt)) = Vg_5S_xy;
				Vg_5N_xy_m2(:,:,cnt2+(1:nt)) = Vg_5N_xy;
				Ua_140E_xy_m2(:,:,cnt2+(1:nt)) = Ua_140E_xy;
				Va_5S_xy_m2(:,:,cnt2+(1:nt)) = Va_5S_xy;
				Va_5N_xy_m2(:,:,cnt2+(1:nt)) = Va_5N_xy;
			end
		end


		if m==1 | m==3 | m==5
			div_h_TP5_m1(cnt1+(1:nt)) = div_h_TP5;
			div_z_TP5_m1(cnt1+(1:nt)) = div_z_TP5;
			surf_TP5_m1(cnt1+(1:nt)) = surf_TP5;
			tend_TP5_m1(cnt1+(1:nt)) = tend_TP5;
			Vnorth_TP5_m1(cnt1+(1:nt)) = Vnorth_TP5;
			Vsouth_TP5_m1(cnt1+(1:nt)) = Vsouth_TP5;
			Uwest_TP5_m1(cnt1+(1:nt)) = Uwest_TP5;
			Vgnorth_TP5_m1(cnt1+(1:nt)) = Vgnorth_TP5;
			Vgsouth_TP5_m1(cnt1+(1:nt)) = Vgsouth_TP5;
			Ugwest_TP5_m1(cnt1+(1:nt)) = Ugwest_TP5;
			Vanorth_TP5_m1(cnt1+(1:nt)) = Vanorth_TP5;
			Vasouth_TP5_m1(cnt1+(1:nt)) = Vasouth_TP5;
			Uawest_TP5_m1(cnt1+(1:nt)) = Uawest_TP5;
			cnt1=cnt1+nt;
		else
			div_h_TP5_m2(cnt2+(1:nt)) = div_h_TP5;
			div_z_TP5_m2(cnt2+(1:nt)) = div_z_TP5;
			surf_TP5_m2(cnt2+(1:nt)) = surf_TP5;
			tend_TP5_m2(cnt2+(1:nt)) = tend_TP5;
			Vnorth_TP5_m2(cnt2+(1:nt)) = Vnorth_TP5;
			Vsouth_TP5_m2(cnt2+(1:nt)) = Vsouth_TP5;
			Uwest_TP5_m2(cnt2+(1:nt)) = Uwest_TP5;
			Vgnorth_TP5_m2(cnt2+(1:nt)) = Vgnorth_TP5;
			Vgsouth_TP5_m2(cnt2+(1:nt)) = Vgsouth_TP5;
			Ugwest_TP5_m2(cnt2+(1:nt)) = Ugwest_TP5;
			Vanorth_TP5_m2(cnt2+(1:nt)) = Vanorth_TP5;
			Vasouth_TP5_m2(cnt2+(1:nt)) = Vasouth_TP5;
			Uawest_TP5_m2(cnt2+(1:nt)) = Uawest_TP5;
			cnt2=cnt2+nt;
		end


	end % for m
	
end % for yr


%----------------------------------------------------------------
% PART 2: patch overlapping state estimates

% weights for patching 

time = datenum([2014 1 1]):datenum([2019 2 28]);
nt = length(time);
nt1 = length(surf_TP5_m1);

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

div_h_TP5_m1(nt1+1:nt) = 0;
div_h_TP5_m2(1:59) = 0;
div_z_TP5_m1(nt1+1:nt) = 0;
div_z_TP5_m2(1:59) = 0;
surf_TP5_m1(nt1+1:nt) = 0;
surf_TP5_m2(1:59) = 0;
tend_TP5_m1(nt1+1:nt) = 0;
tend_TP5_m2(1:59) = 0;
Uwest_TP5_m1(nt1+1:nt) = 0;
Uwest_TP5_m2(1:59) = 0;
Ugwest_TP5_m1(nt1+1:nt) = 0;
Ugwest_TP5_m2(1:59) = 0;
Uawest_TP5_m1(nt1+1:nt) = 0;
Uawest_TP5_m2(1:59) = 0;
Vnorth_TP5_m1(nt1+1:nt) = 0;
Vnorth_TP5_m2(1:59) = 0;
Vgnorth_TP5_m1(nt1+1:nt) = 0;
Vgnorth_TP5_m2(1:59) = 0;
Vanorth_TP5_m1(nt1+1:nt) = 0;
Vanorth_TP5_m2(1:59) = 0;
Vsouth_TP5_m1(nt1+1:nt) = 0;
Vsouth_TP5_m2(1:59) = 0;
Vgsouth_TP5_m1(nt1+1:nt) = 0;
Vgsouth_TP5_m2(1:59) = 0;
Vasouth_TP5_m1(nt1+1:nt) = 0;
Vasouth_TP5_m2(1:59) = 0;


% average

div_h = weight1.*div_h_TP5_m1 + weight2.*div_h_TP5_m2;
div_z = weight1.*div_z_TP5_m1 + weight2.*div_z_TP5_m2;
surf = weight1.*surf_TP5_m1 + weight2.*surf_TP5_m2;
tend = weight1.*tend_TP5_m1 + weight2.*tend_TP5_m2;
Uwest = weight1.*Uwest_TP5_m1 + weight2.*Uwest_TP5_m2;
Ugwest = weight1.*Ugwest_TP5_m1 + weight2.*Ugwest_TP5_m2;
Uawest = weight1.*Uawest_TP5_m1 + weight2.*Uawest_TP5_m2;
Vnorth = weight1.*Vnorth_TP5_m1 + weight2.*Vnorth_TP5_m2;
Vgnorth = weight1.*Vgnorth_TP5_m1 + weight2.*Vgnorth_TP5_m2;
Vanorth = weight1.*Vanorth_TP5_m1 + weight2.*Vanorth_TP5_m2;
Vsouth = weight1.*Vsouth_TP5_m1 + weight2.*Vsouth_TP5_m2;
Vgsouth = weight1.*Vgsouth_TP5_m1 + weight2.*Vgsouth_TP5_m2;
Vasouth = weight1.*Vasouth_TP5_m1 + weight2.*Vasouth_TP5_m2;


% 30-day running mean
div_h_sm = NaN*div_h;
div_z_sm = NaN*div_z;
surf_sm = NaN*surf;
tend_sm = NaN*tend;
Uwest_sm = NaN*tend;
Ugwest_sm = NaN*tend;
Uawest_sm = NaN*tend;
Vnorth_sm = NaN*tend;
Vgnorth_sm = NaN*tend;
Vanorth_sm = NaN*tend;
Vsouth_sm = NaN*tend;
Vgsouth_sm = NaN*tend;
Vasouth_sm = NaN*tend;
for t=16:nt-15
	div_h_sm(t)=mean(div_h(t-15:t+15));
	div_z_sm(t)=mean(div_z(t-15:t+15));
	surf_sm(t)=mean(surf(t-15:t+15));
	tend_sm(t)=mean(tend(t-15:t+15));
	Uwest_sm(t)=mean(Uwest(t-15:t+15));
	Ugwest_sm(t)=mean(Ugwest(t-15:t+15));
	Uawest_sm(t)=mean(Uawest(t-15:t+15));
	Vnorth_sm(t)=mean(Vnorth(t-15:t+15));
	Vgnorth_sm(t)=mean(Vgnorth(t-15:t+15));
	Vanorth_sm(t)=mean(Vanorth(t-15:t+15));
	Vsouth_sm(t)=mean(Vsouth(t-15:t+15));
	Vgsouth_sm(t)=mean(Vgsouth(t-15:t+15));
	Vasouth_sm(t)=mean(Vasouth(t-15:t+15));
end

res_sm = tend_sm-surf_sm+div_h_sm+div_z_sm;


% mean seasonal cycle
% remove leap day, fill to 5 years
tmp=div_h; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
div_h_anom = tmp-repmat(tmpS',[1 5]);

tmp=div_z; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
div_z_anom = tmp-repmat(tmpS',[1 5]);

tmp=surf; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
surf_anom = tmp-repmat(tmpS',[1 5]);

tmp=tend; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
tend_anom = tmp-repmat(tmpS',[1 5]);

tmp=Uwest; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Uwest_anom = tmp-repmat(tmpS',[1 5]);

tmp=Ugwest; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Ugwest_anom = tmp-repmat(tmpS',[1 5]);

tmp=Uawest; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Uawest_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vnorth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vnorth_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vgnorth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vgnorth_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vanorth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vanorth_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vsouth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vsouth_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vgsouth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vgsouth_anom = tmp-repmat(tmpS',[1 5]);

tmp=Vasouth; tmp(365+365+31+29)=[]; 
if length(tmp)<365*5
	tmp(end+1:365*5)=NaN;
else
	tmp=tmp(1:365*5);
end
tmpS = mean(reshape(tmp,365,5),2);
Vasouth_anom = tmp-repmat(tmpS',[1 5]);



% remove leap day
time1=time; time1(365*2+31+29)=[]; 

nt=length(tend_anom);
if length(time1)<nt
	time1(end+1:nt)=NaN;
else
	time1=time1(1:nt);
end


% 30-day running mean
div_h_anom_sm = NaN*tend_anom;
div_z_anom_sm = NaN*tend_anom;
tend_anom_sm = NaN*tend_anom;
surf_anom_sm = NaN*surf_anom;
Uwest_anom_sm = NaN*surf_anom;
Ugwest_anom_sm = NaN*surf_anom;
Uawest_anom_sm = NaN*surf_anom;
Vnorth_anom_sm = NaN*surf_anom;
Vgnorth_anom_sm = NaN*surf_anom;
Vanorth_anom_sm = NaN*surf_anom;
Vsouth_anom_sm = NaN*surf_anom;
Vgsouth_anom_sm = NaN*surf_anom;
Vasouth_anom_sm = NaN*surf_anom;
for t=16:1825-15
	div_h_anom_sm(t)=mean(div_h_anom(t-15:t+15));
	div_z_anom_sm(t)=mean(div_z_anom(t-15:t+15));
	tend_anom_sm(t)=mean(tend_anom(t-15:t+15));
	surf_anom_sm(t)=mean(surf_anom(t-15:t+15));
	Uwest_anom_sm(t)=mean(Uwest_anom(t-15:t+15));
	Ugwest_anom_sm(t)=mean(Ugwest_anom(t-15:t+15));
	Uawest_anom_sm(t)=mean(Uawest_anom(t-15:t+15));
	Vnorth_anom_sm(t)=mean(Vnorth_anom(t-15:t+15));
	Vgnorth_anom_sm(t)=mean(Vgnorth_anom(t-15:t+15));
	Vanorth_anom_sm(t)=mean(Vanorth_anom(t-15:t+15));
	Vsouth_anom_sm(t)=mean(Vsouth_anom(t-15:t+15));
	Vgsouth_anom_sm(t)=mean(Vgsouth_anom(t-15:t+15));
	Vasouth_anom_sm(t)=mean(Vasouth_anom(t-15:t+15));
end

cd /data/averdy/tpose/budgets/
save volume_budget_box.mat *anom *sm div_h div_z tend surf U*west V*north V*south



% spatial patterns

nt = 1885;
weight1x = permute(repmat(weight1',[1 168 51]),[2 3 1]);
weight2x = permute(repmat(weight2',[1 168 51]),[2 3 1]);
weight1y = permute(repmat(weight1',[1 564 51]),[2 3 1]);
weight2y = permute(repmat(weight2',[1 564 51]),[2 3 1]);
weight1z = permute(repmat(weight1',[1 564 168]),[2 3 1]);
weight2z = permute(repmat(weight2',[1 564 168]),[2 3 1]);


U_140E_xy_m1(:,:,nt1+1:nt) = 0;
U_140E_xy_m2(:,:,1:59) = 0;
V_5N_xy_m1(:,:,nt1+1:nt) = 0;
V_5N_xy_m2(:,:,1:59) = 0;
V_5S_xy_m1(:,:,nt1+1:nt) = 0;
V_5S_xy_m2(:,:,1:59) = 0;
Ug_140E_xy_m1(:,:,nt1+1:nt) = 0;
Ug_140E_xy_m2(:,:,1:59) = 0;
Vg_5N_xy_m1(:,:,nt1+1:nt) = 0;
Vg_5N_xy_m2(:,:,1:59) = 0;
Vg_5S_xy_m1(:,:,nt1+1:nt) = 0;
Vg_5S_xy_m2(:,:,1:59) = 0;
W300m_xy_m1(:,:,nt1+1:nt) = 0;
W300m_xy_m2(:,:,1:59) = 0;
surf_xy_m1(:,:,nt1+1:nt) = 0;
surf_xy_m2(:,:,1:59) = 0;
surf_xy_m1(end:end+1,end:end+1,:) = repmat(surf_xy_m1(end,end,:),[2,2,1]);
surf_xy_m2(end:end+1,end:end+1,:) = repmat(surf_xy_m2(end,end,:),[2,2,1]);

U_140E = weight1x.*U_140E_xy_m1 + weight2x.*U_140E_xy_m2;
V_5N = weight1y.*V_5N_xy_m1 + weight2y.*V_5N_xy_m2;
V_5S = weight1y.*V_5S_xy_m1 + weight2y.*V_5S_xy_m2;
Ug_140E = weight1x.*Ug_140E_xy_m1 + weight2x.*Ug_140E_xy_m2;
Vg_5N = weight1y.*Vg_5N_xy_m1 + weight2y.*Vg_5N_xy_m2;
Vg_5S = weight1y.*Vg_5S_xy_m1 + weight2y.*Vg_5S_xy_m2;
W300m = weight1z.*W300m_xy_m1 + weight2z.*W300m_xy_m2;
surf = weight1z.*surf_xy_m1 + weight2z.*surf_xy_m2;


% monthly means
ndays=repmat([31 28 31 30 31 30 31 31 30 31 30 31],[1 5]);
%ndays(26)=29;
%ndays=ndays(1:end-2);

ind=1;
	for n=1:length(ndays)
	U_140E_monthly(:,:,n) = mean(U_140E(:,:,ind:ind+ndays(n)-1),3); % ! "-1"
	V_5N_monthly(:,:,n) = mean(V_5N(:,:,ind:ind+ndays(n)),3);
	V_5S_monthly(:,:,n) = mean(V_5S(:,:,ind:ind+ndays(n)),3);
	Ug_140E_monthly(:,:,n) = mean(Ug_140E(:,:,ind:ind+ndays(n)),3);
	Vg_5N_monthly(:,:,n) = mean(Vg_5N(:,:,ind:ind+ndays(n)),3);
	Vg_5S_monthly(:,:,n) = mean(Vg_5S(:,:,ind:ind+ndays(n)),3);
	W300m_monthly(:,:,n) = mean(W300m(:,:,ind:ind+ndays(n)),3);
	surf_monthly(:,:,n) = mean(surf(:,:,ind:ind+ndays(n)),3);
	ind=ind+ndays(n);
end


% mean seasonal cycle
tmp=U_140E_monthly; 
tmpS = mean(reshape(tmp,168,51,12,5),4);
U_140E_monthly_seas = tmpS;
U_140E_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=V_5N_monthly; 
tmpS = mean(reshape(tmp,564,51,12,5),4);
V_5N_monthly_seas = tmpS;
V_5N_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=V_5S_monthly; 
tmpS = mean(reshape(tmp,564,51,12,5),4);
V_5S_monthly_seas = tmpS;
V_5S_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=Ug_140E_monthly; 
tmpS = mean(reshape(tmp,168,51,12,5),4);
Ug_140E_monthly_seas = tmpS;
Ug_140E_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=Vg_5N_monthly; 
tmpS = mean(reshape(tmp,564,51,12,5),4);
Vg_5N_monthly_seas = tmpS;
Vg_5N_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=Vg_5S_monthly; 
tmpS = mean(reshape(tmp,564,51,12,5),4);
Vg_5S_monthly_seas = tmpS;
Vg_5S_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=W300m_monthly; 
tmpS = mean(reshape(tmp,564,168,12,5),4);
W300m_monthly_seas = tmpS;
W300m_monthly_anom = tmp-repmat(tmpS,[1 1 5]);

tmp=surf_monthly; 
tmpS = mean(reshape(tmp,564,168,12,5),4);
surf_monthly_seas = tmpS;
surf_monthly_anom = tmp-repmat(tmpS,[1 1 5]);


clear *m1 *m2 *sm
save volume_budget_faces.mat *140E* *5N* *5S* W300* surf_monthly*
