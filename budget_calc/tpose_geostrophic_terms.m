% Geostrophic transport for T and S
% adapted from code by V. Tamsitt

% first decompose velocities into geostrophic and ageostrophic components
% then calculate heat/salt transport similar to how the model does it

%----------------------------------------------------------------
% PART 1: read overlapping state estimates

addpath ~/scripts_m

% reference values
Tmean = 18.80;
Smean = 34.96;


% select area
% for averaging
x1 = 109:563; % 140E to eastern boundary
y1 = 64:93; % 5S to 5N
z1 = 1:22; % top 300 m


% select area
%x = 49:528; % 120E to 80W
%x = 109:528; % 140E to 80W
x = 109:563; % 140E to eastern boundary
y = 64:93; % 5S to 5N
z = 1:22; % top 300 m

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

		clear adv*g adv*a adv*tot 
 
		% path to diagnostics
		cd(['/data/SO6/TPOSE_diags/tpose3/' char(mon(m)) num2str(yr) '/last_iter_diags/']);

		diag_budget_file = 'diag_heat_budget';

		% time stepping: 1 day
		nt = length(dir('diag_heat_budget*data'));
		dt = 48;
		ts = dt:dt:(nt*dt);

		% read diagnostics
		% calculate tendencies 
		for t=1:nt

			% THETA field 
			T = rdmds('diag_state',ts(t),'rec',1);

			% SALT field 
			S = rdmds('diag_state',ts(t),'rec',2);

			% velocities: UVEL, VVEL, WVEL
			vel = rdmds('diag_state',ts(t),'rec',3:5);
			U = vel(:,:,:,1); %.*areaWest;
			V = vel(:,:,:,2); %.*areaSouth;
			W = vel(:,:,:,3); %.*areaTop; 

			% geostrophic component
			% phihyd in m2/s2 (P/rho)
			P = rdmds('diag_state',ts(t),'rec',6);
			% ? if eta not included: P = P + g*eta;
			% checked: P includes eta
			P(hFacC==0)=NaN;
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
			Ug(U==0)=0;
			Vg(V==0)=0;

			% make sure velocities are zero where there is land!
			Ug(hFacC==0)=0; 
			Vg(hFacC==0)=0;

			% remaining nans?
			Ug(isnan(Ug))=0;
			Vg(isnan(Vg))=0;

			Ua=U-Ug;
			Va=V-Vg;

			% calculate advection terms using same advection scheme as MITgcm code

			ADVxg=zeros(nx,ny,nz);
			ADVyg=zeros(nx,ny,nz);
			ADVxa=zeros(nx,ny,nz);
			ADVya=zeros(nx,ny,nz);
			ADVxtot=zeros(nx,ny,nz);
			ADVytot=zeros(nx,ny,nz);

			rhoFacC = 1;
			maskC = hFacC; maskC(maskC~=0)=1;
			maskW = hFacW; maskW(maskW~=0)=1;
			maskS = hFacS; maskS(maskS~=0)=1;
			RecipVol = volume;RecipVol(RecipVol==0)=inf;RecipVol = 1./RecipVol;
			recip_dxC = DXC;recip_dxC(recip_dxC==0)=inf;recip_dxC = 1./recip_dxC;
			recip_dyC = DYC;recip_dyC(recip_dyC==0)=inf;recip_dyC = 1./recip_dyC;

			oneSixth = 1/6;
			deltaTloc = 1800; % timestep in seconds

			for k=1:nz
				ugTrans = Ug(:,:,k).*areaWest(:,:,k)*rhoFacC; %u transport
				vgTrans = Vg(:,:,k).*areaSouth(:,:,k)*rhoFacC; %v transport
				uaTrans = (U(:,:,k).*areaWest_noh(:,:,k)-Ug(:,:,k).*areaWest(:,:,k))*rhoFacC; %u transport
				vaTrans = (V(:,:,k).*areaSouth_noh(:,:,k)-Vg(:,:,k).*areaSouth(:,:,k))*rhoFacC; %v transport
				utotTrans = U(:,:,k).*areaWest_noh(:,:,k)*rhoFacC; %u transport
				vtotTrans = V(:,:,k).*areaSouth_noh(:,:,k)*rhoFacC; %v transport
						
				maskLocW = maskW(:,:,k);
				maskLocS = maskS(:,:,k);
				ugT=zeros(nx,ny);
				vgT=zeros(nx,ny);
				uaT=zeros(nx,ny);
				vaT=zeros(nx,ny);
				utotT=zeros(nx,ny);
				vtotT=zeros(nx,ny);


				% TEMPERATURE
				tracer = T(:,:,k);
	 
				for i=3:nx-1
					% weighting terms
					Rjp = (tracer(i+1,:)-tracer(i,:)).*maskLocW(i+1,:);
					Rj  = (tracer(i,:)-tracer(i-1,:)).*maskLocW(i,:);
					Rjm = (tracer(i-1,:)-tracer(i-2,:)).*maskLocW(i-1,:);
				
					ugCFL = Ug(i,:,k);
					ugCFL = abs(Ug(i,:,k)*deltaTloc.*recip_dxC(i,:)); 
					uaCFL = U(i,:,k)-Ug(i,:,k);
					uaCFL = abs((U(i,:,k)-Ug(i,:,k))*deltaTloc.*recip_dxC(i,:));
					utotCFL = U(i,:,k);
					utotCFL = abs(U(i,:,k)*deltaTloc.*recip_dxC(i,:));
					d0g = (2-ugCFL).*(1-ugCFL)*oneSixth;
					d1g = (1-ugCFL.*ugCFL)*oneSixth;
					d0a = (2-uaCFL).*(1-uaCFL)*oneSixth;
					d1a = (1-uaCFL.*uaCFL)*oneSixth;
					d0tot = (2-utotCFL).*(1-utotCFL)*oneSixth;
					d1tot = (1-utotCFL.*utotCFL)*oneSixth;
					ugT(i,:)=0.5*(ugTrans(i,:)+abs(ugTrans(i,:))).*(tracer(i-1,:)+(d0g.*Rj+d1g.*Rjm))...
							 +0.5*(ugTrans(i,:)-abs(ugTrans(i,:))).*(tracer(i,:)-(d0g.*Rj+d1g.*Rjp));
							 uaT(i,:)=0.5*(uaTrans(i,:)+abs(uaTrans(i,:))).*(tracer(i-1,:)+(d0a.*Rj+d1a.*Rjm))...
							 +0.5*(uaTrans(i,:)-abs(uaTrans(i,:))).*(tracer(i,:)-(d0a.*Rj+d1a.*Rjp));
					utotT(i,:)=0.5*(utotTrans(i,:)+abs(utotTrans(i,:))).*(tracer(i-1,:)+(d0tot.*Rj+d1tot.*Rjm))...
							   +0.5*(utotTrans(i,:)-abs(utotTrans(i,:))).*(tracer(i,:)-(d0tot.*Rj+d1tot.*Rjp));
				end
				ADVTxg(:,:,k) = ugT;
				ADVTxa(:,:,k) = uaT;
				ADVTxtot(:,:,k) = utotT;

				for j=3:ny-1
					% weighting terms
					Rip = (tracer(:,j+1)-tracer(:,j)).*maskLocS(:,j+1);
					Ri  = (tracer(:,j)-tracer(:,j-1)).*maskLocS(:,j);
					Rim = (tracer(:,j-1)-tracer(:,j-2)).*maskLocS(:,j-1);
					vgCFL = Vg(:,j,k);
					vgCFL = abs(Vg(:,j,k)*deltaTloc.*recip_dyC(:,j)); 
					vaCFL = V(:,j,k)-Vg(:,j,k);
					vaCFL = abs((V(:,j,k)-Vg(:,j,k))*deltaTloc.*recip_dyC(:,j));
					vtotCFL = V(:,j,k);
					vtotCFL = abs(V(:,j,k)*deltaTloc.*recip_dyC(:,j));
					d0g = (2-vgCFL).*(1-vgCFL)*oneSixth;
					d1g = (1-vgCFL.*vgCFL)*oneSixth;
					d0a = (2-vaCFL).*(1-vaCFL)*oneSixth;
					d1a = (1-vaCFL.*vaCFL)*oneSixth;
					d0tot = (2-vtotCFL).*(1-vtotCFL)*oneSixth;
					d1tot = (1-vtotCFL.*vtotCFL)*oneSixth;
					
					vgT(:,j)=0.5*(vgTrans(:,j)+abs(vgTrans(:,j))).*(tracer(:,j-1)+(d0g.*Ri+d1g.*Rim))...
							  +0.5*(vgTrans(:,j)-abs(vgTrans(:,j))).*(tracer(:,j)-(d0g.*Ri+d1g.*Rip));
					vaT(:,j)=0.5*(vaTrans(:,j)+abs(vaTrans(:,j))).*(tracer(:,j-1)+(d0a.*Ri+d1a.*Rim))...
							  +0.5*(vaTrans(:,j)-abs(vaTrans(:,j))).*(tracer(:,j)-(d0a.*Ri+d1a.*Rip));
					vtotT(:,j)=0.5*(vtotTrans(:,j)+abs(vtotTrans(:,j))).*(tracer(:,j-1)+(d0tot.*Ri+d1tot.*Rim))...
								+0.5*(vtotTrans(:,j)-abs(vtotTrans(:,j))).*(tracer(:,j)-(d0tot.*Ri+d1tot.*Rip));
				end
				ADVTyg(:,:,k) = vgT;
				ADVTya(:,:,k) = vaT;
				ADVTytot(:,:,k) = vtotT;


				% SALT
				tracer = S(:,:,k);
	 
				for i=3:nx-1
					% weighting terms
					Rjp = (tracer(i+1,:)-tracer(i,:)).*maskLocW(i+1,:);
					Rj  = (tracer(i,:)-tracer(i-1,:)).*maskLocW(i,:);
					Rjm = (tracer(i-1,:)-tracer(i-2,:)).*maskLocW(i-1,:);
				
					ugCFL = Ug(i,:,k);
					ugCFL = abs(Ug(i,:,k)*deltaTloc.*recip_dxC(i,:)); 
					%uCFL = | u *deltaT/delta x| (dimensionless), >1 means u is bigger than grid resolution
					uaCFL = U(i,:,k)-Ug(i,:,k);
					uaCFL = abs((U(i,:,k)-Ug(i,:,k))*deltaTloc.*recip_dxC(i,:));
					utotCFL = U(i,:,k);
					utotCFL = abs(U(i,:,k)*deltaTloc.*recip_dxC(i,:));
					d0g = (2-ugCFL).*(1-ugCFL)*oneSixth;
					d1g = (1-ugCFL.*ugCFL)*oneSixth;
					d0a = (2-uaCFL).*(1-uaCFL)*oneSixth;
					d1a = (1-uaCFL.*uaCFL)*oneSixth;
					d0tot = (2-utotCFL).*(1-utotCFL)*oneSixth;
					d1tot = (1-utotCFL.*utotCFL)*oneSixth;

					ugT(i,:)=0.5*(ugTrans(i,:)+abs(ugTrans(i,:))).*(tracer(i-1,:)+(d0g.*Rj+d1g.*Rjm))...
							  +0.5*(ugTrans(i,:)-abs(ugTrans(i,:))).*(tracer(i,:)-(d0g.*Rj+d1g.*Rjp));
					uaT(i,:)=0.5*(uaTrans(i,:)+abs(uaTrans(i,:))).*(tracer(i-1,:)+(d0a.*Rj+d1a.*Rjm))...
							  +0.5*(uaTrans(i,:)-abs(uaTrans(i,:))).*(tracer(i,:)-(d0a.*Rj+d1a.*Rjp));
					utotT(i,:)=0.5*(utotTrans(i,:)+abs(utotTrans(i,:))).*(tracer(i-1,:)+(d0tot.*Rj+d1tot.*Rjm))...
								+0.5*(utotTrans(i,:)-abs(utotTrans(i,:))).*(tracer(i,:)-(d0tot.*Rj+d1tot.*Rjp));
				end
				ADVSxg(:,:,k) = ugT;
				ADVSxa(:,:,k) = uaT;
				ADVSxtot(:,:,k) = utotT;

				for j=3:ny-1
					% weighting terms
					Rip = (tracer(:,j+1)-tracer(:,j)).*maskLocS(:,j+1);
					Ri  = (tracer(:,j)-tracer(:,j-1)).*maskLocS(:,j);
					Rim = (tracer(:,j-1)-tracer(:,j-2)).*maskLocS(:,j-1);
					vgCFL = Vg(:,j,k);
					vgCFL = abs(Vg(:,j,k)*deltaTloc.*recip_dyC(:,j)); 
					vaCFL = V(:,j,k)-Vg(:,j,k);
					vaCFL = abs((V(:,j,k)-Vg(:,j,k))*deltaTloc.*recip_dyC(:,j));
					vtotCFL = V(:,j,k);
					vtotCFL = abs(V(:,j,k)*deltaTloc.*recip_dyC(:,j));
					d0g = (2-vgCFL).*(1-vgCFL)*oneSixth;
					d1g = (1-vgCFL.*vgCFL)*oneSixth;
					d0a = (2-vaCFL).*(1-vaCFL)*oneSixth;
					d1a = (1-vaCFL.*vaCFL)*oneSixth;
					d0tot = (2-vtotCFL).*(1-vtotCFL)*oneSixth;
					d1tot = (1-vtotCFL.*vtotCFL)*oneSixth;
					
					vgT(:,j)=0.5*(vgTrans(:,j)+abs(vgTrans(:,j))).*(tracer(:,j-1)+(d0g.*Ri+d1g.*Rim))...
							  +0.5*(vgTrans(:,j)-abs(vgTrans(:,j))).*(tracer(:,j)-(d0g.*Ri+d1g.*Rip));
					vaT(:,j)=0.5*(vaTrans(:,j)+abs(vaTrans(:,j))).*(tracer(:,j-1)+(d0a.*Ri+d1a.*Rim))...
							  +0.5*(vaTrans(:,j)-abs(vaTrans(:,j))).*(tracer(:,j)-(d0a.*Ri+d1a.*Rip));
					vtotT(:,j)=0.5*(vtotTrans(:,j)+abs(vtotTrans(:,j))).*(tracer(:,j-1)+(d0tot.*Ri+d1tot.*Rim))...
								+0.5*(vtotTrans(:,j)-abs(vtotTrans(:,j))).*(tracer(:,j)-(d0tot.*Ri+d1tot.*Rip));
				end
				ADVSyg(:,:,k) = vgT;
				ADVSya(:,:,k) = vaT;
				ADVSytot(:,:,k) = vtotT;

			end % for k


			% remove mean T
			% adv = u(T-Tref)

			advTWg(:,:,t) = squeeze(ADVTxg(x1(1),y1,z1)) - squeeze(Ug(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Tmean;
			advTEg(:,:,t) = squeeze(ADVTxg(x1(end)+1,y1,z1)) - squeeze(Ug(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Tmean;
			advTSg(:,:,t) = squeeze(ADVTyg(x1,y1(1),z1)) - squeeze(Vg(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Tmean;
			advTNg(:,:,t) = squeeze(ADVTyg(x1,y1(end)+1,z1)) - squeeze(Vg(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Tmean;

			advTWa(:,:,t) = squeeze(ADVTxa(x1(1),y1,z1)) - squeeze(Ua(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Tmean;
			advTEa(:,:,t) = squeeze(ADVTxa(x1(end)+1,y1,z1)) - squeeze(Ua(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Tmean;
			advTSa(:,:,t) = squeeze(ADVTya(x1,y1(1),z1)) - squeeze(Va(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Tmean;
			advTNa(:,:,t) = squeeze(ADVTya(x1,y1(end)+1,z1)) - squeeze(Va(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Tmean;

			advTWtot(:,:,t) = squeeze(ADVTxtot(x1(1),y1,z1)) - squeeze(U(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Tmean;
			advTEtot(:,:,t) = squeeze(ADVTxtot(x1(end)+1,y1,z1)) - squeeze(U(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Tmean;
			advTStot(:,:,t) = squeeze(ADVTytot(x1,y1(1),z1)) - squeeze(V(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Tmean;
			advTNtot(:,:,t) = squeeze(ADVTytot(x1,y1(end)+1,z1)) - squeeze(V(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Tmean;


			advSWg(:,:,t) = squeeze(ADVSxg(x1(1),y1,z1)) - squeeze(Ug(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Smean;
			advSEg(:,:,t) = squeeze(ADVSxg(x1(end)+1,y1,z1)) - squeeze(Ug(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Smean;
			advSSg(:,:,t) = squeeze(ADVSyg(x1,y1(1),z1)) - squeeze(Vg(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Smean;
			advSNg(:,:,t) = squeeze(ADVSyg(x1,y1(end)+1,z1)) - squeeze(Vg(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Smean;

			advSWa(:,:,t) = squeeze(ADVSxa(x1(1),y1,z1)) - squeeze(Ua(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Smean;
			advSEa(:,:,t) = squeeze(ADVSxa(x1(end)+1,y1,z1)) - squeeze(Ua(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Smean;
			advSSa(:,:,t) = squeeze(ADVSya(x1,y1(1),z1)) - squeeze(Va(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Smean;
			advSNa(:,:,t) = squeeze(ADVSya(x1,y1(end)+1,z1)) - squeeze(Va(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Smean;

			advSWtot(:,:,t) = squeeze(ADVSxtot(x1(1),y1,z1)) - squeeze(U(x1(1),y1,z1).*areaWest(x1(1),y1,z1))*Smean;
			advSEtot(:,:,t) = squeeze(ADVSxtot(x1(end)+1,y1,z1)) - squeeze(U(x1(end)+1,y1,z1).*areaWest(x1(end)+1,y1,z1))*Smean;
			advSStot(:,:,t) = squeeze(ADVSytot(x1,y1(1),z1)) - squeeze(V(x1,y1(1),z1).*areaSouth(x1,y1(1),z1))*Smean;
			advSNtot(:,:,t) = squeeze(ADVSytot(x1,y1(end)+1,z1)) - squeeze(V(x1,y1(end)+1,z1).*areaSouth(x1,y1(end)+1,z1))*Smean;


		end % for t


		if m==1 | m==3 | m==5
			 advTWg1(:,:,cnt1+(1:nt)) = advTWg;
			 advTEg1(:,:,cnt1+(1:nt)) = advTEg;
			 advTNg1(:,:,cnt1+(1:nt)) = advTNg;
			 advTSg1(:,:,cnt1+(1:nt)) = advTSg;
			 advTWa1(:,:,cnt1+(1:nt)) = advTWa;
			 advTEa1(:,:,cnt1+(1:nt)) = advTEa;
			 advTNa1(:,:,cnt1+(1:nt)) = advTNa;
			 advTSa1(:,:,cnt1+(1:nt)) = advTSa;
			 advTWtot1(:,:,cnt1+(1:nt)) = advTWtot;
			 advTEtot1(:,:,cnt1+(1:nt)) = advTEtot;
			 advTNtot1(:,:,cnt1+(1:nt)) = advTNtot;
			 advTStot1(:,:,cnt1+(1:nt)) = advTStot;
			 advSWg1(:,:,cnt1+(1:nt)) = advSWg;
			 advSEg1(:,:,cnt1+(1:nt)) = advSEg;
			 advSNg1(:,:,cnt1+(1:nt)) = advSNg;
			 advSSg1(:,:,cnt1+(1:nt)) = advSSg;
			 advSWa1(:,:,cnt1+(1:nt)) = advSWa;
			 advSEa1(:,:,cnt1+(1:nt)) = advSEa;
			 advSNa1(:,:,cnt1+(1:nt)) = advSNa;
			 advSSa1(:,:,cnt1+(1:nt)) = advSSa;
			 advSWtot1(:,:,cnt1+(1:nt)) = advSWtot;
			 advSEtot1(:,:,cnt1+(1:nt)) = advSEtot;
			 advSNtot1(:,:,cnt1+(1:nt)) = advSNtot;
			 advSStot1(:,:,cnt1+(1:nt)) = advSStot;
			 cnt1=cnt1+nt;
		else
			 advTWg2(:,:,cnt2+(1:nt)) = advTWg;
			 advTEg2(:,:,cnt2+(1:nt)) = advTEg;
			 advTNg2(:,:,cnt2+(1:nt)) = advTNg;
			 advTSg2(:,:,cnt2+(1:nt)) = advTSg;
			 advTWa2(:,:,cnt2+(1:nt)) = advTWa;
			 advTEa2(:,:,cnt2+(1:nt)) = advTEa;
			 advTNa2(:,:,cnt2+(1:nt)) = advTNa;
			 advTSa2(:,:,cnt2+(1:nt)) = advTSa;
			 advTWtot2(:,:,cnt2+(1:nt)) = advTWtot;
			 advTEtot2(:,:,cnt2+(1:nt)) = advTEtot;
			 advTNtot2(:,:,cnt2+(1:nt)) = advTNtot;
			 advTStot2(:,:,cnt2+(1:nt)) = advTStot;
			 advSWg2(:,:,cnt2+(1:nt)) = advSWg;
			 advSEg2(:,:,cnt2+(1:nt)) = advSEg;
			 advSNg2(:,:,cnt2+(1:nt)) = advSNg;
			 advSSg2(:,:,cnt2+(1:nt)) = advSSg;
			 advSWa2(:,:,cnt2+(1:nt)) = advSWa;
			 advSEa2(:,:,cnt2+(1:nt)) = advSEa;
			 advSNa2(:,:,cnt2+(1:nt)) = advSNa;
			 advSSa2(:,:,cnt2+(1:nt)) = advSSa;
			 advSWtot2(:,:,cnt2+(1:nt)) = advSWtot;
			 advSEtot2(:,:,cnt2+(1:nt)) = advSEtot;
			 advSNtot2(:,:,cnt2+(1:nt)) = advSNtot;
			 advSStot2(:,:,cnt2+(1:nt)) = advSStot;
			 cnt2=cnt2+nt;
		end

	end % for m

end % for yr


%----------------------------------------------------------------
% PART 2: patch overlapping state estimates

% weights for patching 

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

advTWg1(:,:,nt1+1:nt) = 0;
advTWg2(:,:,1:59) = 0;
advTSg1(:,:,nt1+1:nt) = 0;
advTSg2(:,:,1:59) = 0;
advTNg1(:,:,nt1+1:nt) = 0;
advTNg2(:,:,1:59) = 0;
advTWa1(:,:,nt1+1:nt) = 0;
advTWa2(:,:,1:59) = 0;
advTSa1(:,:,nt1+1:nt) = 0;
advTSa2(:,:,1:59) = 0;
advTNa1(:,:,nt1+1:nt) = 0;
advTNa2(:,:,1:59) = 0;
advSWg1(:,:,nt1+1:nt) = 0;
advSWg2(:,:,1:59) = 0;
advSSg1(:,:,nt1+1:nt) = 0;
advSSg2(:,:,1:59) = 0;
advSNg1(:,:,nt1+1:nt) = 0;
advSNg2(:,:,1:59) = 0;
advSWa1(:,:,nt1+1:nt) = 0;
advSWa2(:,:,1:59) = 0;
advSSa1(:,:,nt1+1:nt) = 0;
advSSa2(:,:,1:59) = 0;
advSNa1(:,:,nt1+1:nt) = 0;
advSNa2(:,:,1:59) = 0;


% average

weight1=weight1(1:nx,1:nz,:);
weight2=weight2(1:nx,1:nz,:);

advTSg = weight1.*advTSg1 + weight2.*advTSg2;
advTNg = weight1.*advTNg1 + weight2.*advTNg2;
advTSa = weight1.*advTSa1 + weight2.*advTSa2;
advTNa = weight1.*advTNa1 + weight2.*advTNa2;
advSSg = weight1.*advSSg1 + weight2.*advSSg2;
advSNg = weight1.*advSNg1 + weight2.*advSNg2;
advSSa = weight1.*advSSa1 + weight2.*advSSa2;
advSNa = weight1.*advSNa1 + weight2.*advSNa2;

weight1=weight1(1:ny,1:nz,:);
weight2=weight2(1:ny,1:nz,:);

advTWg = weight1.*advTWg1 + weight2.*advTWg2;
advTWa = weight1.*advTWa1 + weight2.*advTWa2;
advSWg = weight1.*advSWg1 + weight2.*advSWg2;
advSWa = weight1.*advSWa1 + weight2.*advSWa2;

areaWest140E = squeeze(areaWest(109,:,1:nz));
areaSouth5S = squeeze(areaSouth(:,64,1:nz));
areaSouth5N = squeeze(areaSouth(:,94,1:nz));


% 2d:
advTNg_2d = advTNg;
advTNa_2d = advTNa;
advTSg_2d = advTSg;
advTSa_2d = advTSa;
advTWg_2d = advTWg;
advTWa_2d = advTWa;

advSNg_2d = advSNg;
advSNa_2d = advSNa;
advSSg_2d = advSSg;
advSSa_2d = advSSa;
advSWg_2d = advSWg;
advSWa_2d = advSWa;


% timse series
% scale by volume

advTNg(isnan(advTNg))=0;
advTNa(isnan(advTNa))=0;
advTSg(isnan(advTSg))=0;
advTSa(isnan(advTSa))=0;
advTWg(isnan(advTWg))=0;
advTWa(isnan(advTWa))=0;

advSNg(isnan(advSNg))=0;
advSNa(isnan(advSNa))=0;
advSSg(isnan(advSSg))=0;
advSSa(isnan(advSSa))=0;
advSWg(isnan(advSWg))=0;
advSWa(isnan(advSWa))=0;

advTNg = squeeze(sum(sum(advTNg,1),2));
advTNa = squeeze(sum(sum(advTNa,1),2));
advTSg = squeeze(sum(sum(advTSg,1),2));
advTSa = squeeze(sum(sum(advTSa,1),2));
advTWg = squeeze(sum(sum(advTWg,1),2));
advTWa = squeeze(sum(sum(advTWa,1),2));

advSNg = squeeze(sum(sum(advSNg,1),2));
advSNa = squeeze(sum(sum(advSNa,1),2));
advSSg = squeeze(sum(sum(advSSg,1),2));
advSSa = squeeze(sum(sum(advSSa,1),2));
advSWg = squeeze(sum(sum(advSWg,1),2));
advSWa = squeeze(sum(sum(advSWa,1),2));


advTNg = advTNg'./repmat(vol,[1 nt]);
advTNa = advTNa'./repmat(vol,[1 nt]);
advTSg = advTSg'./repmat(vol,[1 nt]);
advTSa = advTSa'./repmat(vol,[1 nt]);
advTWg = advTWg'./repmat(vol,[1 nt]);
advTWa = advTWa'./repmat(vol,[1 nt]);

advSNg = advSNg'./repmat(vol,[1 nt]);
advSNa = advSNa'./repmat(vol,[1 nt]);
advSSg = advSSg'./repmat(vol,[1 nt]);
advSSa = advSSa'./repmat(vol,[1 nt]);
advSWg = advSWg'./repmat(vol,[1 nt]);
advSWa = advSWa'./repmat(vol,[1 nt]);



% mean seas cycle
tmp=advTNg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTNg_anom = tmp-tmpS;

tmp=advTNa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTNa_anom = tmp-tmpS;

tmp=advTSg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTSg_anom = tmp-tmpS;

tmp=advTSa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTSa_anom = tmp-tmpS;

tmp=advTWg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTWg_anom = tmp-tmpS;

tmp=advTWa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advTWa_anom = tmp-tmpS;

tmp=advSNg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSNg_anom = tmp-tmpS;

tmp=advSNa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSNa_anom = tmp-tmpS;

tmp=advSSg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSSg_anom = tmp-tmpS;

tmp=advSSa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSSa_anom = tmp-tmpS;

tmp=advSWg; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSWg_anom = tmp-tmpS;

tmp=advSWa; tmp(365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(1:365*5),365,5),2);
tmpS=repmat(tmpS',[1 5]); tmpS(end+1:365*5+59)=tmpS(1:59);
advSWa_anom = tmp-tmpS;


% remove leap day
time1=time; time1(365*2+31+29)=[]; 


% monthly means
ndays=repmat([31 28 31 30 31 30 31 31 30 31 30 31],[1 5]);

ind=1;
for n=1:length(ndays)
	 advTWg_anom_monthly(n) = mean(advTWg_anom(ind:ind+ndays(n)));
	 advTWa_anom_monthly(n) = mean(advTWa_anom(ind:ind+ndays(n)));
	 advTNg_anom_monthly(n) = mean(advTNg_anom(ind:ind+ndays(n)));
	 advTNa_anom_monthly(n) = mean(advTNa_anom(ind:ind+ndays(n)));
	 advTSg_anom_monthly(n) = mean(advTSg_anom(ind:ind+ndays(n)));
	 advTSa_anom_monthly(n) = mean(advTSa_anom(ind:ind+ndays(n)));
	 advSWg_anom_monthly(n) = mean(advSWg_anom(ind:ind+ndays(n)));
	 advSWa_anom_monthly(n) = mean(advSWa_anom(ind:ind+ndays(n)));
	 advSNg_anom_monthly(n) = mean(advSNg_anom(ind:ind+ndays(n)));
	 advSNa_anom_monthly(n) = mean(advSNa_anom(ind:ind+ndays(n)));
	 advSSg_anom_monthly(n) = mean(advSSg_anom(ind:ind+ndays(n)));
	 advSSa_anom_monthly(n) = mean(advSSa_anom(ind:ind+ndays(n)));
	 ind=ind+ndays(n);
end


% 30-day running mean
advTSg_sm = NaN*advTSg;
advTSa_sm = NaN*advTSa;
advTNg_sm = NaN*advTNg;
advTNa_sm = NaN*advTNa;
advTWg_sm = NaN*advTWg;
advTWa_sm = NaN*advTWa;
advSSg_sm = NaN*advSSg;
advSSa_sm = NaN*advSSa;
advSNg_sm = NaN*advSNg;
advSNa_sm = NaN*advSNa;
advSWg_sm = NaN*advSWg;
advSWa_sm = NaN*advSWa;
for t=16:nt-15
	 advTSg_sm(t)=mean(advTSg(t-15:t+15));
	 advTSa_sm(t)=mean(advTSa(t-15:t+15));
	 advTNg_sm(t)=mean(advTNg(t-15:t+15));
	 advTNa_sm(t)=mean(advTNa(t-15:t+15));
	 advTWg_sm(t)=mean(advTWg(t-15:t+15));
	 advTWa_sm(t)=mean(advTWa(t-15:t+15));
	 advSSg_sm(t)=mean(advSSg(t-15:t+15));
	 advSSa_sm(t)=mean(advSSa(t-15:t+15));
	 advSNg_sm(t)=mean(advSNg(t-15:t+15));
	 advSNa_sm(t)=mean(advSNa(t-15:t+15));
	 advSWg_sm(t)=mean(advSWg(t-15:t+15));
	 advSWa_sm(t)=mean(advSWa(t-15:t+15));
end


% 30-day running mean
advTSg_anom_sm = NaN*advTSg_anom;
advTSa_anom_sm = NaN*advTSa_anom;
advTNg_anom_sm = NaN*advTNg_anom;
advTNa_anom_sm = NaN*advTNa_anom;
advTWg_anom_sm = NaN*advTWg_anom;
advTWa_anom_sm = NaN*advTWa_anom;
advSSg_anom_sm = NaN*advSSg_anom;
advSSa_anom_sm = NaN*advSSa_anom;
advSNg_anom_sm = NaN*advSNg_anom;
advSNa_anom_sm = NaN*advSNa_anom;
advSWg_anom_sm = NaN*advSWg_anom;
advSWa_anom_sm = NaN*advSWa_anom;
for t=16:nt-1-15
	 advTSg_anom_sm(t)=mean(advTSg_anom(t-15:t+15));
	 advTSa_anom_sm(t)=mean(advTSa_anom(t-15:t+15));
	 advTNg_anom_sm(t)=mean(advTNg_anom(t-15:t+15));
	 advTNa_anom_sm(t)=mean(advTNa_anom(t-15:t+15));
	 advTWg_anom_sm(t)=mean(advTWg_anom(t-15:t+15));
	 advTWa_anom_sm(t)=mean(advTWa_anom(t-15:t+15));
	 advSSg_anom_sm(t)=mean(advSSg_anom(t-15:t+15));
	 advSSa_anom_sm(t)=mean(advSSa_anom(t-15:t+15));
	 advSNg_anom_sm(t)=mean(advSNg_anom(t-15:t+15));
	 advSNa_anom_sm(t)=mean(advSNa_anom(t-15:t+15));
	 advSWg_anom_sm(t)=mean(advSWg_anom(t-15:t+15));
	 advSWa_anom_sm(t)=mean(advSWa_anom(t-15:t+15));
end


cd /data/averdy/tpose/budgets/
save geostr_budgets_box.mat *anom *monthly *a *g *sm 



% spatial patterns
% scale by area

advTWg_2d = advTWg_2d./repmat(areaWest140E(64:93,z),[1 1 nt]);
advTWa_2d = advTWa_2d./repmat(areaWest140E(64:93,z),[1 1 nt]);
advSWg_2d = advSWg_2d./repmat(areaWest140E(64:93,z),[1 1 nt]);
advSWa_2d = advSWa_2d./repmat(areaWest140E(64:93,z),[1 1 nt]);

advTNg_2d = advTNg_2d./repmat(areaSouth5N(109:end-1,z),[1 1 nt]);
advTNa_2d = advTNa_2d./repmat(areaSouth5N(109:end-1,z),[1 1 nt]);
advSNg_2d = advSNg_2d./repmat(areaSouth5N(109:end-1,z),[1 1 nt]);
advSNa_2d = advSNa_2d./repmat(areaSouth5N(109:end-1,z),[1 1 nt]);

advTSg_2d = advTSg_2d./repmat(areaSouth5S(109:end-1,z),[1 1 nt]);
advTSa_2d = advTSa_2d./repmat(areaSouth5S(109:end-1,z),[1 1 nt]);
advSSg_2d = advSSg_2d./repmat(areaSouth5S(109:end-1,z),[1 1 nt]);
advSSa_2d = advSSa_2d./repmat(areaSouth5S(109:end-1,z),[1 1 nt]);


% mean seas cycle
tmp=advTWg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),ny,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTWg_2d_anom = tmp-tmpS;

tmp=advSWg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),ny,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSWg_2d_anom = tmp-tmpS;

tmp=advTSg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTSg_2d_anom = tmp-tmpS;

tmp=advSSg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSSg_2d_anom = tmp-tmpS;

tmp=advTNg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTNg_2d_anom = tmp-tmpS;

tmp=advSNg_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSNg_2d_anom = tmp-tmpS;

tmp=advTWa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),ny,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTWa_2d_anom = tmp-tmpS;

tmp=advSWa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),ny,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSWa_2d_anom = tmp-tmpS;

tmp=advTSa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTSa_2d_anom = tmp-tmpS;

tmp=advSSa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSSa_2d_anom = tmp-tmpS;

tmp=advTNa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advTNa_2d_anom = tmp-tmpS;

tmp=advSNa_2d; tmp(:,:,365+365+31+29)=[]; 
tmpS = mean(reshape(tmp(:,:,1:365*5),nx,nz,365,5),4);
tmpS=repmat(tmpS,[1 1 5]); tmpS(:,:,end+1:365*5+59)=tmpS(:,:,1:59);
advSNa_2d_anom = tmp-tmpS;


% remove leap day
time1=time; time1(365*2+31+29)=[]; 


% monthly means
ndays=repmat([31 28 31 30 31 30 31 31 30 31 30 31],[1 5]);

ind=1;
for n=1:length(ndays)
	 advTWg_2d_monthly(:,:,n) = mean(advTWg_2d(:,:,ind:ind+ndays(n)),3);
	 advTWa_2d_monthly(:,:,n) = mean(advTWa_2d(:,:,ind:ind+ndays(n)),3);
	 advTNg_2d_monthly(:,:,n) = mean(advTNg_2d(:,:,ind:ind+ndays(n)),3);
	 advTNa_2d_monthly(:,:,n) = mean(advTNa_2d(:,:,ind:ind+ndays(n)),3);
	 advTSg_2d_monthly(:,:,n) = mean(advTSg_2d(:,:,ind:ind+ndays(n)),3);
	 advTSa_2d_monthly(:,:,n) = mean(advTSa_2d(:,:,ind:ind+ndays(n)),3);
	 advSWg_2d_monthly(:,:,n) = mean(advSWg_2d(:,:,ind:ind+ndays(n)),3);
	 advSWa_2d_monthly(:,:,n) = mean(advSWa_2d(:,:,ind:ind+ndays(n)),3);
	 advSNg_2d_monthly(:,:,n) = mean(advSNg_2d(:,:,ind:ind+ndays(n)),3);
	 advSNa_2d_monthly(:,:,n) = mean(advSNa_2d(:,:,ind:ind+ndays(n)),3);
	 advSSg_2d_monthly(:,:,n) = mean(advSSg_2d(:,:,ind:ind+ndays(n)),3);
	 advSSa_2d_monthly(:,:,n) = mean(advSSa_2d(:,:,ind:ind+ndays(n)),3);
	 advTWg_2d_anom_monthly(:,:,n) = mean(advTWg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advTWa_2d_anom_monthly(:,:,n) = mean(advTWa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advTNg_2d_anom_monthly(:,:,n) = mean(advTNg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advTNa_2d_anom_monthly(:,:,n) = mean(advTNa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advTSg_2d_anom_monthly(:,:,n) = mean(advTSg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advTSa_2d_anom_monthly(:,:,n) = mean(advTSa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSWg_2d_anom_monthly(:,:,n) = mean(advSWg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSWa_2d_anom_monthly(:,:,n) = mean(advSWa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSNg_2d_anom_monthly(:,:,n) = mean(advSNg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSNa_2d_anom_monthly(:,:,n) = mean(advSNa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSSg_2d_anom_monthly(:,:,n) = mean(advSSg_2d_anom(:,:,ind:ind+ndays(n)),3);
	 advSSa_2d_anom_monthly(:,:,n) = mean(advSSa_2d_anom(:,:,ind:ind+ndays(n)),3);
	 ind=ind+ndays(n);
end


% 30-day running mean
advTSg_2d_sm = NaN*advTSg_2d;
advTSa_2d_sm = NaN*advTSa_2d;
advTNg_2d_sm = NaN*advTNg_2d;
advTNa_2d_sm = NaN*advTNa_2d;
advTWg_2d_sm = NaN*advTWg_2d;
advTWa_2d_sm = NaN*advTWa_2d;
advSSg_2d_sm = NaN*advSSg_2d;
advSSa_2d_sm = NaN*advSSa_2d;
advSNg_2d_sm = NaN*advSNg_2d;
advSNa_2d_sm = NaN*advSNa_2d;
advSWg_2d_sm = NaN*advSWg_2d;
advSWa_2d_sm = NaN*advSWa_2d;
for t=16:nt-15
	 advTSg_2d_sm(:,:,t)=mean(advTSg_2d(:,:,t-15:t+15),3);
	 advTSa_2d_sm(:,:,t)=mean(advTSa_2d(:,:,t-15:t+15),3);
	 advTNg_2d_sm(:,:,t)=mean(advTNg_2d(:,:,t-15:t+15),3);
	 advTNa_2d_sm(:,:,t)=mean(advTNa_2d(:,:,t-15:t+15),3);
	 advTWg_2d_sm(:,:,t)=mean(advTWg_2d(:,:,t-15:t+15),3);
	 advTWa_2d_sm(:,:,t)=mean(advTWa_2d(:,:,t-15:t+15),3);
	 advSSg_2d_sm(:,:,t)=mean(advSSg_2d(:,:,t-15:t+15),3);
	 advSSa_2d_sm(:,:,t)=mean(advSSa_2d(:,:,t-15:t+15),3);
	 advSNg_2d_sm(:,:,t)=mean(advSNg_2d(:,:,t-15:t+15),3);
	 advSNa_2d_sm(:,:,t)=mean(advSNa_2d(:,:,t-15:t+15),3);
	 advSWg_2d_sm(:,:,t)=mean(advSWg_2d(:,:,t-15:t+15),3);
	 advSWa_2d_sm(:,:,t)=mean(advSWa_2d(:,:,t-15:t+15),3);
end


% 30-day running mean
advTSg_2d_anom_sm = NaN*advTSg_2d_anom;
advTSa_2d_anom_sm = NaN*advTSa_2d_anom;
advTNg_2d_anom_sm = NaN*advTNg_2d_anom;
advTNa_2d_anom_sm = NaN*advTNa_2d_anom;
advTWg_2d_anom_sm = NaN*advTWg_2d_anom;
advTWa_2d_anom_sm = NaN*advTWa_2d_anom;
advSSg_2d_anom_sm = NaN*advSSg_2d_anom;
advSSa_2d_anom_sm = NaN*advSSa_2d_anom;
advSNg_2d_anom_sm = NaN*advSNg_2d_anom;
advSNa_2d_anom_sm = NaN*advSNa_2d_anom;
advSWg_2d_anom_sm = NaN*advSWg_2d_anom;
advSWa_2d_anom_sm = NaN*advSWa_2d_anom;
for t=16:nt-1-15
	 advTSg_2d_anom_sm(:,:,t)=mean(advTSg_2d_anom(:,:,t-15:t+15),3);
	 advTSa_2d_anom_sm(:,:,t)=mean(advTSa_2d_anom(:,:,t-15:t+15),3);
	 advTNg_2d_anom_sm(:,:,t)=mean(advTNg_2d_anom(:,:,t-15:t+15),3);
	 advTNa_2d_anom_sm(:,:,t)=mean(advTNa_2d_anom(:,:,t-15:t+15),3);
	 advTWg_2d_anom_sm(:,:,t)=mean(advTWg_2d_anom(:,:,t-15:t+15),3);
	 advTWa_2d_anom_sm(:,:,t)=mean(advTWa_2d_anom(:,:,t-15:t+15),3);
	 advSSg_2d_anom_sm(:,:,t)=mean(advSSg_2d_anom(:,:,t-15:t+15),3);
	 advSSa_2d_anom_sm(:,:,t)=mean(advSSa_2d_anom(:,:,t-15:t+15),3);
	 advSNg_2d_anom_sm(:,:,t)=mean(advSNg_2d_anom(:,:,t-15:t+15),3);
	 advSNa_2d_anom_sm(:,:,t)=mean(advSNa_2d_anom(:,:,t-15:t+15),3);
	 advSWg_2d_anom_sm(:,:,t)=mean(advSWg_2d_anom(:,:,t-15:t+15),3);
	 advSWa_2d_anom_sm(:,:,t)=mean(advSWa_2d_anom(:,:,t-15:t+15),3);
end

cd /data/averdy/tpose/budgets/
save geostr_budgets_faces.mat *2d*


