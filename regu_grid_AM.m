function [RMat] = regu_grid_AM(dipPts,MinMLat,MaxMLat,MinMLon,MaxMLon)

% change obsPts dipole points

%Regularisation matrix grid 
% [X,rho,eta,F] = pcgls(A,L,W,b,k,reorth,sm)

  addpath /Users/annamittelholz/ownCloud/MagRoutines/make_icosahedron

%make a regularisation grid
RMARS = 3930;

% find midpoints of triangles, where f is the frequency of points
f = 20;

figure;
[x,y,z,TRI] = make_icosahedron(f,RMARS,1,1,1);
hold on

%TRI contains the different triangles TRI [index 1, index 2, index 3] give
% the lon and lat for each point --> Calculate midpoint for each set of 3
midpoint = triangulation(TRI, x',y',z');
iC = incenter(midpoint);

    % cols: x y z
figure;
plot3(iC(:,1),iC(:,2),iC(:,3),'ro');

%edge points in spherical coordinates radians
[aze,ele,re] = cart2sph(x,y,z);
%midpoints
[az, el, r] = cart2sph(iC(:,1), iC(:,2), iC(:,3));
%%
%Calculate distance by finding nearest neighour in spherical coordinates
%xs, ys, zs
mid_m = [az, el, r];
az(az<0) = az(az<0)+(2*pi);

%pick area
lonreg = rad2deg(az);
latreg = rad2deg(el);

%subselect area depending on MaxMLat... to deacrease matrix size to
%applicable area
clear  iDipSel
 iDipSel=latreg>=MinMLat & latreg<=MaxMLat & lonreg>=MinMLon & lonreg<=MaxMLon;
% mid_m = mid_m(iDipSel,:);
 
mid_points = [r(iDipSel,:) deg2rad(90-latreg(iDipSel,:)) deg2rad(lonreg(iDipSel,:))]; %[alt colat in radians lon in radians] this is the format needed for building the Gmatrix!!!
iC = iC(iDipSel,:); %midpoints in Cartesian
% only 


for ii = 1:length(mid_points);
        %in cartesian
    [nc,dc] = knnsearch(iC(:,1:3),iC(ii,1:3),'K',4);  %indices of the 3 nearest midpoints 2-4 are the right ones
    nearest = iC(nc(2:4),:)';
    
    nearest_xyz(ii,1:9) = nearest(:)';
    clear n d nearest
     
end

 %these takes the 2 closest points in phi direction and calculates the
 %spheircal coordinates
[azpoint1, elpoint1, rpoint1] = cart2sph(nearest_xyz(:,1), nearest_xyz(:,2), nearest_xyz(:,3));
[azpoint2, elpoint2, rpoint2] = cart2sph(nearest_xyz(:,4), nearest_xyz(:,5),nearest_xyz(:,6));
[azpoint3, elpoint3, rpoint3] = cart2sph(nearest_xyz(:,7), nearest_xyz(:,8),nearest_xyz(:,9));
 
azpoint1(azpoint1<0) = azpoint1(azpoint1<0)+(2*pi);
azpoint2(azpoint2<0) = azpoint2(azpoint2<0)+(2*pi);
azpoint3(azpoint3<0) = azpoint3(azpoint3<0)+(2*pi);
%calculate the distance of these points %midpoint - closest midpoints
mid3 = [rpoint3, pi/2-elpoint3, azpoint3];     %alt, colat in radphi, az (0:360);  %%format needed for building the Gmatrix!!!
mid2 = [rpoint2, pi/2-elpoint2, azpoint2];
mid1 = [rpoint1, pi/2-elpoint1, azpoint1];


%distance between midpoint and midpoints of 3 surrounding triangles

dist3 = sqrt((mid_points(:,2) -  mid3(:,2)).^2 + (mid_points(:,3) -  mid3(:,3)).^2);  % only in phi, theta direction as alt basically constant
dist2 = sqrt((mid_points(:,2) -  mid2(:,2)).^2 + (mid_points(:,3) -  mid2(:,3)).^2);
dist1 = sqrt((mid_points(:,2) -  mid1(:,2)).^2 + (mid_points(:,3) -  mid1(:,3)).^2);

%test
figure;
plot(rad2deg(mid_points(:,2)),rad2deg(mid_points(:,3)),'.')
hold on
plot(rad2deg(mid3(:,2)),rad2deg(mid3(:,3)),'ro')
hold on
plot(rad2deg(mid2(:,2)),rad2deg(mid2(:,3)),'ro')
hold on
plot(rad2deg(mid1(:,2)),rad2deg(mid1(:,3)),'ro')

figure;
plot(rad2deg(dipPts(:,2)),rad2deg(dipPts(:,3)),'.')


 %%
 %for each triangle center point iC and triangles nearest_xyz
 %matrix for each of these points
% ReguMat = zeros(3*length(obsPts), length(mid_m));

thresh = 1500; %800;       
dipDir = 'None';           
z_off = 0;                 
 
%reguMat is final mat deltaB in 3 directions for each midpoint stacked; so
%there are 3 rows per midpoint

%Note: need mid_m in [radius, t (in radians), p (in radians)]
[Reg_MatR] = equivDipolMatOptim2(dipPts,mid_points,1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off); %Gmat for gridpoints

%matrix of nearest mid point distances 
[Reg_MatR_delta1] = equivDipolMatOptim2( dipPts, mid1 ,1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off);
[Reg_MatR_delta2] = equivDipolMatOptim2( dipPts, mid2 ,1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off);
[Reg_MatR_delta3] = equivDipolMatOptim2( dipPts, mid3 ,1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off);

%this produces the delta B part with midpoints - near points      
ReguMat_Br = [Reg_MatR - Reg_MatR_delta1; Reg_MatR - Reg_MatR_delta2; Reg_MatR - Reg_MatR_delta3];

%alpha is the horizotal gradient 
r = [mid_points(:,1); mid_points(:,1); mid_points(:,1)];
alpha = 1./r.*(1./[dist1(:,2); dist2(:,2); dist3(:,2)]) + (1./[sin(dist1(:,2)); sin(dist2(:,2)); sin(dist3(:,2))].*[dist1(:,3); dist2(:,3); dist3(:,3)]);
% 
 size(alpha)
 size(ReguMat_Br)
%horizontal gradirnt of R component
RMat = alpha .* ReguMat_Br;
clear alpha fullalpha r ReguMat_Br Reg_MatR_delta1 Reg_MatR_delta2 Reg_MatR_delta3  Reg_MatR


end
