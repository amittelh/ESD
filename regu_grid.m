function [RMat] = regu_grid(dipPts,MinMLat,MaxMLat,MinMLon,MaxMLon)

% change obsPts dipole points

%Regularisation matrix grid 
% [X,rho,eta,F] = pcgls(A,L,W,b,k,reorth,sm)

  %addpath /Users/annamittelholz/ownCloud/MagRoutines/make_icosahedron

  %make a regularisation grid
RMARS = 3930;

% find midpoints of triangles, where f is the frequency of points
f = 5;

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

%subselect area depending on MaxMLat... to deacrease matrix size to
%applicable area
% iDipSel=el>=MinMLat & el<=MaxMLat & az>=MinMLon & az<=MaxMLon;
% mid_m = mid_m(iDipSel,:);

for ii = 1:length(mid_m);
    
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
 
%calculate the distance of these points %midpoint - closest midpoints
mid3 = [azpoint3,elpoint3,rpoint3];  %phi, theta
mid2 = [azpoint2,elpoint2,rpoint2];
mid1 = [azpoint1,elpoint1,rpoint1];

dist3 = [az, el,r] -  mid3;  %phi, theta
dist2 = [az, el,r] -  mid2;
dist1 = [az, el,r] -  mid1;

 %%
 %for each triangle center point iC and triangles nearest_xyz
 %matrix for each of these points
% ReguMat = zeros(3*length(obsPts), length(mid_m));
ll = 1;

thresh = 1500; %800;        % threshhold great circle distance (1500 for comparisons with MEP code)
dipDir = 'None';             % options are 'None', 'CAD', 'OAD', NB: equivDipolMatOptim.m allows fixed direction, but that is NOT implemented here yet.
z_off = 0;                  % set to zero unless choose OAD option in which case use e.g. 480km  
 

%reguMat is final mat deltaB in 3 directions for each midpoint stacked; so
%there are 3 rows per midpoint

for ii = 1:length(mid_m);
%matrix of midpoints
[Reg_MatR] = equivDipolMatOptim2(dipPts,mid_m(ii,:),1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off);
%matrix of nearest mid point distances 
[Reg_MatR_delta] = equivDipolMatOptim2( dipPts, [mid1(ii,:); mid2(ii,:); mid3(ii,:)] ,1,0,0, 'crit_dist',thresh,'fixed_dir',dipDir,'z_offset',z_off);
%this produces the delta B part with midpoints - near points      
ReguMat_Br = [Reg_MatR - Reg_MatR_delta(:,1:3); Reg_MatR - Reg_MatR_delta(:,4:6); Reg_MatR - Reg_MatR_delta(:,7:9)]  ;

r = [mid_m(ii,1); mid_m(ii,1); mid_m(ii,1)];
% dist# = [az, el,r] 
[alpha] = 1./r.*(1./[dist1(ii,2); dist2(ii,2); dist3(ii,2)]) + (1./[sin(dist1(ii,2)); sin(dist2(ii,2)); sin(dist3(ii,2))].*[dist1(ii,1); dist2(ii,1); dist3(ii,1)]);

fullalpha = [alpha(1)*ones(length(Reg_MatR),3); alpha(2)*ones(length(Reg_MatR),3); alpha(3)*ones(length(Reg_MatR),3)];
ReguMat_final(:,ll:ll+2) = fullalpha.*ReguMat_Br;
clear alpha fullalpha r ReguMat_Br Reg_MatR_delta Reg_MatR


ll = ll+3   ;
end
part1 = length(ReguMat_final)/3;
clear ii 
%%

%RMat = zeros(length(dipPts), max(size(ReguMat_final(1,:))));
tic
for ii = 1:part1
RMat(ii,:) = sum(ReguMat_final(ii:part1:end,:),1);% + ReguMat_final(ii+part1,:) + ReguMat_final(ii+(2*part1),:);
end
toc

%%
end
