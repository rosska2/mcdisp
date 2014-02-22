clear;
% close all
% M = importdata('./ni3p_honeycomb/results/mcdisp.qei', ' ', 19);

M = importdata('./ni3p_honeycomb/results/mcdisp.qei', ' ', 21);
% M = importdata('./pyro2/results/mcdisp.qei', ' ', 21);


% M = importdata('./test_cubic/results/mcdisp.qei', ' ', 18);


% [Msort,is] = sort(M.data(:,6),1);

% ind = find(M.data(:,9)>0);



scale = 4*ones(size(M.data(:,5)));

%%%

% column 5,6,7
% figure;
% scatter(M.data(:,5),M.data(:,9),100,M.data(:,10),'s','filled');
% scatter3(M.data(:,5),M.data(:,6),M.data(:,9),scale,M.data(:,10),'s','filled');
% zlim([0,2])
% caxis([0 0.1])

% figure;
% scatter3(Msmall(:,5),Msmall(:,6),Msmall(:,7),10,Msmall(:,7),'filled');

% %%
% figure;
% scatter(Msmall(:,5),Msmall(:,9),200,Msmall(:,10),'s','filled');
% caxis([0,0.001])



%% binning

dh = 0.1;
dk = 0.1;
dl = 0.1;
dE = 0.1;

maxH = max(M.data(:,5));
minH = min(M.data(:,5));

maxK = max(M.data(:,6));
minK = min(M.data(:,6));

maxL = max(M.data(:,7));
minL = min(M.data(:,7));

maxE = max(M.data(:,9));
minE = min(M.data(:,9));

q = M.data(:,8);
H = M.data(:,5);
K = M.data(:,6);
L = M.data(:,7);
Edata = M.data(:,9);
Idata = M.data(:,10);

Ne = ceil((maxE-minE)/dE)+1; 
Nh = ceil((maxH-minH)/dh)+1;
Nk = ceil((maxK-minK)/dk)+1;
Nl = ceil((maxL-minL)/dl)+1;

I = zeros(Nh,Nk,Nl,Ne);
ns = zeros(Nh,Nk,Nl,Ne);
qs = zeros(Nh,Nk,Nl,Ne);
Iavg = zeros(Nh,Nk,Nl,Ne);


for i = 1:size(H)
    nh = round((H(i)-minH)/dh)+1;
    nk = round((K(i)-minK)/dk)+1;
    nl = round((L(i)-minL)/dl)+1;

    ne = round((Edata(i) -minE)/dE)+1;
        if Idata(i) < 0
%         disp(num2str(ne))
        Idata(i) = 0;
    end
    I(nh,nk,nl,ne) = I(nh, nk, nl, ne) + Idata(i);
    qs(nh,nk,nl,ne) = qs(nh,nk,nl,ne) + q(i);
    ns(nh,nk,nl,ne) = ns(nh, nk, nl ,ne) + 1;
end


% [i,j,k,l] = find(ns);
%%

[i,j,k,l] = ind2sub(size(ns),find(ns));
for ii = 1:size(i,1)
    Iavg(i(ii),j(ii),k(ii),l(ii)) = I(i(ii),j(ii),k(ii),l(ii))./ns(i(ii),j(ii),k(ii),l(ii));
    Qavg(i(ii),j(ii),k(ii),l(ii)) = qs(i(ii),j(ii),k(ii),l(ii))./ns(i(ii),j(ii),k(ii),l(ii));
end




%% 


hh = minH + [0:Nh-1]*dh;
kk = minK + [0:Nk-1]*dk;
ll = minL + [0:Nl-1]*dl;
ee = minE + [0:Ne-1]*dE;

%% convolution with gaussian energy resolution

hhRes = [2*min(hh):mean(diff(hh)):-2*min(hh)];
kkRes = [2*min(kk):mean(diff(kk)):-2*min(kk)];
llRes = [2*min(ll):mean(diff(ll)):-2*min(ll)];
eeRes = [2*min(ee):mean(diff(ee)):2*max(ee)];

 RE = gaussian_area(1, 0.2, 0, eeRes);
 RH = gaussian_area(1, 0.3, 0, hhRes);
 RK = gaussian_area(1, 0.3, 0, kkRes);
 RL = gaussian_area(1, 0.1, 0, llRes);

% tic
% for ii = 1:Nh
%     for jj = 1:Nk
%         for mm = 1:Nl
%            
%                 R(ii,jj,mm,:) = gaussian_area(1,0.2,0,ee);
%            
%         end
%     end
% end
% toc

% tic  
% Iconv = convn(Iavg,R,'same');
% toc

tic
Iconv1 = convn(Iavg,RH','same')*dh;
Iconv2 = convn(Iconv1,RK,'same')*dk;
RL2 = reshape(RL,[1 1 size(RL,2)]);
Iconv3 = convn(Iconv2,RL2,'same')*dl;
RE2 = reshape(RE,[1 1 1 size(RE,2)]);
Iconv = convn(Iconv3,RE2,'same')*dE;
toc


%%

% specify ranges and find the indices appropriate to that range
rE = [0, 0.2];
rH = [-0.2, 0.2];
rK = [-0.2, 0.2];
rL = [-0.1, 0.1];

iE = intersect(find(ee<rE(2)),find(ee>rE(1)));
iH = intersect(find(hh<rH(2)),find(hh>rH(1)));
iK = intersect(find(kk<rK(2)),find(kk>rK(1)));
iL = intersect(find(ll<rL(2)),find(ll>rL(1)));



% add up
I3D = squeeze(sum(Iconv(iH,:,:,:),1));
I3D_H = squeeze(sum(Iconv(:,iK,:,:),2));
Q3D = squeeze(sum(Qavg(iH,:,:,:),1)); 
% I2D = squeeze(sum(I3D(20:25,:,:), 1));
I2D = squeeze(sum(I3D(:,:,iE), 3));
% Q2D = squeeze(Q3D(:,:,round(mean(iE)))); 

I2D_HE = squeeze(sum(I3D_H(:,iL,:),2));

I2D_LE = squeeze(sum(I3D(iK,:,:), 1));
% Q2D_LE = squeeze(Q3D(round(mean(iK)),:,:)); 

I2D_KE = squeeze(sum(I3D(:,iL,:), 2));
% Q2D_KE = squeeze(Q3D(:,round(mean(iL)),:)); 



% [HH,LL] = meshgrid(hh,ll);
[KK,LL] = meshgrid(kk,ll);


figure;
h=pcolor(KK,LL,I2D');
set(h,'EdgeColor','none')
caxis([0 0.5])
xlabel('K');
ylabel('L');
title(num2str(rE));


[HH,EE] = meshgrid(hh,ee);
figure;
h=pcolor(HH,EE,I2D_HE');
set(h,'EdgeColor','none')
caxis([0 0.5])
xlabel('K');


[LL,EE] = meshgrid(ll,ee);
figure;
h=pcolor(LL,EE,I2D_LE');
set(h,'EdgeColor','none')
caxis([0 0.5])
xlabel('L');


[KK,EE] = meshgrid(kk,ee);
figure;
h=pcolor(KK,EE,I2D_KE');
set(h,'EdgeColor','none')
caxis([0 0.5])
xlabel('K');
title(num2str(rL));

%%
% sliceomatic(I3D(:,:,:), hh, ll, ee)
sliceomatic(I3D(:,:,:))

caxis([0 0.5])
% zlim([0 3])
%% try orthogonalslicer

% xvec = sort(hh,2,'descend');
% yvec = sort(kk,2,'descend');
% zvec = sort(ee,2,'descend');
% 
% xvec = sort(hh);
% yvec = sort(kk);
% zvec = sort(ee);
% orthogonalslicer(xvec,yvec,zvec,I3D)

% sliceomatic(M.data(:,10),M.data(:,5),M.data(:,6),M.data(:,7),M.data(:,9))

%%
% tri = delaunay(M.data(:,5),M.data(:,6));
% figure;
% h = trisurf(tri,M.data(:,6),M.data(:,9),M.data(:,10));
% % axis vis3d
% set(h,'EdgeColor', 'none');
% shading flat
% view([0 90])
% 
% caxis([0,0.4])


%% now do Q average
% 
qrange = [0.1:0.05:2];

Ipow = zeros(size(qrange,2)-1, size(Iavg,4));

tic
for i = 1:size(qrange,2)-1
    t1 = find(Qavg>qrange(i));
    t2 = find(Qavg<qrange(i+1));
    linind = intersect(t1,t2);
    [u,v,w,x] = ind2sub(size(Qavg),linind);
    for j = 1:size(u,1)
        Ipow(i,:) = Ipow(i,:) + squeeze(Iavg(u(j),v(j),w(j),:))'; 
    end
    Ipow(i,:) = Ipow(i,:)./size(linind,1);
    Q(i) = (qrange(i) + qrange(i+1))/2;
end
toc


%%
[QQ,EE] = meshgrid(Q,ee);
figure;
h = pcolor(QQ',EE',Ipow);
set(h,'EdgeColor','none');
caxis([0 0.1])

%% cut
cut = intersect(find(ee>0),find(ee<2));
test = sum(Ipow(:,cut),2);
figure;
plot(Q,test)
