clc;clear;close all;format short;
% load('../historical/Nurbs_Nom')
nos=10000; %Number of Simulations
xres=50;
yres=65;
%Baseline Bracket Large Number of Control Points
ku_baseline=3;
nu_baseline=35;
tu_baseline=knot_maker(nu_baseline,ku_baseline);

kv_baseline=3;
nv_baseline=20;
tv_baseline=knot_maker(nv_baseline,kv_baseline);

y_cpoints_baseline=[linspace(0,.9,4) ones(1,10) linspace(1.1,2.9,8) 3*ones(1,10) linspace(3.1,4,4)];
z_cpoints_baseline=[zeros(1,5) linspace(.1,1.9,8) 2*ones(1,10) linspace(1.9,.1,8) zeros(1,5)];
x_cpoints_baseline=zeros(1,nu_baseline+1);

cpoints_baseline=[x_cpoints_baseline;y_cpoints_baseline;z_cpoints_baseline];
cpoints_baseline=repmat(cpoints_baseline,[1,1,nv_baseline+1]);
for i = 1:nv_baseline
    cpoints_baseline(1,:,i+1)=0.5625+cpoints_baseline(1,:,i);
end
cpoints_baseline=cpoints_baseline*25.4;
baseline = nrbmak_lee(cpoints_baseline,{tu_baseline, tv_baseline});
%
nrbeval(baseline,{[0:xres]/xres,[0:yres]/yres});
%

%NURBS for Monitoring
ku=3;
nu_1=30;
tu_1=knot_maker(nu_1,ku);

kv=3;
nv_1=15;
tv_1=knot_maker(nv_1,kv);

nu_2=25;
tu_2=knot_maker(nu_2,ku);

nv_2=10;
tv_2=knot_maker(nv_2,kv);

nu_3=20;
tu_3=knot_maker(nu_3,ku);

nv_3=5;
tv_3=knot_maker(nv_3,kv);

nu_4=10;
tu_4=knot_maker(nu_4,ku);

nv_4=5;
tv_4=knot_maker(nv_4,kv);

Nu_1=basis_fun(ku-1,nu_1,tu_1,linspace(0,1,xres+1));
Nu_1t=Nu_1;
Nu_1=inv(Nu_1'*Nu_1)*Nu_1';

Nv_1=basis_fun(kv-1,nv_1,tv_1,linspace(0,1,yres+1));
Nv_1t=Nv_1;
Nv_1=inv(Nv_1'*Nv_1)*Nv_1';

Nu_2=basis_fun(ku-1,nu_2,tu_2,linspace(0,1,xres+1));
Nu_2t=Nu_2;
Nu_2=inv(Nu_2'*Nu_2)*Nu_2';

Nv_2=basis_fun(kv-1,nv_2,tv_2,linspace(0,1,yres+1));
Nv_2t=Nv_2;
Nv_2=inv(Nv_2'*Nv_2)*Nv_2';

Nu_3=basis_fun(ku-1,nu_3,tu_3,linspace(0,1,xres+1));
Nu_3t=Nu_3;
Nu_3=inv(Nu_3'*Nu_3)*Nu_3';

Nv_3=basis_fun(kv-1,nv_3,tv_3,linspace(0,1,yres+1));
Nv_3t=Nv_3;
Nv_3=inv(Nv_3'*Nv_3)*Nv_3';

Nu_4=basis_fun(ku-1,nu_4,tu_4,linspace(0,1,xres+1));
Nu_4t=Nu_4;
Nu_4=inv(Nu_4'*Nu_4)*Nu_4';

Nv_4=basis_fun(kv-1,nv_4,tv_4,linspace(0,1,yres+1));
Nv_4t=Nv_4;
Nv_4=inv(Nv_4'*Nv_4)*Nv_4';

F1=zeros(3,nu_1+1,yres+1);
F2=zeros(3,nu_2+1,yres+1);
F3=zeros(3,nu_3+1,yres+1);
F4=zeros(3,nu_4+1,yres+1);
R1_1=zeros(1,nos);
R1_2=zeros(1,nos);
R1_3=zeros(1,nos);
R1_4=zeros(1,nos);
E1_Est=zeros(3,xres+1,yres+1);
E2_Est=zeros(3,xres+1,yres+1);
E3_Est=zeros(3,xres+1,yres+1);
E4_Est=zeros(3,xres+1,yres+1);

P1_a=zeros(1,nu_1+1);
P1_b=zeros(1,nu_1+1);
P1_c=zeros(1,nu_1+1);
P1t=zeros(3,nv_1+1,nu_1+1);
G1_Est=zeros(3,nu_1+1,yres+1);
R2_1=zeros(1,nos);

P2_a=zeros(1,nu_2+1);
P2_b=zeros(1,nu_2+1);
P2_c=zeros(1,nu_2+1);
P2t=zeros(3,nv_2+1,nu_2+1);
G2_Est=zeros(3,nu_2+1,yres+1);
R2_2=zeros(1,nos);

P3_a=zeros(1,nu_3+1);
P3_b=zeros(1,nu_3+1);
P3_c=zeros(1,nu_3+1);
P3t=zeros(3,nv_3+1,nu_3+1);
G3_Est=zeros(3,nu_3+1,yres+1);
R2_3=zeros(1,nos);

P4_a=zeros(1,nu_4+1);
P4_b=zeros(1,nu_4+1);
P4_c=zeros(1,nu_4+1);
P4t=zeros(3,nv_4+1,nu_4+1);
G4_Est=zeros(3,nu_4+1,yres+1);
R2_4=zeros(1,nos);
% %variables for nurb evaluation of baseline bracket
pu = zeros(3,1);                                                            
% %variables for nurb evaluation of current bracket
numu = nu_baseline+1; 
numv =nv_baseline+1;
degree = [ku_baseline-1 kv_baseline-1];
tv=tv_baseline;
tu=tu_baseline;

mcv=3*numu;
pv = zeros(mcv,1);                                                            
Nv = zeros(degree(2)+1,1);

Nu = zeros(degree(1)+1,1);
leftv = zeros(degree(2)+1,1);
rightv = zeros(degree(2)+1,1);
leftu = zeros(degree(1)+1,1);
rightu = zeros(degree(1)+1,1);
cpoints = zeros(3*numu,numv);
%variables for nurb evaluation of baseline bracket
pu = zeros(3,1);
% %%%%%
% %Sampling and optimization
%%%%
num_per_file=100;
counter=0;
counter1=1;
cpoints_baseline=reshape(cpoints_baseline,3*numu,numv);
load('../historical/hist1');
ndP1=zeros(nos,(nv_1+1)*(nu_1+1));
ndP2=zeros(nos,(nv_2+1)*(nu_2+1));
ndP3=zeros(nos,(nv_3+1)*(nu_3+1));
ndP4=zeros(nos,(nv_4+1)*(nu_4+1));

for k = 1:nos
    if counter == num_per_file
        counter1=counter1+1;
        counter=0;
        load_name=['../historical/hist' num2str(counter1)];
        load(load_name);
        k
    end
    counter=counter+1;
    cpoints = nrbmak_lee(reshape(cpoints_hist{counter},3,numu,numv),{tu_baseline, tv_baseline});
            p=nrbeval(cpoints,{[0:xres]/xres,[0:yres]/yres});
    for i = 1:nu_1+1
        F1(1,i,:)=Nu_1(i,:)*squeeze(p(1,:,:));
        F1(2,i,:)=Nu_1(i,:)*squeeze(p(2,:,:));
        F1(3,i,:)=Nu_1(i,:)*squeeze(p(3,:,:));
    end
    for i = 1:nu_2+1
        F2(1,i,:)=Nu_2(i,:)*squeeze(p(1,:,:));
        F2(2,i,:)=Nu_2(i,:)*squeeze(p(2,:,:));
        F2(3,i,:)=Nu_2(i,:)*squeeze(p(3,:,:));
    end
    for i = 1:nu_3+1
        F3(1,i,:)=Nu_3(i,:)*squeeze(p(1,:,:));
        F3(2,i,:)=Nu_3(i,:)*squeeze(p(2,:,:));
        F3(3,i,:)=Nu_3(i,:)*squeeze(p(3,:,:));
    end
    for i = 1:nu_4+1
        F4(1,i,:)=Nu_4(i,:)*squeeze(p(1,:,:));
        F4(2,i,:)=Nu_4(i,:)*squeeze(p(2,:,:));
        F4(3,i,:)=Nu_4(i,:)*squeeze(p(3,:,:));
    end
    for i = 1:nv_1+1
        P1t(1,i,:)=Nv_1(i,:)*squeeze(F1(1,:,:))';
        P1t(2,i,:)=Nv_1(i,:)*squeeze(F1(2,:,:))';
        P1t(3,i,:)=Nv_1(i,:)*squeeze(F1(3,:,:))';
    end
    for i = 1:3
        E1_Est(i,:,:)=Nu_1t*squeeze(P1t(i,:,:))'*Nv_1t'-squeeze(p(i,:,:));
        G1_Est(i,:,:)=squeeze(P1t(i,:,:))'*Nv_1t';
    end
    R1_1(k)=sum(sum(sum((E1_Est.^2))));
    R2_1(k)=sum(sum(sum((G1_Est-F1).^2)));
    for i = 1:nv_2+1
        P2t(1,i,:)=Nv_2(i,:)*squeeze(F2(1,:,:))';
        P2t(2,i,:)=Nv_2(i,:)*squeeze(F2(2,:,:))';
        P2t(3,i,:)=Nv_2(i,:)*squeeze(F2(3,:,:))';
    end
    for i = 1:3
        E2_Est(i,:,:)=Nu_2t*squeeze(P2t(i,:,:))'*Nv_2t'-squeeze(p(i,:,:));
        G2_Est(i,:,:)=squeeze(P2t(i,:,:))'*Nv_2t';
    end
    R1_2(k)=sum(sum(sum((E2_Est.^2))));
    R2_2(k)=sum(sum(sum((G2_Est-F2).^2)));
    for i = 1:nv_3+1
        P3t(1,i,:)=Nv_3(i,:)*squeeze(F3(1,:,:))';
        P3t(2,i,:)=Nv_3(i,:)*squeeze(F3(2,:,:))';
        P3t(3,i,:)=Nv_3(i,:)*squeeze(F3(3,:,:))';
    end
    for i = 1:3
        E3_Est(i,:,:)=Nu_3t*squeeze(P3t(i,:,:))'*Nv_3t'-squeeze(p(i,:,:));
        G3_Est(i,:,:)=squeeze(P3t(i,:,:))'*Nv_3t';
    end
    R1_3(k)=sum(sum(sum((E3_Est.^2))));
    R2_3(k)=sum(sum(sum((G3_Est-F3).^2)));
    for i = 1:nv_4+1
        P4t(1,i,:)=Nv_4(i,:)*squeeze(F4(1,:,:))';
        P4t(2,i,:)=Nv_4(i,:)*squeeze(F4(2,:,:))';
        P4t(3,i,:)=Nv_4(i,:)*squeeze(F4(3,:,:))';
    end
    for i = 1:3
        E4_Est(i,:,:)=Nu_4t*squeeze(P4t(i,:,:))'*Nv_4t'-squeeze(p(i,:,:));
        G4_Est(i,:,:)=squeeze(P4t(i,:,:))'*Nv_4t';
    end
    R1_4(k)=sum(sum(sum((E4_Est.^2))));
    R2_4(k)=sum(sum(sum((G4_Est-F4).^2)));
    dP1=P1t-P1t_nom;
    dP2=P2t-P2t_nom;
    dP3=P3t-P3t_nom;
    dP4=P4t-P4t_nom;
    ndP1(k,:)=reshape((squeeze(dP1(1,:,:)).^2+squeeze(dP1(2,:,:)).^2+squeeze(dP1(3,:,:)).^2).^.5,1,(nv_1+1)*(nu_1+1));
    ndP2(k,:)=reshape((squeeze(dP2(1,:,:)).^2+squeeze(dP2(2,:,:)).^2+squeeze(dP2(3,:,:)).^2).^.5,1,(nv_2+1)*(nu_2+1));
    ndP3(k,:)=reshape((squeeze(dP3(1,:,:)).^2+squeeze(dP3(2,:,:)).^2+squeeze(dP3(3,:,:)).^2).^.5,1,(nv_3+1)*(nu_3+1));
    ndP4(k,:)=reshape((squeeze(dP4(1,:,:)).^2+squeeze(dP4(2,:,:)).^2+squeeze(dP4(3,:,:)).^2).^.5,1,(nv_4+1)*(nu_4+1));
    ndP1c(k,:)=reshape(dP1,1,numel(dP1));
    ndP2c(k,:)=reshape(dP2,1,numel(dP2));
    ndP3c(k,:)=reshape(dP3,1,numel(dP3));
    ndP4c(k,:)=reshape(dP4,1,numel(dP4));
end
save('historical_nurbsc','ndP1c','ndP4c');