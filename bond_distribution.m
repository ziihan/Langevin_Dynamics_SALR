clear;
clc;
%%%Smoothing noisy data, using movmean() or movmedian(),
Str=["B","C","D"];
num=[38,20,10];
QQ=[10,14,16,20];
for np=1:length(num)
   someth=[];
  for f=1:length(QQ)
     a=num(np)
     b=QQ(f)
        filedir=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/Epsilon%d/',Str(np),num(np),QQ(f));
        Files=dir(strcat(filedir,'bond_dist_2.dat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A0=load(strcat(filedir,Files(1).name));
        
        BondDistr=[0,0;
    1,0;
    2,0;
    3,0;
    4,0;
    5,0;
    6,0];
FrameCounter=zeros(7,1);
%initialize the bond distribution function
       
     for j=1:length(A0)
    ind=mod(A0(j,1),7)+1;% 
    FrameCounter(ind)=FrameCounter(ind)+1;
    BondDistr(ind,2)=BondDistr(ind,2)+A0(j,2);
     end
     
     CounterMax=max(FrameCounter);
     BondDistr(:,2)=BondDistr(:,2)/CounterMax;
     
      filenameSave=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/Epsilon%d/BDF.dat',Str(np),num(np),QQ(f));
      save(filenameSave,'-ascii','BondDistr');
     A12=BondDistr(:,1).*BondDistr(:,2);
     someth=[someth;b,sum(A12)]
  end
        filenameSave=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/bond_avearge3und5.dat',Str(np),num(np));
      save(filenameSave,'-ascii','someth');
end