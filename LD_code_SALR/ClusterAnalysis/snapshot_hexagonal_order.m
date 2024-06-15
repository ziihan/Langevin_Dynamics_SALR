clear;
clc;
%%%Smoothing noisy data, using movmean() or movmedian(),
NPart=1024;
Str=["A","B","C","D"];
num=[50,38,20,10];
QQ=[10,14,16,20];
%counter=[10,12,14,16];

Lz=200;
Lxy=[200,228.57,320,457.14];

%lxby2 = 0.5*Lx; lyby2 = 0.5*Ly; lzby2 = 0.5*Lz;
skip=10;
sigma2=25.0;
sig_clust_sq=[sigma2*1.103708472^2,sigma2*1.119680565^2,sigma2*1.129086288^2,sigma2*1.135795837^2];

%(1.103708472,1.119680565,1.129086288,1.135795837)
for np=2:length(num)
    Lx=Lxy(np)
    Ly=Lxy(np)
    
    
  for qq=1:length(QQ)   
     q6sqAll=[];  q6ID=[];
     a=num(np)
     b=QQ(qq)
     
     %clust_sq=sig_clust_sq(qq);
     if((qq>2)||(np>2))
        filedir=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/Epsilon%d/',Str(np),num(np),QQ(qq));
        filename=sprintf('Mcom');
        Files=dir(strcat(filedir,filename));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A0=load(strcat(filedir,Files(1).name));
        NuFrame=length(A0)/NPart
%initialize the bond distribution function
 frame=0;
     for f=1:floor(NuFrame/skip)
         frame=frame+1;
         frame
         q6=zeros(NPart,1);
         %csd=A0(((f-1)*skip*NPart+1):((f-1)*skip+1)*NPart,1);
         Pos=A0(((f-1)*skip*NPart+1):((f-1)*skip+1)*NPart,1:2);
    
      for j=1:NPart
          dsq=[];
         for k=1:NPart
             dsx=Pos(k,1)-Pos(j,1);
             dsy=Pos(k,2)-Pos(j,2);
             
             dsx=dsx-Lx*round(dsx/Lx);%PBCs
             dsy=dsy-Ly*round(dsy/Ly);
             
             distsq=dsx*dsx+dsy*dsy;
           
             dsq=[dsq;distsq,k];
         end
         dsqsort=sortrows(dsq,1);
         
         neib_index=dsqsort(1:7,2);
         
         %find the six-nearest neigbours
         for jk=2:7
             dsx=Pos(neib_index(jk),1)-Pos(j,1);
             dsy=Pos(neib_index(jk),2)-Pos(j,2);
      
             dsx=dsx-Lx*round(dsx/Lx);
             dsy=dsy-Ly*round(dsy/Ly);
             
             neighbour_sq=dsx*dsx+dsy*dsy;
            % if(neighbour_sq<sig_clust_sq)
          %   if(neighbour_sq<1000000)
             alpha=atan2(dsy,dsx);
             q6(j)=q6(j)+exp(1i*6.*alpha);
           %  end
         end
             if(abs(q6(j)/6)>1)
                 q6(j)
             end       
      end
      q6=q6/6.;
      q6sq=(abs(q6)).^2;%%%local order parameter
     % q6sqAll=[q6sqAll;q6sq];%%
      %  q6ID=[q6ID;q6sq.'];%%%%output for VMD visualization
        
        %if f==499
        %   shot_csd=[Pos,csd];
           shot_q6sq=[Pos,q6sq]; 
        %end
    % filenameCsd=sprintf('/Volumes/LaCie/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/Epsilon%d/q6_movie/csdIDPos_shot%d.dat',Str(np),num(np),QQ(qq),f);
    % save(filenameCsd,'-ascii','shot_csd');    
     
          filenameSave=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/MoreAttractive/Epsilon%d/q6IDPos_shot%d.dat',Str(np),num(np),QQ(qq),f);
     save(filenameSave,'-ascii','shot_q6sq');   
     end
     end
%       filenameSave=sprintf('/Volumes/LaCie/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/Epsilon%d/q6Sq.dat',Str(np),num(np),QQ(qq));
%      save(filenameSave,'-ascii','q6sqAll');
%            filenameSave=sprintf('/Volumes/LaCie/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/Epsilon%d/q6_colorID.dat',Str(np),num(np),QQ(qq));
%      save(filenameSave,'-ascii','q6ID');

  end
end
