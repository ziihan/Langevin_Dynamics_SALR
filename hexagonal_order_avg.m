clear;
clc;
%%%Smoothing noisy data, using movmean() or movmedian(),
NPart=1024;
Str="A";%["A","B","C","D"];
num=50;%[50,38,20,10];
QQ=[3,5];
%counter=[10,12,14,16];


%200,228.57,320,457.14


%(1.103708472,1.119680565,1.129086288,1.135795837)
for np=1:length(num)
        q6sq_avg=[];  
  for qq=1:length(QQ)   

     a=num(np)
     b=QQ(qq)
         
        filedir=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/Epsilon%d/',Str(np),num(np),QQ(qq));
        filename=sprintf('q6Sq.dat');
        Files=dir(strcat(filedir,filename));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A0=load(strcat(filedir,Files(1).name));
        AA0=A0(2001:length(A0),1);
        q6sq_avg=[q6sq_avg;b,mean(AA0)];
%initialize the bond distribution function
 frame=0;  
  end
 % q6sq_avg=[0.0,0.0;q6sq_avg];
        filenameSave=sprintf('/Volumes/IBI4-ZTanA/ProteinDiffusion/DATA/Langevin_Q2D_SLAR_1024/%s_phi_dot%d/q6sq_avg_more.dat',Str(np),num(np));
     save(filenameSave,'-ascii','q6sq_avg');
end
