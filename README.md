# MATPOWER-FACTS-PROGRAMS
I hereby attach the list matpower software which also includes several programs developed by me during final year project to calculate the optimal size and location of FACTS device in power system network.
APPENDIX A

It contains the bus data of the buses used in programming.

Bus Data Format :

1.	bus number (positive integer)
2.	bus type
a.	PQ bus          = 1
b.	PV bus          = 2
c.	reference bus   = 3
d.	isolated bus    = 4
3.	Pd, real power demand (MW)
4.	Qd, reactive power demand (MVAr)
5.	Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
6.	Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
7.	Vm, voltage magnitude (p.u.)
8.	Va, voltage angle (degrees)
9.	baseKV, base voltage (kV)
10.	maxVm, maximum voltage magnitude (p.u.)
11.	minVm, minimum voltage magnitude (p.u.)
12.	base MVA= 100.








1. IEEE 5 bus data

Bus no	Type	Pd	Qd	Gs	Bs	Vm	Va	baseKV	Vmax	Vmin
1	2	0	0	0	0	1	0	230	1.1	0.9
2	1	300	98.61	0	0	1	0	230	1.1	0.9
3	2	300	98.61	0	0	1	0	230	1.1	0.9
4	2	400	131.47	0	0	1	0	230	1.1	0.9
5	2	0	0	0	0	1	0	230	1.1	0.9


2. IEEE 9 bus data

Bus no	Type	Pd	Qd	Gs	Bs	Vm	Va	baseKV	Vmax	Vmin
1	3	0	0	0	0	1	0	345	1.1	0.9
2	2	0	0	0	0	1	0	345	1.1	0.9
3	2	0	0	0	0	1	0	345	1.1	0.9
4	1	0	0	0	0	1	0	345	1.1	0.9
5	1	90	30	0	0	1	0	345	1.1	0.9
6	1	0	0	0	0	1	0	345	1.1	0.9
7	1	100	35	0	0	1	0	345	1.1	0.9
8	1	0	0	0	0	1	0	345	1.1	0.9
9	1	125	50	0	0	1	0	345	1.1	0.9









APPENDIX B

1. Standard 5 bus system ’CASE 5’

function mpc = case5
%CASE5  Power flow data for modified 5 bus, 5 gen case based on PJM 5-bus system
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
%     Proceedings of the 2010 IEEE Power & Energy Society General Meeting
 
%   Created by Rui Bo in 2006, modified in 2010, 2014.
%   Distributed with permission.
 
%   MATPOWER
 
%% MATPOWER Case Format : Version 2
mpc.version = '2';
 
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;
 
%% bus data
%   bus_i   type    Pd  Qd  Gs  Bs  area    Vm  Va  baseKV  zone    Vmax    Vmin
mpc.bus = [
    1   2   0   0   0   0   1   1   0   230 1   1.1 0.9;
    2   1   300 98.61   0   0   1   1   0   230 1   1.1 0.9;
    3   2   300 98.61   0   0   1   1   0   230 1   1.1 0.9;
    4   3   400 131.47  0   0   1   1   0   230 1   1.1 0.9;
    5   2   0   0   0   0   1   1   0   230 1   1.1 0.9;
];
 
%% generator data
%   bus Pg  Qg  Qmax    Qmin    Vg  mBase   status  Pmax    Pmin    Pc1 Pc2 Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc    ramp_10 ramp_30 ramp_q  apf
mpc.gen = [
    1   40  0   30  -30 1   100 1   40  0   0   0   0   0   0   0   0   0   0   0   0;
    1   170 0   127.5   -127.5  1   100 1   170 0   0   0   0   0   0   0   0   0   0   0   0;
    3   323.49  0   390 -390    1   100 1   520 0   0   0   0   0   0   0   0   0   0   0   0;
    4   0   0   150 -150    1   100 1   200 0   0   0   0   0   0   0   0   0   0   0   0;
    5   466.51  0   450 -450    1   100 1   600 0   0   0   0   0   0   0   0   0   0   0   0;
];
 
%% branch data
%   fbus    tbus    r   x   b   rateA   rateB   rateC   ratio   angle   status  angmin  angmax
mpc.branch = [
    1   2   0.00281 0.0281  0.00712 400 400 400 0   0   1   -360    360;
    1   4   0.00304 0.0304  0.00658 0   0   0   0   0   1   -360    360;
    1   5   0.00064 0.0064  0.03126 0   0   0   0   0   1   -360    360;
    2   3   0.00108 0.0108  0.01852 0   0   0   0   0   1   -360    360;
    3   4   0.00297 0.0297  0.00674 0   0   0   0   0   1   -360    360;
    4   5   0.00297 0.0297  0.00674 240 240 240 0   0   1   -360    360;
];
 
%%-----  OPF Data  -----%%
%% generator cost data
%   1   startup shutdown    n   x1  y1  ... xn  yn
%   2   startup shutdown    n   c(n-1)  ... c0
mpc.gencost = [
    2   0   0   2   14  0;
    2   0   0   2   15  0;
    2   0   0   2   30  0;
    2   0   0   2   40  0;
    2   0   0   2   10  0;
];












2. Standard 5 bus system ’CASE 9’

function mpc = case9
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.
 
%   MATPOWER
 
%% MATPOWER Case Format : Version 2
mpc.version = '2';
 
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;
 
%% bus data
%   bus_i   type    Pd  Qd  Gs  Bs  area    Vm  Va  baseKV  zone    Vmax    Vmin
mpc.bus = [
    1   3   0   0   0   0   1   1   0   345 1   1.1 0.9;
    2   2   0   0   0   0   1   1   0   345 1   1.1 0.9;
    3   2   0   0   0   0   1   1   0   345 1   1.1 0.9;
    4   1   0   0   0   0   1   1   0   345 1   1.1 0.9;
    5   1   90  30  0   0   1   1   0   345 1   1.1 0.9;
    6   1   0   0   0   0   1   1   0   345 1   1.1 0.9;
    7   1   100 35  0   0   1   1   0   345 1   1.1 0.9;
    8   1   0   0   0   0   1   1   0   345 1   1.1 0.9;
    9   1   125 50  0   0   1   1   0   345 1   1.1 0.9;
];
 
%% generator data
%   bus Pg  Qg  Qmax    Qmin    Vg  mBase   status  Pmax    Pmin    Pc1 Pc2 Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc    ramp_10 ramp_30 ramp_q  apf
mpc.gen = [
    1   0   0   300 -300    1   100 1   250 10  0   0   0   0   0   0   0   0   0   0   0;
    2   163 0   300 -300    1   100 1   300 10  0   0   0   0   0   0   0   0   0   0   0;
    3   85  0   300 -300    1   100 1   270 10  0   0   0   0   0   0   0   0   0   0   0;
];
 
%% branch data
%   fbus    tbus    r   x   b   rateA   rateB   rateC   ratio   angle   status  angmin  angmax
mpc.branch = [
    1   4   0   0.0576  0   250 250 250 0   0   1   -360    360;
    4   5   0.017   0.092   0.158   250 250 250 0   0   1   -360    360;
    5   6   0.039   0.17    0.358   150 150 150 0   0   1   -360    360;
    3   6   0   0.0586  0   300 300 300 0   0   1   -360    360;
    6   7   0.0119  0.1008  0.209   150 150 150 0   0   1   -360    360;
    7   8   0.0085  0.072   0.149   250 250 250 0   0   1   -360    360;
    8   2   0   0.0625  0   250 250 250 0   0   1   -360    360;
    8   9   0.032   0.161   0.306   250 250 250 0   0   1   -360    360;
    9   4   0.01    0.085   0.176   250 250 250 0   0   1   -360    360;
];
 
%%-----  OPF Data  -----%%
%% generator cost data
%   1   startup shutdown    n   x1  y1  ... xn  yn
%   2   startup shutdown    n   c(n-1)  ... c0
mpc.gencost = [
    2   1500    0   3   0.11    5   150;
    2   2000    0   3   0.085   1.2 600;
    2   3000    0   3   0.1225  1   335;
];

































3. To calculate best rated SVC at a particular bus


function [minLoss,rat]=SVC(cs,busno,mxrat)
%% function to automate best gen SVC rating for min losses at particular bus
%syntax [minLoss,rat]=SVC('case file',busno,maximum rating of SVC available)
 
st=loadcase(cs);
r=runopf(st);         
loss=get_losses(r);         %losses with no SVC
m=size(loss);
tLoss=0;
p=0;
 
for i=1:m                 %adding total real losses with no SVC
    tLoss=tLoss+real(loss(i));  
end
W=tLoss;
st1=loadcase(st);
for j=1:mxrat                              %SVC rated upto mxrat MVAr.
    st1.bus(busno,4)=st.bus(busno,4)-j;     % changing reactive power supplied at bus due to SVC
                                                                %of rating j
    r=runopf(st1);
    loss=get_losses(r);                     % loss is a matrix of losses
    tLoss=0;
    for i=1:m           %adding total real losses with SVC of j rating
        tLoss=tLoss+real(loss(i));
    end
    if W>tLoss                                            % comparing to get value of minimal losses
        W=tLoss;
        p=j;
    end
end
minLoss=W;
rat=p;
%fprintf('\n w=%f     rating=%d \n',W,p);
end










4. To calculate overall best SVC



function [minLoss,rat,pos]=SVC1(cs,mxrat)
%% function to find overall best  SVC the network
%   syntax  [minLoss,bestRating,position]=SVC1('casefile',max rat. of SVC available)
 
st=loadcase(cs);
[W,R]=SVC(st,1,mxrat);  
ar(1,1)=W;
ar(1,2)=R;
 
n=size(st.bus);    % size of bus system
for i=2:n                     % making array of losses of all buses
    [l,r]=SVC(st,i,mxrat);
    ar(i,1)=l;
    ar(i,2)=r;
end
 
mn=ar(1,1);
pos=1;
 
for i=2:n                       %finding minimum losses in whole array
    if ar(i,1)<mn
        mn=ar(i,1);
        pos=i;
    end
end
     
    minLoss=mn;           %returning values to function parameters   
    rat=ar(pos,2);
   
   % printing whole loss data array
  fprintf(' Bus no          Line Loss          SVC used')
   for i=1:n
       fprintf('\n   %d           %fMW       %fMvar \n',i,ar(i,1),ar(i,2));
   end
   
    fprintf(' \n best SVC for min branch losses is of %d MVAr at bus %d with line loss %d MW \n',rat,pos,minLoss);
end





5. To find best voltage profile

function [v]=bestVoltageProfileSVC(cs,mxrat,mode)
%% bus voltage profile for best used SVC for both modes
%syntax [matrix]=bestVoltageProfileSVC(casefile,max rate SVC,mode(g/m))
if mode=='g'
    [minLoss,rat,pos]=SVC1(cs,mxrat); % calculating best SVC for gen mode
else 
    [minLoss,rat,pos]=mSVC1(cs,mxrat); % calculating best SVC for mot mode
end
    
st=loadcase(cs);
[m,n]=size(st.bus); % calculating bus size
 
%% changing bus data at position calculated for best SVC
if mode=='g'
    st.bus(pos,4)=st.bus(pos,4)-rat;   % generating mode
else
    st.bus(pos,4)=st.bus(pos,4)+rat;   % motoring mode
end
    
results=runopf(st); % running opf for best rated SVC
 
 
for i=1:m
    v(i,1)=i;     %storing bus no
    v(i,2)=results.bus(i,8); % storing bus voltage mag of ith bus .V Magnitude data stored in 8th column of result matrix
    v(i,3)=results.bus(i,9);% storing bus voltageangle  of ith bus .V Magnitude data stored in 9th column of result matrix
end
%% printing the matrix m
fprintf('\nbest SVC for min branch losses is of %d MVAr at bus %d with line loss %.6f MW \n',rat,pos,minLoss);
 
fprintf('bus no      Voltage Mag    Voltage Angle(degree)\n');
for i=1:m
    fprintf('%.0f              %f          %f\n',v(i,1),v(i,2),v(i,3));
end
%% calling function to plot the voltage profile
pl(v);
end









6. Best voltage profile for double SVC system in generating mode


function [v]=bestVoltageProfileSVC2(cs,mRat)
%% syntax [voltage matrix]=bestVoltageProfileSVC2(cs,mxrat1,mxrat2,mode)
% plots voltage profile with best double SVC system in generating mode
[ml,rat1,rat2,ps1,ps2]=dSVC2(cs,mRat); % calling function for doube SVC
st=loadcase(cs);
[m,n]=size(st.bus);
% changing ratings of SVC
st.bus(ps1,4)=st.bus(ps1,4)-rat1;
st.bus(ps2,4)=st.bus(ps2,4)-rat2;
 
results=runopf(st); %running opf with SVC installed system
for i=1:m
    v(i,1)=i;     % storing bus no
    v(i,2)=results.bus(i,8); % storing bus voltage mag of ith bus .V Magnitude data stored in 8th column of result matrix
    v(i,3)=results.bus(i,9);% storing bus voltageangle  of ith bus .V Magnitude data stored in 9th column of result matrix
end
fprintf('\n best installation of double SVC system if \n \n');
fprintf('SVC of rating %f MVAr at bus  %f and of %f MVAr at bus %f \n',rat1,ps1,rat2,ps2);
fprintf('minimal losses=%.6f MW',ml);
fprintf('bus no      Voltage Mag    Voltage Angle(degree)\n');
 
%% printing voltage profile matrix
for i=1:m
    fprintf('%.0f              %f          %f\n',v(i,1),v(i,2),v(i,3));
end
%% calling function to plot the voltage profile
pl(v);
end










7. Calculate best voltage profile and plot graph

function p=pl(v)
%% function to plot the voltage profile stored in matrix v
[m,n]=size(v);
for i=1:m
    plot(v(i,2),v(i,3));
    hold on;
    text(v(i,2),v(i,3),num2str(i));
end
plot(v(1:m,2),v(1:m,3),'rx--');
xlabel('bus voltage angle in degrees');
ylabel('bus voltage magnitude in pu');
title('Graph of voltage profile with best SVC used');
hold off;
grid;
end
















8. SVC in motoring mode

function [minLoss,rat]=mSVC(cs,busno,mxrat)
%% function to automate best motoring SVC rating for min losses at particular bus
%syntax [minLoss,rat]=SVC('case file',busno,maximum rating of SVC available)
 
st=loadcase(cs);
r=runopf(st);         
loss=get_losses(r);         %losses with no SVC
m=size(loss);
tLoss=0;
p=0;
 
for i=1:m                 %adding total real losses with no SVC
    tLoss=tLoss+real(loss(i));  
end
W=tLoss;
st1=loadcase(st);
for j=1:mxrat                              %SVC rated upto mxrat MVAr.
    st1.bus(busno,4)=st.bus(busno,4)+j;     % changing reactive power supplied at bus due to SVC in motoring mode
    r=runopf(st1);
    loss=get_losses(r);                     % loss is a matrix of losses
    tLoss=0;
    for i=1:m           %adding total real losses with SVC of j rating
        tLoss=tLoss+real(loss(i));
    end
    if W>tLoss                                            % comparing to get value of minimal losses
        W=tLoss;
        p=j;
    end
end
minLoss=W;
rat=p;
%fprintf('\n w=%f     rating=%d \n',W,p);
end






9. Overall best motoring SVC

function [minLoss,rat,pos]=mSVC1(cs,mxrat)
%% function to find overall best motoring SVC the network
%   syntax  [minLoss,bestRating,position]=SVC1('casefile',max rat. of SVC available)
 
st=loadcase(cs);
[W,R]=mSVC(st,1,mxrat);  
ar(1,1)=W;
ar(1,2)=R;
 
n=size(st.bus);    % size of bus system
for i=2:n                     % making array of losses of all buses
    [l,r]=mSVC(st,i,mxrat);
    ar(i,1)=l;
    ar(i,2)=r;
end
 
mn=ar(1,1);
pos=1;
 
for i=2:n                       %finding minimum losses in whole array
    if ar(i,1)<mn
        mn=ar(i,1);
        pos=i;
    end
end
     
    minLoss=mn;           %returning values to function parameters   
    rat=ar(pos,2);
   
   % printing whole loss data array
  fprintf(' Bus no          Line Loss          SVC used')
   for i=1:n
       fprintf('\n   %d           %fMW       %fMvar \n',i,ar(i,1),ar(i,2));
   end
   
    fprintf(' \n best SVC for min branch losses is of %d MVAr at bus %d with line loss %d MW \n',rat,pos,minLoss);
end





10. Calculating best SVC for a double SVC system

function [ml,rat1,rat2]=dSVC(cs,pos1,pos2,mRat)
%% function for double SVC installation in gen mode
%Syntax [minloss1,minloss2,rating1,rating2]=(busno1,busno2,mxarating)
 
st=loadcase(cs);
p=runopf(st);
tL=get_losses(p);        %losses matrix without SVC
TL=real(sum(tL));        %total real losses
st1=loadcase(cs);
ml=TL;
rat1=0;rat2=0;
for i=1:mRat
    st1.bus(pos1,4)=st.bus(pos1,4)-i;
    for j=1:mRat
        st1.bus(pos2,4)=st.bus(pos2,4)-j;
        tL1=real(sum(get_losses(runopf(st1))));   % getting real losses with SVC at pos1 of ith rating 
                                                                        %and at pos2 of jth rating
        if ml>tL1             %  updating minimum value of losses in branches with SVC of i and j rating in line
            ml=tL1;
            rat1=i;
            rat2=j;
        end
    end
end
%fprintf('\n minimal losses of %f MW with SVC of rating %f MVAr and %f MVAr \n',ml,rat1,rat2);
end
        
    








11. Overall best SVC in generating mode for double bus system

function [ml,rat1,rat2,ps1,ps2]=dSVC2(cs,mRat)
%% Double SVC system in gen mode
%Syntax [mnLoss,ratOfSvc1,ratOfSvc2,svc1Pos,svc2Pos]=dSVC2(CASEfile,maxRatning)
%requires a lot of time
%complexity = NC2 * mRat * mRat where N= number of buses
st=loadcase(cs);
[ml,rat1,rat2]=dSVC(st,1,2,mRat); %storing initial variables for best SVC at pos 1 and pos 2
n=size(st.bus);    % bus size
ps1=1;
ps2=2;
for i=1:n
    for j=i+1:n
        [mL,raT1,raT2]=dSVC(st,i,j,mRat); %calculating minimal losses for SVC at i,j
        
        if ml>mL
            ml=mL;
            rat1=raT1;
            rat2=raT2;
            ps1=i;
            ps2=j;
        end
    end
end
 
fprintf('\n best installation of double SVC system if \n \n');
fprintf('SVC of rating %f MVAr at bus  %f and of %f MVAr at bus %f \n',rat1,ps1,rat2,ps2);
fprintf('minimal losses=%.6f MW',ml);
end
        








12. Double SVC in motoring mode

function [ml,rat1,rat2]=mdSVC(cs,pos1,pos2,mRat)
%% function for double SVC installation in motoring mode
%Syntax [minloss1,minloss2,rating1,rating2]=(busno1,busno2,mxarating)
 
st=loadcase(cs);
p=runopf(st);
tL=get_losses(p);        %losses matrix without SVC
TL=real(sum(tL));        %total real losses
st1=loadcase(cs);
ml=TL;
rat1=0;rat2=0;
for i=1:mRat
    st1.bus(pos1,4)=st.bus(pos1,4)+i;
    for j=1:mRat
        st1.bus(pos2,4)=st.bus(pos2,4)+j;
        tL1=real(sum(get_losses(runopf(st1))));   % getting real losses with SVC at pos1 of ith rating 
                                                                        %and at pos2 of jth rating
        if ml>tL1             %  updating minimum value of losses in branches with SVC of i and j rating in line
            ml=tL1;
            rat1=i;
            rat2=j;
        end
    end
end
%fprintf('\n minimal losses of %f MW with SVC of rating %f MVAr and %f MVAr \n',ml,rat1,rat2);
end
        
    








13. Overall best SVC in double SVC system in motoring mode

function [ml,rat1,rat2,ps1,ps2]=mdSVC2(cs,mRat)
%% Double SVC system in motoring mode
%Syntax [mnLoss,ratOfSvc1,ratOfSvc2,svc1Pos,svc2Pos]=dSVC2(CASEfile,maxRatning)
%requires a lot of time
%complexity = NC2 * mRat * mRat where N= number of buses
st=loadcase(cs);
[ml,rat1,rat2]=mdSVC(st,1,2,mRat); %storing initial variables for best SVC at pos 1 and pos 2
n=size(st.bus);    % bus size
ps1=1;
ps2=2;
for i=1:n
    for j=i+1:n
        [mL,raT1,raT2]=mdSVC(st,i,j,mRat); %calculating minimal losses for SVC at i,j
        
        if ml>mL
            ml=mL;
            rat1=raT1;
            rat2=raT2;
            ps1=i;
            ps2=j;
        end
    end
end
 
fprintf('\n best installation of double SVC system if \n \n');
fprintf('SVC of rating %f MVAr at bus  %f and of %f MVAr at bus %f \n',rat1,ps1,rat2,ps2);
fprintf('minimal losses=%.6f MW',ml);
end
        








14. Plotting all loss graphs in generating mode with SVC

function a=plotAllGraphs(cs,mxrat)
%% syntax plotAllGraphs(casefile,max rating svc)
% plots graphs in generating mode
st=loadcase(cs);
st1=loadcase(st); 
 m=size(st1.bus); % m is a matrix os 2*1
 p=m(1);
% calculating dimenstions of subplot
 if p==9      % p is the number of buses
    n1=3;n2=3;
elseif p==14
    n1=4;n2=4;
elseif p==5
    n1=3;n2=2;
elseif p==6
    n1=3;n2=2;
elseif p==4
    n1=2;n2=2;
elseif p==24
    n1=5;n2=5;
elseif p==30
    n1=6;n2=5;
elseif p==33
    n1=6;n2=6;
elseif p==39
    n1=7;n2=6;
end
        
 for i=1:p
     subplot(n1,n2,i);       %making a square subplots
     plotGraph(st,i,mxrat);     % plotting all graphs
     xlabel(['Reactive Power Injected by SVC MVAr in bus  ' num2str(i) ' ']);
     ylabel('Real Power Losses in MW');
 end 
 hold off;
end














15. Plotting loss graph for single bus


function ar=plotGraph(cs,busno,mxrat)
%% plotting graph
st=loadcase(cs);
 st1=loadcase(st); 
 m=size(st1.bus);
for j=1:mxrat                              %SVC rated upto mxrat MVAr.
    st1.bus(busno,4)=st.bus(busno,4)-j;     % changing reactive power supplied at bus due to SVC
                                                                %of rating j
    r=runopf(st1);
    loss=get_losses(r);                     % loss is a matrix of losses
    tLoss=0;
    for i=1:m           %adding total real losses with SVC of j rating
        tLoss=tLoss+real(loss(i));
    end
    ar(j)=tLoss;
end
%plotting graph
plot(ar);
grid;
xlabel('Reactive Power Injected by SVC MVAr');
ylabel('Real Power Loss');
end













16. Plotting loss graphs for Standard Bus system in motoring mode

function a=mplotAllGraphs(cs,mxrat)
%% syntax plotAllGraphs(casefile,max rating svc)
% plots graphs in motoring mode
st=loadcase(cs);
st1=loadcase(st); 
 m=size(st1.bus); % m is a matrix os 2*1
 p=m(1);
% calculating dimenstions of subplot
 if p==9      % p is the number of buses
    n1=3;n2=3;
elseif p==14
    n1=4;n2=4;
elseif p==5
    n1=3;n2=2;
elseif p==6
    n1=3;n2=2;
elseif p==4
    n1=2;n2=2;
elseif p==24
    n1=5;n2=5;
elseif p==30
    n1=6;n2=5;
elseif p==33
    n1=6;n2=6;
elseif p==39
    n1=7;n2=6;
end
        
 for i=1:p
     subplot(n1,n2,i);       %making a square subplots
     mplotGraph(st,i,mxrat);     % plotting all graphs
     xlabel(['Reactive Power taken by SVC MVAr from bus  ' num2str(i) ' ']);
     ylabel('Real Power Loss in MW');
 end 
 hold off;
end







17. Plotting graph of losses at given bus in motoring mode

function ar=mplotGraph(cs,busno,mxrat)
%% plotting graph
st=loadcase(cs);
 st1=loadcase(st); 
 m=size(st1.bus);
for j=1:mxrat                              %SVC rated upto mxrat MVAr.
    st1.bus(busno,4)=st.bus(busno,4)+j;     % changing reactive power taken at bus due to SVC
                                                                %of rating j
    r=runopf(st1);
    loss=get_losses(r);                     % loss is a matrix of losses
    tLoss=0;
    for i=1:m           %adding total real losses with SVC of j rating
        tLoss=tLoss+real(loss(i));
    end
    ar(j)=tLoss;
end
%plotting graph
plot(ar);
grid;
xlabel('Reactive Power Injected by SVC MVAr');
ylabel('Real Power Loss');
end












18. Calculating bus power Factor

function pf=powerFactor(P,Q)
%% program to calculate PF using power triangle
a=P*P+Q*Q;
b=sqrt(a);
pf=P/b;
end




















19. Calculating bus power factor at best SVC

function [pf]=powerFactorAtBestSVC(cs,mxrat,mode)
%% load overall pf for best used SVC for both modes
%syntax [overall power factor]=powerFactorAtBestSVC(casefile,max rate SVC,mode(g/m))
if mode=='g'
    [minLoss,rat,pos]=SVC1(cs,mxrat); % calculating best SVC for gen mode
else 
    [minLoss,rat,pos]=mSVC1(cs,mxrat); % calculating best SVC for mot mode
end
    
st=loadcase(cs);
[m,n]=size(st.bus); % calculating bus size
 
%% changing bus data at position calculated for best SVC
if mode=='g'
    st.bus(pos,4)=st.bus(pos,4)-rat;   % generating mode
else
    st.bus(pos,4)=st.bus(pos,4)+rat;   % motoring mode
end
    
results=runopf(st); % running opf for best rated SVC
 
results1=runopf(cs); %without the SVC
%% claculating total active and reactive power at load and Overall Power Factor
P=0;P1=0;
Q=0;Q1=0;
for i=1:m
    P=P+results.bus(i,3);  % active power stored in 3rd column of load side
    Q=Q+results.bus(i,4); % reactive power stored in 4th column of load side
     P1=P1+results1.bus(i,3); 
    Q1=Q1+results1.bus(i,4);
end
pf=powerFactor(P,Q);
pf1=powerFactor(P1,Q1); %pf without SVC
 
fprintf('\nbest SVC for min branch losses is of %d MVAr at bus %d with line loss %.6f MW \n',rat,pos,minLoss);
fprintf('\n\n power factor without SVC is =%.4f and with SVC is=%.4f\n',pf1,pf);
end
