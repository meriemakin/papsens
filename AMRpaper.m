%Modeling of AMR on network of paper fibers sputtered with permalloy
%Author: Meriem Akin
%Started: July 14th 2015
clear
%tic
%number of random samples
nrs=1;
sav=[];
for countcount=1:nrs
%generate pathes
population=1000;
all_magR=[];
%distribution
dist = 1;
%0 if it is uniform
%1 if it is normal
%2 if it is exponential
%3 if it is rayleigh
%standard deviation of normal distribution
devnor = 0.3; %>0.28 optimal
%mean of exponential distribution
meanexp = 1; %>0.36 optimal value
% scale parameter of rayleigh
scrgh = 0.43015; %0.43015 optimal value
%hold angle or length to uniform distribution
%and vary the other only with a different distribution
unif=1; %0 dont do anything,
%1 angle should be uniform, 2 length should be uniform
figure
for count=1:population
    %tic
    %Define geometry shape
    %1: square, 2: stripe.
    geometry = 1;
    %1:square
    %    a
    %--------
    %|      |
    %|      |a
    %|      |
    %--------
    %2:stripe
    %                a
    %------------------------------
    %|                            |
    %|                            |b
    %------------------------------
    %Define geometry dimensions
    a = 1; % a is square side length, stripe long side length
    b = 0.01*a; % b is the stripe small side length

    %Define bonding box
    startX=0.0;
    startY=0.0;
    endX=a;
    endY=a;
    if geometry == 2
        endY=b;
    end

    %Measurement path
    measStartX = startX;   
    measStartY= startY+0.5*a;
    if geometry == 2
        measStartY = startY + 0.5*b;
    end
    measEndX = a;
    measEndY = startY + 0.5*a;
    if geometry == 2
        measEndY = startY + 0.5*b;
    end

    %Electrical resistance
    %R = rho * L/A, where L is the length of the fiber, A cross sectional area
    %of the permalloy coating, rho is resistivity
    %Assume that rho/A = 1 (a constant) along the fibers, and R a function of length only.

    %Generate a path and calculate individual electrical resistance
    randdis = 0.3*a; %can be varied to study effect of fiber length
    randangle1 = pi;
    randangle2 = 0.5*pi;

    %store angles and magnetoresistances
    angles=[];
    magRpaper=[];
    magRsmooth=[];
    listR = [];
    listanglecorr = [];
    %generate random path per sample
    generated=0;
    %Loop over angle between current and magnetic field
    for angle=-pi:pi/12:pi
        R_mag_total = 0;
        R_mag_smooth = 0;
        listR_mag = [];
        if (generated == 0)
            x_temp = measStartX;
            y_temp = measStartY;
            listXY = [x_temp, y_temp];
            lastpoint = false;
            changeangle = 0; %0-case: -pi/2 to pi/2, 1-case:-pi/2 to 0; 2-case: 0 to pi/2
            while (abs(x_temp-measEndX)> 0.1*a)
                   factor1 = rand+(-0.5);%shifting since we want angles between (-pi/4,pi/4)
                   factor2 = rand;
                   if (dist == 1)
                       if(unif ~= 1)
                           factor1 = (devnor*randn+0.5);
                           while (factor1 < 0 || factor1 > 1)
                               factor1 = (devnor*randn+0.5);
                           end
                            factor1 = factor1+(-0.5);%shifting since we want angles between (-pi/4,pi/4)
                       end
                       if(unif ~= 2)
                           factor2 = (devnor*randn+0.5);
                           while (factor2 < 0 || factor2 > 1)
                               factor2 = (devnor*randn+0.5);
                           end
                       end
                   elseif (dist == 2)
                       interval1 = rand;
                       interval2 = rand;
                       if (interval1 <= 0.5)
                           if(unif ~= 1)
                               factor1 = exprnd(meanexp)+0.5;
                               while (factor1 < 0.5 || factor1 > 1)
                                  factor1 = exprnd(meanexp)+0.5;
                               end
                               factor1 = factor1+(-0.5);%shifting since we want angles between (-pi/4,pi/4)
                           end
                       else
                           if(unif ~= 1)
                               factor1 = -1.*exprnd(meanexp)+0.5;
                               while (factor1 < 0 || factor1 > 0.5)
                                  factor1 = -1.*exprnd(meanexp)+0.5;
                               end
                               factor1 = factor1+(-0.5);%shifting since we want angles between (-pi/4,pi/4)
                           end
                       end
                       if (interval2 <= 0.5)
                           if(unif ~= 2)
                               factor2 = exprnd(meanexp)+0.5;
                               while (factor2 < 0.5 || factor2 > 1)
                                  factor2 = exprnd(meanexp)+0.5;
                               end
                           end
                       else
                           if(unif ~= 2)
                               factor2 = -1.*exprnd(meanexp)+0.5;
                               while (factor2 < 0 || factor2 > 0.5)
                                  factor2 = -1.*exprnd(meanexp)+0.5;
                               end
                           end
                       end
                   elseif (dist == 3)
                       if(unif ~= 1)
                               factor1 = raylrnd(scrgh);
                               while (factor1 < 0 || factor1 > 1)
                                  factor1 = raylrnd(scrgh);
                               end
                               factor1 = factor1+(-0.5);%shifting since we want angles between (-pi/4,pi/4)
                       end
                       if(unif ~= 2)
                               factor2 = raylrnd(scrgh);
                               while (factor2 < 0 || factor2 > 1)
                                  factor2 = raylrnd(scrgh);
                               end
                       end
                   end
                   x_temp = listXY(end,1)+cos(factor1*randangle1)*factor2*randdis;
                   y_temp = listXY(end,2)+sin(factor1*randangle1)*factor2*randdis;
                   rhopar_temp = factor2*randdis + 0.02*(factor2*randdis);
                   rhoper_temp = factor2*randdis - 0.02*(factor2*randdis);
                   rho_temp = factor2*randdis;
                   if(changeangle == 1)
                      factor1 = (0.5)*rand+(-0.5);
                      if ((dist == 1) && (unif ~= 1))
                           factor1 = (devnor*randn+0.5);
                           while (factor1 < 0 || factor1 > 1)
                               factor1 = (devnor*randn+0.5);
                           end
                           factor1 = 0.5*factor1+(-0.5);
                      elseif ((dist == 2) && (unif ~= 1))
                           interval = rand;
                           if (interval <= 0.5)
                               factor1 = exprnd(meanexp)+0.5;
                               while (factor1 < 0.5 || factor1 > 1)
                                  factor1 = exprnd(meanexp)+0.5;
                               end
                               factor1 = 0.5*factor1+(-0.5);
                           else
                               factor1 = -1.*exprnd(meanexp)+0.5;
                               while (factor1 < 0 || factor1 > 0.5)
                                  factor1 = -1.*exprnd(meanexp)+0.5;
                               end
                               factor1 = 0.5*factor1+(-0.5);
                           end
                      elseif ((dist == 3) && (unif ~= 1))
                          factor1 = raylrnd(scrgh);
                           while (factor1 < 0 || factor1 > 1)
                               factor1 = raylrnd(scrgh);
                           end
                           factor1 = 0.5*factor1+(-0.5);
                      end
                      x_temp = listXY(end,1)+cos(factor1*randangle1)*factor2*randdis;
                      y_temp = listXY(end,2)+sin(factor1*randangle1)*factor2*randdis;
                      changeangle = 0;
                   elseif (changeangle == 2)
                      factor1 = 0.5*rand;   
                      if ((dist == 1) && (unif ~= 1))
                           factor1 = (devnor*randn+0.5);
                           while (factor1 < 0 || factor1 > 1)
                               factor1 = (devnor*randn+0.5);
                           end
                           factor1 = 0.5*factor1;
                      elseif ((dist == 2) && (unif ~= 1))
                           interval = rand;
                           if (interval <= 0.5)
                               factor1 = exprnd(meanexp)+0.5;
                               while (factor1 < 0.5 || factor1 > 1)
                                  factor1 = exprnd(meanexp)+0.5;
                               end
                               factor1 = 0.5*factor1;
                           else
                               factor1 = -1.*exprnd(meanexp)+0.5;
                               while (factor1 < 0 || factor1 > 0.5)
                                  factor1 = -1.*exprnd(meanexp)+0.5;
                               end
                               factor1 = 0.5*factor1;
                           end
                      elseif ((dist == 3) && (unif ~= 1))
                           factor1 = raylrnd(scrgh);
                           while (factor1 < 0 || factor1 > 1)
                               factor1 = raylrnd(scrgh);
                           end
                           factor1 = 0.5*factor1;
                      end
                      x_temp = listXY(end,1)+cos(factor1*randangle1)*factor2*randdis;
                      y_temp = listXY(end,2)+sin(factor1*randangle1)*factor2*randdis;
                      changeangle = 0;
                   end
                   if((x_temp < endX) && (x_temp > startX) && (y_temp> startY) && (y_temp < endY))
                       listXY = [listXY; [x_temp, y_temp]];
                   elseif (x_temp > endX)
                       x_temp = measEndX;
                       y_temp = measEndY;
                       listXY = [listXY; [x_temp, y_temp]];
                       lastpoint = true;
                   elseif (y_temp > endY)
                       %find the intersection with the (y=endY)-line
                       slope = (y_temp-listXY(end,2))/(x_temp-listXY(end,1));
                       offset = y_temp - slope*x_temp;
                       x_temp = (endY-offset)/slope;    
                       y_temp = endY;
                       %recalculate distance
                       dis=sqrt((x_temp-listXY(end,1))*(x_temp-listXY(end,1))+(y_temp-listXY(end,2))*(y_temp-listXY(end,2)));
                       rhopar_temp = dis + 0.02*dis;
                       rhoper_temp = dis - 0.02*dis;
                       rho_temp = dis;
                       listXY = [listXY; [x_temp, y_temp]];
                       changeangle = 1;
                   elseif (y_temp < startY)
                       %find the intersection with the (y=startY)-line
                       slope = (y_temp-listXY(end,2))/(x_temp-listXY(end,1));
                       offset = y_temp - slope*x_temp;
                       x_temp = (startY-offset)/slope;
                       y_temp = startY;
                       %recalculate distance
                       dis=sqrt((x_temp-listXY(end,1))*(x_temp-listXY(end,1))+(y_temp-listXY(end,2))*(y_temp-listXY(end,2)));
                       rhopar_temp = dis + 0.02*dis;
                       rhoper_temp = dis - 0.02*dis;
                       rho_temp = dis;
                       listXY = [listXY; [x_temp, y_temp]];
                       changeangle = 2;
                   end
                   %calculate angle
                   angle_corr = atan(abs(listXY(end,2)-listXY(end-1,2))/abs(listXY(end,1)-listXY(end-1,1)));
                   angle_corr = angle-angle_corr;
                   R_temp = rhoper_temp + (rhopar_temp - rhoper_temp)*cos(angle_corr)*cos(angle_corr);
                   listR_mag = [listR_mag; R_temp];
                   listanglecorr = [listanglecorr; angle_corr];
                   listR = [listR; rho_temp];
            end
            generated = 1;
            %connect last point to measurement end point
            if (lastpoint == false)
                  dis=sqrt((measEndX-listXY(end,1))*(measEndX-listXY(end,1))+(measEndY-listXY(end,2))*(measEndY-listXY(end,2)));
                  rhopar_temp = dis + 0.02*dis;
                  rhoper_temp = dis - 0.02*dis;
                  rho_temp = dis;
                  listXY = [listXY; [measEndX, measEndY]];
                  angle_corr = atan(abs(listXY(end,2)-listXY(end-1,2))/abs(listXY(end,1)-listXY(end-1,1)));
                  angle_corr2 = angle-angle_corr; 
                  R_temp = rhoper_temp + (rhopar_temp - rhoper_temp)*cos(angle_corr2)*cos(angle_corr2);
                  listR_mag = [listR_mag; R_temp];
                  listR = [listR; rho_temp];
                  listanglecorr = [listanglecorr; angle_corr];
            end
        elseif (generated == 1)
            for m=1:length(listR)
                dis = listR(m);
                rhopar_temp = dis + 0.02*dis;
                rhoper_temp = dis - 0.02*dis;
                angle_corr = angle-listanglecorr(m);
                R_temp = rhoper_temp + (rhopar_temp - rhoper_temp)*cos(angle_corr)*cos(angle_corr);
                listR_mag = [listR_mag; R_temp];
            end
        end
        %Calculate total electrical resistance
        R_mag_total = sum(listR_mag);
        R_mag_smooth = ((a-0.02*a) + ((a+0.02*a) - (a-0.02*a))*cos(angle)*cos(angle));
        angles = [angles;angle];
        magRpaper = [magRpaper;(R_mag_total)];
        magRsmooth = [magRsmooth;(R_mag_smooth)];
    end
    %plotyy(angles, magRpaper, angles, magRsmooth)
    %Plot path
    %figure
    plot(listXY(:,1),listXY(:,2),'o-')
    hold on
    %toc
    all_magR=[all_magR,magRpaper];
    if (count==population)
        %find the maximum
        all_magRmax = max(all_magR,[],2);
        %find the minimum
        all_magRmin = min(all_magR,[],2);
        %caluclate the parallel resistance of all pathes
        all_magR_temp = 1./all_magR;
        totpar=sum(all_magR_temp,2);
        totpar=1./totpar;
        %calculate the average
        all_magR = mean(all_magR,2);
        %calculate the errors
        all_magRmax = all_magRmax - all_magR;
        all_magRmin = all_magR - all_magRmin;
        %add the reference solution
        all_magR = [all_magR,angles,magRsmooth,all_magRmax,all_magRmin,totpar];
    end  
end

dRRpaper = (max(all_magR(:,1))-min(all_magR(:,1)))/min(all_magR(:,1));
sav=[sav;dRRpaper];
dRRglass = (max(all_magR(:,3))-min(all_magR(:,3)))/min(all_magR(:,3));

%figure
%[ax h1 h2] = plotyy(all_magR(:,2),all_magR(:,3),all_magR(:,2),all_magR(:,1))
%set(h1,'LineStyle','--','Color','b','Marker','o') 
%set(h2,'LineStyle','--','Color','r','Marker','*') 
%xlabel('angles [rad]')
%axes(ax(1));ylabel('smooth');
%axes(ax(2));ylabel('paper');

%figure
%[ax h1 h2] = plotyy(all_magR(:,2),all_magR(:,3),all_magR(:,2),all_magR(:,6))
%set(h1,'LineStyle','--','Color','b','Marker','o') 
%set(h2,'LineStyle','--','Color','r','Marker','*') 
%xlabel('angles [rad]')
%axes(ax(1));ylabel('smooth');
%axes(ax(2));ylabel('paper');
 
%plot(all_magR(:,2),all_magR(:,1));
%errorbar(all_magR(:,2), all_magR(:,1), all_magR(:,4), all_magR(:,5));
%xlabel('angles [rad]');
%ylabel('R paper');

%toc
end
