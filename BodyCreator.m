function [ Poz,Vel ] = BodyCreator(Rd,Rink,Count,Geom)
% The function generates the initial distributions of particles 
% inside the volume of the droplet.
%
% The input parameters:
%      Rd - initial radius of the droplet; 
%      Rink -  radius of inclusions; 
%      Count -  number of inclusions in the droplet's volume, 
%      Geom - the geometry of the task, 2 - for 2D geometry; 3 for 3D geometry 
%
% Output parametes:
%     Poz - matrix of the position of inclusions [particles id, coordinates]
%     Vel - matrix of the initial velocity of inclusions [particles id,
%     speed components]
% =========================================================================
function [ P,V ] = body2D(rd,rink,count);
    Sd = ( pi*rd^2 )/2;
    Sink = ( pi*rink^2 )/2;
    
    if ( Sd * 0.7 ) <= ( Sink * count )
           count = 0.7 * Sd / Sink;
    end
    
    P = zeros( count,2 );
    V = P;
    
    n = 1;
    while n ~= count + 1
%  while n ~= count + 1 - trzeba brac count+1 bo petla najpierw sprawdza
%  warunek a potem licze 
        rho = rd * rand;
        theta  = 2 * pi * rand;
        [ X,Y ] = pol2cart( theta,rho );
        Votknul = 1;
        
        for ii = 1:count
            if norm( [P(ii,1) - X, P(ii,2) - Y ]) <= 2*rink
                Votknul = 0;
            end
        end
        
        if Votknul
            P(n,:) = [X,Y];
            n = n + 1;
        end
    end
% -------------------------------------------------------------------------
    draw = 1;
    if draw 
        figure;
        xx = rd * sin(0:0.01:2*pi);
        yy = rd * cos(0:0.01:2*pi);
        plot(xx,yy);
        grid on;
        hold on;
        for ii = 1 : size(P,1)
            xx = rink * sin(0:0.01:2*pi)+P(ii,1);
            yy = rink * cos(0:0.01:2*pi)+P(ii,2);
            plot(xx,yy);
        end
    end
% -------------------------------------------------------------------------   
end
% =========================================================================
function [ P,V ] = bodyGR(rd,rink);
   
    count = floor( pi*rd/rink );
    P = zeros(count,2);
    V = P;
    X = 0;
    Y = 0;
    dtheta = 2 * asin( rink/rd );
    for ii = 1: count 
        [ X,Y ] = pol2cart(dtheta*(ii-1),rd);
        P( ii,1 ) = X;
        P( ii,2 ) = Y;
    end
%---------------------------------------------------------------
      draw = 1;
    if draw 
        figure;
        xx = rd * sin(0:0.01:2*pi);
        yy = rd * cos(0:0.01:2*pi);
        plot(xx,yy);
        grid on;
        hold on;
        for ii = 1 : size(P,1)
            xx = rink * sin(0:0.01:2*pi)+P(ii,1);
            yy = rink * cos(0:0.01:2*pi)+P(ii,2);
            plot(xx,yy);
        end
    end
end
% =========================================================================
function [ P,V ] = body3D(rd,rink,count);
    Sd = ( pi*rd^2 )/2;
    Sink = ( pi*rink^2 )/2;
    
    if ( Sd * 0.7 ) <= ( Sink * count )
           count = 0.7 * Sd / Sink;
    end
    
    P = zeros( count,3 );
    V = P;
    
    n = 1;
    while n ~= count + 1
%  while n ~= count + 1 - trzeba brac count+1 bo petla najpierw sprawdza
%  warunek a potem licze 
        R = rd * rand;
        THETA  =  pi * rand;
        PHI = 2*pi * rand;
        [X,Y,Z] = sph2cart(THETA,PHI,R);
        Votknul = 1;
        
        for ii = 1:count
            if norm( [P(ii,1) - X, P(ii,2) - Y,P(ii,3) - Z ]) <= 2*rink
                Votknul = 0;
            end
        end
        
        if Votknul
            P(n,:) = [X,Y,Z];
            n = n + 1;
        end
    end
% -------------------------------------------------------------------------
    draw = 1;
     figure;
        xx = rd * sin(0:0.01:2*pi);
        yy = rd * cos(0:0.01:2*pi);
        zz = zeros(1,length(xx));
        plot3(xx,yy,zz);
        grid on;
        hold on;
        plot3(xx,zz,yy);
        plot3(zz,xx,yy);
%         u = (0:0.1*pi:2*pi)';
%     v = [0:0.1*pi:2*pi];
%     X1 = rd*sin(u)*cos(v) + P(ii,1);
%     Y1 = rd*sin(u)*sin(v) + P(ii,2);
%     Z1 = rd*cos(u)*ones(size(v)) + P(ii,3);
    u = (0:0.1*pi:2*pi)';
    v = [0:0.1*pi:2*pi];
for ii = 1 : size(P,1)
    X1 = rink*sin(u)*cos(v) + P(ii,1);
    Y1 = rink*sin(u)*sin(v) + P(ii,2);
    Z1 = rink*cos(u)*ones(size(v)) + P(ii,3);
    hS = mesh(X1,Y1,Z1);
end
end
% =========================================================================

switch Geom
    case 2
        [ Poz,Vel ] = body2D(Rd,Rink,Count);
    case 3
        [ Poz,Vel ] = body3D(Rd,Rink,Count);
    case 4
        [ Poz,Vel ] = bodyGR(Rd,Rink);
    otherwise
        disp('Unknown metod')
end



end