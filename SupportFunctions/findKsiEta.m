function array=findKsiEta(data)
    [m,n]=size(data);
    
    for i=1:m
        [P,Q]=convert2PQ(data(i,:));
    
        geom100 = distances_to_set(P,Q,'100'); % Generate geometry parameters
        geom110 = distances_to_set(P,Q,'110');  
        geom111 = distances_to_set(P,Q,'111');
        
        d=[geom100(1,:) geom110(1,:) geom111(1,:)];
        k=[geom100(2,:) geom110(2,:) geom111(2,:)];
        e=[geom100(3,:) geom110(3,:) geom111(3,:)];
        p=[geom100(4,:) geom110(4,:) geom111(4,:)]/(pi/2);
        s=[ones(1,length(geom100(1,:)))*100 ones(1,length(geom110(1,:)))*110 ones(1,length(geom111(1,:)))*111];
        
        lowdist=d<=0.0001;
        highphi=p>=0.9999;
        select=lowdist & highphi;
        
        dist=d(select);
        ksi=k(select);
        eta=e(select);
        phi=p(select);
        set=s(select);
        
        [val, ind]=min(dist);
        
        if ~isempty(dist)
            array(i,1)=set(ind);
            array(i,2)=ksi(ind)*180/pi;
            array(i,3)=eta(ind)*180/pi;
            array(i,4)=phi(ind);
            array(i,5)=data(i,19);
        else
            array(i,1)=0;
            array(i,2)=0;
            array(i,3)=0;
            array(i,4)=0;
            array(i,5)=0;
        end
    end
end

function geom = distances_to_set(P,Q,whichaxes,dismax)
% geom = distances_to_set(P,Q,whichaxes,dismax)
%
% Calculates the geometry parameters for a given grain boundary relative to
% a given set of axes.
%
% P and Q are rotation matrices giving the orientations of the two grains.
% The grain boundary normal is fixed at [1,0,0].
%
% whichaxes is one of '100', '110', or '111'
%
% dismax is an optional parameter specifying the maximum distance to
% include. It defaults to slightly less than 1, which is the largest
% allowable distance before some anomalies start to appear.
%
% Result geom is a 4xn matrix, where the rows are distance, ksi, eta, and
% phi. It keeps all of the hits up to a distance of dismax.
%
% distance is 2*sin(delta/2) where delta is the angle of closest approach
% between a misorientation axis and one of the axes.  Note there are 24
% representations of the rotations and 3, 6, or 4 equivalent high-symmetry
% axes, so it calculates as many as 144 distances.  But only ones below
% the cutoff dismax are kept.
%
% Once it's picked the closest approximation to the boundary for a given
% axis and coset element, it finds the parameters ksi, eta, phi defining
% that idealized boundary (since the axis is defined, it's a 3-space).
%
% These are:
% phi, the angle between the rotation axis and the boundary plane normal
% (taken as the mean of the normals represented in the two actual grain
% orientations, which works when dismax is less than 1)
%
% ksi, the misorientation angle
%
% eta, a parameter giving the second axis of the boundary plane normal in
% terms of specified directions ('dirs') perpendicular to each
% high-symmetry axis.


    if ~exist('dismax','var') || isempty(dismax),
        dismax = 0.999999 ;  % Force the distance to be strictly less than one, allowing for roundoff
    end % Note if dismax >= 1, you're likely to get a warning about m1 being singular.

    switch whichaxes
        case {'110'}
            % Define 110 axes, normalize
            axes = [1  1  1  1  0  0 ;
                    1 -1  0  0  1  1 ;
                    0  0  1 -1  1 -1 ]/sqrt(2) ;

            % Define a crystal direction perpendicular to each rotation axis.
            % The formalism demands that this be an axis of at least two-fold
            % symmetry.
            dirs = [0  0  0  0  1  1 ;
                    0  0  1  1  0  0 ;
                    1  1  0  0  0  0 ] ;
        case {'111'}

            % Define 111 axes, normalize
            axes = [1  1 -1 -1 ;
                    1 -1  1 -1 ;
                    1 -1 -1  1 ]/sqrt(3) ;

            dirs = [1  1  1  1 ;
                   -1  1  1 -1 ;
                    0  0  0  0 ]/sqrt(2) ;

        case {'100'}
            % Define 100 axes, normalize
            axes = [1  0  0 ;
                    0  1  0 ;
                    0  0  1 ] ;

            dirs = [0  0  1 ;
                    1  0  0 ;
                    0  1  0 ] ;
        otherwise
            error('Undefined axis set')
    end

    naxes  = size(axes,2) ;
    period = pi*naxes/6 ;

    %  Define the symmetry operators

    rotX90  = [ 1  0  0 ;     %  Rotation by +90 degrees around X axis
                0  0 -1 ;
                0  1  0 ] ;

    rotY90  = [ 0  0  1 ;     %  Rotation by +90 degrees around Y axis 
                0  1  0 ;
               -1  0  0 ] ;   

    rotZ90  = [ 0 -1  0 ;      %  Rotation by +90 degrees around Z axis
                1  0  0 ;
                0  0  1 ] ;

    rotZ90m = [ 0  1  0 ;      %  Rotation by -90 degrees around Z axis
               -1  0  0 ;
                0  0  1 ] ;

    % Create 24 symmetry equivalent variants of Q
    % This is the coset appropriate for the rotation convention where Q'*P
    % is the misorientation represented in the grain frame.  If you're
    % getting odd results, e.g. misorientations that you know are CSL are
    % coming out entirely wrong, you may be using the opposite convention;
    % try replacing P and Q with P' and Q'.
    V = cell(24,1);

    V{1}  = Q;    
    V{2}  = V{1}*rotX90 ;        % Rotate the vectors three times around X by +90 degrees
    V{3}  = V{2}*rotX90 ;
    V{4}  = V{3}*rotX90 ;

    for j = 1:12                   % Rotate three times around Y by +90 degrees
      V{j+4} = V{j}*rotY90 ;
    end

    for j = 1:4
      V{j+16} = V{j} *rotZ90;    % Rotate three times around Z by +90 degrees
      V{j+20} = V{j} *rotZ90m;    % Rotate three times around Z by -90 degrees   
    end


    % Preallocate all parameter lists at their maximum possible sizes.
    % Redundant representations will be removed at the end.
    distances = zeros(1,24*naxes);
    phis      = zeros(1,24*naxes);
    ksis      = zeros(1,24*naxes);
    etas      = zeros(1,24*naxes);

    thisindex = 0;  % Number of hits found so far

    % Step through all combinations of symmetrically-equivalent axes and coset
    % elements V{j}.
    for i = 1:naxes

        ax   = axes(:,i) ;    % ax is a high-symmetry axis   

        dir = dirs(:,i) ;       %  This is the pivot vector used to partition 
                                %  the rotation around axis "i"
        dir2 = cross(ax,dir);   %  Completing the orthonormal coordinate set.
                                %  theta1 and theta2 are defined in the plane
                                %  spanned by (dir,dir2).


        for j = 1:24    % For each symmetry-related variant of the second grain

            Q     = V{j} ; 
            R     = Q'*P ;   %  This rotates any vector in cube P into a vector in cube Q

            q = mat2quat(R) ;   % Calculation from here on out is much easier with quaternions.
            axi = q(2:4)'/sqrt(sum(q(2:4).^2)); % Normalized rotation axis
            psi = 2*acos(q(1)); % Rotation angle

            dotp  = axi*ax ;

            % Compute rotational distance from boundary P/Q to the rotation set "i" 
            % This formula produces 2*sin(delta/2), where delta is the angle of
            % closest approach.
            dis   = 2*sqrt(abs(1 - dotp*dotp))*sin(psi/2) ;

            if dis < dismax
                thisindex = thisindex + 1;

                theta = 2*atan(dotp*tan(psi/2)) ; % angle of rotation about ax that most closely approximates R

                % Compute the normal of the best-fitting GB in grain 1
                n1    = P(1,:)' ;
                n2    = Q(1,:)' ;

                RA = quat2mat([cos(theta/2);sin(theta/2)*ax]);
                % RA is the rotation about ax that most closely approximates R

                % From this point on we're dealing with the idealized rotation RA, not
                % the original rotation R.
                m1    = n1 + RA'*n2 ;

                % The next problem comes up only for very large distances,
                % which are normally cut off
                if norm(m1) < 0.000001
                    disp('m1 is singular!!!')
                end

                m1    = m1/norm(m1)  ;  % Halfway between the two normal vectors from the two grains
                m2    = RA*m1 ;   % And the same represented in the other grain

                % Compute the inclination angle for the common rotation axis
                phi   = real(acos(abs(m1'*ax))) ; % "real" because of numerical problems when they're exactly parallel

                % Partition the total rotation angle "theta"
 
                if abs(ax'*m1) > 0.9999      % Check if the best-fitting GB is pure twist
                    theta1 = - theta/2 ;   % eta is meaningless for a twist boundary.
                    theta2 =   theta/2 ;
                else

                    theta1 = atan2(dir2'*m1,dir'*m1);
                    theta2 = atan2(dir2'*m2,dir'*m2);
                    % It's projecting m1 and m2 into the plane normal to ax and
                    % then determining the rotation angles of them relative to
                    % dir.

                end        

                % Reduce both angles to interval (-period/2,period/2],
                % semi-open with a small numerical error.
                theta2  = theta2 - round(theta2/period)*period ;
                theta1  = theta1 - round(theta1/period)*period ;

                % This implements the semi-open interval in order to avoid an
                % annoying numerical problem where certain representations are
                % double-counted.
                if abs(theta2+period/2)<0.000001,
                    theta2 = theta2 + period;
                end
                if abs(theta1+period/2)<0.000001,
                    theta1 = theta1 + period;
                end

                % Since this is only being run on fcc elements, which are
                % centrosymmetric, and all dir vectors are 2-fold axes, then
                % the operations of swapping theta1 and theta2, and of
                % multilying both by -1, are symmetries for the energy
                % function. This lets us fold everything into a small right
                % triangle in (ksi,eta) space:
                ksi     = abs(theta2 - theta1) ;
                eta     = abs(theta2 + theta1) ;

                % And store them in the vectors
                distances(thisindex) = dis;
                ksis(thisindex)      = ksi;
                etas(thisindex)      = eta;
                phis(thisindex)      = phi;
            end
        end      
    end   

    % Dump the excess pre-allocated ones and sort the rest in order of distance
    [distances,sortindex] = sort(distances(1:thisindex));
    ksis = ksis(sortindex);
    etas = etas(sortindex);
    phis = phis(sortindex);

    % Clean up redundancy.  Double-counting the same representation of one
    % boundary messes up the weighting functions in weightedmeanenergy.m

    % First round everything to 1e-6, so that negligible numerical
    % differences are dropped
    distances = 1e-6*round(distances*1e6);
    ksis = 1e-6*round(ksis*1e6);
    etas = 1e-6*round(etas*1e6);
    phis = 1e-6*round(phis*1e6);

    % And finally create the 4 x thisindex array of geometrical parameters
    geom = unique([distances',ksis',etas',phis'],'rows')';

end

function q = mat2quat(m)
% q = mat2quat(m)
%
% Auxiliary function converts a rotation matrix, assumed orthonormal, into
% a unit quaternion.
    t = m(1,1)+m(2,2)+m(3,3);
    e0 = sqrt(1+t)/2;
    if t > -0.999999999
        e = [m(2,3)-m(3,2);m(3,1)-m(1,3);m(1,2)-m(2,1)]/(4*e0);
    else
        e0 = 0;
        e3 = sqrt(-(m(1,1)+m(2,2))/2);
        if abs(e3) > 2e-8   % Check for singularity, allowing numerical error
            e = [m(1,3)/(2*e3) ; m(2,3)/(2*e3) ; e3];
        else
            e1 = sqrt((m(1,1)+1)/2);
            if e1 ~= 0
                e = [e1;m(2,1)/(2*e1);0];
            else
                e = [0;1;0];
            end
        end
    end
    
    q = [e0;-e];
end


function m = quat2mat(q)
% m = quat2mat(q)
%
% Auxiliary function converts a quaternion into a rotation matrix with no
% assumption about normalization.
    e0 = q(1);
    e1 = q(2);
    e2 = q(3);
    e3 = q(4);

    m = [e0^2+e1^2-e2^2-e3^2 , 2*(e1*e2-e0*e3) , 2*(e1*e3+e0*e2); ...
          2*(e1*e2+e0*e3) , e0^2-e1^2+e2^2-e3^2 , 2*(e2*e3-e0*e1); ...
          2*(e1*e3-e0*e2) , 2*(e2*e3+e0*e1) , e0^2-e1^2-e2^2+e3^2 ]...
          /(e0^2+e1^2+e2^2+e3^2);
end
