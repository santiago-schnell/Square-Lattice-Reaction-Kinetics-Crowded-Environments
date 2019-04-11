function CA2(runs,N,timesteps,f,r,g,p,filename,disp);

% CA2 uses the LATTICE GAS model of molecular dynamics on a membrane to simulate 
% the interaction of Substrate, Enzyme, Complex and Product molecules involved in the 
% Michaelis Menten reaction E + S <-> C -> E + P.

% Paramaters:

% runs:      number of times to repeat the simulation before averaging results
% N:         dimension of (N x N x N) cubic grid
% timesteps: number of timesteps simulations runs over before termination
% f:         E + S -> C reaction probability == k1
% r:         C -> E + S reaction probability == kNEG1
% g:         C -> E + P reaction probability == k2
% p(i)       VECTOR of initial molecular densities: [E S P C O]
% filename:  name of output .csv file of molecular density data over time
% disp:      BOOLEAN: display output during execution? 1=YES, 0=NO.

filedump=zeros(timesteps,6);    %temporary virtual file store

% ---------
% RUNS LOOP
% ---------

for run=1:runs  % counter for number of iterations of the simulation

% -----------------------
% VARIABLE INITIALISATION
% -----------------------
    
ndata=zeros(timesteps,4);   % array recording no's of each molecule over time    

movieloop=0;                % counter for recording animations
map = [0,0,1;               % colormap of output display
       0,1,0;
       1,0,0;
       1,1,0;
       0,0,0;
       1,1,1];

error = 0;                  %error flag for debugging purposes

% A is an N x N matrix representing the 2D membrane.
% each cell contains one of the following entries:
%   6 - the cell is empty
%   1 - the cell is occupied by an ENZYME (E) molecule
%   2 - the cell is occupied by a SUBSTRATE (S) molecule
%   3 - the cell is occupied by a PRODUCT (P) molecule
%   4 - the cell is occupied by an enzyme-substrate COMPLEX (C) molecule
%   5 - OBSTRUCTION (O) - no molecule can occupy this cell.

A = 6*ones(N,N,N);    %An initially empty membrane


% 3D Matrix M contains 5 lists of the (x,y) coordinates of each molecule type
% (E, S, P, C or O) in 5 adjacent 2D lists.
% e.g. M(3,7,1) contains the X co-ord of the 7th molecule of type 3 (P)
%      M(3,7,2) contains the Y co-ord of the 7th molecule of type 3 (P)

M = zeros(1,5,3);


% Initial densities (p(X)) dictate the no of initial molecules of each type required
% Density = fraction of grid sites occupied by molecule type.
% (Initial densities specified in initialisation paramaters)

% Translate to absolute no's of each molecule type (n(i))

for i = 1:5
    n(i) = ceil(N*N*N*p(i));
end;

% Gamma is a counter that monitors the number of E + S -> C reactions that have 
% occurred throughout each simulation

gamma = 0;

% Randomly assign the initial molecules to grid sites...

for i = 1:5
    for j = 1:n(i)
        placed=0;
        while placed == 0 
            x = ceil(N*rand(1));
            y = ceil(N*rand(1));
            z = ceil(N*rand(1));
            if A(x,y,z) == 6
                A(x,y,z) = i;
                M(j,i,1) = x;
                M(j,i,2) = y;
                M(j,i,3) = z;
                placed = 1;
            end
        end
    end
end

% --------
% MAINLOOP
% --------

% Increment thru time randomly selecting molecules and updating their
% behaviour based on pre-determined rules.

for mainloop = 1:timesteps
    
time = 0;   % time is a counter that ensures each timestep is of length 1/(total no.
            % of molecules) such that in a given timestep on average each molecule 
            % will move once.

%nTOT = norm(n,1)-n(5);  %nTOT keeps track of the total no of mols in simulation
            
% ----------
% INNER LOOP
% ----------
            
while time < 1    
        
    nTOT = norm(n,1)-n(5);  %nTOT keeps track of the total no of mols in simulation
    
    %select a molecule at random
    mol=ceil(nTOT*rand(1));
    for type=1:4
        if mol > n(type), mol = mol-n(type);
            else break;
        end
    end
    x = M(mol,type,1);
    y = M(mol,type,2);
    z = M(mol,type,3);
    
    % randomly choose direction of movement
    % 1=UP, 2=DOWN, 3=LEFT, 4=RIGHT
    i = ceil(6*rand(1));
       
    switch i
        case {1}
            xdest = x;          %xdest, ydest contain the co-ords a molecule wants
            ydest = y - 1;      %to move to.
            zdest = z;
        case {2}
            xdest = x;
            ydest = y + 1;
            zdest = z;
        case {3}
            xdest = x - 1;
            ydest = y;
            zdest = z;
        case (4) 
            xdest = x + 1;
            ydest = y;
            zdest = z;
        case (5) 
            xdest = x;
            ydest = y;
            zdest = z - 1;
        otherwise 
            xdest = x;
            ydest = y;
            zdest = z + 1;
    end
    
    
    % the membrane has cyclic boundary conditions (go off one end and appear on the other)
    % boundarycheck is a routine that checks if the chosen destination site is off the 
    % edge of the membrane and if so redirects the molecule to the opposite edge
    [xdest,ydest,zdest] = boundarycheck3D(xdest,ydest,zdest,N);   
    
    
    % For E,S,P: If destination site is empty molecule moves there and matrices 
    % A & M are updated...
    if type < 4             %(E,S or P)
        
        if A(xdest,ydest,zdest) == 6      % if destination site is empty...
            A(x,y,z) = 6;             % ...move molecule to destination site.
            A(xdest,ydest,zdest) = type;
            M(mol,type,1) = xdest;
            M(mol,type,2) = ydest;
            M(mol,type,3) = zdest;
            
        else                        % if not...
            
        % For E,S: if destination site occupied by S or E respectively, reaction 
        % occurs with probability f producing C in destination site and leaving 
        % original site empty
     
            switch type
                case {1}                        % if chosen molecule is an E...
                    if A(xdest,ydest,zdest) == 2      % ...and destination contains an S...
                        if rand(1) < f          % ...and with probability f...
                            A(xdest,ydest,zdest) = 4; % ...place C in destination site...
                            [M,n(4)] = addmol3D(M,4,xdest,ydest,zdest,n(4));
                            A(x,y,z) = 6;         % ...and clear original site
                            [M,n(1)] = removemol3D(M,1,mol,n(1));
                            row = findmol3D(M,2,xdest,ydest,zdest,n(2));
                            if row == 0, error = 1
                            else [M,n(2)] = removemol3D(M,2,row,n(2)); end
                            gamma = gamma + 1;  % update reaction counter
                        end
                    end
                case {2}                        % if chosen molecule is an S...
                    if A(xdest,ydest,zdest) == 1      % ...and destination contains an E...
                        if rand(1) < f          % ...and with probability f...
                            A(xdest,ydest,zdest) = 4; % ...place C in destination site...
                            [M,n(4)] = addmol3D(M,4,xdest,ydest,zdest,n(4));
                            A(x,y,z) = 6;         % ...and clear original site
                            [M,n(2)] = removemol3D(M,2,mol,n(2));
                            row = findmol3D(M,1,xdest,ydest,zdest,n(1));
                            if row == 0, error = 1
                            else [M,n(1)] = removemol3D(M,1,row,n(1)); end
                            gamma = gamma + 1;  % update reaction counter
                        end
                    end
            end
        end
    end
    
    % For C...
    
    if type == 4
        
        %check availibility of nearest neighbours. store in near(i)
        near = zeros(1,6);
        if y>1, if A(x,y-1,z) == 6, near(1)=1; end    % UP
        else if A(x,N,z) == 6, near(1)=1; end, end
        if y<N, if A(x,y+1,z) == 6, near(2)=1; end    % DOWN
        else if A(x,1,z) == 6, near(2)=1; end, end
        if x>1, if A(x-1,y,z) == 6, near(3)=1; end    % LEFT
        else if A(N,y,z) == 6, near(3)=1; end, end
        if x<N, if A(x+1,y,z) == 6, near(4)=1; end    % RIGHT
        else if A(1,y,z) == 6, near(4)=1; end, end
        if z>1, if A(x,y,z-1) == 6, near(5)=1; end    % BELOW
        else if A(x,y,N) == 6, near(5)=1; end, end
        if z<N, if A(x,y,z+1) == 6, near(6)=1; end    % ABOVE
        else if A(x,y,1) == 6, near(6)=1; end, end
        
        
        % choose a random destination site from vacant neighbours and take a random 
        % no, q... 
        % if q < r:        Place S at destn site, replace C at initial site with E 
        % else if q < r+g: Place P at destn site, replace C at initial site with E
        % else:            Move C to destination site
                
        
        if norm(near,1) > 0     %if there is an empty neighbour site
            
            % Pick random direction from available sites
            dir = ceil(norm(near,1)*rand(1));
            for i = 1:5
                if near(i) == 0 & dir >= i, dir = dir + 1; end
            end
            
            switch dir
                case {1}
                    xdest = x;          %xdest, ydest contain the co-ords a molecule wants
                    ydest = y - 1;      %to move to.
                    zdest = z;
                case {2}
                    xdest = x;
                    ydest = y + 1;
                    zdest = z;
                case {3}
                    xdest = x - 1;
                    ydest = y;
                    zdest = z;
                case (4) 
                    xdest = x + 1;
                    ydest = y;
                    zdest = z;
                case (5) 
                    xdest = x;
                    ydest = y;
                    zdest = z - 1;
                otherwise 
                    xdest = x;
                    ydest = y;
                    zdest = z + 1;
            end
            
            [xdest,ydest,zdest] = boundarycheck3D(xdest,ydest,zdest,N);  %check not going off edge
            
            % Remove C from (x,y)
            [M,n(4)] = removemol3D(M,4,mol,n(4));
            
            q = rand(1);
            
            if q < r + g                         % if q < r+g...
                A(x,y,z) = 1;                      % ...Place E at (x,y)
                [M,n(1)] = addmol3D(M,1,x,y,z,n(1));
            else 
                A(x,y,z) = 6;                      % ...else (x,y) is made empty.
            end
            
            if q < r, newtype = 2;               % if q < r new molecule is an S
            else if q < r + g, newtype = 3;      % else if q < r+g new molecule is a P
                else newtype = 4;                % else molecule stays as a C
                end       
            end
            
            A(xdest,ydest,zdest) = newtype;            % Place S, P or C at (xdest,ydest)
            [M,n(newtype)] = addmol3D(M,newtype,xdest,ydest,zdest,n(newtype));
        end
    end
    
    time = time + 1/nTOT;                        % increment time
    
end

% --------------
% END INNER LOOP
% --------------

ndata(mainloop,1:4)=n(1:4);                      % record no's of each molecule
ndata(mainloop,5)=gamma;                         % record gamma

if mod(mainloop,5) == 0 & disp==1                % output current state to screen...
                                                 % ...every 5 iterations of MAINLOOP
    figure(1);
    %subplot(3,4,[1,2,5,6])                       % plot membrane
    %image(A); colormap(map);
    %axis equal;
    %axis tight;
    
    subplot(3,4,[3,4,7,8])                       % plot bar chart of no's of each...
    bar(n(1:4));                                 % ...molecule
    axis square;
    axis([0,5,0,1000]);
    axis off;
    text(1,-20,'E');
    text(2,-20,'S');
    text(3,-20,'P');
    text(4,-20,'C');
    
    subplot(3,4,[9,10,11,12])                    % plot progress bars
    barh([(run-1+mainloop/timesteps)/runs,mainloop/timesteps]);
    axis([0,1,0,3]);
    axis off;
    pause(0);
    
    % movieloop = movieloop+1; 
    % MOV(movieloop) = getframe(figure(1));
end;

end;

% ------------
% END MAINLOOP
% ------------

check3D     % a routine that verifies the data in matrices A and M agree by the end of 
            % the simulation

% Add data from each run in order to compute an overall average
filedump(:,1:5) = filedump(:,1:5)+ndata;       

if disp == 0                % if graphical output is off, display progress bar only
    figure(1);
    barh((run-1)/runs);
    axis([0,1,0,2]);
    pause(0);
end;

figure(2);
plot(filedump(:,1:5));      % plot averaged data of runs so far

end;

% -------------
% END RUNS LOOP
% -------------

filedump(1:6,6) = [runs,N,timesteps,f,r,g]';      % add config data to virtual file
filedump(7:11,6) = p';

filedump(:,1:5) = filedump(:,1:5)/runs;           % convert sum to average
csvwrite(filename,filedump);                      % write virtual file to disk

%movie2avi(MOV,'movie.avi')