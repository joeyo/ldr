%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Q , fQ , info ] = run_xxx( X , X_indices , r , method , Q_0 )

    %%%%%%%%%
    % input checking
    %%%%%%%%%
    if nargin < 5 || isempty(Q_0)
        % initial random guess
        Q_0 = project_stiefel(randn(size(X,1),r));
    end
    if nargin < 4 || isempty(method)
        % default stiefel
        method = 'stiefel';
    end
    if nargin < 3 || isempty(r)
        % default the standard choice
        r = size(X,1);
    end
    if nargin < 2 || isempty(X_indices) || isempty(X)
        % unacceptable
        fprintf('ERROR: Bad X and/or X_indices inputs to run_xxx.')
        keyboard
    end
    
    %%%%%%%%%
    % form SA and SB
    %%%%%%%%%
    SA = X(:,X_indices)*X(:,X_indices)';
    SB = X*X';
   
    
    %%%%%%%%%
    % Project in the specified manner
    %%%%%%%%%
    switch(lower(method))
        
        case 'heuristic'
            % the heuristic manner
            % start timer
            t0 = tic;

            % do svd to get the bases in the d space
            [ U , ~ , ~ ] = svd(SB\SA);
            % the result is the first r columns of U
            Q = U(:,1:r);
            % and the objective
            fQ = f_xxx( Q , [] , [] , SA , SB );
            %
            % info 
            info(1).time = toc(t0);
            info(1).iter = 1;

        case 'stiefel'
            
            % get initial value
            [ Q_0 , ~ ] = run_xxx( X , X_indices, r , 'heuristic');
            
            % run the method
            [ Q , fQ , info ] = minimize_stiefel_sd( 'f_xxx' , Q_0 , [] , [] , [] , SA , SB );
            
        case 'stiefel_trust'
            
            [ Q , fQ , info ] = minimize_stiefel_trust( 'f_xxx' ,  Q_0 , [] , [] , [] , SA , SB );
            
        case 'grassmann_trust'
            
            [ Q , fQ , info ] = minimize_grassmann_trust( 'f_xxx' ,  Q_0 , [] , [] , [] , SA , SB );

        case 'stiefel_mosd'
            
            [ Q , fQ , info ] = minimize_stiefel_mosd( 'f_xxx' ,  Q_0 , [] , [] , [] , SA , SB );
            
        case 'grassmann_mosd'
            
            [ Q , fQ , info ] = minimize_grassmann_mosd( 'f_xxx' ,  Q_0 , [] , [] , [] , SA , SB );

            
        otherwise
            fprintf('Your method is not implemented.  Try stiefel or heuristic.');
            error()
    end
    
            
            