%%%%%%%%%%%%%%%%%%%

function [ f , gradf ] = f_xxx( Q , X , X_indices , SA , SB )

    % ideally, the summary matrices SA and SB will have been computed
    % externally.  if not, we do it here, but it is inefficient
    if nargin < 5 || isempty(SA) || isempty(SA)
        %%%%%%%%%
        % form SA and SB
        %%%%%%%%%
        SA = X(:,X_indices)*X(:,X_indices)';
        SB = X*X';
    end
        
    %%%%%%%%%%
    % now eval
    %%%%%%%%%%
    % both singletons
    trA = trace(Q'*SA*Q);
    trB = trace(Q'*SB*Q);
    
    % the function to *minimize*
    f = - trA / trB ;
    
    % the grad in Q ... usual quotient rule
    gradf =  - ( (2*trB)*SA*Q - (2*trA)*SB*Q )./(trB^2);
    
end 