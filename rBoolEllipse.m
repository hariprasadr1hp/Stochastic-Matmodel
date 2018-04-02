
function info_Ellipse = rBoolEllipse(lambda,a1,a2,p,b,W)
% rBoolEllipse  Generates a random configuration of ellipses in 2D.
%
% The output is a matrix containing ellipse centres(x,y coordinates),
% long semi-axes, short semi-axes and Orientation angles
%
% matrix = rBoolEllipse(lambda,b,a1,a2,p,W)
%
% lambda - intensity
% b      - short semi-axis
% a1     - a random number such that 'a1' >= 'b' with a probabilty 'p'
% a2     - a random number such that 'a2' >= 'a1' with a probabilty '1-p'
% p      - probability
% W      - Window Area
x_min = W(1,1);
x_max = W(1,2);
y_min = W(2,1);
y_max = W(2,2);

% Generatre position distribution with window area and intensity lambda
x_space = (x_max - x_min + 2*a2);  % samplespace in x-direction
y_space = (y_max - y_min + 2*a2);  % samplespace in y-direction
sample_space = x_space*y_space;
total_points = poissrnd(lambda*sample_space); 
info_Ellipse = zeros(total_points,5);    
                           
% Store values of diffrent quantities in Matrix
% Each row represents one point
for i = 1:total_points
    info_Ellipse(i,1) = (rand * ((x_max - x_min) + 2*a2)) + (x_min - a2) ;   % X coordinate of the Ellipse centre
    info_Ellipse(i,2) = (rand * ((y_max - y_min) + 2*a2)) + (y_min - a2) ;   % Y coordinate of the Ellipse centre   
    major_axis = rand;
    % setting the semi major axis to 'a1' or 'a2' depending on 
    % the probablity 'p' of 'a1'
    if major_axis <= p
        info_Ellipse(i,3) = a1;              % longer Semi-axis with a1 size
    else
        info_Ellipse(i,3) = a2;          % longer Semi-axis with a2 size 
    end 
    info_Ellipse(i,4) = b;                 % minor Semi-axis
    info_Ellipse(i,5) = unifrnd(0,pi); % Random angle between x-axis and long semi-axis    
end

end
