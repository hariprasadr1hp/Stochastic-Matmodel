function Envelope = calcEnvelope(lambda,a1,a2,b,p,W,spacing,m,k_max)

% Envelope  Calculates the Monte carlo envolopes
%       95%-envelope of the three quadratic Minkowski functions
%
% lamda = intensity
% a1,a2  = to determine a probabilitic value for the Semi major axis 'a'
% b = semi minor axis
% p = probablity of a1 or a2
% W = window size specification
% spacing = pixel Spacing
% m  = number of sample realisations
% k_max= maximum range of k in pixel length

minkovFnTable = zeros(k_max+1,4,m);
lowIndex = floor(m*.025) +1;
highIndex = floor(m*.975) +1;
envol_data = zeros(k_max+1,6); % table column [L1,U1,L2,U2,L3,U3]

for i = 1:m
    infoEllipse = rBoolEllipse(lambda,b,a1,a2,p,W);
    rbool =digitizeEllSys(infoEllipse,W,spacing);  % digitization of the generated Realisation
    minkovFnTable(:,:,i) =estQMinkowskiFcts(rbool,k_max,spacing); %  spacing
end


% Sorting the minkovFnTable for each r
sortedMinkovFnTable = sort(minkovFnTable,3); 
envol_data(:,1)= sortedMinkovFnTable(:,1,1);
for j= 1:k_max+1
    envol_data(j,1)= sortedMinkovFnTable(j,2,lowIndex);
    envol_data(j,2)= sortedMinkovFnTable(j,2,highIndex);
    envol_data(j,3)= sortedMinkovFnTable(j,3,lowIndex);
    envol_data(j,4)= sortedMinkovFnTable(j,3,highIndex);
    envol_data(j,5)= sortedMinkovFnTable(j,4,lowIndex);
    envol_data(j,6)= sortedMinkovFnTable(j,4,highIndex);
end
Envelope = envol_data;
end