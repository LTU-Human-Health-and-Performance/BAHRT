function [angle] = find_angle(vector1, vector2)
%find_angle finds the angle between two vectors in an array of vectors (each column in the matrix being a 3D-vector)
%   Detailed explanation goes here
vector1 = vector1 ./ sqrt(vector1(1,:).^2 + vector1(2,:).^2 + vector1(3,:).^2);
vector2 = vector2 ./ sqrt(vector2(1,:).^2 + vector2(2,:).^2 + vector2(3,:).^2);

angle = acosd(dot(vector1, vector2));

end

