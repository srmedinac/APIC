
function A = rotationTranslationMatrix(theta, c)
% rotationTranslationMatrix: Returns a rotation ans translation matrix for
% a given angle and offset.

if ~exist('theta', 'var') || isempty(theta)
    theta = 0;
end

if ~exist('c', 'var') || isempty(c)
    c = zeros(1, 2);
end

A = [ cos(theta) -sin(theta) 0  ;
    sin(theta)  cos(theta) 0  ;
    c(1)        c(2)       1 ];

end