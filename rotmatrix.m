function [matrix] = rotmatrix(basis2,basis1)
%rotmatrix This function finds the rotation matrix between two input basis'

matrix(1,1) = dot(basis1(:,1),basis2(:,1));
matrix(1,2) = dot(basis1(:,2),basis2(:,1));
matrix(1,3) = dot(basis1(:,3),basis2(:,1));
matrix(2,1) = dot(basis1(:,2),basis2(:,1));
matrix(2,2) = dot(basis1(:,2),basis2(:,2));
matrix(2,3) = dot(basis1(:,2),basis2(:,3));
matrix(3,1) = dot(basis1(:,3),basis2(:,1));
matrix(3,2) = dot(basis1(:,3),basis2(:,2));
matrix(3,3) = dot(basis1(:,3),basis2(:,3));

end