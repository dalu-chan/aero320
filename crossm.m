function crossProductMatrix = crossm(vector)
    % Check if the input vector is of the correct size
    if numel(vector) ~= 3
        error('Input vector must have exactly 3 elements');
    end

    % Extract components of the vector
    x = vector(1);
    y = vector(2);
    z = vector(3);

    % Create the cross product matrix
    crossProductMatrix = [0, -z, y;
                          z, 0, -x;
                         -y, x, 0];
end