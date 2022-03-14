function output = odd_int(input)
    % highest odd integer <= input
    % input - integer scalar
    output = input - 1*(abs(mod(input,2)-1));
end

