function structA = mergestruct(structA,structB)
    % copies fields of struct B to struct A
    f = fieldnames(structB);
    for i = 1:length(f)
        structA.(f{i}) = structB.(f{i});
    end
end

