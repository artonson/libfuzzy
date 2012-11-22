function column = ensure_column(vector)
%ENSURE_COLUMN Return column given row or column
%   Use to make sure a vector is a column vector.

vector_dimenstions = size(vector);
if length(vector_dimenstions) > 2 || ...
        (length(vector_dimenstions) == 2 && sum(vector_dimenstions > 1) == 2)
    error('libfuzzy:ensure_column:too_many_dims', ...
        'given is not a vector')
end

if vector_dimenstions(1) > vector_dimenstions(2)
    column = vector;
else
    column = vector';
end

end

