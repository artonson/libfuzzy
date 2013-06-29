function projection = projection_operator(elements)
%PROJECTION_OPERATOR Compute projection to a linear subspace spanned over
%                    elements

projection = elements * elements';

end
