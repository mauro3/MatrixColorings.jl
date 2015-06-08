export isequal_sparse

function isequal_sparse(A_1::SparseMatrixCSC, A_2::SparseMatrixCSC)
    (   size(A_1)==size(A_2)
    && nonzeros(A_1)==nonzeros(A_2)
    && rowvals(A_1)==rowvals(A_2)
    && A_1.colptr==A_2.colptr)
end

if VERSION < v"0.4.0-dev"
    # TODO: remove these two for 0.4
    rowvals(S::SparseMatrixCSC) = S.rowval
    nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)
end
