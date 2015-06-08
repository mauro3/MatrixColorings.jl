export isequal_sparse

function isequal_sparse(A_1::SparseMatrixCSC, A_2::SparseMatrixCSC)
    (   size(A_1)==size(A_2)
    && nonzeros(A_1)==nonzeros(A_2)
    && rowvals(A_1)==rowvals(A_2)
    && A_1.colptr==A_2.colptr)
end

function ==(A1::SparseMatrixCSC, A2::SparseMatrixCSC)
    size(A1)!=size(A2) && return false
    rows1 = rowvals(A1)
    vals1 = nonzeros(A1)
    vals2 = nonzeros(A2)

    m, n = size(A1)
    for i = 1:n
        nz1 = nzrange(A1, i)
        nz2 = nzrange(A2, i)
        if nz1==nz2
            for j in nz1
                vals1[j]!=vals2[j] && return false
            end
        else # do slow indexing on A2
            for j in nz1
                jj = rows1[j]
                vals1[j]!=A2[jj,i] && return false
            end
        end
    end
    return true
end

if VERSION < v"0.4.0-dev"
    # TODO: remove these two for 0.4
    rowvals(S::SparseMatrixCSC) = S.rowval
    nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)
end
