using LinearAlgebra

function upper_triangulate(A)

    for j = 1:N
        pivot!(A,j)
        eliminate!(A, j)
    end

    return M
end

function eliminate!(A, j)

    pivot = A[j, j]
    pivot_row = A[j, :]

    for k = j+1:M
        fac = A[k, j]/pivot
        A[k, :] = A[k, :] - fac*pivot_row
    end

end

function pivot!(A,j)
    pivot = A[j,j]
    pivot_row = j
    new_pivot_row = j 

    for k = j+1:M
        if abs(A[k,j]) > abs(pivot)
            pivot = abs(A[k,j])
            new_pivot_row = k
        end
    end
    if(pivot_row != new_pivot_row)
        pivotMatrix = Matrix{Float64}(I, M, N)
        swap!(pivotMatrix, pivot_row, new_pivot_row)
        push!(pivotMatricies, pivotMatrix)
        swap!(A,pivot_row,new_pivot_row)
    end
    # temp = pivotMatrix[j,:]
    # pivotMatrix[j,:] = pivotMatrix[pivot_row,:]
    # pivotMatrix[pivot_row,:] = temp
    # push!(pivotMatricies, pivotMatrix)
    # temp = A[j,:]
    # A[j,:] = A[pivot_row,:]
    # A[pivot_row,:] = temp
end

function swap!(A,j,k)
    A[j,:], A[k,:] = A[k,:], A[j,:]
end

function pivotMultiply()
    PM = pop!(pivotMatricies)
    for i = 1:size(pivotMatricies)[1]
        PM = PM * pop!(pivotMatricies)
    end
    return PM
end

function computeLUP(A)
    UTM = copy(A)
    upper_triangulate(UTM)
    PM = pivotMultiply()
    LTM = (PM*A)/UTM
    return UTM, LTM, PM
end

function LUPsolve(A,b) 
end

A = Float64[-2 2 -1; 6 -6 7; 3 -8 4]
# A = Float64[1 2 3; 4 5 6; 7 8 0]
M = size(A,1)
N = size(A,2)
pivotMatricies = []

UTM, LTM, PM = computeLUP(A)

