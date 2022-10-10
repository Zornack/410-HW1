using LinearAlgebra

function upper_triangulate!(A, LTM)

    for j = 1:N
        pivot!(A,j)
        LTM = eliminate!(A, j, LTM)
    end

    return LTM
end

function eliminate!(A, j, LTM)
    pivot = A[j, j]
    pivot_row = A[j, :]

    if(size(pivotMatricies)[1] > 1)
        LTM = pivotMatricies[size(pivotMatricies)[1]] * LTM
    end

    for k = j+1:M
        fac = A[k, j]/pivot
        LTM[k,j] = fac
        A[k, :] = A[k, :] - fac*pivot_row
    end

    return LTM

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
    else
        pivotMatrix = Matrix{Float64}(I, M, N)
        push!(pivotMatricies, pivotMatrix)
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
    LTM = zeros(M,N)
    LTM = upper_triangulate!(UTM,LTM)
    PM = pivotMultiply()
    for i in 1:M
        LTM[i,i] = 1
    end
    return UTM, LTM, PM
end

function backwardsSub(A, b)
    N = size(A,2)
    x = zeros(N)
    for i in N:-1:1
        x[i] = b[i]
        for j in i+1:1:N
            x[i] = x[i]-A[i,j]*x[j]
        end
        x[i] = x[i]/A[i,i]
    end
    return x
end

function forwardsSub(A, b)
    N = size(A,2)
    x = zeros(N)
    for i in 1:N
        x[i] = b[i]
        for j in i-1:-1:1
            x[i] = x[i]-A[i,j]*x[j]
        end
        x[i] = x[i]/A[i,i]
    end
    return x
end

function LUPsolve(A,b) 
    UTM, LTM, PM = computeLUP(A)

end




# A = Float64[-2 2 -1; 6 -6 7; 3 -8 4]
# A = Float64[1 2 3; 4 5 6; 7 8 0]
# A = Float64[0 1 2 1; 1 0 0 0; 2 1 2 1; 1 2 4 3]
# A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
# A = Float64[2 4 2; 4 -10 2; 1 2 4]
A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
M = size(A,1)
N = size(A,2)
pivotMatricies = []
# L = zeros(M,N)

UTM, LTM, PM = computeLUP(A)

