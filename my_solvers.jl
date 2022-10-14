
using Plots

"""
    upper_triangulate!(A, LTM)
Takes in matrix A and zero matrix of the same dimensions of A. Computes LU factorization,
modifying A into the upper triangle matrix and returning the lower triangle matrix. 
"""
function upper_and_lower_triangulate!(A, LTM)
    M = size(A,1)
    N = size(A,2)

    for j = 1:N
        pivot!(A,j)
        LTM = eliminate!(A, j, LTM)
    end

    return LTM
end

"""
    pivot!(A,j)
Partially pivots matrix A, putting the largest value below A[j,j]
into pivot position. Also stores the pivot operation applied to A
in order to pivot the LTM. 
"""
function pivot!(A,j)
    M = size(A,1)
    N = size(A,2)
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
        pivotMatrix = makeIdentity(M,N)
        swap!(pivotMatrix, pivot_row, new_pivot_row)
        push!(pivotMatricies, pivotMatrix)
        swap!(A,pivot_row,new_pivot_row)
    else
        pivotMatrix = makeIdentity(M,N)
        push!(pivotMatricies, pivotMatrix)
    end
end

"""
    eliminate!(A, j, LTM)
Eliminates the elements in the pivot row below the pivot element. 
Also pivots the lower triangular matrix based on the past elimination's
pivot. Returns the pivoted lower triangular matrix. 
"""
function eliminate!(A, j, LTM)
    M = size(A,1)
    N = size(A,2)
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

"""
    makeIdentity(m,n)
Creates and returns an identity matrix of size M,N.
"""
function makeIdentity(m,n)
    iden = zeros(m,n)
    for i = 1:m
        iden[i,i] = 1
    end
    return iden
end

"""
    swap!(A,j,k)
Swaps rows j and k of matrix A. 
"""
function swap!(A,j,k)
    A[j,:], A[k,:] = A[k,:], A[j,:]
end

"""
    pivotMultiply()
Multiplies through all the pivot matricies used during elimination
to calculate and return the final pivot matrix of a LUP factorization. 
"""
function pivotMultiply()
    PM = pop!(pivotMatricies)
    for i = 1:size(pivotMatricies)[1]
        PM = PM * pop!(pivotMatricies)
    end
    return PM
end

"""
    computeLUP(A)
Computes the LUP facotirzation of matrix A. Returns the upper teriangular,
lower triangluar and pivot matrix. 
"""
function computeLUP(A)
    M = size(A,1)
    N = size(A,2)
    UTM = copy(A)
    LTM = zeros(M,N)
    LTM = upper_and_lower_triangulate!(UTM,LTM)
    PM = pivotMultiply()
    for i in 1:M
        LTM[i,i] = 1
    end
    return UTM, LTM, PM
end

"""
    backwardsSub(A, b)
Performs backwards substition on matrix [A b], returning the solution as 
the vector x.
"""
function backwardsSub(A, b)
    M = size(A,1)
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

"""
    forwardsSub(A, b)
Performs forward substition on matrix [A b], returning the solution as 
the vector x.
"""
function forwardsSub(A, b)
    M = size(A,1)
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

"""
    LUPsolve(A,b) 
Solves A*x=b for the matrix [A b] by first computing the LUP factorization of A
and then solving L*ƀ=y and U*y = x. Returns the vector x. 
"""
function LUPsolve(A,b) 
    N = size(A,1)
    UTM, LTM, PM = computeLUP(A)
    y = forwardsSub(LTM, PM*b)
    x = backwardsSub(UTM, y)
    return x
end

"""
    testMatrix(N)
Creates an N by N text matrix A and random vector b for
testing purposes. Returns A and b. 
"""
function testMatrix(N)
    B = rand(N,N)
    A = transpose(B)*B+makeIdentity(N,N)
    b = rand(N,1)
    return A,b
end

"""
    confirmAccuracy(N)
Creates a random N by N matrix A and random N by 1 vector b. 
Performs LUP factorization on A and solves [A b], then confirms their 
accuracy by asserting that LTM * UTM ≈ PM*A and that A*x-b ≈ 0. 
"""
function confirmAccuracy(N)
    B = rand(N,N)
    A = transpose(B)*B+makeIdentity(N,N)
    b = rand(N,1)
    UTM, LTM, PM = computeLUP(A)
    x = LUPsolve(A,b)
    @assert LTM*UTM ≈ PM*A
    @assert ((A*x-b).+1) ≈ (zeros(N,1).+1)
end

"""
    timeCompute(N)
Computes the time to compute the LUP factorization of 
5 N by N matricies and stores the matrix sizes in xCompute
and the time to solve in yCompute.
"""
function timeCompute(N)
    A,b = testMatrix(N)
    computeLUP(A)
    for i in 1:5
        push!(xCompute,N)
        A,b = testMatrix(N)
        push!(yCompute, @elapsed computeLUP(A))
    end
end

"""
    timeCompute(N)
Computes the time to solve 5 N by N matricies 
and stores the matrix sizes in xSolves and the time to solve inySolves.
"""
function timeSolves(N)
    A,b = testMatrix(N)
    LUPsolve(A,b)
    for i in 1:5
        push!(xSolves, N)
        A,b = testMatrix(N)
        push!(ySolves, @elapsed LUPsolve(A,b))
    end
end

pivotMatricies = []
xCompute = Float64[]
yCompute = Float64[]
xSolves = Float64[]
ySolves = Float64[]