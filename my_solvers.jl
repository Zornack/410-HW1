

"""
upper_triangulate!(A, LTM)
Compute the product of matrices ‘a‘ and ‘b‘ and store in ‘c‘.
"""
function upper_triangulate!(A, LTM)
    M = size(A,1)
    N = size(A,2)

    for j = 1:N
        pivot!(A,j)
        LTM = eliminate!(A, j, LTM)
    end

    return LTM
end


"""
eliminate!(A, j, LTM)
Compute the product of matrices ‘a‘ and ‘b‘ and store in ‘c‘.
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

function makeIdentity(m,n)
    iden = zeros(m,n)
    for i = 1:m
        iden[i,i] = 1
    end
    return iden
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
    M = size(A,1)
    N = size(A,2)
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

function LUPsolve(A,b) 
    N = size(A,1)
    UTM, LTM, PM = computeLUP(A)
    y = forwardsSub(LTM, PM*b)
    x = backwardsSub(UTM, y)
    return x
    @assert LTM*UTM ≈ PM*A
    (@assert (A*x-b).+1) ≈ (zeros(N,1)+.1)
end

function testMatrix()
    B = rand(N,N)
    A = transpose(B)*B+makeIdentity(N,N)
    b = rand(N,1)
    return B,A,B
end
function confirmAccuracy(N)
    B = rand(N,N)
    A = transpose(B)*B+makeIdentity(N,N)
    b = rand(N,1)
    UTM, LTM, PM = computeLUP(A)
    x = LUPsolve(A,b)
    @assert LTM*UTM ≈ PM*A
    @assert ((A*x-b).+1) ≈ (zeros(N,1).+1)
    # result = A*x-b
    # for i = 1:N
    #     if result[i] > 0
    #         if result[i] < 1e-10
    #             result[i] = 0
    #         end
    #     else
    #         if result[i] > -1e-10
    #             result[i] = 0
    #         end
    #     end
    # end
    # @assert result == zeros(N,1)
    # return LTM, UTM, PM, A, b, x
end



# A = Float64[-2 2 -1; 6 -6 7; 3 -8 4]
# A = Float64[1 2 3; 4 5 6; 7 8 0]
# A = Float64[0 1 2 1; 1 0 0 0; 2 1 2 1; 1 2 4 3]
# A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
# A = Float64[2 4 2; 4 -10 2; 1 2 4]
# A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
# A = Float64[6 -2 2; 12 -8 6; 3 -13 3]
# A = Float64[1 2 4; 2 1 3; 3 2 4]
# A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
# A = Float64[2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
# A = Float64[10 -7 0; -3 2 6; 5 -1 5]
# A = Float64[1 1 1; 0 2 5; 2 5 -1]
# b = Float64[6; -4; 27]
pivotMatricies = []
# L = zeros(M,N)

# UTM, LTM, PM = computeLUP(A)

# x = LUPsolve(A,b)

# global time1 = 0
# for i=1:5
#     global time1 = time1 + @elapsed confirmAccuracy(10)
# end
# time1 = time1/5

# global time2 = 0
# for i=1:5
#     global time2 = time2 + @elapsed confirmAccuracy(100)
# end
# time2 = time2/5

# global time3 = 0
# for i=1:5
#     global time3 = time3 + @elapsed confirmAccuracy(1000)
# end
# time3 = time3/5

# A*x-b ≈ u

# @assert ((A*x-b).+1) ≈ (u.+1)