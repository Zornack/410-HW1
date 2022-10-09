using LinearAlgebra
# matrix = Float64[2 1 2; -2 2 -1; 4 1 2]
# matrix = Float64[2 1; 3 -1]
# matrix = Float64[6 -2 2 16; 0 -4 2 -6; 0 -12 2 -27]
# matrix = Float64[-8 1 -2; -3 -1 7; 2 -6 -1]
matrix = Float64[1 2 -1 1 6; -1 1 2 -1 3; 2 -1 2 2 14; 1 1 -1 2 8]
matrix2 = Float64[1 2 -1 1 6; -1 1 2 -1 3; 2 -1 2 2 14; 1 1 -1 2 8]
M = size(matrix,1)
N = size(matrix,2)
UTM = zeros(M,M)
# lElements =  Vector{Float64}()
offset = 0
pivot = 1

for i = 1:N-1
    print(i)
end
# for row in range(2,size(matrix,1))
#     if (matrix[row+size(matrix,1)*offset] > matrix[pivot+size(matrix,1)*offset]) 
#        global pivot = row
#     end
# end
# c = matrix[1,:]
# matrix[1,:] = matrix[pivot,:]
# matrix[pivot,:] = c
# matrix
# pivot = 1
# for x in range(2,size(matrix,1))
#     for row in range(2,size(matrix,1))
#         matrix[row+pivot-1,:] = matrix[row+pivot-1,:] - matrix[row+size(matrix,1)*(pivot-1)]/matrix[pivot+size(matrix,1)*(pivot-1)]*matrix[pivot,:]
#     end
#     global pivot = pivot + 1
# end

# for row in range(3,size(matrix,1))
#     matrix[row,:] = matrix[row,:] - matrix[row+size(matrix,1)*(pivot-1)]/matrix[pivot+size(matrix,1)*(pivot-1)]*matrix[pivot,:]
# end

# for row in range(2,size(matrix,1))
#     matrix[row,:] = matrix[row,:] - matrix[row]/matrix[pivot]*matrix[pivot,:]
# end

# for rowCount in range(2,size(matrix,1))
#     for row in range(2,size(matrix,1)-offset)
#         if (abs(matrix[row+size(matrix,1)*offset]) > abs(matrix[pivot+size(matrix,1)*offset])) 
#            global pivot = row
#         end
#     end
#     c = matrix[1+offset,:]
#     matrix[1+offset,:] = matrix[pivot,:]
#     matrix[pivot,:] = c
#     matrix
#     oldpivot = pivot
#     global pivot = 1
#     for row in range(2,size(matrix,1)-offset)
#         matrix[row+offset,:] = matrix[row+offset,:] - matrix[row+offset+size(matrix,1)*offset]/matrix[pivot+offset+size(matrix,1)*offset]*matrix[pivot+offset,:]
#     end
#     global offset = offset + 1
# end

# for row in range(2, size(matrix,1))
#     matrix[row,:] = matrix[row,:]-matrix[row]/matrix[pivot]*matrix[pivot,:]
# end
# matrix

# for rowCount in range(1,size(matrix,1))
#     matrix[rowCount,:] = matrix[rowCount,:] / matrix[rowCount,rowCount]
# end

LTM = Matrix{Float64}(I, M, M)
for rowCount in range(1,size(matrix,2)-1)
    # for rows in range(2+offset,size(matrix,1))
    #     if (abs(matrix[rows+size(matrix,1)*offset]) > abs(matrix[pivot+size(matrix,1)*offset])) 
    #         global pivot = rows
    #     end
    # end
    # if(pivot != 1)
    #     c = matrix[1+offset,:]
    #     matrix[1+offset,:] = matrix[pivot,:]
    #     matrix[pivot,:] = c
    #     global pivot = 1
    # end
    # for rows in range(2+offset, size(matrix,1))
    #     if (abs(matrix[rows, offset+1] > abs(matrix[pivot+offset,offset+1])))
    #         global pivot = rows
    #     end
    # end
    # if(pivot != 1)
    #     global c = matrix[1+offset, :]
    #     matrix[1+offset,:] = matrix[pivot,:]
    #     matrix[pivot,:] = c
    #     global pivot = 1
    # end
    for row in range(2+offset, size(matrix,1))
        LTM[row+size(matrix,1)*offset] = matrix[row+size(matrix,1)*offset]/matrix[pivot+offset+size(matrix,1)*offset]
        matrix[row,:] = matrix[row,:]-matrix[row+size(matrix,1)*offset]/matrix[pivot+offset+size(matrix,1)*offset]*matrix[pivot+offset,:]
    end
    global offset = offset + 1
end

# variables = ones(N-1)
# offset = 0
# for unknown in range(1,1)
#     variables[N-offset] = matrix[size(matrix,1),size(matrix,2)]/matrix[size(matrix,1),size(matrix,2)]
#     offset = offset + 1
# end

variables = zeros(M)
b = view(matrix, :, 5)
u = view(matrix, :, 1:4)
for i in N-1:-1:1
    variables[i] = b[i]
    for j in i+1:1:N-1
        variables[i] = variables[i]-u[i,j]*variables[j]
    end
    variables[i] = variables[i]/u[i,i]
end
matrix

