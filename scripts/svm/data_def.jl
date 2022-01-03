using LinearAlgebra, DelimitedFiles, FileIO, Random

# read data
path = raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\classif"
ion_dat = readdlm(string(path, "/ionosphere.data"), ',')

function partitionTrainTest(data, at = 0.7)
  n = size(data, 1)
  idx = shuffle(1:n)
  train_idx = view(idx, 1:floor(Int, at*n))
  test_idx = view(idx, (floor(Int, at*n)+1):n)
  return data[train_idx,:], data[test_idx,:]
end
ion_dat_train, ion_dat_test = partitionTrainTest(ion_dat, 0.7)

X_train, X_test = Matrix{Float64}(ion_dat_train[:, 1:end-1]), Matrix{Float64}(ion_dat_test[:, 1:end-1])
y_train, y_test = Float64[class == "g" ? 1 : -1 for class in ion_dat_train[:, end]], Float64[class == "g" ? 1 : -1 for class in ion_dat_test[:, end]]

n_dat_train, n_feat = size(X_train)
Krbf(x1, x2; γ = 1/n_feat) = exp(-γ * norm(x1 .- x2, 2)^2)
A = zeros(1, n_dat_train)
A[1, :] .= y_train
b = [0.]
nvar = n_dat_train
Cp = 10.
c = .-ones(nvar)
lvar, uvar = zeros(nvar), Cp .* ones(nvar) 
