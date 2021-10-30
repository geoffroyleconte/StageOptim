using MAT, LinearAlgebra
file = matopen(raw"C:\Users\Geoffroy Leconte\Documents\doctorat\code\datasets\ball-topple\17MFBallManLog.mat")
varnames = names(file)
d = read(file, "keaLog")