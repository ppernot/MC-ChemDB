# Generate HR photon flux

lambda = seq(50,200,by=0.1)
flux = lambda * 0
flux[lambda > 73.5 & lambda <= 74.5] = 1

write.csv(cbind(lambda,flux),
          file='surf73_HR.txt',
          sep = " ",
          row.names = FALSE)
