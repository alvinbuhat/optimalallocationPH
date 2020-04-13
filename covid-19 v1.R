setwd("~/Desktop/Optimal/r code")

start_time <- Sys.time()
library(readxl)
library(REdaS)
ARV = read_excel("GPS.xlsx", sheet="CitiesMunicipalities")
HCF = read_excel("GPS.xlsx", sheet = "Testing Centers")

Lat = na.omit(ARV$LATITUDE)
Lon = na.omit(ARV$LONGITUDE)
Mun = matrix(0, nrow = length(Lat), 2)
Mun[,1] = Lat
Mun[,2] = Lon
Mun=deg2rad(Mun)

Lat_t = na.omit(HCF$Latitude)
Lon_t = na.omit(HCF$Longitude)
h = matrix(0, length(Lat_t), 2)
Ll=na.omit(HCF$Limit)


h[,1] = Lat_t
h[,2] = Lon_t
h = deg2rad(h)

# no. of HCF
n=length(Lat_t)
# no. of ARV
m=length(Lat)

L=matrix(0, n,1)
for (j in 1:n){
  L[j,1] = Ll[j]
}

#k=(mean(I[1:17])^2-sd(I[1:17])^2/length(I[1:17]))/(sd(I[1:17])^2-mean(I[1:17]))
k=1.561449648
#k=0.084869853
A=100

I=na.omit(ARV$INFECTED)
M=sum(I)/sum(L)

for (i in 1:m){
  if (I[i]==0){
    I[i]=1e-5
  }
}

d = matrix(0, m, n)
for (i in 1:m){
  for (j in 1:n){
    d[i,j] = 6731*acos(cos(Mun[i,1])*cos(h[j,1])*cos(h[j,2]-Mun[i,2])+sin(Mun[i,1])*sin(h[j,1]))
  }
}
#write.csv(d, "d.csv")
fd = exp(-k*d^2)

D=rep(0,n)
for (j in 1:n){
  D[j] = sum(fd[,j]*I)
}

T = matrix(0, m, n)

e <- function(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10){
  e=0
  for (i in 1:m){
    T[i,1] = s1*fd[i,1]*I[i]/D[1]
    T[i,2] = s2*fd[i,2]*I[i]/D[2]
    T[i,3] = s3*fd[i,3]*I[i]/D[3]
    T[i,4] = s4*fd[i,4]*I[i]/D[4]
    T[i,5] = s5*fd[i,5]*I[i]/D[5]
    T[i,6] = s6*fd[i,6]*I[i]/D[6]
    T[i,7] = s7*fd[i,7]*I[i]/D[7]
    T[i,8] = s8*fd[i,8]*I[i]/D[8]
    T[i,9] = s9*fd[i,9]*I[i]/D[9]
    T[i,10] = s10*fd[i,10]*I[i]/D[10]
    T[i,11] = (A-s1-s2-s3-s4-s5-s6-s7-s8-s9-s10)*fd[i,11]*I[i]/D[11]
    e = e + (sum(T[i,])/I[i] - A/sum(I))^2
  }
  e=1e5*e
  return(e)
}

ee <- function(s) e(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10])

AA=matrix(0,nrow=(m+2*n),ncol=n-1)
b=rep(0,(m+2*n))

for (i in 1:m){
  for (j in 1:(n-1)){
    AA[i,j] = (fd[i,11]/D[11])-(fd[i,j]/D[j])
  }
  b[i] = A*fd[i,11]/D[11]-1
}
AA[(m+1):(m+n-1),]=diag(n-1)
AA[(m+n),]=rep(0,n-1)-1
b[m+n] = -A
p1=m+n+1
p2=m+2*n-1
AA[p1:p2,]=-diag(n-1)
b[p1:p2] = -M*L[1:10]
AA[(m+2*n),]=rep(1,10)
b[m+2*n]=A-M*L[11]
result <- constrOptim(rep(9,10), ee, NULL, AA, b)
s=matrix(0, n,1)
s[1:(n-1),] = result$par
s[n,1]=A-sum(result$par)
#write.csv(t(s), "sj.csv")
#print(s)

A=300000
s=A/100*s
#print(s)

tc = c("RITM", "San Lazaro", "UP-NIH", "Lung Center", "Baguio Gen", "TMC, Pasig", "St. Luke's, QC", "Bicol PHL", "WVMC", "VSMMC", "SPMC")
s_j = data.frame( testing_center = tc, S_j = s)
print(s_j)

end_time <- Sys.time()
print(end_time - start_time)

