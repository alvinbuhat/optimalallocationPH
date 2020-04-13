setwd("~/Desktop/Optimal/r code")

start_time <- Sys.time()

library(readxl)
library(REdaS)
ARV = read_excel("GPS.xlsx", sheet="CitiesMunicipalities")
HCF = read_excel("GPS.xlsx", sheet = "withtop10")

Lat = na.omit(ARV$LATITUDE)
Lon = na.omit(ARV$LONGITUDE)
Mun = matrix(0, nrow = length(Lat), 2)
Mun[,1] = Lat
Mun[,2] = Lon
Mun=deg2rad(Mun)

Lat_t = HCF$Latitude
Lon_t = HCF$Longitude
h = matrix(0, length(Lat_t), 2)
h[,1] = Lat_t
h[,2] = Lon_t
h = deg2rad(h)

L = HCF$lim

# no. of HCF
n=length(Lat_t)
# no. of ARV
m=length(Lat)

#k=(mean(I)^2-sd(I)^2/length(I))/(sd(I)^2-mean(I))
#k=(mean(I[1:17])^2-sd(I[1:17])^2/length(I[1:17]))/(sd(I[1:17])^2-mean(I[1:17]))
#k=0.08
k=1.561449648
#k=0.084869853
#k=1.1
A=100

I=na.omit(ARV$INFECTED)

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
d[14,30]=0
#write.csv(d, "d.csv")

fd = exp(-k*d^2)

D=rep(0,n)
for (j in 1:n){
  D[j] = sum(fd[,j]*I)
}

T = matrix(0, m, n)

e <- function(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29){
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
    T[i,11] = s10*fd[i,11]*I[i]/D[11]
    T[i,12] = s10*fd[i,12]*I[i]/D[12]
    T[i,13] = s10*fd[i,13]*I[i]/D[13]
    T[i,14] = s10*fd[i,14]*I[i]/D[14]
    T[i,15] = s10*fd[i,15]*I[i]/D[15]
    T[i,16] = s10*fd[i,16]*I[i]/D[16]
    T[i,17] = s10*fd[i,17]*I[i]/D[17]
    T[i,18] = s10*fd[i,18]*I[i]/D[18]
    T[i,19] = s10*fd[i,19]*I[i]/D[19]
    T[i,20] = s10*fd[i,20]*I[i]/D[20]
    T[i,21] = s10*fd[i,21]*I[i]/D[21]
    T[i,22] = s10*fd[i,22]*I[i]/D[22]
    T[i,23] = s10*fd[i,23]*I[i]/D[23]
    T[i,24] = s10*fd[i,24]*I[i]/D[24]
    T[i,25] = s10*fd[i,25]*I[i]/D[25]
    T[i,26] = s10*fd[i,26]*I[i]/D[26]
    T[i,27] = s10*fd[i,27]*I[i]/D[27]
    T[i,28] = s10*fd[i,28]*I[i]/D[28]
    T[i,29] = s10*fd[i,29]*I[i]/D[29]
    T[i,30] = (A-s1-s2-s3-s4-s5-s6-s7-s8-s9-s10-s11-s12-s13-s14-s15-s16-s17-s18-s19-s20-s21-s22-s23-s24-s25-s26-s27-s28-s29)*fd[i,30]*I[i]/D[30]
    e = e + (sum(T[i,])/I[i] - A/sum(I))^2
  }
  e=1e5*e
  return(e)
}

ee <- function(s) e(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15],s[16],s[17],s[18],s[19],s[20],s[21],s[22],s[23],s[24],s[25],s[26],s[27],s[28],s[29])

AA=matrix(0,nrow=(m+2*n),ncol=n-1)
b=rep(0,(m+2*n))

for (i in 1:m){
  for (j in 1:(n-1)){
    AA[i,j] = (fd[i,30]/D[30])-(fd[i,j]/D[j])
  }
  b[i] = A*fd[i,30]/D[30]-1
}
AA[(m+1):(m+n-1),]=diag(n-1)
AA[(m+n),]=rep(0,n-1)-1
b[m+n] = -A
p1=m+n+1
p2=m+2*n-1
AA[p1:p2,]=-diag(n-1)
b[p1:p2] = -M*L[1:29]
AA[(m+2*n),]=rep(1,29)
b[m+2*n]=A-M*L[30]
result <- constrOptim(rep(floor(80/n),n-1), ee, NULL, AA, b)
s=matrix(0, n,1)
s[1:(n-1),] = result$par
s[n,1]=A-sum(result$par)
#write.csv(t(s), "sj.csv")
#print(s)

A=300000
s=A/100*s
#print(s)

s_j = data.frame( testing_center = HCF$`Testing Center`, S_j = s)
print(s_j)

end_time <- Sys.time()
print(end_time - start_time)

