library(ggplot2)

x <- rnorm(500, 0, 1)
z <- rnorm(500, 0, 1)

m1 <- cbind(x,z)
m2 <- scale(m1)
c1 <- chol(var(m2))
m3 <- m2 %*% solve(c1)

# la matrice de corr?lation que tu souhaites obtenir
cor1 <- matrix(c(1,0.8,0.8,1),2)

m4 <- m3 %*% chol(cor1)
m5 <- sweep( m4, 2, attr(m2, 'scaled:scale'), '*')
m5 <- sweep( m5, 2, attr(m2, 'scaled:center'), '+')

cor(m5)   # coefficient de corrélation

x1 <- m5[,1]
x2 <- m5[,2]

coef(lm(x2 ~ x1)) #pour avoir les coefficients de x2 = a + bx1 

qplot(x1,x2,xlab="X1",ylab="X2") + geom_abline(intercept = 0.08506075, slope = 0.82355114, color = 'red') 

v1 <- pnorm(x1,0,1)
v2 <- pnorm(x2,0,1)

coef(lm(v2 ~ v1)) #pour avoir les coefficients de v2 = a + b.v1 

qplot(v1,v2,xlab="V1",ylab="V2") + geom_abline(intercept = 0.1221782, slope = 0.8102443, color = 'red') 
# les points obtenus sont plus dispersés que pour les points (X1,X2) mais on a quasiment la meme correlation.


qlapl <- function(x) {
  #y <- 1:500
  for (i in 1:length(x)){
    if (x[i]<=1/2) {y <- log(2*x)}
    else {y <- -log(2*(1-x))}
  }
  return(y)
}

y1 <- qcauchy(v1,0,1)
y2 <- qlapl(v2)

#nuage des points
qplot(y1, y2, colour = 'red')
#vu la forme du nuage de points on remarque que la variable y1 ne 'varie' pas lorsque y2 varie, on en deduit qu'il n y a pas d'interractions entre les deux variables, d'ou correlation nulle
lm(y2~y1)
#ceci est aussi confirme par le faible coefficient de correlation (de pearson) qu'on obtient alors un coefficient de correlation de l'ordre de 1.645729e-05 qui est tres faible, 
#theoriquement, les variables V1 et V2 sont correles, du coup on peut rien dire concernant la correlation entre les deux lorsque on compose V1 et V2 par les fonctions G1 inverse et G2 inverse!.
#partie 3
data1 = read.csv('/Users/taha/Desktop/EDU/Stat1 Projet/^DJI.csv')
data2 = read.csv('/Users/taha/Desktop/EDU/Stat1 Projet/^FCHI.csv')
#on fait un merge sur date pour que les valeurs NULL apparaissent
data = merge(data1, data2, by = 'Date')
library(dplyr)
data = data%>%
  filter(Open.y != 'null')
# pour Dow jones on code la variation relative
data$var_cac <- NA
data$var_cac[1] = 0
for (i in (2:length(data$Open.x)))
{data$var_cac[i] = (data$Open.x[i]-data$Open.x[i-1])/data$Open.x[i-1]}

#la var pour CAC 40

data$var_cac1 <- NA
data$var_cac1[1] = 0
for (j in (2:length(data$Open.y)))
  {data$var_cac1[j] = (as.double(data$Open.y[j])-as.double(data$Open.y[j-1]))/as.double(data$Open.y[j-1])}


qplot(data$var_cac1, data$var_cac, xlab="Y1",ylab="Y2")  

lm(data$var_cac~data$var_cac1)
# a partir du nuage des points, on peut dire que les deux variables sont tres faiblement correles au sens de pearson


#Fonction de repartition 

dis1 = ecdf(data$var_cac1)
dis2 = ecdf(data$var_cac)
V1 = dis1(data$var_cac1)
V2 = dis2(data$var_cac)

qplot(V1,V2) + geom_abline(intercept = 0.1939, slope = 0.1939, color = 'red') 
lm(V1~V2)
#la correlation entre V1 et V2 (au sens de pearson) est beaucoup plus importantes que pour Y1 et Y2 

library(VineCopula)
ro <- BiCopEst(V1, V2, family = 1, method = "mle")
summary(ro)
#le parametre obtenue ro est 0.63, tres proche de la valeur de la regression de V2 en V1  
# pour le parametre ro = 0.63 on ecrit:
plot <- BiCop(family = 1, par = 0.63) #on simule une Copule gaussienne
plot(plot)

contour(plot)

max(V1)
max(V2)
min(V1)
min(V2)

# on une valeur de V1 ou de egale 1 

for (i in 1:length(V1))
if (V1[i] == 1) print(i) # pour i = 2001 on V1 = 1

# on elimine ces observations parce que dans le cas d'une distribution normale, X a des valeurs infinies car F(X) = 1
fus = cbind(V1,V2)
U = subset(fus, V1 !=1 & V1 != 0 &  V2 !=1 & V2 != 0)

BiCopKDE(U[,1], U[,2], type = "contour")
#les lignes de niveaux pour la densité de la copule estimée sur l’échantillon (V1, V2) est tres proche en forme du contour de la densite d'un couple gaussien de correlation ro = 0.63. On puet alors conclure que cette une bonne approximation 




