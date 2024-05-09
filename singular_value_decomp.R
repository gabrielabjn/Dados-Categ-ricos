# SINGULAR VALUE DECOMPOSITION --------------------------------------------------

A<-matrix(c(395,147,694,2456,152,327,1758,916,1347),3)
U<-rep(1,3)
U<-as.matrix(rep(1,3),3)

n<-t(U)%*%A%*%U
#ou
n<-sum(A)

P <- (1/n[[1]])*A #matriz de proporcoes

r<-P%*%U # definicao da massa de linhas
c<-t(P)%*%U # definicao da massa de colunas

Dr<-diag(as.vector(r))
Dc<-diag(as.vector(c)) 

S<- solve(Dr^(1/2))%*%(P-r%*%t(c))%*%solve(Dc^(1/2))
S

S2<-sum(diag(S%*%t(S)))

Vt<-t(U)%*%U

dec <- svd(S)

# ------------------------------------------------------------------------------

R<-solve(Dr)%*% P # PERIFL DE LINHA
C<-solve(Dc%*%t(P)) # PERFIL DE COLUNA

centroidedeC<-t(C)%*%c
# verifique que eh igual a massa de colunas

phi <- solve(Dr^(1/2))%*%dec$u 
tau <- solve(Dc^(1/2))%*%dec$v

Fik<-phi%*%diag(dec$d)
Gjk<-tau%*%diag(dec$d)

plot(x=Fik[,1],y = Fik[,2], pch = 3, ylim = c(-0.4,0.4), xlim = c(-1,1))
points(x=Gjk[,1],y = Gjk[,2], pch = 19)

# fraude acontece mais em grande oslo
# fraude apresenta uma associação maior com grande oslo
