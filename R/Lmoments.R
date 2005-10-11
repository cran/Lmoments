#File:		Lmoments.R
#Version:	1.0-1
#Author:	Juha Karvanen 
#Date:		11 October 2005


Lmoments<-function(data,rmax=4)
{
if (!is.matrix(data)) data<-as.matrix(data);
n<-dim(data)[1];
p<-dim(data)[2]
x<-array(,c(p,n));
L<-array(,c(p,rmax));
for (i in 1:p)
{
x[i,]<-sort(data[,i]);
}

bcoef<-array(,c(rmax,n));
bcoefm<-array(,c(rmax,p,n));
b<-array(,c(p,rmax));

bcoef[1,]<-seq(0,1,by=(1/(n-1)));
bcoefm[1,,]<-t(array(rep(bcoef[1,],p),c(n,p))); 
b[,1]<-rowMeans(x);
b[,2]<-rowMeans(bcoefm[1,,]*x);

L[,1]=b[,1];

for (r in 2:(rmax-1))
{
    rr<-r+1;
	bcoef[r,]<-bcoef[r-1,]*seq((-(r-1)/(n-r)),1,by=(1/(n-r)));
	bcoefm[r,,]<-t(array(rep(bcoef[r,],p),c(n,p))); 
	b[,rr]<-rowMeans(bcoefm[r,,]*x);
}

for (r in 1:(rmax-1))
{
	L[,r+1]<-0;
	for (k in 0:r)
	{
	    kk<-k+1;
	    L[,r+1]<-L[,r+1]+(-1)^(r-k)*gamma(r+k+1)/(gamma(k+1)^2)/gamma(r-k+1)*b[,kk];
	}
}

L
}

Lmomcov<-function(data,rmax=4)
{
C<-t(array(c(1,0,0,0,-1,2,0,0,1,-6,6,0,-1,12,-30,20),c(4,4)))
if (!is.matrix(data)) data<-as.matrix(data);

n<-dim(data)[1];
p<-dim(data)[2]
x<-array(,c(p,n));
for (i in 1:p)
{
x[i,]<-sort(data[,i]);
}

bcoef<-array(,c(rmax,n));
b1coef<-array(,c(rmax,n));
b2coef<-array(,c(rmax,rmax,n));
b20coef<-array(,c(rmax,n));
bcoefm<-array(,c(rmax,p,n));
b<-array(,c(p,rmax));

bcoef[1,]<-seq(0,1,by=(1/(n-1)));
bcoefm[1,,]<-t(array(rep(bcoef[1,],p),c(n,p))); 
b[,1]<-rowMeans(x);
b[,2]<-rowMeans(bcoefm[1,,]*x);

for (r in 2:(rmax-1))
{
    rr<-r+1;
	bcoef[r,]<-bcoef[r-1,]*seq((-(r-1)/(n-r)),1,by=(1/(n-r)));
	bcoefm[r,,]<-t(array(rep(bcoef[r,],p),c(n,p))); 
	b[,rr]<-rowMeans(bcoefm[r,,]*x);
}

b20coef[1,]<-seq((-1)/(n-2),(n-2)/(n-2),by=1/(n-2))
b1coef[1,]<-seq(0,(n-1)/(n-2),by=1/(n-2));
for (k in 1:(rmax-1))
{
if (k>1) 
{
b20coef[k,]<-b20coef[k-1,]*seq((-k)/(n-1-k),(n-1-k)/(n-1-k),by=1/(n-1-k));
b1coef[k,]<-b1coef[k-1,]*seq((-k+1)/(n-1-k),(n-k)/(n-1-k),by=1/(n-1-k));
}
b2coef[k,1,]<-seq((-k-1)/(n-k-2),(n-k-2)/(n-k-2),by=1/(n-k-2));
for (l in 2:(rmax-1))
{
b2coef[k,l,]<-b2coef[k,l-1,]*seq((-k-l)/(n-k-l-1),(n-k-1-l)/(n-k-l-1),by=1/(n-k-l-1))
}}

theta<-array(,c(4,4))
xx<-outer(as.vector(x),as.vector(x));
xx[!upper.tri(xx)]<-0;
for (k in 0:(rmax-1))
{
for (l in 0:(rmax-1))
{
if (k>0 & l>0)
{
term1<-outer(b1coef[k,],b2coef[k,l,])
term2<-outer(b1coef[l,],b2coef[l,k,])
} 
if (k==0 & l>0)
{
term1<-outer(rep(1,n),b20coef[l,])
term2<-outer(b1coef[l,],rep(1,n))
} 
if (k>0 & l==0)
{
term1<-outer(b1coef[k,],rep(1,n))
term2<-outer(rep(1,n),b20coef[k,])
} 
if (k==0 & l==0)
{
term1<-outer(rep(1,n),rep(1,n))
term2<-outer(rep(1,n),rep(1,n))
}
term1[!upper.tri(term1)]<-0;
term2[!upper.tri(term2)]<-0;

jointbb<-sum((term1+term2)*xx)/(n*(n-1));

theta[k+1,l+1]<-b[,k+1]*b[,l+1]-jointbb;
}}
covmatrix<-C%*%theta%*%t(C);
covmatrix
}



Lcoefs<-function(data)
{
Lmom<-Lmoments(data,rmax=4);
Lmom[,4]<-Lmom[,4]/Lmom[,2];
Lmom[,3]<-Lmom[,3]/Lmom[,2];
Lmom
}



t1lmoments<-function(data)
{
rmax<-4;
if (!is.matrix(data)) data<-as.matrix(data);

n<-dim(data)[1];
p<-dim(data)[2]
x<-array(,c(p,n));
L<-array(0,c(p,rmax));
for (j in 1:p)
{
x[j,]<-sort(data[,j]);
}

i<-3:n;
s11<-c(1,(i-1)/(i-2)*(n-i)/(n-i+1));
s12<-c(1,(i-1)/(i-2)*(n-i-1)/(n-i+1));
s13<-c(1,(i-1)/(i-2)*(n-i-2)/(n-i+1));
s14<-c(1,(i-1)/(i-2)*(n-i-3)/(n-i+1));
s21<-c(1,(i-1)/(i-3)*(n-i)/(n-i+1));
s22<-c(1,(i-1)/(i-3)*(n-i-1)/(n-i+1));
s23<-c(1,(i-1)/(i-3)*(n-i-2)/(n-i+1));
s31<-c(1,(i-1)/(i-4)*(n-i)/(n-i+1));
s32<-c(1,(i-1)/(i-4)*(n-i-1)/(n-i+1));
s41<-c(1,(i-1)/(i-5)*(n-i)/(n-i+1));
s21[2]<-1;
s22[2]<-1;
s23[2]<-1;
s31[2]<-1;
s31[3]<-1;
s32[2]<-1;
s32[3]<-1;
s41[2]<-1;
s41[3]<-1;
s41[4]<-1;

c11<-choose(n-2,1)/choose(n,3)*cumprod(s11);
c12<-choose(n-2,2)/choose(n,4)*cumprod(s12);
c13<-choose(n-2,3)/choose(n,5)*cumprod(s13);
c14<-choose(n-2,4)/choose(n,6)*cumprod(s14);
c21<-choose(n-3,1)/choose(n,4)*cumprod(s21);
c22<-choose(n-3,2)/choose(n,5)*cumprod(s22);
c23<-choose(n-3,3)/choose(n,6)*cumprod(s23);
c31<-choose(n-4,1)/choose(n,5)*cumprod(s31);
c32<-choose(n-4,2)/choose(n,6)*cumprod(s32);
c41<-choose(n-5,1)/choose(n,6)*cumprod(s41);
c21[1]<-0;
c22[1]<-0;
c23[1]<-0;
c31[1]<-0;
c31[2]<-0;
c32[1]<-0;
c32[2]<-0;
c41[1]<-0;
c41[2]<-0;
c41[3]<-0;

for (j in 1:p)
{
   	L[j,1]<-sum(c11[1:(n-2)]*x[j,2:(n-1)]);
   	L[j,2]<-sum((c21[1:(n-2)]-c12[1:(n-2)])*x[j,2:(n-1)])/2;
	L[j,3]<-sum((c31[1:(n-2)]-2*c22[1:(n-2)]+c13[1:(n-2)])*x[j,2:(n-1)])/3;
	L[j,4]<-sum((c41[1:(n-2)]-3*c32[1:(n-2)]+3*c23[1:(n-2)]-c14[1:(n-2)])*x[j,2:(n-1)])/4;
}

L
}


