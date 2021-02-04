u27 = ffgen(('x^3-'x+1)*Mod(1,3),'u);
codf27(s) = [if(x==32,0,u27^(x-97))|x<-Vec(Vecsmall(s)),x==32||x>=97&&x<=122];

data = readvec("input.txt")[1];
c = data[2];
c1 = data[1];
n = data[3];

cc = codf27(c);
cc1 = codf27(c1);

\\ associe à un élément de F_{3^3} son nombre "sympathique": espace=0 -> 32, a=1 -> 0, ... z=u^25 -> 25
arr = [1, u27, u27^2, u27 + 2, u27^2 + 2*u27, 2*u27^2 + u27 + 2, u27^2 + u27 + 1, u27^2 + 2*u27 + 2, 2*u27^2 + 2, u27 + 1, u27^2 + u27, u27^2 + u27 + 2, u27^2 + 2, 2, 2*u27, 2*u27^2, 2*u27 + 1, 2*u27^2 + u27, u27^2 + 2*u27 + 1, 2*u27^2 + 2*u27 + 2, 2*u27^2 + u27 + 1, u27^2 + 1, 2*u27 + 2, 2*u27^2 + 2*u27, 2*u27^2 + 2*u27 + 1, 2*u27^2 + 1];

\\ associe à un élément de F_{3^3] le bon caractère ASCII
dec_char(c) = {
  my(x);
  if(c==0, x=32, x=select(f(e)=(e == c), arr, 1)[1] - 1 + 97);
};

\\ fonction de décodage
decf27(s) = Strchr(vector(#s, i, dec_char(s[i])));

\\ algo d'exponentiation rapide itératif
fast_exp(m, n) =  {
	if (n == 0, m);
	acc = matid(matsize(m)[1]);
	while (n > 1,
		if (n%2 == 0, m=m ^2; n = n / 2, acc=m*acc; m = m^2;n = (n - 1) / 2;
		);
	);
	acc = m * acc;
	acc;
};

\\ calcule les n premiers termes d'un lfsr d'état initial u0 et polynôme feedback c
\\ inefficace
linrec(u0, c, n) = {
  v = u0;
  if(type(c)=="t_POL",c=Vecrev(c));
  print("k=0 -> " , decf27(v));
\\  for(k=1,n, if(k%41==0,print("k=", k ," -> ", decf27(v)),); s = -sum(i=1,#u0,v[i]*c[i]); for(i=1,#u0-1, v[i]=v[i+1]); v[#u0] = s;);
  for(k=1,n, s = -sum(i=1,#u0,v[i]*c[i]); for(i=1,#u0-1, v[i]=v[i+1]); v[#u0] = s;);
  print(decf27(v));
}

\\print(decf27( Mod(x, x^40 +x + u27)^n * cc2));
\\linrec(cc, x^40+x+u27, ;;2^40-1);

\\ initialise la matrice encodant la relation de récurrence
m=matrix(40,40,i,j,if(j==i+1,1,0));
m[40,1]=-u27;
m[40,2]=-1;

\\ m0=m^(-n) * m_n, comme on souhaite calculer les termes en sens inverse on inverse m
m=m^(-1);
\\ d'après le cours, méthode en O(log n) pour calculer les termes suivants au bout de n itérations du lfsr
print(decf27(fast_exp(m, n) * cc~));
