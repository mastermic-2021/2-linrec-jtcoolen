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

dec_char(c) = {
  my(x);
  if(c==0, x=32, x=select(f(e)=(e == c), arr, 1)[1] - 1 + 97);
};

decf27(s) = Strchr(vector(#s, i, dec_char(s[i])));

\\ calcule les n premiers termes d'un lfsr d'état initial u0 et polynôme feedback c
linrec(u0, c, n) = {
  v = u0;
  if(type(c)=="t_POL",c=Vecrev(c));
  print("k=0 -> " , decf27(v));
  for(k=1,n, if(k%41==0,print("k=", k ," -> ", decf27(v)),); s = -sum(i=1,#u0,v[i]*c[i]); for(i=1,#u0-1, v[i]=v[i+1]); v[#u0] = s;);
  v;
}

\\print(decf27( Mod(x, x^40 +x + u27)^n * cc2));
linrec(cc, x^40+x+u27, 2^40-1);
