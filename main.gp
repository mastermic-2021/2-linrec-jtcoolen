u27 = ffgen(('x^3-'x+1)*Mod(1,3),'u);
codf27(s) = [if(x==32,0,u27^(x-97))|x<-Vec(Vecsmall(s)),x==32||x>=97&&x<=122];


c = "snpmmgzhfsqrmevdkajq hvxxnixufdeiheiurgq";
c1 = " abcdefghijklmnopqrstuvwxyz";
c2 = "lkttoafqagecglgkflpanbgkacwudnvnxreaspmxsnpmmgzhfsqrmevdkajq hvxxnixufdeiheiurgq";

c3 = "snpmmgzhfsqrmevdkajq hvxxnixufdeiheiurgqlkttoafqagecglgkflpanbgkacwudnvnxreaspmx";
n = 929583887302112;
\\ n = 50;

cc = Vecrev(codf27(c));
cc1 = codf27(c1);
cc2 = codf27(c2);
cc3 = codf27(c3);

\\ le polynôme x^40+x+u, inversé
pol = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 1,u];
\\print(#cc, " ", #pol);
print("#cc=", #cc);
print("cc=", cc);
\\print(pol);


\\ associe à un élément de F_{3^3} son nombre "sympathique": espace=0 -> 32, a=1 -> 0, ... z=u^25 -> 25
dec_char(c) = {
my(x);
if(c== 0, x=32,
if(c== 1, x=0,
if(c== u27, x=1,
if(c== u27^2, x=2,  
if(c== u27 + 2, x=3,
if(c== u27^2 + 2*u27, x=4,
if(c== 2*u27^2 + u27 + 2, x=5,
if(c== u27^2 + u27 + 1, x=6,
if(c== u27^2 + 2*u27 + 2, x=7,
if(c== 2*u27^2 + 2, x=8,
if(c== u27 + 1, x=9,
if(c== u27^2 + u27, x=10,
if(c== u27^2 + u27 + 2, x=11,
if(c== u27^2 + 2, x=12,
if(c== 2, x=13,
if(c== 2*u27,x= 14,
if(c== 2*u27^2, x=15,
if(c== 2*u27 + 1, x=16,
if(c== 2*u27^2 + u27, x=17,
if(c== u27^2 + 2*u27 + 1, x=18,
if(c== 2*u27^2 + 2*u27 + 2, x=19,
if(c== 2*u27^2 + u27 + 1, x=20,
if(c== u27^2 + 1, x=21,
if(c== 2*u27 + 2, x=22,
if(c== 2*u27^2 + 2*u27, x=23,
if(c== 2*u27^2 + 2*u27 + 1, x=24,
if(c== 2*u27^2 + 1, x=25, 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0);x; };

decf27(s) = Strchr(vector(#s, i, if(dec_char(s[i] )==32, 32, dec_char(s[i]) + 97 )));

\\ calcule les n premiers termes d'un lfsr d'état initial u0 et polynôme feedback c
linrec(u0, c, n) = {
  v = u0;
  print("v=", #v);
  print("c=",#c);
  print("dec v=" , decf27(v));
  for(k=1,n, print("k=", k ," -> ", decf27(v)); s = -sum(i=1,#u0,v[i]*c[i]); for(i=1,#u0-1, v[i]=v[i+1]); v[#u0] = s;);
  v;
}

\\ fonction du cours, bien trop gourmand en mémoire
linrec2(u0,c,n) = {
  if(type(c)=="t_POL",c=Vecrev(c));
  d = #c - 1;
  if(type(u0)=="t_POL",u0=Vecrev(u,d));
  v = Vec(u0,n);
  for(k=d+1,n,  print("k=", k ," -> ", decf27(v));  v[k] = -sum(i=0,d-1,v[k-d+i]*c[i+1]));
  v;
}

\\print(decf27(cc1));
linrec(cc, pol, 496561832237/40);
