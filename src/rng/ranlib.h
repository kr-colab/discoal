/* Prototypes for all user accessible RANLIB routines */

extern void advnst(long k);
extern double genbet(double aa,double bb);
extern double genchi(double df);
extern double genexp(double av);
extern double genf(double dfn, double dfd);
extern double gengam(double a,double r);
extern void genmn(double *parm,double *x,double *work);
extern void genmul(long n,double *p,long ncat,long *ix);
extern double gennch(double df,double xnonc);
extern double gennf(double dfn, double dfd, double xnonc);
extern double gennor(double av,double sd);
extern void genprm(long *iarray,int larray);
extern double genunf(double low,double high);
extern void getsd(long *iseed1,long *iseed2);
extern void gscgn(long getset,long *g);
extern long ignbin(long n,double pp);
extern long ignnbn(long n,double p);
extern long ignlgi(void);
extern long ignpoi(double mu);
extern long ignuin(long low,long high);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern void phrtsd(char* phrase,long* seed1,long* seed2);
extern double ranf(void);
extern void setall(long iseed1,long iseed2);
extern void setant(long qvalue);
extern void setgmn(double *meanv,double *covm,long p,double *parm);
extern void setsd(long iseed1,long iseed2);
extern double sexpo(void);
extern double sgamma(double a);
extern double snorm(void);

