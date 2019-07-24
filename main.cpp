/*

	Lyapunov images based on the article by Mario Markus from Spektrum der Wissenschaft in 1995

	Implementation by Marc Meidlinger, july 2019

*/

#include "stdio.h"
#include "time.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "stdint.h"


// const definitions

// supported Lyapunov functions
const int32_t ID_FKT_I=1;
const int32_t ID_FKT_II=2;
const int32_t ID_FKT_SICO=3;
const int32_t ID_FKT_III=7;
const int32_t ID_FKT_VII=11;
const int32_t ID_FKT_IX=13;
const int32_t ID_FKT_X=14;
const int32_t ID_FKT_METADET=16;
const int32_t ID_FKT_METAABSC=17;
const int32_t ID_FKT_LSIN=18;
const int32_t ID_FKT_ATAN=22;

const char COLORCOLLECTIONDIR[]="COLORCOLLECTION\\";

const double PI05=0.5*M_PI;

enum { WAS_F=1, WAS_ABL };

enum {
	FKTTYP_NORMAL=1,
	FKTTYP_ABSCHNITTSWEISE,
	FKTTYP_DETACHED,
	FKTTYP_METADET
};

const int32_t MAXRGBITERS=64;
const int32_t MAXINTANZ=32;
const int32_t ID_FAERBUNG_INTERVALL=2;


// struct definitions

struct Bitmap {
	int32_t xlen,ylen,bytes,ybytes;
	uint8_t* bmp;

	Bitmap();
	virtual ~Bitmap();
	int32_t setlenxy(const int32_t,const int32_t);
	void save(const char*);
	void disp(void);
};

struct IterDouble {
	int32_t anz,nr;
	double wert,delta,l,r;

	IterDouble(const double,const double,const int32_t);
	virtual int32_t iterStart(void);
	virtual int32_t iterWeiter(void);
};

struct Function {
	int32_t id;
	int32_t typ;
	IterDouble* iterb;

	virtual void eval(const double,const double,double&) { };
	virtual void eval(const double,const double,double&,double&) { };
	virtual void evalabl(const double,const double,double&) { }

	virtual void save(FILE *) { };
	virtual int32_t load(const int32_t,FILE *) { return 0; };
	virtual int32_t iterStart(void) { return 0; };
	virtual int32_t iterWeiter(void) { return 0; }
	virtual char* fktStr(char*) { return 0; }
	virtual char* ablStr(char*) { return 0; }
	virtual void set_iterb(IterDouble* adr) { iterb=adr; }
	virtual void set_b(const double) { return; }
};

struct FunctionI : public Function {
	// f=r*x*(1-x)
	// g=f'
	FunctionI();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual int32_t load(const int32_t,FILE *);
	virtual int32_t iterStart(void) { return 1; }; // einmal geht
	virtual int32_t iterWeitloer(void) { return 0; } // aber nicht weiter
	virtual char* fktStr(char* s) { sprintf(s,"N(%i) r*x*(1-x)",id); return s; }
	virtual char* ablStr(char* s) { sprintf(s,"N(%i) [=f'] r-2rx",id); return s; }
};

struct FunctionII : public Function {
	// f=b*sin²(x+r)
	// g=f'
	double b,b2;

	FunctionII();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual int32_t load(const int32_t aid,FILE *);
	virtual void set_b(const double d) { b=d; b2=d+d; }
	virtual int32_t iterStart(void);
	virtual int32_t iterWeiter(void);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionVII : public FunctionII {
	// f=b*sin^2(x+r)
	// g=r-2rx

	FunctionVII();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double,const double,double&);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionMetaABSC : public FunctionII {
	// sectionally defined function
	// f=fint or fext depending on x and i0
	// g=derivative of f depending on x and i1
	double I0MIN,I1MIN;
	double I0MAX,I1MAX;
	Function *fint,*fext;

	FunctionMetaABSC();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	void setfint(Function* p) { fint=p; }
	void setfext(Function* p) { fext=p; }
	void setsections(const double a,const double b,const double c,const double d) { I0MIN=a; I0MAX=b; I1MIN=c; I1MAX=d; }
	virtual void set_b(const double d);
	virtual char* ablStr(char* s);
	virtual int32_t load(const int32_t aid,FILE *);
	virtual int32_t iterStart(void);
	virtual int32_t iterWeiter(void);
};

struct FunctionX : public FunctionII {
	// f=r*sin^2(x-r)+b*sin^2(x+2r)
	// g=rx-b*fin^2(r*x-b)
	FunctionX();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionMetaDet : public FunctionII {
	// f=pointer to f
	// g=pointer to abl
	Function *f, *abl;
	int32_t fwas,ablwas;

	FunctionMetaDet();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	void setF(Function* p,const int32_t a) { f=p; fwas=a; }
	void setAbl(Function* p,const int32_t a) { abl=p; ablwas=a; }
	virtual void set_b(const double d);
	virtual char* ablStr(char* s);
	virtual int32_t load(const int32_t aid,FILE *);
	virtual int32_t iterStart(void);
	virtual int32_t iterWeiter(void);
};

struct FunctionIX : public FunctionII {
	// f=b*sin(x+r)+b*sin^2(b*x+r)
	// g=sin^2(x+r*b)-r*x
	FunctionIX();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionSICO : public FunctionII {
	// f=b*sin(x+r*cos(x+r))
	// g=f'
	FunctionSICO();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double,const double,double&);
	virtual void save(FILE *);
	virtual int32_t load(const int32_t aid,FILE *);
	virtual int32_t iterStart(void);
	virtual int32_t iterWeiter(void);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionLSIN : public FunctionII {
	// f=r*sin(x)*(1-b*sin(x+r))
	// g=f'
	FunctionLSIN();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionATAN : public FunctionII {
	// f=b*atan((x+r)*sin(x+r))
	// g=f'
	FunctionATAN();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double,const double,double&);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct FunctionIII : public FunctionII {
	// f=b*sin(x+r)*sin(x-r);
	FunctionIII();
	virtual void eval(const double,const double,double&);
	virtual void eval(const double,const double,double&,double&);
	virtual void evalabl(const double x,const double r,double& abl);
	virtual void save(FILE *);
	virtual char* fktStr(char* s);
	virtual char* ablStr(char* s);
};

struct ColIntv {
	double gl,gr,breite;
    int32_t lr,lg,lb,rr,rg,rb,dr,dg,db; 

	ColIntv() { gr=gl=breite=0; };
	void setgrenzel(const double w) { gl=w; breite=gr-gl; }
    void setgrenzer(const double w) { gr=w; breite=gr-gl; }
    void setfarbel(const int32_t,const int32_t,const int32_t);
    void setfarber(const int32_t,const int32_t,const int32_t);
    void precalc(void);
    virtual int32_t farbe(const double,int32_t&,int32_t&,int32_t&);
    virtual void save(FILE *f);
    virtual int32_t load(FILE *f);
};

struct IntervalColoring {
	int32_t id;

    ColIntv *ints[MAXINTANZ];
	double mingl,maxgl;
    int32_t intanz;
    int32_t lr,lg,lb,rr,rg,rb;
    double minwert,maxwert;

	IntervalColoring() { intanz=0; mingl=maxgl=0; id=ID_FAERBUNG_INTERVALL; };
    virtual ~IntervalColoring();
	virtual void save(FILE*);
    virtual int32_t load(const int32_t,FILE*);
    virtual int32_t farbe(const double,int32_t&,int32_t&,int32_t&);
	int32_t Addintervall(ColIntv*);
	void clear();
    void setfarbel(const int32_t ar,const int32_t ag,const int32_t ab) { lr=ar; lg=ag; lb=ab;  }
    void setfarber(const int32_t ar,const int32_t ag,const int32_t ab) { rr=ar; rg=ag; rb=ab;  }
};

struct IntRect {
	int32_t top,left,bottom,right;
};

struct Point32_t {
	double x,y;
};

struct Ljapunow {
	Function *fkt;
	IntervalColoring *farbe;
    int32_t lenx,leny,iter0,iter1,seqlen;
    int32_t iter0h,iter1h,iter1d;
    double INViter1d;
    char fn[1024];
    char sequence[256];
    double* exps;
    double x0;
    Point32_t upperleft,lowerleft,lowerright;
	IterDouble* iterC; 

    Ljapunow();
    virtual ~Ljapunow();

    int32_t calc(const int32_t start,const int32_t ende);
	int32_t iterStart(void);
	int32_t iterWeiter(void);

    int32_t loadpar(char *fn);
    int32_t loadexp(char *fn);
    int32_t loadcolor(char *fn);
	char* getSequence(char* s);
    void createBmp(Bitmap*);
    void setfarbe(IntervalColoring*);
    void setFunction(Function *f) { fkt=f; }
    void setSequence(char *s);
    void setPosition(const double,const double,const double,const double,const double,const double);
    void setlen(const int32_t xl,const int32_t yl);
    void setiter(const int32_t i0,const int32_t i1);
    // saving values
	void savepar(char *fn);
    void saveexp(char *fn);
    void savebmp(char *fn,Bitmap*);
	void savedescr(const char* fn);
	// rect manipulations
	void crop(const int32_t,const int32_t,const int32_t,const int32_t);
	void tile(char*,const int32_t,const int32_t);
	void centerPixel(const int32_t,const int32_t);
	void rot(const int32_t);
	void stretch(const double,const double);
};


// forward declarations

Function* loadFunction(FILE*);
Function* getNewFunction(const int32_t);
IntervalColoring* loadfaerbung(FILE*);
inline double fastsin(double);
inline double fastcos(double);

char dez(const char c);
char* stripext(char*);
char* upper(char*);
void writehex(FILE*,const char*);
inline double maximumD(const double,const double);
inline int32_t maximumI(const int32_t,const int32_t);
char* chomp(char *);
char* upper(char*);
char* removeStr(const char*,const char*,char*);
int32_t getFirstColorFile(char*);
int32_t getNextColorFile(char*);


// defines as small functions

#define SEQPOSINC(POS)\
	POS++;\
    if (POS>=seqlen) POS=0;


// globals

Ljapunow* ljap=NULL;
FILE *ffarbe=NULL;


// function definitions

// load functions

IntervalColoring* loadfaerbung(FILE *f) {
	int32_t aid;
	if (!f) return 0;
    fscanf(f,"ID\n%i\n",&aid);
    IntervalColoring *p=NULL;
    if (aid==ID_FAERBUNG_INTERVALL) p=new IntervalColoring;
	if (p) p->load(aid,f);

	return p;
};

Function* getNewFunction(const int32_t aid) {
	Function* p=NULL;
    switch (aid) {
        case ID_FKT_I: p=new FunctionI(); break;
        case ID_FKT_II: p=new FunctionII(); break;
        case ID_FKT_SICO: p=new FunctionSICO(); break;
        case ID_FKT_III: p=new FunctionIII(); break;
        case ID_FKT_VII: p=new FunctionVII(); break;
        case ID_FKT_IX: p=new FunctionIX(); break;
        case ID_FKT_X: p=new FunctionX(); break;
        case ID_FKT_LSIN: p=new FunctionLSIN(); break;
        case ID_FKT_ATAN: p=new FunctionATAN(); break;
        case ID_FKT_METADET: p=new FunctionMetaDet(); break;
        case ID_FKT_METAABSC: p=new FunctionMetaABSC(); break;
        default: printf("unknown function\n",aid); p=NULL; break;
    }
    return p;
}

Function* loadFunction(FILE *f) {
	int32_t aid;
	char tmp[1024];
    while (1) {
		fgets(tmp,1000,f); chomp(tmp);
		upper(tmp);
		if (tmp[0]=='#') continue;
		if (strcmp(tmp,"ID") == 0) break; else return NULL;
	}

	fscanf(f,"%i\n",&aid);
    Function *p = getNewFunction(aid);
	if (p) p->load(aid,f);

	return p;
};

int32_t getFirstColorFile(char* ff) {
	// construct a list of color file names
	// in subdirectory COLORCOLLECTIONDIR
	if (ffarbe) fclose(ffarbe);
	char tmp[1024];
	sprintf(tmp,"dir %s*.par /x /b >_ljap.tmp.farbe 2>nul",COLORCOLLECTIONDIR);
	system(tmp);
	ffarbe=fopen("_ljap.tmp.farbe","rt");
	if (!ffarbe) return 0;
	if (fgets(tmp,1024,ffarbe) == NULL) { fclose(ffarbe);	ffarbe=NULL; return 0; }
	chomp(tmp);
	sprintf(ff,"%s%s",COLORCOLLECTIONDIR,tmp);
	
	return 1;
}

int32_t getNextColorFile(char* ff) {
	if (!ffarbe) return 0;
	if (feof(ffarbe)) {	
		fclose(ffarbe);	
		ffarbe=NULL; 
		return 0; 
	}
	char tmp[1024];
	if (fgets(tmp,1024,ffarbe) == NULL) { 
		fclose(ffarbe); 
		ffarbe=NULL; 
		return 0; 
	}
	chomp(tmp);
	sprintf(ff,"%s%s",COLORCOLLECTIONDIR,tmp);
	return 1;
}

char dez(const char c) {
	if ((c>='A')&&(c<='F')) { return c-'A'+10; }
	if ((c>='a')&&(c<='f')) { return c-'a'+10; }
	if ((c>='0')&&(c<='9')) { return c-'0'; }
	return 0;
}

void writehex(FILE *f,const char* s) {
	for(uint32_t i=0;i<strlen(s);i+=2) {
		uint8_t c = 16*dez(s[i]) + dez(s[i+1]);
		fwrite(&c,sizeof(c),1,f);
	}
}

char* stripext(char* s) {
	if (!s) return NULL;

	int32_t l=strlen(s);
	for(int32_t i=(l-1);i>=0;i--) if (s[i]=='.') { s[i]=0; return s; }

	return s;
}

char* upper(char* s) {
	if (!s) return 0;
	for(uint32_t i=0;i<strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
		if (s[i]=='ö') s[i]='Ö';
		if (s[i]=='ü') s[i]='Ü';
		if (s[i]=='ä') s[i]='Ä';
	}

	return s;
}

// fastsin, fastcos von github: fasttrig.as

inline double fastsin(double x) {
	/*
		based on:
		Fast Polynomial Approximations to Sine and Cosine
		Charles K Garrett, 2012
	*/ 
	
	// range -pi..pi
	if (x < -3.14159265) {
		double d=(3.14159265-x) / 6.28318531;
		x += floor(d)*6.28318531;
	} else if (x > 3.14159265) {
		double d=(x+3.14159265) / 6.28318531;
		x -= floor(d)*6.28318531;
	}
	
	double x2 = x * x;
	return
	(((((-2.05342856289746600727e-08*x2+2.70405218307799040084e-06)*x2
	-1.98125763417806681909e-04)*x2+8.33255814755188010464e-03)*x2
	-1.66665772196961623983e-01)*x2+9.99999707044156546685e-01)*x;

}

inline double fastcos(double x) {
	return fastsin(x+PI05);
}

char* removeStr(const char* q,const char* was,char* erg) {
	char* p=strstr(q,was);
	if (p) {
		for(int32_t i=0;i<(p-q);i++) erg[i]=q[i];
		strcpy(&erg[p-q],p+strlen(was));
	} else {
		strcpy(erg,q);
	}

	return erg;
}

inline double maximumD(const double a,const double b) {
	if (a > b) return a;

	return b;
}

inline int32_t minimumI(const int32_t a,const int32_t b) {
	if (a < b) return a;
	return b;
}

char* chomp(char* s) {
	if (!s) return 0;
	for(int32_t i=strlen(s);i>=0;i--) if (s[i]<32) s[i]=0; else break;
	return s;
}


// struct IterDouble

IterDouble::IterDouble(const double a,const double b,const int32_t an) {
	l=a; r=b;
	anz=an;
	if (an==1) delta=b-a; else delta=(b-a)/(an-1);
	wert=l;
	nr=0;
}

int32_t IterDouble::iterStart(void) {
	wert=l;
	nr=1;
	return 1;
}

int32_t IterDouble::iterWeiter(void) {
	wert+=delta;
	nr++;
	if (nr>anz) return 0;
	return 1;
}


// Function LSIN

char* FunctionLSIN::fktStr(char* s) {
	sprintf(s,"N(%i) r*sin(x)*(1-%le*sin(x+r))",id,b);
	return s;
}

char* FunctionLSIN::ablStr(char* s) {
	sprintf(s,"N(%i) [==f'(x)] -r*(%le*sin(2x+r)-cos(x))",id,b);
	return s;
}

FunctionLSIN::FunctionLSIN() {
	id = ID_FKT_LSIN;
	b=2.7;
	typ=FKTTYP_NORMAL;
}

void FunctionLSIN::eval(const double x,const double r,double& fx) {
	fx=r*fastsin(x)*(1-b*fastsin(x+r));
}

void FunctionLSIN::save(FILE *f) {
	char tmp[1024];
	fprintf(f,"ID\n%i\n#FUNCTION LSIN\nB\n%le\n",id,b);
}

void FunctionLSIN::eval(const double x,const double r,double& fx,double& abl) {
	fx=r*fastsin(x)*(1-b*fastsin(x+r));
	abl=-r*(b*fastsin(x+x+r)-fastcos(x));
}

void FunctionLSIN::evalabl(const double x,const double r,double& abl) {
    abl=-r*(b*fastsin(x+x+r)-fastcos(x));
}


// Function ATAN

char* FunctionATAN::fktStr(char* s) {
	sprintf(s,"N(%i) %le*atan((x+r)*sin(x+r))",id,b);
	return s;
}

char* FunctionATAN::ablStr(char* s) {
	sprintf(s,"N(%i) [==f'(x)] %le*( sin(x+r) + (x+r)*cos(x+r) ) / ( 1+(x+r)^2*sin^2(x+r)",id,b);
	return s;
}

FunctionATAN::FunctionATAN() {
    id = ID_FKT_ATAN;
    b=2.7;
    typ=FKTTYP_NORMAL;
}

void FunctionATAN::eval(const double x,const double r,double& fx) {
	const double xr=x+r;
	fx=b*atan(xr*fastsin(xr));
}

void FunctionATAN::save(FILE *f) {
	fprintf(f,"ID\n%i\n#FUNCTION ATAN\nB\n%le\n",id,b);
}

void FunctionATAN::eval(const double x,const double r,double& fx,double& abl) {
	const double xr=x+r;
	const double si=fastsin(xr);
	const double xsi=xr*si;
	fx=b*atan(xsi);
	abl=b*(si+xr*fastcos(xr))/(1+xsi*xsi);
}

void FunctionATAN::evalabl(const double x,const double r,double& abl) {
	const double xr=x+r;
	const double si=fastsin(xr);
	const double xsi=xr*si;
	abl=b*(si+xr*fastcos(xr))/(1+xsi*xsi);
}


// Function II

char* FunctionII::fktStr(char* s) {
	sprintf(s,"N(%i) %le*sin^2(x+r)",id,b);
	return s;
}

char* FunctionII::ablStr(char* s) {
	sprintf(s,"N(%i) [==f'(x)] %lf*sin(x+r)*cos(x+r)",id,b2);
	return s;
}

int32_t FunctionII::iterStart(void) {
	if (!iterb) return 0; // gibt keinen
	if (iterb->iterStart()<=0) return 0;
	set_b(iterb->wert);
	return 1;
}

int32_t FunctionII::iterWeiter(void) {
	if (!iterb) return 0;
	double w;
	if (iterb->iterWeiter()) {
		set_b(iterb->wert);
		return 1;
	}

	return 0;
}

FunctionII::FunctionII() {
    id = ID_FKT_II;
    b=2.7; b2=5.4;
    typ=FKTTYP_NORMAL;
}

void FunctionII::eval(const double x,const double r,double& fx) {
	const double si=fastsin(x+r);
	fx=b*si*si;
}

void FunctionII::eval(const double x,const double r,double& fx,double& abl) {
	const double xr=x+r;
	const double si=fastsin(xr);
	fx=b*si*si;
    abl=b2*si*fastcos(xr);
}

void FunctionII::evalabl(const double x,const double r,double& abl) {
	const double xr=x+r;
    abl=b2*fastsin(xr)*fastcos(xr);
}

void FunctionII::save(FILE *f) {
	fprintf(f,"ID\n%i\n#FUNCTION II\nB\n%le\n",id,b);
}

int32_t FunctionII::load(const int32_t aid,FILE *f) {
	if (aid!=id) { return 0; }

	int32_t pnotw=1,param=0;
    char puffer[1000];
    double w; 
	int32_t i=0;
    while (i<pnotw) {
		fgets(puffer,1000,f); chomp(puffer);
		upper(puffer);
		if (puffer[0]=='#') continue; // Bemerkung
        else i++;
        if (strcmp(puffer,"B")==0) { fscanf(f,"%le\n",&w); param++; set_b(w); }
	}
    if (param!=pnotw) return -1;

	return 1;
}


// Function III

char* FunctionIII::fktStr(char* s) {
	sprintf(s,"N(%i) %le*sin(x+r)*sin(x-r)",id,b);
	return s;
}

char* FunctionIII::ablStr(char* s) {
	sprintf(s,"N(%i) [==f'(x)] %le*sin(2*x)",id,b);
	return s;
}

FunctionIII::FunctionIII() {
	id = ID_FKT_III;
	b=2.7; 
	typ=FKTTYP_NORMAL;
}

void FunctionIII::eval(const double x,const double r,double& fx) {
	fx=b*fastsin(x+r)*fastsin(x-r);
}

void FunctionIII::eval(const double x,const double r,double& fx,double& abl) {
	fx=b*fastsin(x+r)*fastsin(x-r);
	abl=b*fastsin(x+x);
}

void FunctionIII::evalabl(const double x,const double r,double& abl) {
	abl=b*fastsin(x+x);
}

void FunctionIII::save(FILE *f) {
	fprintf(f,"ID\n%i\n#FUNCTION III\nB\n%le\n",id,b);
}


// Function VII

char* FunctionVII::fktStr(char* s) {
	sprintf(s,"D(%i) %le*sin^2(x+r)",id,b);
	return s;
}

FunctionVII::FunctionVII() {
	id = ID_FKT_VII;
	b=2.7; b2=b+b;
	typ=FKTTYP_DETACHED;
}

void FunctionVII::eval(const double x,const double r,double& fx) {
	const double si=fastsin(x+r);
	fx=b*si*si;
}

void FunctionVII::eval(const double x,const double r,double& fx,double& abl) {
	const double si=fastsin(x+r);
	fx=b*si*si;
	const double rx=r*x;
	abl=r-rx-rx;
}

void FunctionVII::evalabl(const double x,const double r,double& abl) {
	const double si=fastsin(x+r);
	const double rx=r*x;
	abl=r-rx-rx;
}

char* FunctionVII::ablStr(char* s) {
	sprintf(s,"DET(%i) r-2rx",id);
	return s;
}

void FunctionVII::save(FILE *f) {
	fprintf(f,"ID\n%i\n#DETACHED FUNCTION VII\nB\n%le\n",id,b);
}


// Function IX

char* FunctionIX::fktStr(char* s) {
	sprintf(s,"DET(%i) %le*sin(x+r)%+le*sin^2(%le*x+r)",id,b,b,b,b);
	return s;
}

FunctionIX::FunctionIX() {
	id = ID_FKT_IX;
	b=2.7; b2=b+b;
	typ=FKTTYP_DETACHED;
}

void FunctionIX::eval(const double x,const double r,double& fx) {
	const double si=fastsin(b*x+r);
	fx=b*fastsin(x+r)+b*si*si;
}

void FunctionIX::eval(const double x,const double r,double& fx,double& abl) {
	const double si=fastsin(b*x+r);
	fx=b*fastsin(x+r)+b*si*si;
	const double si2=fastsin(x+r*b);
	abl=si2*si2-r*x;
}

void FunctionIX::evalabl(const double x,const double r,double& abl) {
	const double si=fastsin(b*x+r);
	const double si2=fastsin(x+r*b);
	abl=si2*si2-r*x;
}

char* FunctionIX::ablStr(char* s) {
	sprintf(s,"DET(%i) sin^2(x%+le*r)-r*x",id,b);
	return s;
}

void FunctionIX::save(FILE *f) {
	fprintf(f,"ID\n%i\n#DETACHED FUNCTION IX\nB\n%le\n",id,b);
}


// Function X

char* FunctionX::fktStr(char* s) {
	sprintf(s,"DET(%i) r*sin^2(x-r)%+le*sin^3(x+2*r)",id,b,b,b);
	return s;
}

FunctionX::FunctionX() {
	id = ID_FKT_X;
	b=2.7; b2=b+b;
	typ=FKTTYP_DETACHED;
}

void FunctionX::eval(const double x,const double r,double& fx) {
	const double si=fastsin(x-r);
	const double si2=fastsin(x+r+r);
	fx=r*si*si+b*si2*si2*si2;
}

void FunctionX::eval(const double x,const double r,double& fx,double& abl) {
	const double si=fastsin(x-r);
	const double si2=fastsin(x+r+r);
	fx=r*si*si+b*si2*si2*si2;
	const double rx=r*x;
	const double si3=fastsin(rx-b);
	const double si4=si3*si3;
	abl=rx-b*si4*si4;
}

void FunctionX::evalabl(const double x,const double r,double& abl) {
	const double rx=r*x;
	const double si3=fastsin(rx-b);
	const double si4=si3*si3;
	abl=rx-b*si4*si4;
}

char* FunctionX::ablStr(char* s) {
	sprintf(s,"DET(%i) rx-%le*sin^4(rx-%le)",id,b,b);
	return s;
}

void FunctionX::save(FILE *f) {
	fprintf(f,"ID\n%i\n#DETACHED FUNCTION X\nB\n%le\n",id,b);
}


// Function MetaDet

char* FunctionMetaDet::fktStr(char* s) {
	char t1[1024],t2[1024];
	if (fwas==WAS_F) f->fktStr(t1); else f->ablStr(t1);

	removeStr(t1,"[==f'(x)]",t2);
	sprintf(s,"METADET(%i) f(x)=%s",id,t2);

	return s;
}

char* FunctionMetaDet::ablStr(char* s) {
	char t1[1024],t2[1024];
	if (ablwas==WAS_F) abl->fktStr(t2); else abl->ablStr(t2);

	sprintf(s,"METADET(%i) g(x)=%s",id,t2);

	return s;
}

void FunctionMetaDet::set_b(const double d) {
	if (f) f->set_b(d);
	if (abl) abl->set_b(d);
}

int32_t FunctionMetaDet::load(const int32_t aid,FILE* ff) {
	if (aid!=id) { return 0; }

    int32_t pnotw=2,param=0;
    char puffer[1000];
    double w;
    int32_t i=0;
    while (i<pnotw) {
		fgets(puffer,1000,ff); chomp(puffer);
        if (puffer[0]=='#') continue; // Bemerkung
        else i++;
        upper(puffer);
		if (strcmp(puffer,"FWAS")==0) { fscanf(ff,"%i\n",&fwas); param++; }
        else if (strcmp(puffer,"ABLWAS")==0) { fscanf(ff,"%i\n",&ablwas); param++; }
	}
	
    if (param!=pnotw) return 0;
	if ((f=loadFunction(ff)) == NULL) return 0;
	if ((abl=loadFunction(ff)) == NULL) return 0;

	return 1;
}

int32_t FunctionMetaDet::iterStart(void) {
	if (!iterb) return 0; // gibt keinen
	if (iterb->iterStart()<=0) return 0;
	if (f) f->set_b(iterb->wert);
	if (abl) abl->set_b(iterb->wert);
	return 1;
}

int32_t FunctionMetaDet::iterWeiter(void) {
	if (!iterb) return 0;
	double w;
	if (iterb->iterWeiter()) {
		if (f) f->set_b(iterb->wert);
		if (abl) abl->set_b(iterb->wert);
		return 1;
	}

	return 0;
}

FunctionMetaDet::FunctionMetaDet() {
	id = ID_FKT_METADET;
    b=2.7; b2=b+b;
    f=abl=NULL;
    fwas=ablwas=WAS_F;
    typ=FKTTYP_METADET;
}

void FunctionMetaDet::eval(const double x,const double r,double& fx) {
	double tmp;

	if (fwas==WAS_F) f->eval(x,r,fx); else f->eval(x,r,tmp,fx);
}

void FunctionMetaDet::eval(const double  x,const double  r,double & fx,double & ab) {
	double tmp;
	if (fwas==WAS_F) f->eval(x,r,fx); else f->evalabl(x,r,fx);
	if (ablwas==WAS_F) abl->eval(x,r,ab); else abl->evalabl(x,r,ab);
}

void FunctionMetaDet::save(FILE *ff) {
	fprintf(ff,"ID\n%i\n#METADET\nFWAS\n%i\nABLWAS\n%i\n",id,fwas,ablwas);
	fprintf(ff,"#FKT\n");
	f->save(ff);
	fprintf(ff,"#ABL\n");
	abl->save(ff);
}


// Function MetaAbsc

char* FunctionMetaABSC::fktStr(char* s) {
	char tmp2[1024],t3[1024],tmp3[1024],tmp4[1024];
	sprintf(s,"METAABSC(%i)\n",id);
	sprintf(&s[strlen(s)],"Initial iterations if %.5lf <= x <= %.5lf: f(x)=%s else f(x)=%s\n",
		I0MIN,I0MAX,fint->fktStr(tmp2),fext->fktStr(t3));
	return s;
}

char* FunctionMetaABSC::ablStr(char* s) {
	char tmp2[1024],t3[1024],tmp3[1024],tmp4[1024];
	sprintf(s,"METAABSC(%i)\n",id);
	sprintf(&s[strlen(s)],"Computing iterations: if %.5lf <= x <= %.5lf: g(x)=%s else g(x)=%s",
		I1MIN,I1MAX,
		fint->ablStr(tmp3),
		fext->ablStr(tmp4)
	);

	return s;
}

void FunctionMetaABSC::set_b(const double d) {
	fint->set_b(d);
	fext->set_b(d);
}

int32_t FunctionMetaABSC::load(const int32_t aid,FILE* ff) {
	if (aid!=id) { return 0; }

	int32_t pnotw=6,param=0;
    char puffer[1000];
    double w;
    int32_t i=0;
    while (i<pnotw) {
		fgets(puffer,1000,ff); chomp(puffer);
		upper(puffer);
        if (puffer[0]=='#') continue; // Bemerkung
        else i++;

        if (strcmp(puffer,"I0MAX")==0) { fscanf(ff,"%le\n",&w); I0MAX=w; param++; }
        else if (strcmp(puffer,"I0MIN")==0) { fscanf(ff,"%le\n",&w); I0MIN=w; param++; }
        else if (strcmp(puffer,"I1MAX")==0) { fscanf(ff,"%le\n",&w); I1MAX=w; param++; }
        else if (strcmp(puffer,"I1MIN")==0) { fscanf(ff,"%le\n",&w); I1MIN=w; param++; }
        if (strcmp(puffer,"FINT")==0) {
			if ((fint=loadFunction(ff)) == NULL) return -1;
			param++;
		} else if (strcmp(puffer,"FEXT")==0) {
			if ((fext=loadFunction(ff)) == NULL) return -1;
			param++;
		}
	}
    
    if (param!=pnotw) { printf("Parameters missing\n"); return 0; }

	return 1;
}

int32_t FunctionMetaABSC::iterStart(void) {
	if (!iterb) return 0;
	if (iterb->iterStart()<=0) return 0;
	if (fint) fint->set_b(iterb->wert);
	if (fext) fext->set_b(iterb->wert);
	
	return 1;
}

int32_t FunctionMetaABSC::iterWeiter(void) {
	if (!iterb) return 0;

	double w;
	if (iterb->iterWeiter()) {
		if (fint) fint->set_b(iterb->wert);
		if (fext) fext->set_b(iterb->wert);
		return 1;
	}

	return 0;
}

FunctionMetaABSC::FunctionMetaABSC() {
	id = ID_FKT_METAABSC;
    b=2.7; b2=b+b;
    fint=fext=NULL;
    typ=FKTTYP_ABSCHNITTSWEISE;
}

void FunctionMetaABSC::eval(const double x,const double r,double& fx) {
	double tmp;

	if ((x > I0MAX) || (x < I0MIN)) fext->eval(x,r,fx);
	else fint->eval(x,r,fx);
}

void FunctionMetaABSC::eval(const double x,const double r,double& fx,double& abl) {
	double tmp;

	if ((x > I1MAX) || (x < I1MIN)) fext->eval(x,r,fx,abl);
	else fint->eval(x,r,fx,abl);
}

void FunctionMetaABSC::save(FILE *ff) {
	fprintf(ff,"ID\n%i\n#METAABSC\n",id);
	fprintf(ff,"I0MIN\n%le\n",I0MIN);
	fprintf(ff,"I0MAX\n%le\n",I0MAX);
	fprintf(ff,"I1MIN\n%le\n",I1MIN);
	fprintf(ff,"I1MAX\n%le\n",I1MAX);
	fprintf(ff,"FINT\n");
	fint->save(ff);
	fprintf(ff,"FEXT\n");
	fext->save(ff);
}


// Function SICO

char* FunctionSICO::fktStr(char* s) {
	sprintf(s,"N(%i) %le*sin(x+r*cos(x+r)",id,b);
	return s;
}

int32_t FunctionSICO::iterStart(void) {
	if (!iterb) return 0;
	if (iterb->iterStart()<=0) return 0;
	set_b(iterb->wert);
	return 1;
}

int32_t FunctionSICO::iterWeiter(void) {
	if (!iterb) return 0;
	double w;
	if (iterb->iterWeiter()) {
		set_b(iterb->wert);
		return 1;
	}

	return 0;
}

FunctionSICO::FunctionSICO() {
    id = ID_FKT_SICO;
    b=2.7;
    typ=FKTTYP_NORMAL;
}

void FunctionSICO::evalabl(const double x,const double r,double& abl) {
	const double xr=x+r;
	const double xrc=x+r*fastcos(xr);
	abl=b*(1-r*fastsin(xr))*fastcos(xrc);
}

void FunctionSICO::eval(const double x,const double r,double& fx) {
	fx=b*fastsin(x+r*fastcos(x+r));
}

void FunctionSICO::eval(const double x,const double r,double& fx,double& abl) {
	const double xr=x+r;
	const double xrc=x+r*fastcos(xr);
	fx=b*fastsin(xrc);
	abl=b*(1-r*fastsin(xr))*fastcos(xrc);
}

char* FunctionSICO::ablStr(char* s) {
	sprintf(s,"N(%i) [==f'(x)] %le*(1-r*sin(x+r))*cos(x+r*cos(x+r))",id,b);
	return s;
}

void FunctionSICO::save(FILE *f) {
	fprintf(f,"ID\n%i\n#FUNCTION SICO\nB\n%le\n",id,b);
}

int32_t FunctionSICO::load(const int32_t aid,FILE *f) {
	if (aid!=id) { return 0; }

    int32_t pnotw=1,param=0;
    char puffer[1000];
    double w;
	int32_t i=0;
    while (i<pnotw) {
		fgets(puffer,1000,f); chomp(puffer);
        if (puffer[0]=='#') continue; // Bemerkung
        else i++;
        upper(puffer);
        if (strcmp(puffer,"B")==0) { fscanf(f,"%le\n",&w); param++; set_b(w); }
	}
    if (param!=pnotw) return 0;

	return 1;
}


// Function I

FunctionI::FunctionI() {
	id = ID_FKT_I;
    typ=FKTTYP_NORMAL;
}

void FunctionI::eval(const double x,const double r,double &fx) { 
	fx=r*x*(1-x); 
}

void FunctionI::eval(const double x,const double r,double &fx,double &abl) {
	const double rx=r*x;
	fx=rx*(1-x);
    abl=r-rx-rx;
}

void FunctionI::evalabl(const double x,const double r,double &abl) {
	const double rx=r*x;
	abl=r-rx-rx;
}

void FunctionI::save(FILE *f) {
	fprintf(f,"ID\n%i\n#FUNCTION I\n",id);
}

int32_t FunctionI::load(const int32_t aid,FILE *) {
	if (aid!=id) return 0;
	return 1;
}


// struct ColIntv

int32_t ColIntv::load(FILE *f) {
	int32_t pnotw=4,param=0;
    char puffer[1000];
    double w;
	int32_t r,g,b;
    int32_t i=0;
	while (i<pnotw) {
		fgets(puffer,1000,f); chomp(puffer);
        if (puffer[0]=='#') continue; 
		else i++;
		upper(puffer);
		if (strcmp(puffer,"GRENZEL")==0) { fscanf(f,"%le\n",&w); param++; setgrenzel(w); }
        else if (strcmp(puffer,"GRENZER")==0) { fscanf(f,"%le\n",&w); param++; setgrenzer(w); }
        else if (strcmp(puffer,"FARBEL")==0) { fscanf(f,"%i\n%i\n%i\n",&r,&g,&b); param++; setfarbel(r,g,b); }
        else if (strcmp(puffer,"FARBER")==0) { fscanf(f,"%i\n%i\n%i\n",&r,&g,&b); param++; setfarber(r,g,b); }
	}
	
    if (param!=pnotw) return 0;

	return 1;
};

void ColIntv::save(FILE *f) {
    fprintf(f,"GRENZEL\n%le\n",gl);
    fprintf(f,"GRENZER\n%le\n",gr);
    fprintf(f,"FARBEL\n%i\n%i\n%i\n",lr,lg,lb);
    fprintf(f,"FARBER\n%i\n%i\n%i\n",rr,rg,rb);
}

void ColIntv::precalc(void) {
	dr=rr-lr; dg= rg-lg; db=rb-lb;
}

int32_t ColIntv::farbe(const double w,int32_t&r,int32_t&g,int32_t&b) {
	if ((w>=gl)&&(w<gr)) { 
		const double wg=(w-gl)/breite;
	    r=lr+int32_t(wg*dr);
	    g=lg+int32_t(wg*dg);
		b=lb+int32_t(wg*db);
	    return 1;
	}

    return 0;
}

void ColIntv::setfarbel(const int32_t ar,const int32_t ag,const int32_t ab) {
	lr=ar; lg=ag; lb=ab; dr=rr-lr; 
	precalc();
}

void ColIntv::setfarber(const int32_t ar,const int32_t ag,const int32_t ab) {
	rr=ar; rg=ag; rb=ab; 
	precalc();
}


// struct IntervalColoring

void IntervalColoring::save(FILE *f) {
	fprintf(f,"ID\n%i\n",id);
    fprintf(f,"FARBEL\n%i\n%i\n%i\n",lr,lg,lb);
    fprintf(f,"FARBER\n%i\n%i\n%i\n",rr,rg,rb);
    fprintf(f,"INTANZ\n%i\n",intanz);
    for(int32_t i=0;i<intanz;i++) ints[i]->save(f);
}

int32_t IntervalColoring::load(const int32_t aid,FILE *f) {
	if (aid!=id) return 0;
        
	int32_t r,g,b,i;
    char puffer[1000];
    int32_t lesanz=0,param=0;
    i=0;
    while (i<3) {
		fgets(puffer,1000,f); chomp(puffer);
        if (puffer[0]=='#') continue;
        else i++;
        upper(puffer);
        
        if (strcmp(puffer,"FARBEL")==0) { 
			fscanf(f,"%i\n%i\n%i\n",&r,&g,&b); 
			setfarbel(r,g,b); 
			param++; 
		} else if (strcmp(puffer,"FARBER")==0) { 
			fscanf(f,"%i\n%i\n%i\n",&r,&g,&b); 
			setfarber(r,g,b); 
			param++; 
		} else if (strcmp(puffer,"INTANZ")==0) { 
			fscanf(f,"%i\n",&lesanz); 
			param++; 
		}
	}
	
	if (lesanz > MAXINTANZ) return 0;

    clear();
    for(i=0;i<lesanz;i++) {
		ColIntv* ints=new ColIntv; 
		ints->load(f);
        Addintervall(ints);
	}

    return 1;
};

int32_t IntervalColoring::Addintervall(ColIntv* p) {
	if (intanz<MAXINTANZ) {
		ints[intanz]=p;
		intanz++;
        if ((p->gl<mingl)||(intanz==1)) mingl=p->gl;
        if ((p->gr>maxgl)||(intanz==1)) maxgl=p->gr;
        return 1;
	}
	
	return 0;
}

int32_t IntervalColoring::farbe(const double w,int32_t&r,int32_t&g,int32_t&b) {
	if (w<mingl) { 
		r=lr; g=lg; b=lb; 
		return 1;
	}
	if (w>maxgl) { 
		r=rr; g=rg; b=rb; 
		return 1;
	}
    for(int32_t i=0;i<intanz;i++) {
		if (ints[i]->farbe(w,r,g,b) > 0) return 1;
	}
	
    return 0;
}

void IntervalColoring::clear(void) {
	if (ints) {
		for(int32_t i=0;i<intanz;i++) delete ints[i];
		intanz=0;
	}
}

IntervalColoring::~IntervalColoring() {
	clear();
}


// struct ljapunow

int32_t Ljapunow::loadpar(char *fn) {
	FILE *f=fopen(fn,"rt");
	if (!f) return 0;

    char puffer[1000],puffer2[1000];
    int32_t pnotw=11,param=0,tlenx=0,tleny=0;
    int32_t tmpiter0=100,tmpiter1=200;
    double w1,w2; 
	while (!feof(f)) {
		fgets(puffer,1000,f); chomp(puffer);
        if ((puffer[0]=='#')||(puffer[0]=='.')) continue;
        
		upper(puffer);
		if (strcmp(puffer,"FUNKTION")==0) {
			Function *p = loadFunction(f);
            if (p) { 
				setFunction(p); 
				param++; 
			} else return -1;
		} else 
		if (strcmp(puffer,"FAERBUNG")==0) {
			IntervalColoring *p=loadfaerbung(f);
			if (p) { 
				setfarbe(p); 
				param++; 
			} else return 0;
		} else if (strcmp(puffer,"LENX")==0) { 
			fscanf(f,"%i\n",&tlenx); 
			param++; 
		} else if (strcmp(puffer,"LENY")==0) { 
			fscanf(f,"%i\n",&tleny); 
			param++; 
		} else if (strcmp(puffer,"ITER0")==0) { 
			fscanf(f,"%i\n",&tmpiter0); 
			param++; 
		} else if (strcmp(puffer,"ITER1")==0) { 
			fscanf(f,"%i\n",&tmpiter1); 
			param++; 
		} else if (strcmp(puffer,"X0")==0) { 
			fscanf(f,"%le\n",&w1); 
			x0=w1; 
			param++; 
		} else if (strcmp(puffer,"OL")==0) { 
			fscanf(f,"%le\n%le\n",&w1,&w2); 
			upperleft.x=w1; upperleft.y=w2; 
			param++; 
		} else if (strcmp(puffer,"UL")==0) { 
			fscanf(f,"%le\n%le\n",&w1,&w2); 
			lowerleft.x=w1; 
			lowerleft.y=w2; 
			param++; 
		}
		else if (strcmp(puffer,"UR")==0) { 
			fscanf(f,"%le\n%le\n",&w1,&w2); 
			lowerright.x=w1; lowerright.y=w2; 
			param++; 
		} else if (strcmp(puffer,"SEQUENZ")==0) {
			param++;
            fgets(puffer2,1000,f); chomp(puffer2);
            setSequence(puffer2);
		} else {
			printf("Unknown parameter %s\n",puffer);
			return 0;
		}
	}
    if (param!=pnotw) { 
		printf("Not enough parameters,\n");
		return 0;
	}
    setlen(tlenx,tleny);
    fclose(f);

    setiter(tmpiter0,tmpiter1); 

	return 1;
};

int32_t Ljapunow::loadcolor(char *fn) {
	FILE *f=fopen(fn,"rt");
    if (!f) return 0;
    char puffer[1000],puffer2[1000];
    int32_t pnotw=1,param=0;
    while (!feof(f)) {
		fgets(puffer,1000,f); chomp(puffer);
        if ((puffer[0]=='#')||(puffer[0]=='.')) continue; // Weiter
        upper(puffer);
        if (strcmp(puffer,"FAERBUNG")==0) {
			IntervalColoring *p = loadfaerbung(f);
            if (p) { 
				setfarbe(p); 
				param++; 
				break;
			} else return 0;
        }
	}
	if (param!=pnotw) { printf("Ljapunow::loadfarbe. Zuwenig Parameter (%i/%i)\n",param,pnotw); return -1; }
    fclose(f);

	return 1;
};

Ljapunow::Ljapunow() {
    fkt=NULL;
    farbe=NULL;
    lenx=leny=600;
    iter0=50; iter1=100;
    iter1d=100; iter1h=50; iter0h=50;
    INViter1d=1.0; INViter1d /= iter1d;
    x0=0.5;
    seqlen=0; 
	exps=0;
};

Ljapunow::~Ljapunow() {
	if (exps) delete[] exps;
	if (fkt) delete fkt;
    if (farbe) delete farbe;
};

int32_t Ljapunow::calc(const int32_t astart,const int32_t aende) {
	double AB[16]; 
    double px,abl1,abl2,tmp;
	double lambda;
    double xmin=9999999,xmax=-999999;
    int32_t start=astart;
    if (start<0) start=0;
    if (start>=leny) start=leny-1;
    int32_t ende=aende;
    if (ende<0) ende=0;
    if (ende>=leny) ende=leny-1;

	uint32_t offset;
    Point32_t vx,vy;
    vx.x=(lowerright.x-lowerleft.x)/lenx; vx.y=(lowerright.y-lowerleft.y)/lenx;
    vy.x=(upperleft.x-lowerleft.x)/leny; vy.y=(upperleft.y-lowerleft.y)/leny;

    time_t t0=time(NULL);
    const int32_t NOCH0=128;
	int32_t noch=NOCH0;

	offset=start*lenx;
	for(uint32_t y=start;y<=ende;y++) {
		if (--noch==0) {
			noch=NOCH0;
			time_t b=time(NULL);
			double d=difftime(b,t0);
			d /= (y-start);
			d *= (ende-y);
			printf("row %i --- %.0lf sec to go ---\n",y,d);
		}

		AB[0]=lowerleft.x+y*vy.x;
		AB[1]=lowerleft.y+y*vy.y;

		for(int32_t x=0;x<lenx;x++) {
			px=x0;
            int32_t seqpos=0;

			// initial iterations to settle a bit
			for(int32_t i=0;i<iter0h;i++) {
				fkt->eval(px,AB[sequence[seqpos]],tmp); 
				SEQPOSINC(seqpos);
                fkt->eval(tmp,AB[sequence[seqpos]],px); 
                SEQPOSINC(seqpos);
			} // i
            
            lambda = 0.0;
            // lyapunov value computing iterations
			for(uint32_t i=0;i<iter1h;i++) {
				fkt->eval(px,AB[sequence[seqpos]],tmp,abl1); 
				SEQPOSINC(seqpos);
                fkt->eval(tmp,AB[sequence[seqpos]],px,abl2); 
                SEQPOSINC(seqpos);
				const double ab=fabs(abl1*abl2);
				if (ab > 1E-300) lambda += log(ab);
			} // i

            exps[offset]=lambda * INViter1d;
            offset++;
			AB[0]+=vx.x;
			AB[1]+=vx.y;
		} // x
	} // y

	return 1;
};

void Ljapunow::setfarbe(IntervalColoring* f) {
	if (farbe) delete farbe;
	farbe=f;
}

void Ljapunow::setiter(const int32_t i0,const int32_t i1) {
	iter0=i0; iter0h=(i0 >> 1);
	iter1=i1; iter1h=(i1 >> 1); 
	iter1d=(iter1h << 1);
	INViter1d=1.0; INViter1d /= iter1d;
}

void Ljapunow::savedescr(const char* fn) {
	FILE *fff=fopen(fn,"wt");
	char tmp2[1024],tmp3[1024];

	fprintf(fff,"x0=%lf, %i initial and %i computing iterations\n",x0,iter0,iter1);
	fprintf(fff,"Trajectory function f(x)=%s\n",fkt->fktStr(tmp2));
	fprintf(fff,"Computing function g(x)=%s\n",fkt->ablStr(tmp2));
	fprintf(fff,"Sequence %s. Center (%.2lf/%.2lf) size=%.10lf\n",
		getSequence(tmp3),
		(lowerleft.x+lowerright.x)*0.5,
		0.5*(lowerleft.y+upperleft.y),
		sqrt( maximumD(
			( (lowerleft.x-lowerright.x)*(lowerleft.x-lowerright.x) + (lowerleft.y-lowerright.y)*(lowerleft.y-lowerright.y) ),
			( (lowerleft.x-upperleft.x)*(lowerleft.x-upperleft.x) + (lowerleft.y-upperleft.y)*(lowerleft.y-upperleft.y) )
		) )
	);

	if (farbe->id == ID_FAERBUNG_INTERVALL) {
		IntervalColoring* previ=(IntervalColoring*)farbe;
		fprintf(fff,"Coloring with linear RGB value interpolation in several intervals\nLyapunov exponents less than %lf: RGB(%i,%i,%i)\n",previ->ints[0]->gl,previ->lr,previ->lg,previ->lb);
		for(int32_t i=0;i<previ->intanz;i++) {
			fprintf(fff,"in [%.2lf..%.2lf] (%i,%i,%i)..(%i,%i,%i)\n",previ->ints[i]->gl,previ->ints[i]->gr,
			previ->ints[i]->lr,previ->ints[i]->lg,previ->ints[i]->lb,previ->ints[i]->rr,previ->ints[i]->rg,previ->ints[i]->rb);
		}
		fprintf(fff,"greater than %.2lf: (%i,%i,%i)\n",previ->ints[previ->intanz-1]->gr,previ->rr,previ->rg,previ->rb);
	}
	fclose(fff);
}

int32_t Ljapunow::iterStart(void) {
	if ((!fkt) || (!farbe) ) return 0;
	if (fkt->iterStart()<=0) return 0;
	return 1;
}

int32_t Ljapunow::iterWeiter(void) {
	if ((!fkt) || (!farbe) ) return 0;
	if (fkt->iterWeiter()<=0) return 0;
	return 1;
}

void Ljapunow::setlen(const int32_t xl,const int32_t yl) {
	if (exps) delete[] exps;
	lenx=((xl >> 2) << 2);
	leny=((yl >> 2) << 2);
    exps=new double[lenx*leny];
}

void Ljapunow::saveexp(char *fn) {
	FILE *f=fopen(fn,"wb");
    fwrite(&lenx,sizeof(lenx),1,f);
    fwrite(&leny,sizeof(leny),1,f);
    fwrite(exps,sizeof(double),lenx*leny,f);
    fclose(f);
}

int32_t Ljapunow::loadexp(char *fn) {
    FILE *f=fopen(fn,"rb");
	if (!f) return 0;

    int32_t wx,wy;
    fread(&wx,sizeof(lenx),1,f);
	fread(&wy,sizeof(leny),1,f);
	if ((wx!=lenx)||(wy!=leny)) {
		fclose(f);
		return 0;
	}
	double *ex=new double[lenx];
	int32_t off=0;
	for(int32_t y=0;y<leny;y++) {
		fread(ex,sizeof(double),lenx,f);
		for(int32_t i=0;i<lenx;i++) exps[off++]=ex[i];
	}

	delete[] ex;
	fclose(f);
	
	return 1;
}

void Ljapunow::createBmp(Bitmap* bmp) {
	bmp->setlenxy(lenx,leny);
    int32_t offsetbmp=0,offsetexp=0;
    int32_t r,g,b;
    for(int32_t y=0;y<leny;y++) {
		for(int32_t x=0;x<lenx;x++) {
			farbe->farbe(exps[offsetexp],r,g,b);
            bmp->bmp[offsetbmp]=b;
            bmp->bmp[offsetbmp+1]=g;
            bmp->bmp[offsetbmp+2]=r;
            offsetbmp+=3; offsetexp++;
		} // x
    } // y
}

void Ljapunow::savebmp(char *fn,Bitmap* bmp) {
	int32_t neu=0;
	if (!bmp) {
		bmp=new Bitmap;
		neu=1;
	}

	createBmp(bmp);
	bmp->save(fn);
	if (neu) delete bmp;
};

char* Ljapunow::getSequence(char* s) {
	s[0]=0;
	for(int32_t i=0;i<seqlen;i++) s[i]='A'+sequence[i];
	s[seqlen]=0;
	return s;
}

void Ljapunow::savepar(char *fn) {
	FILE *f=fopen(fn,"wt");
    fprintf(f,"FUNKTION\n");
    if (fkt) fkt->save(f);
    fprintf(f,"FAERBUNG\n");
    if (farbe) farbe->save(f);
    fprintf(f,"LENX\n%i\n",lenx);
    fprintf(f,"LENY\n%i\n",leny);
    fprintf(f,"ITER0\n%i\n",iter0);
    fprintf(f,"ITER1\n%i\n",iter1);
    fprintf(f,"X0\n%le\n",x0);
    fprintf(f,"SEQUENZ\n");
    for(int32_t i=0;i<seqlen;i++) fprintf(f,"%c",'A'+sequence[i]);
    fprintf(f,"\n");
    fprintf(f,"OL\n%le\n%le\n",upperleft.x,upperleft.y);
    fprintf(f,"UL\n%le\n%le\n",lowerleft.x,lowerleft.y);
    fprintf(f,"UR\n%le\n%le\n",lowerright.x,lowerright.y);

    fclose(f);
}

void Ljapunow::centerPixel(const int32_t px,const int32_t py) {
	Point32_t M,PX;
	M.x=lowerleft.x+0.5*(lowerright.x-lowerleft.x)+0.5*(upperleft.x-lowerleft.x);
	M.y=lowerleft.y+0.5*(lowerright.y-lowerleft.y)+0.5*(upperleft.y-lowerleft.y);

	Point32_t vx,vy;
    vx.x=(lowerright.x-lowerleft.x)/lenx; vx.y=(lowerright.y-lowerleft.y)/lenx;
    vy.x=(upperleft.x-lowerleft.x)/leny; vy.y=(upperleft.y-lowerleft.y)/leny;

	PX.x = lowerleft.x + px*vx.x + (leny-py)*vy.x;
	PX.y = lowerleft.y + px*vx.y + (leny-py)*vy.y;

	lowerleft.x=lowerleft.x+PX.x-M.x;
	lowerleft.y=lowerleft.y+PX.y-M.y;
	lowerright.x=lowerright.x+PX.x-M.x;
	lowerright.y=lowerright.y+PX.y-M.y;
	upperleft.x=upperleft.x+PX.x-M.x;
	upperleft.y=upperleft.y+PX.y-M.y;
}

void Ljapunow::tile(char* fnprefix,const int32_t anzx,const int32_t anzy) {
	Point32_t siclowerleft,siclowerright,sicupperleft;
	
	siclowerleft.x=lowerleft.x;
	siclowerleft.y=lowerleft.y;
	siclowerright.x=lowerright.x;
	siclowerright.y=lowerright.y;
	sicupperleft.x=upperleft.x;
	sicupperleft.y=upperleft.y;

    Point32_t vx,vy;
    vx.x=(lowerright.x-lowerleft.x)/anzx; vx.y=(lowerright.y-lowerleft.y)/anzx;
    vy.x=(upperleft.x-lowerleft.x)/anzy; vy.y=(upperleft.y-lowerleft.y)/anzy;

	int32_t ctr=1, anz=anzx*anzy;
	char tmp[1000];
	for(int32_t x=0;x<anzx;x++) {
		for(int32_t y=0;y<anzy;y++) {
			printf("\ntile %i/%i ... ",ctr,anz);
			lowerleft.x=siclowerleft.x+x*vx.x+y*vy.x;
			lowerleft.y=siclowerleft.y+x*vx.y+y*vy.y;
			lowerright.x=siclowerleft.x+(x+1)*vx.x+y*vy.x;
			lowerright.y=siclowerleft.y+(x+1)*vx.y+y*vy.y;
			upperleft.x=siclowerleft.x+x*vx.x+(y+1)*vy.x;
			upperleft.y=siclowerleft.y+x*vx.y+(y+1)*vy.y;

			calc(0,leny-1);
			sprintf(tmp,"_walktile_%s_%06i.par",fnprefix,ctr); 
			savepar(tmp);
			sprintf(tmp,"_walktile_%s_%06i.bmp",fnprefix,ctr); 
			savebmp(tmp,NULL);
			ctr++;
		} // y
	} // x
}

void Ljapunow::crop(const int32_t pulneux,const int32_t pulneuy,const int32_t porneux,const int32_t porneuy) {
	Point32_t lowerleftn,lowerrightn,upperleftn;

	Point32_t vx,vy;
	vx.x=(lowerright.x-lowerleft.x)/lenx; vx.y=(lowerright.y-lowerleft.y)/lenx;
    vy.x=(upperleft.x-lowerleft.x)/leny; vy.y=(upperleft.y-lowerleft.y)/leny;

	lowerleftn.x = lowerleft.x + pulneux*vx.x + (leny-pulneuy)*vy.x;
	lowerleftn.y = lowerleft.y + pulneux*vx.y + (leny-pulneuy)*vy.y;

	lowerrightn.x = lowerleft.x + porneux*vx.x + (leny-pulneuy)*vy.x;
	lowerrightn.y = lowerleft.y + porneux*vx.y + (leny-pulneuy)*vy.y;

	upperleftn.x = lowerleft.x + pulneux*vx.x + (leny-porneuy)*vy.x;
	upperleftn.y = lowerleft.y + pulneux*vx.y + (leny-porneuy)*vy.y;

	lowerleft.x=lowerleftn.x;
	lowerleft.y=lowerleftn.y;
	lowerright.x=lowerrightn.x;
	lowerright.y=lowerrightn.y;
	upperleft.x=upperleftn.x;
	upperleft.y=upperleftn.y;
}

void Ljapunow::stretch(const double fx,const double fy) {
	Point32_t lowerleftn,upperleftn,lowerrightn,M;

	M.x = 0.5*(lowerright.x+upperleft.x);
	M.y = 0.5*(lowerright.y+upperleft.y);

	lowerleftn.x = M.x + (lowerleft.x-M.x)*fx;
	lowerleftn.y = M.y + (lowerleft.y-M.y)*fy;

	lowerrightn.x = M.x + (lowerright.x-M.x)*fx;
	lowerrightn.y = M.y + (lowerright.y-M.y)*fy;

	upperleftn.x = M.x + (upperleft.x-M.x)*fx;
	upperleftn.y = M.y + (upperleft.y-M.y)*fy;

	lowerleft.x=lowerleftn.x;
	lowerleft.y=lowerleftn.y;
	lowerright.x=lowerrightn.x;
	lowerright.y=lowerrightn.y;
	upperleft.x=upperleftn.x;
	upperleft.y=upperleftn.y;

}

void Ljapunow::rot(const int32_t deg) {
	// rotation by deg degree
	double w=deg*2.0*M_PI/360.0;
	double cosa=fastcos(w),sina=fastsin(w);

	Point32_t upperleftn,lowerleftn,lowerrightn,M,tmp;

	M.x = 0.5*(lowerright.x+upperleft.x);
	M.y = 0.5*(lowerright.y+upperleft.y);

	lowerleftn.x = (lowerleft.x-M.x)*cosa-sina*(lowerleft.y-M.y)+M.x;
	lowerleftn.y = (lowerleft.x-M.x)*sina+cosa*(lowerleft.y-M.y)+M.y;

	lowerrightn.x = (lowerright.x-M.x)*cosa-sina*(lowerright.y-M.y)+M.x;
	lowerrightn.y = (lowerright.x-M.x)*sina+cosa*(lowerright.y-M.y)+M.y;

	upperleftn.x = (upperleft.x-M.x)*cosa-sina*(upperleft.y-M.y)+M.x;
	upperleftn.y = (upperleft.x-M.x)*sina+cosa*(upperleft.y-M.y)+M.y;

	lowerleft.x=lowerleftn.x;
	lowerleft.y=lowerleftn.y;
	lowerright.x=lowerrightn.x;
	lowerright.y=lowerrightn.y;
	upperleft.x=upperleftn.x;
	upperleft.y=upperleftn.y;
}

void Ljapunow::setPosition(
	const double lolex,const double loley,
	const double lorix,const double loriy,
	const double uplex,const double upley
) {
	lowerleft.x=lolex;
	lowerleft.y=loley;
	lowerright.x=lorix;
	lowerright.y=loriy;
	upperleft.x=uplex;
	upperleft.y=upley;
}

void Ljapunow::setSequence(char *s) {
	s[255]=0;
	upper(s);
    seqlen=strlen(s);
    for(int32_t i=0;i<seqlen;i++) {
		if (s[i]=='A') sequence[i]=0;
		else if (s[i]=='B') sequence[i]=1;
		else {
			printf("Error in sequence.\n");
			seqlen=0;
			return;
		}
	}
}


// struct Bitmap

int32_t Bitmap::setlenxy(const int32_t xl,const int32_t yl) {
	if ((!bmp)||(xlen!=xl)||(ylen!=yl)) {
		disp();
		xlen=xl; ylen=yl;
		ybytes=3*xlen;
		bytes=ylen*ybytes;
		bmp=new uint8_t[bytes];
		if (!bmp) return 0;
	}
	
	return 1;
}

void Bitmap::save(const char* fn) {
	FILE *fbmp = fopen(fn,"wb");
	writehex(fbmp,"424DF6C62D00000000003600000028000000");
	uint32_t w = xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w = ylen; 
	fwrite(&w,sizeof(w),1,fbmp);
	writehex(fbmp,"0100180000000000C0C62D00C40E0000C40E00000000000000000000");
	fwrite(bmp,bytes,1,fbmp);

	fclose(fbmp);
}

Bitmap::Bitmap(void) {
	xlen=ylen=bytes=ybytes=0;
	bmp=0;
}

void Bitmap::disp(void) {
	if (bmp) { delete[] bmp; bmp=0; bytes=xlen=ylen=0; }
}

Bitmap::~Bitmap(void) {
	disp();
}


// main routine

int32_t main(int32_t argc,char** argv) {
	srand(time(NULL));
	ljap=new Ljapunow;
	ffarbe=NULL;
	int32_t iterfilecount=1;
	int32_t tilefilenr=1;

	char tmp[1000],utmp[1000];
	int32_t defect=0;

	while (1) {
		printf("\n\n\nLjapunow\n");
		if (defect>0) {
			printf("\n\nError in function. Load anew recommended\n\n");
		} else {
			if (ljap->fn[0]>0) printf("File %s\n",ljap->fn);
			if (ljap->fkt) printf("function %s\n",ljap->fkt->fktStr(tmp)); else printf("Function undefiniert\n");
			printf("upper left(%le|%le)\nlower left(%le|%le)\nlower right(%le|%le)\n",ljap->upperleft.x,ljap->upperleft.y,ljap->lowerleft.x,ljap->lowerleft.y,ljap->lowerright.x,ljap->lowerright.y);
			printf("Image size (%i|%i)\n",ljap->lenx,ljap->leny);
			printf("sequence %s\n",ljap->getSequence(tmp));
			printf("iterations (%i|%i)\n",ljap->iter0,ljap->iter1);
			printf("================================\n\n");
		}

		printf("\n> ");
		gets(tmp);
		defect=0; 

		chomp(tmp);	sprintf(utmp,"%s",tmp);
		for(int32_t i=(strlen(tmp)-1);i>=0;i--) if (tmp[i]==')') { tmp[i]=0; break; }
		upper(utmp);

		if (strcmp(utmp,"E")==0) break;
		else if (strstr(utmp,"LOAD(")==utmp) {
			char d[1024],fn[500],b[500];
			sprintf(d,"%s",&tmp[5]);
			stripext(d);
			sprintf(fn,"%s.par",d);
			if (ljap->loadpar(fn) > 0) {
				printf("Parameters loaded\n"); 
				sprintf(fn,"%s.ljd",d);
				// if not existent, no problem
				if (ljap->loadexp(fn) > 0) printf("Lyapunov values loaded\n");
				strcpy(ljap->fn,d);
			} else { 
				defect=1; 
				printf("Error loading %s\n",fn);
				ljap->fn[0]=0;
			}
		} else if (strstr(utmp,"LOADCOLOR(")==utmp) {
			char fn[1024];
			if (strstr(utmp,".PAR")==NULL) {
				sprintf(fn,"%s.par",&utmp[10]);
			} else strcpy(fn,&utmp[10]);

			if (ljap->loadcolor(fn) <= 0) defect=1;
		} else if (strstr(utmp,"WALKSEQ(")==utmp) {
			int32_t anz,slen;
			
			if (sscanf(&utmp[8],"%i,%i",&anz,&slen) != 2) { 
				printf("Error\n");
				continue; 
			}

			char ts[500];
			if (slen>64) slen=64;
			ts[slen]=0;
			
			for(int n=0;n<anz;n++) {
				for(int32_t i=0;i<slen;i++) ts[i]='A'+rand()%2;
				ljap->setSequence(ts);
				printf("%s ",ts);
				ljap->calc(0,ljap->leny-1);

				char fn[1024],orig[1024];
				sprintf(orig,"_walkseq_%04i_%s",n+1,ts); 
				sprintf(fn,"%s.bmp",orig);
				ljap->savebmp(fn,NULL);
				sprintf(fn,"%s.par",orig);
				ljap->savepar(fn);
				sprintf(fn,"%s.ljd",orig);
				ljap->saveexp(fn);
			} // n
		} else if (strstr(utmp,"SETSIZE(")==utmp) {
			int32_t xl,yl;
			if (sscanf(&utmp[8],"%i,%i",&xl,&yl) != 2) { printf("Error\n");continue; }
			ljap->setlen(xl,yl);
		} else if (strstr(utmp,"SETITER(")==utmp) {
			int32_t a,b;
			if (sscanf(&utmp[8],"%i,%i",&a,&b) != 2) { printf("Error\n");continue; }
			ljap->setiter(a,b);
		} else if (strstr(utmp,"ROTATEDEG(")==utmp) {
			int32_t w;
			if (sscanf(&utmp[10],"%i",&w) != 1) { printf("Error\n"); continue; }
			ljap->rot(w);
		} else if (strstr(utmp,"STRETCH(")==utmp) {
			double fx;
			if (sscanf(&tmp[8],"%le",&fx) != 1) { printf("Error\n");continue; }
			ljap->stretch(fx,fx);
		} else if (strstr(utmp,"SETSEQUENCE(")==utmp) {
			ljap->setSequence(&utmp[12]);
		} else if (strstr(utmp,"SETPOSITION(")==utmp) {
			double a,b,c,d,e,f;
			if (
				sscanf(&tmp[12],"%lf,%lf,%lf,%lf,%lf,%lf",&a,&b,&c,&d,&e,&f)
			== 6) ljap->setPosition(a,b,c,d,e,f);
			else printf("Errror\n");
		} else if (strstr(utmp,"WALKB(")==utmp) {
			double a,b;
			int32_t n;
			if (sscanf(&tmp[6],"%le,%le,%i",&a,&b,&n) != 3) {
				printf("Unknown parameters\n");
				continue;
			}
			IterDouble itd(a,b,n);

			if (ljap->fkt->id != ID_FKT_I) {
				ljap->fkt->set_iterb(&itd); 

				Bitmap bmp;
				ljap->iterStart();
				do {
					printf("b=%.10lf ",itd.wert);
					ljap->fkt->set_b(itd.wert);
					ljap->calc(0,ljap->leny-1);
					char fn[1000];
					sprintf(tmp,"_walkb%04i_b_%+.10lf",iterfilecount,itd.wert);
					sprintf(fn,"%s.bmp",tmp); 
					ljap->savebmp(fn,NULL);
					sprintf(fn,"%s.par",tmp); 
					ljap->savepar(fn);
					sprintf(fn,"%s.ljd",tmp); 
					ljap->saveexp(fn);
					iterfilecount++; 
				} while (ljap->iterWeiter());
			}
		} else if (!strcmp(utmp,"WALKSECTION")) {
			if (ljap->fkt->typ != FKTTYP_ABSCHNITTSWEISE) continue;

			FunctionMetaABSC *fvi=(FunctionMetaABSC*)ljap->fkt;

			int32_t ctr=1;

			char tmp[1024];
			const double i0START=-1.0;
			const double bis=1.0;
			const double delta=0.5;

			for(double i0min=i0START;i0min < bis;i0min += delta) {
				printf("i0=%f to %f\n",i0min,bis);
				for(double i0max=(i0min+delta);i0max < bis;i0max += delta) {
					for(double i1min=i0START;i1min < bis;i1min += delta) {
						printf("i1=%f to %f\n",i1min,bis);
						for(double i1max=(i1min+delta);i1max < bis;i1max += delta) {
							fvi->setsections(i0min,i0max,i1min,i1max);
							ljap->calc(0,ljap->leny-1);
							sprintf(tmp,"_walksection%04i.bmp",ctr); 
							ljap->savebmp(tmp,NULL);
							sprintf(tmp,"_walksection%04i.par",ctr); 
							ljap->savepar(tmp);
							ctr++;
						}
					}
				} 
			} 
		} else if (strcmp(utmp,"WALKRGB")==NULL) {
			Bitmap bmp;
			srand(time(NULL));
			int r1,g1,b1,r2,g2,b2,idx;
			for(int i=0;i<MAXRGBITERS;i++) {
				idx=rand()%ljap->farbe->intanz;
				#define RNDF rand()%256;
				r1=RNDF; g1=RNDF; b1=RNDF;
				r2=RNDF; g2=RNDF; b2=RNDF;
				
				ljap->farbe->ints[idx]->setfarbel(r1,g1,b1);
				ljap->farbe->ints[idx]->setfarber(r2,g2,b2);
				char fn[1024];
				sprintf(fn,"_walkrgb_%04i.bmp",i+1);
				ljap->savebmp(fn,&bmp);
				printf(".");
				sprintf(fn,"_walkrgb_%04i.par",i+1); 
				ljap->savepar(fn);

				char ff[1024]; 
				int32_t c=1; 
				if (getFirstColorFile(ff) <= 0) continue; 
				do { 
					if (ljap->loadcolor(ff) > 0) { 
						sprintf(fn,"_walkcolordir_%04i.par",c);  
						ljap->savepar(fn); 
						printf("."); 
						sprintf(fn,"_walkcolordir_%04i.bmp",c); 
						ljap->savebmp(fn,NULL);
						c++; 
					} 
				} while (getNextColorFile(ff) > 0); 
			} // i
		} else if (strstr(utmp,"WALKDET(")==utmp) {
			double b0,b1;
			int32_t fktid,abl0,abl1,n;
			if (sscanf(&tmp[8],"%i,%i,%i,%lf,%lf,%i",
				&fktid,&abl0,&abl1,&b0,&b1,&n) != 6) 
			{
				printf("Error parameters\n");
				continue;
			}

			char fn[1024];
			Function* fktp=getNewFunction(fktid);
			if (!fktp) {
				printf("Error. Function not recognized.\n");
				continue;
			}
			FunctionMetaDet* hierp=(FunctionMetaDet*)getNewFunction(ID_FKT_METADET);
			Function *sicp=ljap->fkt;
			ljap->fkt=hierp;
			hierp->f=fktp;

			for(int32_t abl=abl0;abl<=abl1;abl++) {
				if (
					(abl == ID_FKT_METADET) ||
					(abl == ID_FKT_METAABSC)
				) continue;
				
				Function* ablp=getNewFunction(abl);
				if (!ablp) continue; // not existent
				printf("derivative %i\n",abl);
				hierp->abl=ablp;

				IterDouble itd(b0,b1,n);
				ljap->fkt->set_iterb(&itd); // auch bei eigentlich nicht unterstützenden Functionen

				for(int32_t fwas=1;fwas<=2;fwas++) for(int32_t ablwas=1;ablwas<=2;ablwas++) {
					ljap->iterStart();
					hierp->fwas=fwas;
					hierp->ablwas=ablwas;

					do {
						ljap->fkt->set_b(itd.wert);
						printf("%lf",itd.wert);
						ljap->calc(0,ljap->leny-1);
						char fn[1000];
						sprintf(tmp,"_walkdet%02i_%02i_%04i_b_%+.10lf",fktid,abl,iterfilecount,itd.wert);
						sprintf(fn,"%s.bmp",tmp); 
						ljap->savebmp(fn,NULL);
						sprintf(fn,"%s.par",tmp); 
						ljap->savepar(fn);
						sprintf(fn,"%s.ljd",tmp); 
						ljap->saveexp(fn);

					} while (ljap->iterWeiter());
					iterfilecount++;
				} 

				delete ablp;
			}

			delete fktp;
			delete hierp;
			ljap->fkt=sicp;
		} else if ( (strstr(utmp,"CROP(")==utmp) || (strstr(utmp,"C(")==utmp) ) {
			int32_t pulneux,pulneuy,porneux,porneuy;

			if (sscanf(&utmp[5],"%i,%i,%i,%i",&pulneux,&pulneuy,&porneux,&porneuy) != 4) { printf("Error\n");continue; }

			ljap->crop(pulneux,pulneuy,porneux,porneuy);
		} else if (strstr(utmp,"WALKTILE(")==utmp) {
			int32_t anzx,anzy;

			if (sscanf(&utmp[9],"%i,%i",&anzx,&anzy) != 2) { printf("Error\n");continue; }
			sprintf(tmp,"%04i",tilefilenr++);

			ljap->tile(tmp,anzx,anzy);
		} else if (strstr(utmp,"CENTER(")==utmp) {
			int32_t x,y;
			if (sscanf(&utmp[7],"%i,%i",&x,&y) != 2) { printf("Error\n");continue; }
			ljap->centerPixel(x,y);
		} else if ( (strstr(utmp,"SAVE(")==utmp) || (strstr(utmp,"WR(")==utmp) ) {
			char fn[1000];
			char fn2[1024]; strcpy(fn2,&tmp[5]);
			int32_t lp=strlen(fn2)-1;
			while (lp>=0) if (fn2[lp] == '.') { fn2[lp]=0; break; } else lp--;
			sprintf(fn,"%s.bmp",fn2); ljap->savebmp(fn,0);
			sprintf(fn,"%s.par",fn2); ljap->savepar(fn);
			sprintf(fn,"%s.ljd",fn2); ljap->saveexp(fn);
			sprintf(fn,"%s.descr",&tmp[5]);
			if (ljap->lenx >= 600) ljap->savedescr(fn);
		} else if (strstr(utmp,"RUN")==utmp) {
			int32_t start,ende;
			if (strcmp(utmp,"RUN")==0) {
				start=0;
				ende=ljap->leny-1;
			} else {
				if (sscanf(&tmp[4],"%i,%i",&start,&ende) != 2) {
					start=0;
					ende=ljap->leny-1;
				}
			}
			time_t a,b;
			a=time(NULL);
			ljap->calc(start,ende);
			b=time(NULL);
			double d=difftime(b,a);
			printf("Time used %.2lf sec\n",d);
			sprintf(tmp,"tmpljap.bmp"); ljap->savebmp(tmp,NULL);
			sprintf(tmp,"tmpljap.par"); ljap->savepar(tmp);
			sprintf(tmp,"tmpljap.ljd"); ljap->saveexp(tmp);
		}
	} // while

	delete ljap;

	return 0;
}
