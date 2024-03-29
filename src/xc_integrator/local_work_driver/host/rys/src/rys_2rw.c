#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define MAX(a,b)    ((a) < (b) ? (b) : (a))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))

void rys_2rw(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]) {
  int jump2[41] =
    { 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6,
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8
    };

  double e;
  int m, n;
  double t, x, y, f1, r1, r2, w1, w2;
  int tcase;

  m = 0;
  for (n = 0; n < nt; ++n) {
    t = tval[n];
    if (t <= 3e-7) {
      r1 = .130693606237085 - t * .0290430236082028;
      r2 = 2.86930639376291 - t * .637623643058102;
      wts[m] *= .652145154862545 - t * .122713621927067;
      wts[m + 1] *= .347854845137453 - t * .210619711404725;
      rts[m] = r1 / (r1 + 1.);
      rts[m + 1] = r2 / (r2 + 1.);
      m += 2;
      goto L200;
    }

    tcase = (int) MIN ((t + 1.0), 41.);
    switch (jump2[tcase - 1])
      {
      case 1:
	goto L2100;
      case 2:
	goto L2200;
      case 3:
	goto L2300;
      case 4:
	goto L2400;
      case 5:
	goto L2500;
      case 6:
	goto L2600;
      case 7:
	goto L2700;
      case 8:
	goto L2800;
      }

  L2100:
    f1 = ((((((((t * -8.36313918003957e-8 + 1.21222603512827e-6) * t -
		1.15662609053481e-5) * t + 9.25197374512647e-5) * t -
	      6.40994113129432e-4) * t + .00378787044215009) * t -
	    .0185185172458485) * t + .0714285713298222) * t -
	  .199999999997023) * t + .333333333333318;
    w1 = (t + t) * f1 + exp(-t);
    r1 = (((((((t * -2.35234358048491e-9 + 2.49173650389842e-8) * t -
	       4.558315364581e-8) * t - 2.447252174587e-6) * t +
	     4.743292959463e-5) * t - 5.33184749432408e-4) * t +
	   .00444654947116579) * t - .0290430236084697) * t +
      .130693606237085;
    r2 = (((((((t * -2.4740490232917e-8 + 2.36809910635906e-7) * t +
	       1.83536773631e-6) * t - 2.066168802076e-5) * t -
	     1.345693393936e-4) * t - 5.88154362858038e-5) * t +
	   .0532735082098139) * t - .637623643056745) * t +
      2.86930639376289;
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2200:
    x = t - 2.;
    f1 = ((((((((((x * -1.61702782425558e-10 + 1.96215250865776e-9) * x -
		  2.14234468198419e-8) * x + 2.17216556336318e-7) * x -
		1.98850171329371e-6) * x + 1.62429321438911e-5) * x -
	      1.16740298039895e-4) * x + 7.24888732052332e-4) * x -
	    .00379490003707156) * x + .0161723488664661) * x -
	  .0529428148329736) * x + .115702180856167;
    w1 = (t + t) * f1 + exp(-t);
    r1 = (((((((((x * -6.36859636616415e-12 + 8.4741706477627e-11) * x -
		 5.152207846962e-10) * x - 3.846389873308e-10) * x +
	       8.47225338838e-8) * x - 1.85306035634293e-6) * x +
	     2.47191693238413e-5) * x - 2.49018321709815e-4) * x +
	   .00219173220020161) * x - .0163329339286794) * x +
      .0868085688285261;
    r2 = (((((((((x * 1.45331350488343e-10 + 2.07111465297976e-9) * x -
		 1.878920917404e-8) * x - 1.725838516261e-7) * x +
	       2.247389642339e-6) * x + 9.76783813082564e-6) * x -
	     1.93160765581969e-4) * x - .00158064140671893) * x +
	   .0485928174507904) * x - .430761584997596) * x +
      1.8040097453795;
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2300:
    x = t - 4.;
    f1 = ((((((((((x * -2.62453564772299e-11 + 3.24031041623823e-10) * x
		  - 3.614965656163e-9) * x + 3.760256799971e-8) * x -
		3.553558319675e-7) * x + 3.022556449731e-6) * x -
	      2.290098979647e-5) * x + 1.526537461148e-4) * x -
	    8.81947375894379e-4) * x + .00433207949514611) * x -
	  .0175257821619926) * x + .0528406320615584;
    w1 = (t + t) * f1 + exp(-t);
    r1 = ((((((((x * -4.11560117487296e-12 + 7.10910223886747e-11) * x -
		1.73508862390291e-9) * x + 5.93066856324744e-8) * x -
	      9.76085576741771e-7) * x + 1.08484384385679e-5) * x -
	    1.12608004981982e-4) * x + .00116210907653515) * x -
	  .00989572595720351) * x + .0612589701086408;
    r2 = (((((((((x * -1.80555625241001e-10 + 5.44072475994123e-10) * x +
		 1.60349804524e-8) * x - 1.497986283037e-7) * x -
	       7.017002532106e-7) * x + 1.85882653064034e-5) * x -
	     2.04685420150802e-5) * x - .00249327728643089) * x +
	   .0356550690684281) * x - .260417417692375) * x +
      1.12155283108289;
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2400:
    e = exp(-t);
    x = 1. / t;
    y = t - 7.5;
    w1 = ((((((x * .46897511375022 - .69955602298985) * x +
	      .53689283271887) * x - .32883030418398) * x +
	    .24645596956002) * x - .49984072848436) * x -
	  3.1501078774085e-6) * e + sqrt(x * .785398163397448);
    f1 = (w1 - e) / (t + t);
    r1 = (((((((((((((y * -1.43632730148572e-16 + 2.38198922570405e-16) *
		     y + 1.3583196188e-14) * y - 7.064522786879e-14) * y -
		   7.719300212748e-13) * y + 7.802544789997e-12) * y +
		 6.628721099436e-11) * y - 1.775564159743e-9) * y +
	       1.71382882399e-8) * y - 1.497500187053e-7) * y +
	     2.283485114279e-6) * y - 3.76953869614706e-5) * y +
	   4.74791204651451e-4) * y - .00460448960876139) * y +
      .0372458587837249;
    r2 = ((((((((((((y * 2.487916227989e-14 - 1.36113510175724e-13) * y -
		    2.224334349799e-12) * y + 4.190559455515e-11) * y -
		  2.222722579924e-10) * y - 2.624183464275e-9) * y +
		6.128153450169e-8) * y - 4.383376014528e-7) * y -
	      2.4995220023291e-6) * y + 1.0323664788832e-4) * y -
	    .00144614664924989) * y + .0135094294917224) * y -
	  .0953478510453887) * y + .54476524568679;
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2500:
    e = exp(-t);
    x = 1. / t;
    w1 = (((x * -.18784686463512 + .22991849164985) * x - .49893752514047)
	  * x - 2.1916512131607e-5) * e + sqrt(x * .785398163397448);
    f1 = (w1 - e) / (t + t);
    r1 = ((((t * -1.01041157064226e-5 + .00119483054115173) * t -
	    .0673760231824074) * t + 1.25705571069895) * t + (((x *
								-8576.09422987199
								+
								5910.05939591842)
							       * x -
							       1708.07677109425)
							      * x +
							      264.536689959503)
	  * x - 23.8570496490846) * e + .275255128608411 / (t -
							    .275255128608411);
    r2 = (((t * 3.39024225137123e-4 - .0934976436343509) * t -
	   4.2221648330632) * t +
	  (((x * -2084.57050986847 - 1049.99071905664) * x +
	    339.891508992661) * x - 156.184800325063) * x +
	  8.00839033297501) * e + 2.72474487139158 / (t -
						      2.72474487139158);
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2600:
    e = exp(-t);
    x = 1. / t;
    w1 = ((x * .1962326414943 - .4969524146449) * x - 6.0156581186481e-5)
      * e + sqrt(x * .785398163397448);
    f1 = (w1 - e) / (t + t);
    r1 = ((((t * -1.14906395546354e-6 + 1.76003409708332e-4) * t -
	    .0171984023644904) * t - .137292644149838) * t + (x *
							      -47.5742064274859
							      +
							      9.21005186542857)
	  * x - .0231080873898939) * e + .275255128608411 / (t -
							     .275255128608411);
    r2 = (((t * 3.64921633404158e-4 - .0971850973831558) * t -
	   4.02886174850252) * t + (x * -135.831002139173 -
				    86.6891724287962) * x +
	  2.98011277766958) * e + 2.72474487139158 / (t -
						      2.72474487139158);
    w2 = ((f1 - w1) * r1 + f1) * (r2 + 1.) / (r2 - r1);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2700:
    e = exp(-t);
    w1 = sqrt(.785398163397448 / t);
    w2 = (t * 4.468573893084 - 77.9250653461045) * e + w1 *
      .0917517095361369;
    r1 = (t * -.87894730749888 + 10.9243702330261) * e + .275255128608411
      / (t - .275255128608411);
    r2 = (t * -9.28903924275977 + 81.0642367843811) * e +
      2.72474487139158 / (t - 2.72474487139158);
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = r1 / (r1 + 1.);
    rts[m + 1] = r2 / (r2 + 1.);
    m += 2;
    goto L200;

  L2800:
    w1 = sqrt(.785398163397448 / t);
    w2 = w1 * .0917517095361369;
    wts[m] *= w1 - w2;
    wts[m + 1] *= w2;
    rts[m] = .275255128608411 / t;
    rts[m + 1] = 2.72474487139158 / t;
    m += 2;
  L200:
    ;
  }
}
