/*
 * Copyright (c) 2012-2020 MIRACL UK Ltd.
 *
 * This file is part of MIRACL Core
 * (see https://github.com/miracl/core).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* CORE Elliptic Curve Functions */
/* SU=m, SU is Stack Usage (Weierstrass Curves) */

//#define HAS_MAIN

#include "ecp_ED25519.h"

/* test for P=O point-at-infinity */
int ECP_ED25519_isinf(ECP_ED25519 *P)
{

#if CURVETYPE_ED25519==EDWARDS
    return (FP_F25519_iszilch(&(P->x)) && FP_F25519_equals(&(P->y), &(P->z)));
#endif
#if CURVETYPE_ED25519==WEIERSTRASS
    return (FP_F25519_iszilch(&(P->x)) && FP_F25519_iszilch(&(P->z)));
#endif
#if CURVETYPE_ED25519==MONTGOMERY
    return FP_F25519_iszilch(&(P->z));
#endif

}

/* Conditional swap of P and Q dependant on d */
static void ECP_ED25519_cswap(ECP_ED25519 *P, ECP_ED25519 *Q, int d)
{
    FP_F25519_cswap(&(P->x), &(Q->x), d);
#if CURVETYPE_ED25519!=MONTGOMERY
    FP_F25519_cswap(&(P->y), &(Q->y), d);
#endif
    FP_F25519_cswap(&(P->z), &(Q->z), d);
}

#if CURVETYPE_ED25519!=MONTGOMERY
/* Conditional move Q to P dependant on d */
static void ECP_ED25519_cmove(ECP_ED25519 *P, ECP_ED25519 *Q, int d)
{
    FP_F25519_cmove(&(P->x), &(Q->x), d);
#if CURVETYPE_ED25519!=MONTGOMERY
    FP_F25519_cmove(&(P->y), &(Q->y), d);
#endif
    FP_F25519_cmove(&(P->z), &(Q->z), d);
}

/* return 1 if b==c, no branching */
static int teq(sign32 b, sign32 c)
{
    sign32 x = b ^ c;
    x -= 1; // if x=0, x now -1
    return (int)((x >> 31) & 1);
}
#endif // CURVETYPE_ED25519!=MONTGOMERY

#if CURVETYPE_ED25519!=MONTGOMERY
/* Constant time select from pre-computed table */
static void ECP_ED25519_select(ECP_ED25519 *P, ECP_ED25519 W[], sign32 b)
{
    ECP_ED25519 MP;
    sign32 m = b >> 31;
    sign32 babs = (b ^ m) - m;

    babs = (babs - 1) / 2;

    ECP_ED25519_cmove(P, &W[0], teq(babs, 0)); // conditional move
    ECP_ED25519_cmove(P, &W[1], teq(babs, 1));
    ECP_ED25519_cmove(P, &W[2], teq(babs, 2));
    ECP_ED25519_cmove(P, &W[3], teq(babs, 3));
    ECP_ED25519_cmove(P, &W[4], teq(babs, 4));
    ECP_ED25519_cmove(P, &W[5], teq(babs, 5));
    ECP_ED25519_cmove(P, &W[6], teq(babs, 6));
    ECP_ED25519_cmove(P, &W[7], teq(babs, 7));

    ECP_ED25519_copy(&MP, P);
    ECP_ED25519_neg(&MP);  // minus P
    ECP_ED25519_cmove(P, &MP, (int)(m & 1));
}
#endif

/* Test P == Q */
/* SU=168 */
int ECP_ED25519_equals(ECP_ED25519 *P, ECP_ED25519 *Q)
{
    FP_F25519 a, b;

    FP_F25519_mul(&a, &(P->x), &(Q->z));
    FP_F25519_mul(&b, &(Q->x), &(P->z));
    if (!FP_F25519_equals(&a, &b)) return 0;

#if CURVETYPE_ED25519!=MONTGOMERY
    FP_F25519_mul(&a, &(P->y), &(Q->z));
    FP_F25519_mul(&b, &(Q->y), &(P->z));
    if (!FP_F25519_equals(&a, &b)) return 0;
#endif

    return 1;

}

/* Set P=Q */
/* SU=16 */
void ECP_ED25519_copy(ECP_ED25519 *P, ECP_ED25519 *Q)
{
    FP_F25519_copy(&(P->x), &(Q->x));
#if CURVETYPE_ED25519!=MONTGOMERY
    FP_F25519_copy(&(P->y), &(Q->y));
#endif
    FP_F25519_copy(&(P->z), &(Q->z));
}

/* Set P=-Q */
#if CURVETYPE_ED25519!=MONTGOMERY
/* SU=8 */
void ECP_ED25519_neg(ECP_ED25519 *P)
{
#if CURVETYPE_ED25519==WEIERSTRASS
    FP_F25519_neg(&(P->y), &(P->y));
    FP_F25519_norm(&(P->y));
#else
    FP_F25519_neg(&(P->x), &(P->x));
    FP_F25519_norm(&(P->x));
#endif

}
#endif

/* Set P=O */
void ECP_ED25519_inf(ECP_ED25519 *P)
{
    FP_F25519_zero(&(P->x));
#if CURVETYPE_ED25519!=MONTGOMERY
    FP_F25519_one(&(P->y));
#endif
#if CURVETYPE_ED25519!=EDWARDS
    FP_F25519_zero(&(P->z));
#else
    FP_F25519_one(&(P->z));
#endif
}

/* Calculate right Hand Side of curve equation y^2=RHS */
/* SU=56 */
void ECP_ED25519_rhs(FP_F25519 *v, FP_F25519 *x)
{
#if CURVETYPE_ED25519==WEIERSTRASS
    /* x^3+Ax+B */
    FP_F25519 t;
    FP_F25519_sqr(&t, x);
    FP_F25519_mul(&t, &t, x);

#if CURVE_A_ED25519 == -3

        FP_F25519_neg(v, x);
        FP_F25519_norm(v);
        FP_F25519_imul(v, v, -CURVE_A_ED25519);
        FP_F25519_norm(v);
        FP_F25519_add(v, &t, v);
#else 
        FP_F25519_copy(v, &t);
#endif

    FP_F25519_rcopy(&t, CURVE_B_ED25519);

    FP_F25519_add(v, &t, v);
    FP_F25519_reduce(v);
#endif

#if CURVETYPE_ED25519==EDWARDS
    /* (Ax^2-1)/(Bx^2-1) */
    FP_F25519 t, one;
    FP_F25519_sqr(v, x);
    FP_F25519_one(&one);
    FP_F25519_rcopy(&t, CURVE_B_ED25519);

    FP_F25519_mul(&t, v, &t);
    FP_F25519_sub(&t, &t, &one);
    FP_F25519_norm(&t);
#if CURVE_A_ED25519 == 1
        FP_F25519_sub(v, v, &one);
#endif

#if CURVE_A_ED25519 == -1
        FP_F25519_add(v, v, &one);
        FP_F25519_norm(v);
        FP_F25519_neg(v, v);
#endif
    FP_F25519_norm(v);
    FP_F25519_inv(&t, &t, NULL);
    FP_F25519_mul(v, v, &t);
    FP_F25519_reduce(v);
#endif

#if CURVETYPE_ED25519==MONTGOMERY
    /* x^3+Ax^2+x */
    FP_F25519 x2, x3;
    FP_F25519_sqr(&x2, x);
    FP_F25519_mul(&x3, &x2, x);
    FP_F25519_copy(v, x);
    FP_F25519_imul(&x2, &x2, CURVE_A_ED25519);
    FP_F25519_add(v, v, &x2);
    FP_F25519_add(v, v, &x3);
    FP_F25519_reduce(v);
#endif
}

#if CURVETYPE_ED25519==MONTGOMERY

/* Set P=(x,{y}) */

int ECP_ED25519_set(ECP_ED25519 *P, BIG_256_56 x)
{
    //BIG_256_56 m, b;
    FP_F25519 rhs;
    //BIG_256_56_rcopy(m, Modulus_F25519);

    FP_F25519_nres(&rhs, x);

    ECP_ED25519_rhs(&rhs, &rhs);

    //FP_F25519_redc(b, &rhs);
    //if (BIG_256_56_jacobi(b, m) != 1)
    if (!FP_F25519_qr(&rhs,NULL))
    {
        ECP_ED25519_inf(P);
        return 0;
    }

    FP_F25519_nres(&(P->x), x);
    FP_F25519_one(&(P->z));
    return 1;
}

/* Extract x coordinate as BIG */
int ECP_ED25519_get(BIG_256_56 x, ECP_ED25519 *P)
{
    ECP_ED25519 W;
    ECP_ED25519_copy(&W, P);
    ECP_ED25519_affine(&W);
    if (ECP_ED25519_isinf(&W)) return -1;
    FP_F25519_redc(x, &(W.x));
    return 0;
}


#else
/* Extract (x,y) and return sign of y. If x and y are the same return only x */
/* SU=16 */
int ECP_ED25519_get(BIG_256_56 x, BIG_256_56 y, ECP_ED25519 *P)
{
    ECP_ED25519 W;
    ECP_ED25519_copy(&W, P);
    ECP_ED25519_affine(&W);

    if (ECP_ED25519_isinf(&W)) return -1;

    FP_F25519_redc(y, &(W.y));
    FP_F25519_redc(x, &(W.x));

    return FP_F25519_sign(&(W.x));
}

/* Set P=(x,{y}) */
/* SU=96 */
int ECP_ED25519_set(ECP_ED25519 *P, BIG_256_56 x, BIG_256_56 y)
{
    FP_F25519 rhs, y2;

    FP_F25519_nres(&y2, y);
    FP_F25519_sqr(&y2, &y2);
    FP_F25519_reduce(&y2);

    FP_F25519_nres(&rhs, x);
    ECP_ED25519_rhs(&rhs, &rhs);

    if (!FP_F25519_equals(&y2, &rhs))
    {
        ECP_ED25519_inf(P);
        return 0;
    }

    FP_F25519_nres(&(P->x), x);
    FP_F25519_nres(&(P->y), y);
    FP_F25519_one(&(P->z));
    return 1;
}

/* Set P=(x,y), where y is calculated from x with sign s */
/* SU=136 */
int ECP_ED25519_setx(ECP_ED25519 *P, BIG_256_56 x, int s)
{
    FP_F25519 rhs,hint;
    FP_F25519_nres(&rhs, x);

    ECP_ED25519_rhs(&rhs, &rhs);

    if (!FP_F25519_qr(&rhs,&hint))
    {
        ECP_ED25519_inf(P);
        return 0;
    }

    FP_F25519_nres(&(P->x), x);
    FP_F25519_sqrt(&(P->y), &rhs,&hint);

    if (FP_F25519_sign(&(P->y))!=s)
        FP_F25519_neg(&(P->y), &(P->y));
    FP_F25519_reduce(&(P->y));
    FP_F25519_one(&(P->z));
    return 1;
}

#endif

void ECP_ED25519_cfp(ECP_ED25519 *P)
{   /* multiply point by curves cofactor */
    BIG_256_56 c;
    int cf = CURVE_Cof_I_ED25519;
    if (cf == 1) return;
    if (cf == 4)
    {
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        return;
    }
    if (cf == 8)
    {
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        return;
    }
    BIG_256_56_rcopy(c, CURVE_Cof_ED25519);
    ECP_ED25519_mul(P, c);
    return;
}

/* Hunt and Peck a BIG to a curve point */
/*
void ECP_ED25519_hap2point(ECP_ED25519 *P,BIG_256_56 h)
{
    BIG_256_56 x;
    BIG_256_56_copy(x,h);

	for (;;)
	{
#if CURVETYPE_ED25519!=MONTGOMERY
		ECP_ED25519_setx(P,x,0);
#else
		ECP_ED25519_set(P,x);
#endif
		BIG_256_56_inc(x,1); BIG_256_56_norm(x);
		if (!ECP_ED25519_isinf(P)) break;
	}
}
*/
/* Constant time Map to Point */
void ECP_ED25519_map2point(ECP_ED25519 *P,FP_F25519 *h)
{
#if CURVETYPE_ED25519==MONTGOMERY
// Elligator 2
    int qres;
    BIG_256_56 a;
    FP_F25519 X1,X2,w,t,one,A,N,D,hint;
    //BIG_256_56_zero(a); BIG_256_56_inc(a,CURVE_A_ED25519); BIG_256_56_norm(a); FP_F25519_nres(&A,a);
    FP_F25519_from_int(&A,CURVE_A_ED25519);
    FP_F25519_copy(&t,h);
    FP_F25519_sqr(&t,&t);   // t^2

    if (PM1D2_F25519 == 2)
        FP_F25519_add(&t,&t,&t);  // 2t^2
    if (PM1D2_F25519 == 1)
        FP_F25519_neg(&t,&t);      // -t^2
    if (PM1D2_F25519 > 2)
        FP_F25519_imul(&t,&t,QNRI_F25519); // precomputed QNR
    FP_F25519_norm(&t);  // z.t^2

    FP_F25519_one(&one);
    FP_F25519_add(&D,&t,&one); FP_F25519_norm(&D);  // Denominator D=1+z.t^2

    FP_F25519_copy(&X1,&A);
    FP_F25519_neg(&X1,&X1);  FP_F25519_norm(&X1);  // X1 = -A/D
    FP_F25519_copy(&X2,&X1);
    FP_F25519_mul(&X2,&X2,&t);             // X2 = -At/D

    FP_F25519_sqr(&w,&X1); FP_F25519_mul(&N,&w,&X1);  // w=X1^2, N=X1^3
    FP_F25519_mul(&w,&w,&A); FP_F25519_mul(&w,&w,&D); FP_F25519_add(&N,&N,&w);  // N = X1^3+ADX1^2
    FP_F25519_sqr(&t,&D);
    FP_F25519_mul(&t,&t,&X1);  
    FP_F25519_add(&N,&N,&t);  // N=X1^3+ADX1^2+D^2X1  // Numerator of gx =  N/D^3
    FP_F25519_norm(&N);

    FP_F25519_mul(&t,&N,&D);                   // N.D
    qres=FP_F25519_qr(&t,&hint);  // *** exp
    FP_F25519_inv(&w,&t,&hint);
    FP_F25519_mul(&D,&w,&N);     // 1/D
    FP_F25519_mul(&X1,&X1,&D);    // get X1
    FP_F25519_mul(&X2,&X2,&D);    // get X2
    FP_F25519_cmove(&X1,&X2,1-qres);
    FP_F25519_redc(a,&X1);

    ECP_ED25519_set(P,a);
    return;
#endif

#if CURVETYPE_ED25519==EDWARDS
// Elligator 2 - map to Montgomery, place point, map back
    int qres,ne,rfc,qnr;
    BIG_256_56 x,y;
    FP_F25519 X1,X2,t,w,one,A,w1,w2,B,Y,K,D,hint;
    FP_F25519_one(&one);

#if MODTYPE_F25519 != GENERALISED_MERSENNE
// its NOT goldilocks!
// Figure out the Montgomery curve parameters

    FP_F25519_rcopy(&B,CURVE_B_ED25519);
#if CURVE_A_ED25519 ==  1
        FP_F25519_add(&A,&B,&one);  // A=B+1
        FP_F25519_sub(&B,&B,&one);   // B=B-1
#else
        FP_F25519_sub(&A,&B,&one);  // A=B-1
        FP_F25519_add(&B,&B,&one);  // B=B+1
#endif
    FP_F25519_norm(&A); FP_F25519_norm(&B);

    FP_F25519_div2(&A,&A);    // (A+B)/2
    FP_F25519_div2(&B,&B);    // (B-A)/2
    FP_F25519_div2(&B,&B);    // (B-A)/4

    FP_F25519_neg(&K,&B); FP_F25519_norm(&K);
    //FP_F25519_inv(&K,&K,NULL);    // K
    FP_F25519_invsqrt(&K,&w1,&K);

    rfc=RIADZ_F25519;
    if (rfc)
    { // RFC7748 method applies
        FP_F25519_mul(&A,&A,&K);   // = J
        FP_F25519_mul(&K,&K,&w1);
        //FP_F25519_sqrt(&K,&K,NULL);
    } else { // generic method
        FP_F25519_sqr(&B,&B);
    }
#else
    FP_F25519_from_int(&A,156326);
    rfc=1;
#endif
// Map to this Montgomery curve X^2=X^3+AX^2+BX

    FP_F25519_copy(&t,h); 
    FP_F25519_sqr(&t,&t);   // t^2

    if (PM1D2_F25519 == 2)
    {
        FP_F25519_add(&t,&t,&t);  // 2t^2
        qnr=2;
    }
    if (PM1D2_F25519 == 1)
    {
        FP_F25519_neg(&t,&t);      // -t^2
        qnr=-1;
    }
    if (PM1D2_F25519 > 2)
    {
        FP_F25519_imul(&t,&t,QNRI_F25519);  // precomputed QNR
        qnr=QNRI_F25519;
    }
    FP_F25519_norm(&t);
    FP_F25519_add(&D,&t,&one); FP_F25519_norm(&D);  // Denominator=(1+z.u^2)

    FP_F25519_copy(&X1,&A);
    FP_F25519_neg(&X1,&X1);  FP_F25519_norm(&X1);    // X1=-(J/K).inv(1+z.u^2)
    FP_F25519_mul(&X2,&X1,&t);  // X2=X1*z.u^2

// Figure out RHS of Montgomery curve in rational form gx1/d^3

    FP_F25519_sqr(&w,&X1); FP_F25519_mul(&w1,&w,&X1);  // w=X1^2, w1=X1^3
    FP_F25519_mul(&w,&w,&A); FP_F25519_mul(&w,&w,&D); FP_F25519_add(&w1,&w1,&w);  // w1 = X1^3+ADX1^2
    FP_F25519_sqr(&w2,&D);
    if (!rfc)
    {
        FP_F25519_mul(&w,&X1,&B);
        FP_F25519_mul(&w2,&w2,&w); //
        FP_F25519_add(&w1,&w1,&w2);   // w1=X1^3+ADX1^2+BD^2X1
    } else {
        FP_F25519_mul(&w2,&w2,&X1);  
        FP_F25519_add(&w1,&w1,&w2);  // w1=X1^3+ADX1^2+D^2X1  // was &X1
    }
    FP_F25519_norm(&w1);

    FP_F25519_mul(&B,&w1,&D);     // gx1=num/den^3 - is_qr num*den (same as num/den, same as num/den^3)
    qres=FP_F25519_qr(&B,&hint);  // ***
    FP_F25519_inv(&w,&B,&hint);
    FP_F25519_mul(&D,&w,&w1);     // 1/D
    FP_F25519_mul(&X1,&X1,&D);    // get X1
    FP_F25519_mul(&X2,&X2,&D);    // get X2
    FP_F25519_sqr(&D,&D);

    FP_F25519_imul(&w1,&B,qnr);       // now for gx2 = Z.u^2.gx1
    FP_F25519_rcopy(&w,CURVE_HTPC_ED25519);   // qnr^C3  
    FP_F25519_mul(&w,&w,&hint);    // modify hint for gx2
    FP_F25519_mul(&w2,&D,h);

    FP_F25519_cmove(&X1,&X2,1-qres);  // pick correct one
    FP_F25519_cmove(&B,&w1,1-qres);
    FP_F25519_cmove(&hint,&w,1-qres);
    FP_F25519_cmove(&D,&w2,1-qres);
     
    FP_F25519_sqrt(&Y,&B,&hint);   // sqrt(num*den)
    FP_F25519_mul(&Y,&Y,&D);       // sqrt(num/den^3)

/*
    FP_F25519_sqrt(&Y,&B,&hint);   // sqrt(num*den)
    FP_F25519_mul(&Y,&Y,&D);       // sqrt(num/den^3)

    FP_F25519_imul(&B,&B,qnr);     // now for gx2 = Z.u^2.gx1
    FP_F25519_rcopy(&w,CURVE_HTPC_ED25519);   // qnr^C3  
    FP_F25519_mul(&hint,&hint,&w);    // modify hint for gx2

    FP_F25519_sqrt(&Y3,&B,&hint);  // second candidate
    FP_F25519_mul(&D,&D,h);
    FP_F25519_mul(&Y3,&Y3,&D);

    FP_F25519_cmove(&X1,&X2,1-qres);  // pick correct one
    FP_F25519_cmove(&Y,&Y3,1-qres);
*/
// correct sign of Y
    FP_F25519_neg(&w,&Y); FP_F25519_norm(&w);
    FP_F25519_cmove(&Y,&w,qres^FP_F25519_sign(&Y));

    if (!rfc)
    {
        FP_F25519_mul(&X1,&X1,&K);
        FP_F25519_mul(&Y,&Y,&K);
    }

#if MODTYPE_F25519 == GENERALISED_MERSENNE
// GOLDILOCKS isogeny
    FP_F25519_sqr(&t,&X1);  // t=u^2
    FP_F25519_add(&w,&t,&one); // w=u^2+1
    FP_F25519_norm(&w);
    FP_F25519_sub(&t,&t,&one); // t=u^2-1
    FP_F25519_norm(&t);
    FP_F25519_mul(&w1,&t,&Y);  // w1=v(u^2-1)
    FP_F25519_add(&w1,&w1,&w1);
    FP_F25519_add(&X2,&w1,&w1);
    FP_F25519_norm(&X2);       // w1=4v(u^2-1)
    FP_F25519_sqr(&t,&t);      // t=(u^2-1)^2
    FP_F25519_sqr(&Y,&Y);      // v^2
    FP_F25519_add(&Y,&Y,&Y);
    FP_F25519_add(&Y,&Y,&Y);
    FP_F25519_norm(&Y);        // 4v^2
    FP_F25519_add(&B,&t,&Y);  // w2=(u^2-1)^2+4v^2
    FP_F25519_norm(&B);                                   // X1=w1/w2 - X2=w1, B=w2

    FP_F25519_sub(&w2,&Y,&t);   // w2=4v^2-(u^2-1)^2
    FP_F25519_norm(&w2);
    FP_F25519_mul(&w2,&w2,&X1); // w2=u(4v^2-(u^2-1)^2)
    FP_F25519_mul(&t,&t,&X1);   // t=u(u^2-1)^2
    FP_F25519_div2(&Y,&Y);      // 2v^2
    FP_F25519_mul(&w1,&Y,&w);  // w1=2v^2(u^2+1)
    FP_F25519_sub(&w1,&t,&w1);  // w1=u(u^2-1)^2 - 2v^2(u^2+1)
    FP_F25519_norm(&w1);                                   // Y=w2/w1

    FP_F25519_mul(&t,&X2,&w1);    // output in projective to avoid inversion
    FP_F25519_copy(&(P->x),&t);
    FP_F25519_mul(&t,&w2,&B);
    FP_F25519_copy(&(P->y),&t);
    FP_F25519_mul(&t,&w1,&B);
    FP_F25519_copy(&(P->z),&t);

    return;

#else
    FP_F25519_add(&w1,&X1,&one); FP_F25519_norm(&w1); // (s+1)
    FP_F25519_sub(&w2,&X1,&one); FP_F25519_norm(&w2); // (s-1)
    FP_F25519_mul(&t,&w1,&Y);
    FP_F25519_mul(&X1,&X1,&w1);

    if (rfc)
        FP_F25519_mul(&X1,&X1,&K);

    FP_F25519_mul(&Y,&Y,&w2);      // output in projective to avoid inversion
    FP_F25519_copy(&(P->x),&X1);
    FP_F25519_copy(&(P->y),&Y);
    FP_F25519_copy(&(P->z),&t);
    return;
#endif

#endif

#if CURVETYPE_ED25519==WEIERSTRASS
    int sgn,ne;
    BIG_256_56 a,x,y;
    FP_F25519 X1,X2,X3,t,w,one,A,B,Y,D;
    FP_F25519 D2,hint,GX1;

#if HTC_ISO_ED25519 != 0
// Map to point on isogenous curve
    int i,k,isox,isoy,iso=HTC_ISO_ED25519;
    FP_F25519 xnum,xden,ynum,yden;
    BIG_256_56 z;
    FP_F25519_rcopy(&A,CURVE_Ad_ED25519);
    FP_F25519_rcopy(&B,CURVE_Bd_ED25519);
#else
    FP_F25519_from_int(&A,CURVE_A_ED25519);
    FP_F25519_rcopy(&B,CURVE_B_ED25519);
#endif

    FP_F25519_one(&one);
    FP_F25519_copy(&t,h);
    sgn=FP_F25519_sign(&t);

#if CURVE_A_ED25519 != 0 || HTC_ISO_ED25519 != 0

        FP_F25519_sqr(&t,&t);
        FP_F25519_imul(&t,&t,RIADZ_F25519);  // Z from hash-to-point draft standard
        FP_F25519_add(&w,&t,&one);     // w=Zt^2+1
        FP_F25519_norm(&w);

        FP_F25519_mul(&w,&w,&t);      // w=Z^2*t^4+Zt^2
        FP_F25519_mul(&D,&A,&w);      // A=Aw
                               
        FP_F25519_add(&w,&w,&one); FP_F25519_norm(&w);
        FP_F25519_mul(&w,&w,&B);
        FP_F25519_neg(&w,&w);          // -B(w+1)
        FP_F25519_norm(&w);

        FP_F25519_copy(&X2,&w);        // Numerators
        FP_F25519_mul(&X3,&t,&X2);

// x^3+Ad^2x+Bd^3
        FP_F25519_sqr(&GX1,&X2);
        FP_F25519_sqr(&D2,&D); FP_F25519_mul(&w,&A,&D2); FP_F25519_add(&GX1,&GX1,&w); FP_F25519_norm(&GX1); FP_F25519_mul(&GX1,&GX1,&X2); FP_F25519_mul(&D2,&D2,&D); FP_F25519_mul(&w,&B,&D2); FP_F25519_add(&GX1,&GX1,&w); FP_F25519_norm(&GX1);

        FP_F25519_mul(&w,&GX1,&D);
        int qr=FP_F25519_qr(&w,&hint);   // qr(ad) - only exp happens here
        FP_F25519_inv(&D,&w,&hint);     // d=1/(ad)
        FP_F25519_mul(&D,&D,&GX1);     // 1/d
        FP_F25519_mul(&X2,&X2,&D);     // X2/=D
        FP_F25519_mul(&X3,&X3,&D);     // X3/=D
        FP_F25519_mul(&t,&t,h);        // t=Z.u^3
        FP_F25519_sqr(&D2,&D);

        FP_F25519_mul(&D,&D2,&t);
        FP_F25519_imul(&t,&w,RIADZ_F25519);
        FP_F25519_rcopy(&X1,CURVE_HTPC_ED25519);     
        FP_F25519_mul(&X1,&X1,&hint); // modify hint

        FP_F25519_cmove(&X2,&X3,1-qr); 
        FP_F25519_cmove(&D2,&D,1-qr);
        FP_F25519_cmove(&w,&t,1-qr);
        FP_F25519_cmove(&hint,&X1,1-qr);

        FP_F25519_sqrt(&Y,&w,&hint);  // first candidate if X2 is correct
        FP_F25519_mul(&Y,&Y,&D2);

/*
        FP_F25519_sqrt(&Y,&w,&hint);  // first candidate if X2 is correct
        FP_F25519_mul(&Y,&Y,&D2);

        FP_F25519_mul(&D2,&D2,&t);     // second candidate if X3 is correct
        FP_F25519_imul(&w,&w,RIADZ_F25519); 

        FP_F25519_rcopy(&X1,CURVE_HTPC_ED25519);     
        FP_F25519_mul(&hint,&hint,&X1); // modify hint

        FP_F25519_sqrt(&Y3,&w,&hint);
        FP_F25519_mul(&Y3,&Y3,&D2);

        FP_F25519_cmove(&X2,&X3,!qr); 
        FP_F25519_cmove(&Y,&Y3,!qr); 
*/
        ne=FP_F25519_sign(&Y)^sgn;
        FP_F25519_neg(&w,&Y); FP_F25519_norm(&w);
        FP_F25519_cmove(&Y,&w,ne);

#if HTC_ISO_ED25519 != 0

// (X2,Y) is on isogenous curve
        k=0;
        isox=iso;
        isoy=3*(iso-1)/2;

// xnum
        FP_F25519_rcopy(&xnum,PC_ED25519[k++]);
        for (i=0;i<isox;i++)
        {
            FP_F25519_mul(&xnum,&xnum,&X2); 
            FP_F25519_rcopy(&w,PC_ED25519[k++]);
            FP_F25519_add(&xnum,&xnum,&w); FP_F25519_norm(&xnum);
        }

// xden
        FP_F25519_copy(&xden,&X2);
        FP_F25519_rcopy(&w,PC_ED25519[k++]);
        FP_F25519_add(&xden,&xden,&w); FP_F25519_norm(&xden);
 
        for (i=0;i<isox-2;i++)
        {
            FP_F25519_mul(&xden,&xden,&X2);
            FP_F25519_rcopy(&w,PC_ED25519[k++]);
            FP_F25519_add(&xden,&xden,&w); FP_F25519_norm(&xden);
        }

// ynum
        FP_F25519_rcopy(&ynum,PC_ED25519[k++]);
        for (i=0;i<isoy;i++)
        {
            FP_F25519_mul(&ynum,&ynum,&X2); 
            FP_F25519_rcopy(&w,PC_ED25519[k++]);
            FP_F25519_add(&ynum,&ynum,&w); FP_F25519_norm(&ynum);
        }

// yden 
        FP_F25519_copy(&yden,&X2);
        FP_F25519_rcopy(&w,PC_ED25519[k++]);
        FP_F25519_add(&yden,&yden,&w); FP_F25519_norm(&yden);
      
        for (i=0;i<isoy-1;i++)
        {
            FP_F25519_mul(&yden,&yden,&X2); 
            FP_F25519_rcopy(&w,PC_ED25519[k++]);
            FP_F25519_add(&yden,&yden,&w); FP_F25519_norm(&yden);
        }

        FP_F25519_mul(&ynum,&ynum,&Y);

// return in projective coordinates
        FP_F25519_mul(&t,&xnum,&yden);
        FP_F25519_copy(&(P->x),&t);

        FP_F25519_mul(&t,&ynum,&xden);
        FP_F25519_copy(&(P->y),&t);

        FP_F25519_mul(&t,&xden,&yden);
        FP_F25519_copy(&(P->z),&t);
        return;
#else

        FP_F25519_redc(x,&X2);
        FP_F25519_redc(y,&Y);
        ECP_ED25519_set(P,x,y);
        return;
#endif
#else 
// SVDW - Shallue and van de Woestijne
        FP_F25519_from_int(&Y,RIADZ_F25519);
        ECP_ED25519_rhs(&A,&Y);  // A=g(Z)
        FP_F25519_rcopy(&B,SQRTm3_F25519);
        FP_F25519_imul(&B,&B,RIADZ_F25519); // Z*sqrt(-3)

        FP_F25519_sqr(&t,&t);
        FP_F25519_mul(&Y,&A,&t);   // tv1=u^2*g(Z)
        FP_F25519_add(&t,&one,&Y); FP_F25519_norm(&t); // tv2=1+tv1
        FP_F25519_sub(&Y,&one,&Y); FP_F25519_norm(&Y); // tv1=1-tv1
        FP_F25519_mul(&D,&t,&Y);
        FP_F25519_mul(&D,&D,&B);

        FP_F25519_copy(&w,&A);
        FP_F25519_tpo(&D,&w);   // tv3=inv0(tv1*tv2*z*sqrt(-3)) and sqrt(g(Z))   // ***

        FP_F25519_mul(&w,&w,&B);  // tv4=Z.sqrt(-3).sqrt(g(Z))
        if (FP_F25519_sign(&w)==1)
        { // depends only on sign of constant RIADZ
            FP_F25519_neg(&w,&w);
            FP_F25519_norm(&w);
        }
        FP_F25519_mul(&w,&w,&B);  // Z.sqrt(-3)
        FP_F25519_mul(&w,&w,h);    // u
        FP_F25519_mul(&w,&w,&Y);   // tv1
        FP_F25519_mul(&w,&w,&D);  // tv3   // tv5=u*tv1*tv3*tv4*Z*sqrt(-3)

        FP_F25519_from_int(&X1,RIADZ_F25519);
        FP_F25519_copy(&X3,&X1);
        FP_F25519_neg(&X1,&X1); FP_F25519_norm(&X1); FP_F25519_div2(&X1,&X1); // -Z/2
        FP_F25519_copy(&X2,&X1);
        FP_F25519_sub(&X1,&X1,&w); FP_F25519_norm(&X1);
        FP_F25519_add(&X2,&X2,&w); FP_F25519_norm(&X2);
        FP_F25519_add(&A,&A,&A);
        FP_F25519_add(&A,&A,&A);
        FP_F25519_norm(&A);      // 4*g(Z)
        FP_F25519_sqr(&t,&t);
        FP_F25519_mul(&t,&t,&D);
        FP_F25519_sqr(&t,&t);    // (tv2^2*tv3)^2
        FP_F25519_mul(&A,&A,&t); // 4*g(Z)*(tv2^2*tv3)^2

        FP_F25519_add(&X3,&X3,&A); FP_F25519_norm(&X3);

        ECP_ED25519_rhs(&w,&X2);
        FP_F25519_cmove(&X3,&X2,FP_F25519_qr(&w,NULL));                           // ***
        ECP_ED25519_rhs(&w,&X1);
        FP_F25519_cmove(&X3,&X1,FP_F25519_qr(&w,NULL));                           // ***
        ECP_ED25519_rhs(&w,&X3);
        FP_F25519_sqrt(&Y,&w,NULL);                                        // ***
        FP_F25519_redc(x,&X3);

        ne=FP_F25519_sign(&Y)^sgn;
        FP_F25519_neg(&w,&Y); FP_F25519_norm(&w);
        FP_F25519_cmove(&Y,&w,ne);

        FP_F25519_redc(y,&Y);
        ECP_ED25519_set(P,x,y);
        return;
#endif

#endif
}


/* Map octet to point */
/*
void ECP_ED25519_mapit(ECP_ED25519 *P, octet *W)
{
    BIG_256_56 q, x;
    DBIG_256_56 dx;
    BIG_256_56_rcopy(q, Modulus_F25519);
    BIG_256_56_dfromBytesLen(dx, W->val,W->len);
    BIG_256_56_dmod(x, dx, q);
    ECP_ED25519_hap2point(P,x);
    ECP_ED25519_cfp(P);
}
*/
/* Convert P to Affine, from (x,y,z) to (x,y) */
/* SU=160 */
void ECP_ED25519_affine(ECP_ED25519 *P)
{
    FP_F25519 one, iz;
    if (ECP_ED25519_isinf(P)) return;
    FP_F25519_one(&one);
    if (FP_F25519_equals(&(P->z), &one)) return;

    FP_F25519_inv(&iz, &(P->z),NULL);
    FP_F25519_mul(&(P->x), &(P->x), &iz);

#if CURVETYPE_ED25519==EDWARDS || CURVETYPE_ED25519==WEIERSTRASS

    FP_F25519_mul(&(P->y), &(P->y), &iz);
    FP_F25519_reduce(&(P->y));

#endif

    FP_F25519_reduce(&(P->x));
    FP_F25519_copy(&(P->z), &one);
}

/* SU=120 */
void ECP_ED25519_outputxyz(ECP_ED25519 *P)
{
    BIG_256_56 x, z;
    if (ECP_ED25519_isinf(P))
    {
        printf("Infinity\n");
        return;
    }
    FP_F25519_reduce(&(P->x));
    FP_F25519_redc(x, &(P->x));
    FP_F25519_reduce(&(P->z));
    FP_F25519_redc(z, &(P->z));

#if CURVETYPE_ED25519!=MONTGOMERY
    BIG_256_56 y;
    FP_F25519_reduce(&(P->y));
    FP_F25519_redc(y, &(P->y));
    printf("(");
    BIG_256_56_output(x);
    printf(",");
    BIG_256_56_output(y);
    printf(",");
    BIG_256_56_output(z);
    printf(")\n");

#else
    printf("(");
    BIG_256_56_output(x);
    printf(",");
    BIG_256_56_output(z);
    printf(")\n");
#endif
}

/* SU=16 */
/* Output point P */
void ECP_ED25519_output(ECP_ED25519 *P)
{
    BIG_256_56 x;
    if (ECP_ED25519_isinf(P))
    {
        printf("Infinity\n");
        return;
    }
    ECP_ED25519_affine(P);
#if CURVETYPE_ED25519!=MONTGOMERY
    BIG_256_56 y;
    FP_F25519_redc(x, &(P->x));
    FP_F25519_redc(y, &(P->y));
    printf("(");
    BIG_256_56_output(x);
    printf(",");
    BIG_256_56_output(y);
    printf(")\n");
#else
    FP_F25519_redc(x, &(P->x));
    printf("(");
    BIG_256_56_output(x);
    printf(")\n");
#endif
}

/* SU=16 */
/* Output point P */
void ECP_ED25519_rawoutput(ECP_ED25519 *P)
{
    BIG_256_56 x, z;

#if CURVETYPE_ED25519!=MONTGOMERY
    BIG_256_56 y;
    FP_F25519_redc(x, &(P->x));
    FP_F25519_redc(y, &(P->y));
    FP_F25519_redc(z, &(P->z));
    printf("(");
    BIG_256_56_output(x);
    printf(",");
    BIG_256_56_output(y);
    printf(",");
    BIG_256_56_output(z);
    printf(")\n");
#else
    FP_F25519_redc(x, &(P->x));
    FP_F25519_redc(z, &(P->z));
    printf("(");
    BIG_256_56_output(x);
    printf(",");
    BIG_256_56_output(z);
    printf(")\n");
#endif
}

/* SU=88 */
/* Convert P to octet string */
void ECP_ED25519_toOctet(octet *W, ECP_ED25519 *P, bool compress)
{
#if CURVETYPE_ED25519==MONTGOMERY
    BIG_256_56 x;
    ECP_ED25519_get(x, P);
    W->len = MODBYTES_256_56;// + 1;
    //W->val[0] = 6;
    BIG_256_56_toBytes(&(W->val[0]), x);
#else
    BIG_256_56 x, y;
    bool alt=false;
    ECP_ED25519_affine(P);
    ECP_ED25519_get(x, y, P);

#if (MBITS-1)%8 <= 4
#ifdef ALLOW_ALT_COMPRESS_ED25519
    alt=true;
#endif
#endif
    if (alt)
    {
        BIG_256_56_toBytes(&(W->val[0]), x);
        if (compress)
        {
            W->len = MODBYTES_256_56;
            W->val[0]|=0x80;
            if (FP_F25519_islarger(&(P->y))==1) W->val[0]|=0x20;
        } else {
            W->len = 2 * MODBYTES_256_56;
            BIG_256_56_toBytes(&(W->val[MODBYTES_256_56]), y);
        }
    } else {
        BIG_256_56_toBytes(&(W->val[1]), x);
        if (compress)
        {
            W->val[0] = 0x02;
            if (FP_F25519_sign(&(P->y)) == 1) W->val[0] = 0x03;
            W->len = MODBYTES_256_56 + 1;
        } else {
            W->val[0] = 0x04;
            W->len = 2 * MODBYTES_256_56 + 1;
            BIG_256_56_toBytes(&(W->val[MODBYTES_256_56 + 1]), y);
        }
    }
#endif
}

/* SU=88 */
/* Restore P from octet string */
int ECP_ED25519_fromOctet(ECP_ED25519 *P, octet *W)
{
#if CURVETYPE_ED25519==MONTGOMERY
    BIG_256_56 x;
    BIG_256_56_fromBytes(x, &(W->val[0]));
    if (ECP_ED25519_set(P, x)) return 1;
    return 0;
#else
    BIG_256_56 x, y;
    bool alt=false;
    int sgn,cmp,typ = W->val[0];

#if (MBITS-1)%8 <= 4
#ifdef ALLOW_ALT_COMPRESS_ED25519
    alt=true;
#endif
#endif

    if (alt)
    {
        W->val[0]&=0x1f;
        BIG_256_56_fromBytes(x, &(W->val[0]));
        W->val[0]=typ;
        if ((typ&0x80)==0)
        {
            BIG_256_56_fromBytes(y, &(W->val[MODBYTES_256_56]));
            if (ECP_ED25519_set(P, x, y)) return 1;
            return 0;
        } else {
            if (!ECP_ED25519_setx(P,x,0)) return 0;
            sgn=(typ&0x20)>>5;
            cmp=FP_F25519_islarger(&(P->y));
            if ((sgn==1 && cmp!=1) || (sgn==0 && cmp==1)) ECP_ED25519_neg(P);
            return 1;
        }

    } else {
        BIG_256_56_fromBytes(x, &(W->val[1]));
        if (typ == 0x04)
        {
            BIG_256_56_fromBytes(y, &(W->val[MODBYTES_256_56 + 1]));
            if (ECP_ED25519_set(P, x, y)) return 1;
        }
        if (typ == 0x02 || typ == 0x03)
        {
            if (ECP_ED25519_setx(P, x, typ & 1)) return 1;
        }
    }
    return 0;
#endif
}


/* Set P=2P */
/* SU=272 */
void ECP_ED25519_dbl(ECP_ED25519 *P)
{
#if CURVETYPE_ED25519==WEIERSTRASS
    FP_F25519 t0, t1, t2, t3, x3, y3, z3, b;

#if CURVE_A_ED25519 == 0
        FP_F25519_sqr(&t0, &(P->y));                   //t0.sqr();
        FP_F25519_mul(&t1, &(P->y), &(P->z));          //t1.mul(z);

        FP_F25519_sqr(&t2, &(P->z));                   //t2.sqr();
        FP_F25519_add(&(P->z), &t0, &t0);      //z.add(t0);
        FP_F25519_norm(&(P->z));                   //z.norm();
        FP_F25519_add(&(P->z), &(P->z), &(P->z));  //z.add(z);
        FP_F25519_add(&(P->z), &(P->z), &(P->z));  //z.add(z);
        FP_F25519_norm(&(P->z));                   //z.norm();

        FP_F25519_imul(&t2, &t2, 3 * CURVE_B_I_ED25519);   //t2.imul(3*ROM.CURVE_B_I);
        FP_F25519_mul(&x3, &t2, &(P->z));          //x3.mul(z);

        FP_F25519_add(&y3, &t0, &t2);              //y3.add(t2);
        FP_F25519_norm(&y3);                       //y3.norm();
        FP_F25519_mul(&(P->z), &(P->z), &t1);      //z.mul(t1);

        FP_F25519_add(&t1, &t2, &t2);              //t1.add(t2);
        FP_F25519_add(&t2, &t2, &t1);              //t2.add(t1);
        FP_F25519_sub(&t0, &t0, &t2);              //t0.sub(t2);
        FP_F25519_norm(&t0);                       //t0.norm();
        FP_F25519_mul(&y3, &y3, &t0);              //y3.mul(t0);
        FP_F25519_add(&y3, &y3, &x3);              //y3.add(x3);
        FP_F25519_mul(&t1, &(P->x), &(P->y));          //t1.mul(y);

        FP_F25519_norm(&t0);                   //x.norm();
        FP_F25519_mul(&(P->x), &t0, &t1);      //x.mul(t1);
        FP_F25519_add(&(P->x), &(P->x), &(P->x));  //x.add(x);
        FP_F25519_norm(&(P->x));                   //x.norm();
        FP_F25519_copy(&(P->y), &y3);              //y.copy(y3);
        FP_F25519_norm(&(P->y));                   //y.norm();
#else

        if (CURVE_B_I_ED25519 == 0)                 //if (ROM.CURVE_B_I==0)
            FP_F25519_rcopy(&b, CURVE_B_ED25519);      //b.copy(new FP(new BIG(ROM.CURVE_B)));

        FP_F25519_sqr(&t0, &(P->x));                   //t0.sqr();  //1    x^2
        FP_F25519_sqr(&t1, &(P->y));                   //t1.sqr();  //2    y^2
        FP_F25519_sqr(&t2, &(P->z));                   //t2.sqr();  //3

        FP_F25519_mul(&t3, &(P->x), &(P->y));          //t3.mul(y); //4
        FP_F25519_add(&t3, &t3, &t3);              //t3.add(t3);
        FP_F25519_norm(&t3);                       //t3.norm();//5

        FP_F25519_mul(&z3, &(P->z), &(P->x));          //z3.mul(x);   //6
        FP_F25519_add(&z3, &z3, &z3);              //z3.add(z3);
        FP_F25519_norm(&z3);                       //z3.norm();//7

        if (CURVE_B_I_ED25519 == 0)                     //if (ROM.CURVE_B_I==0)
            FP_F25519_mul(&y3, &t2, &b);           //y3.mul(b); //8
        else
            FP_F25519_imul(&y3, &t2, CURVE_B_I_ED25519); //y3.imul(ROM.CURVE_B_I);

        FP_F25519_sub(&y3, &y3, &z3);              //y3.sub(z3); //y3.norm(); //9  ***
        FP_F25519_add(&x3, &y3, &y3);              //x3.add(y3);
        FP_F25519_norm(&x3);                       //x3.norm();//10

        FP_F25519_add(&y3, &y3, &x3);              //y3.add(x3); //y3.norm();//11
        FP_F25519_sub(&x3, &t1, &y3);              //x3.sub(y3);
        FP_F25519_norm(&x3);                       //x3.norm();//12
        FP_F25519_add(&y3, &y3, &t1);              //y3.add(t1);
        FP_F25519_norm(&y3);                       //y3.norm();//13
        FP_F25519_mul(&y3, &y3, &x3);              //y3.mul(x3); //14
        FP_F25519_mul(&x3, &x3, &t3);              //x3.mul(t3); //15
        FP_F25519_add(&t3, &t2, &t2);              //t3.add(t2);  //16
        FP_F25519_add(&t2, &t2, &t3);              //t2.add(t3); //17

        if (CURVE_B_I_ED25519 == 0)                 //if (ROM.CURVE_B_I==0)
            FP_F25519_mul(&z3, &z3, &b);           //z3.mul(b); //18
        else
            FP_F25519_imul(&z3, &z3, CURVE_B_I_ED25519); //z3.imul(ROM.CURVE_B_I);

        FP_F25519_sub(&z3, &z3, &t2);              //z3.sub(t2); //z3.norm();//19
        FP_F25519_sub(&z3, &z3, &t0);              //z3.sub(t0);
        FP_F25519_norm(&z3);                       //z3.norm();//20  ***
        FP_F25519_add(&t3, &z3, &z3);              //t3.add(z3); //t3.norm();//21

        FP_F25519_add(&z3, &z3, &t3);              //z3.add(t3);
        FP_F25519_norm(&z3);                       //z3.norm(); //22
        FP_F25519_add(&t3, &t0, &t0);              //t3.add(t0); //t3.norm(); //23
        FP_F25519_add(&t0, &t0, &t3);              //t0.add(t3); //t0.norm();//24
        FP_F25519_sub(&t0, &t0, &t2);              //t0.sub(t2);
        FP_F25519_norm(&t0);                       //t0.norm();//25

        FP_F25519_mul(&t0, &t0, &z3);              //t0.mul(z3);//26
        FP_F25519_add(&y3, &y3, &t0);              //y3.add(t0); //y3.norm();//27
        FP_F25519_mul(&t0, &(P->y), &(P->z));          //t0.mul(z);//28
        FP_F25519_add(&t0, &t0, &t0);              //t0.add(t0);
        FP_F25519_norm(&t0);                       //t0.norm(); //29
        FP_F25519_mul(&z3, &z3, &t0);              //z3.mul(t0);//30
        FP_F25519_sub(&(P->x), &x3, &z3);              //x3.sub(z3); //x3.norm();//31
        FP_F25519_add(&t0, &t0, &t0);              //t0.add(t0);
        FP_F25519_norm(&t0);                       //t0.norm();//32
        FP_F25519_add(&t1, &t1, &t1);              //t1.add(t1);
        FP_F25519_norm(&t1);                       //t1.norm();//33
        FP_F25519_mul(&(P->z), &t0, &t1);              //z3.mul(t1);//34

        FP_F25519_norm(&(P->x));                   //x.norm();
        FP_F25519_copy(&(P->y), &y3);              //y.copy(y3);
        FP_F25519_norm(&(P->y));                   //y.norm();
        FP_F25519_norm(&(P->z));                   //z.norm();
#endif
#endif

#if CURVETYPE_ED25519==EDWARDS
    /* Not using square for multiplication swap, as (1) it needs more adds, and (2) it triggers more reductions */

    FP_F25519 C, D, H, J;

    FP_F25519_sqr(&C, &(P->x));                        //C.sqr();

    FP_F25519_mul(&(P->x), &(P->x), &(P->y));      //x.mul(y);
    FP_F25519_add(&(P->x), &(P->x), &(P->x));      //x.add(x);
    FP_F25519_norm(&(P->x));                       //x.norm();

    FP_F25519_sqr(&D, &(P->y));                        //D.sqr();

#if CURVE_A_ED25519 == -1
        FP_F25519_neg(&C, &C);             //  C.neg();
#endif
    FP_F25519_add(&(P->y), &C, &D);    //y.add(D);
    FP_F25519_norm(&(P->y));               //y.norm();
    FP_F25519_sqr(&H, &(P->z));            //H.sqr();
    FP_F25519_add(&H, &H, &H);             //H.add(H);

    FP_F25519_sub(&J, &(P->y), &H);        //J.sub(H);
    FP_F25519_norm(&J);                    //J.norm();

    FP_F25519_mul(&(P->x), &(P->x), &J);   //x.mul(J);
    FP_F25519_sub(&C, &C, &D);             //C.sub(D);
    FP_F25519_norm(&C);                    //C.norm();
    FP_F25519_mul(&(P->z), &(P->y), &J);   //z.mul(J);
    FP_F25519_mul(&(P->y), &(P->y), &C);   //y.mul(C);


#endif

#if CURVETYPE_ED25519==MONTGOMERY
    FP_F25519 A, B, AA, BB, C;

    FP_F25519_add(&A, &(P->x), &(P->z));       //A.add(z);
    FP_F25519_norm(&A);                    //A.norm();
    FP_F25519_sqr(&AA, &A);            //AA.sqr();
    FP_F25519_sub(&B, &(P->x), &(P->z));       //B.sub(z);
    FP_F25519_norm(&B);                    //B.norm();

    FP_F25519_sqr(&BB, &B);            //BB.sqr();
    FP_F25519_sub(&C, &AA, &BB);           //C.sub(BB);
    FP_F25519_norm(&C);                    //C.norm();
    FP_F25519_mul(&(P->x), &AA, &BB);  //x.mul(BB);
    FP_F25519_imul(&A, &C, (CURVE_A_ED25519 + 2) / 4); //A.imul((ROM.CURVE_A+2)/4);

    FP_F25519_add(&BB, &BB, &A);           //BB.add(A);
    FP_F25519_norm(&BB);                   //BB.norm();
    FP_F25519_mul(&(P->z), &BB, &C);   //z.mul(C);

#endif
}

#if CURVETYPE_ED25519==MONTGOMERY

/* Set P+=Q. W is difference between P and Q and is affine */
void ECP_ED25519_add(ECP_ED25519 *P, ECP_ED25519 *Q, ECP_ED25519 *W)
{
    FP_F25519 A, B, C, D, DA, CB;

    FP_F25519_add(&A, &(P->x), &(P->z)); //A.add(z);
    FP_F25519_sub(&B, &(P->x), &(P->z)); //B.sub(z);

    FP_F25519_add(&C, &(Q->x), &(Q->z)); //C.add(Q.z);

    FP_F25519_sub(&D, &(Q->x), &(Q->z)); //D.sub(Q.z);
    FP_F25519_norm(&A);            //A.norm();

    FP_F25519_norm(&D);            //D.norm();
    FP_F25519_mul(&DA, &D, &A);        //DA.mul(A);

    FP_F25519_norm(&C);            //C.norm();
    FP_F25519_norm(&B);            //B.norm();
    FP_F25519_mul(&CB, &C, &B);    //CB.mul(B);

    FP_F25519_add(&A, &DA, &CB);   //A.add(CB);
    FP_F25519_norm(&A);            //A.norm();
    FP_F25519_sqr(&(P->x), &A);        //A.sqr();
    FP_F25519_sub(&B, &DA, &CB);   //B.sub(CB);
    FP_F25519_norm(&B);            //B.norm();
    FP_F25519_sqr(&B, &B);         //B.sqr();

    FP_F25519_mul(&(P->z), &(W->x), &B); //z.mul(B);

}

#else

/* Set P+=Q */
/* SU=248 */
void ECP_ED25519_add(ECP_ED25519 *P, ECP_ED25519 *Q)
{
#if CURVETYPE_ED25519==WEIERSTRASS

    int b3;
    FP_F25519 t0, t1, t2, t3, t4, x3, y3, z3, b;

#if CURVE_A_ED25519 == 0

        b3 = 3 * CURVE_B_I_ED25519;             //int b=3*ROM.CURVE_B_I;
        FP_F25519_mul(&t0, &(P->x), &(Q->x));      //t0.mul(Q.x);
        FP_F25519_mul(&t1, &(P->y), &(Q->y));      //t1.mul(Q.y);
        FP_F25519_mul(&t2, &(P->z), &(Q->z));      //t2.mul(Q.z);
        FP_F25519_add(&t3, &(P->x), &(P->y));      //t3.add(y);
        FP_F25519_norm(&t3);                   //t3.norm();

        FP_F25519_add(&t4, &(Q->x), &(Q->y));      //t4.add(Q.y);
        FP_F25519_norm(&t4);                   //t4.norm();
        FP_F25519_mul(&t3, &t3, &t4);          //t3.mul(t4);
        FP_F25519_add(&t4, &t0, &t1);          //t4.add(t1);

        FP_F25519_sub(&t3, &t3, &t4);          //t3.sub(t4);
        FP_F25519_norm(&t3);                   //t3.norm();
        FP_F25519_add(&t4, &(P->y), &(P->z));      //t4.add(z);
        FP_F25519_norm(&t4);                   //t4.norm();
        FP_F25519_add(&x3, &(Q->y), &(Q->z));      //x3.add(Q.z);
        FP_F25519_norm(&x3);                   //x3.norm();

        FP_F25519_mul(&t4, &t4, &x3);          //t4.mul(x3);
        FP_F25519_add(&x3, &t1, &t2);          //x3.add(t2);

        FP_F25519_sub(&t4, &t4, &x3);          //t4.sub(x3);
        FP_F25519_norm(&t4);                   //t4.norm();
        FP_F25519_add(&x3, &(P->x), &(P->z));      //x3.add(z);
        FP_F25519_norm(&x3);                   //x3.norm();
        FP_F25519_add(&y3, &(Q->x), &(Q->z));      //y3.add(Q.z);
        FP_F25519_norm(&y3);                   //y3.norm();
        FP_F25519_mul(&x3, &x3, &y3);          //x3.mul(y3);

        FP_F25519_add(&y3, &t0, &t2);          //y3.add(t2);
        FP_F25519_sub(&y3, &x3, &y3);          //y3.rsub(x3);
        FP_F25519_norm(&y3);                   //y3.norm();
        FP_F25519_add(&x3, &t0, &t0);          //x3.add(t0);
        FP_F25519_add(&t0, &t0, &x3);          //t0.add(x3);
        FP_F25519_norm(&t0);                   //t0.norm();
        FP_F25519_imul(&t2, &t2, b3);              //t2.imul(b);

        FP_F25519_add(&z3, &t1, &t2);          //z3.add(t2);
        FP_F25519_norm(&z3);                   //z3.norm();
        FP_F25519_sub(&t1, &t1, &t2);          //t1.sub(t2);
        FP_F25519_norm(&t1);                   //t1.norm();
        FP_F25519_imul(&y3, &y3, b3);              //y3.imul(b);

        FP_F25519_mul(&x3, &y3, &t4);          //x3.mul(t4);
        FP_F25519_mul(&t2, &t3, &t1);          //t2.mul(t1);
        FP_F25519_sub(&(P->x), &t2, &x3);          //x3.rsub(t2);
        FP_F25519_mul(&y3, &y3, &t0);          //y3.mul(t0);
        FP_F25519_mul(&t1, &t1, &z3);          //t1.mul(z3);
        FP_F25519_add(&(P->y), &y3, &t1);          //y3.add(t1);
        FP_F25519_mul(&t0, &t0, &t3);          //t0.mul(t3);
        FP_F25519_mul(&z3, &z3, &t4);          //z3.mul(t4);
        FP_F25519_add(&(P->z), &z3, &t0);          //z3.add(t0);

        FP_F25519_norm(&(P->x));               //x.norm();
        FP_F25519_norm(&(P->y));               //y.norm();
        FP_F25519_norm(&(P->z));               //z.norm();
#else

        if (CURVE_B_I_ED25519 == 0)             //if (ROM.CURVE_B_I==0)
            FP_F25519_rcopy(&b, CURVE_B_ED25519);  //b.copy(new FP(new BIG(ROM.CURVE_B)));

        FP_F25519_mul(&t0, &(P->x), &(Q->x));      //t0.mul(Q.x); //1
        FP_F25519_mul(&t1, &(P->y), &(Q->y));      //t1.mul(Q.y); //2
        FP_F25519_mul(&t2, &(P->z), &(Q->z));      //t2.mul(Q.z); //3

        FP_F25519_add(&t3, &(P->x), &(P->y));      //t3.add(y);
        FP_F25519_norm(&t3);                   //t3.norm(); //4
        FP_F25519_add(&t4, &(Q->x), &(Q->y));      //t4.add(Q.y);
        FP_F25519_norm(&t4);                   //t4.norm();//5
        FP_F25519_mul(&t3, &t3, &t4);          //t3.mul(t4);//6

        FP_F25519_add(&t4, &t0, &t1);          //t4.add(t1); //t4.norm(); //7
        FP_F25519_sub(&t3, &t3, &t4);          //t3.sub(t4);
        FP_F25519_norm(&t3);                   //t3.norm(); //8
        FP_F25519_add(&t4, &(P->y), &(P->z));      //t4.add(z);
        FP_F25519_norm(&t4);                   //t4.norm();//9
        FP_F25519_add(&x3, &(Q->y), &(Q->z));      //x3.add(Q.z);
        FP_F25519_norm(&x3);                   //x3.norm();//10
        FP_F25519_mul(&t4, &t4, &x3);          //t4.mul(x3); //11
        FP_F25519_add(&x3, &t1, &t2);          //x3.add(t2); //x3.norm();//12

        FP_F25519_sub(&t4, &t4, &x3);          //t4.sub(x3);
        FP_F25519_norm(&t4);                   //t4.norm();//13
        FP_F25519_add(&x3, &(P->x), &(P->z));      //x3.add(z);
        FP_F25519_norm(&x3);                   //x3.norm(); //14
        FP_F25519_add(&y3, &(Q->x), &(Q->z));      //y3.add(Q.z);
        FP_F25519_norm(&y3);                   //y3.norm();//15

        FP_F25519_mul(&x3, &x3, &y3);          //x3.mul(y3); //16
        FP_F25519_add(&y3, &t0, &t2);          //y3.add(t2); //y3.norm();//17

        FP_F25519_sub(&y3, &x3, &y3);          //y3.rsub(x3);
        FP_F25519_norm(&y3);                   //y3.norm(); //18

        if (CURVE_B_I_ED25519 == 0)             //if (ROM.CURVE_B_I==0)
            FP_F25519_mul(&z3, &t2, &b);       //z3.mul(b); //18
        else
            FP_F25519_imul(&z3, &t2, CURVE_B_I_ED25519); //z3.imul(ROM.CURVE_B_I);

        FP_F25519_sub(&x3, &y3, &z3);          //x3.sub(z3);
        FP_F25519_norm(&x3);                   //x3.norm(); //20
        FP_F25519_add(&z3, &x3, &x3);          //z3.add(x3); //z3.norm(); //21

        FP_F25519_add(&x3, &x3, &z3);          //x3.add(z3); //x3.norm(); //22
        FP_F25519_sub(&z3, &t1, &x3);          //z3.sub(x3);
        FP_F25519_norm(&z3);                   //z3.norm(); //23
        FP_F25519_add(&x3, &x3, &t1);          //x3.add(t1);
        FP_F25519_norm(&x3);                   //x3.norm(); //24

        if (CURVE_B_I_ED25519 == 0)             //if (ROM.CURVE_B_I==0)
            FP_F25519_mul(&y3, &y3, &b);       //y3.mul(b); //18
        else
            FP_F25519_imul(&y3, &y3, CURVE_B_I_ED25519); //y3.imul(ROM.CURVE_B_I);

        FP_F25519_add(&t1, &t2, &t2);          //t1.add(t2); //t1.norm();//26
        FP_F25519_add(&t2, &t2, &t1);          //t2.add(t1); //t2.norm();//27

        FP_F25519_sub(&y3, &y3, &t2);          //y3.sub(t2); //y3.norm(); //28

        FP_F25519_sub(&y3, &y3, &t0);          //y3.sub(t0);
        FP_F25519_norm(&y3);                   //y3.norm(); //29
        FP_F25519_add(&t1, &y3, &y3);          //t1.add(y3); //t1.norm();//30
        FP_F25519_add(&y3, &y3, &t1);          //y3.add(t1);
        FP_F25519_norm(&y3);                   //y3.norm(); //31

        FP_F25519_add(&t1, &t0, &t0);          //t1.add(t0); //t1.norm(); //32
        FP_F25519_add(&t0, &t0, &t1);          //t0.add(t1); //t0.norm();//33
        FP_F25519_sub(&t0, &t0, &t2);          //t0.sub(t2);
        FP_F25519_norm(&t0);                   //t0.norm();//34

        FP_F25519_mul(&t1, &t4, &y3);          //t1.mul(y3);//35
        FP_F25519_mul(&t2, &t0, &y3);          //t2.mul(y3);//36
        FP_F25519_mul(&y3, &x3, &z3);          //y3.mul(z3);//37
        FP_F25519_add(&(P->y), &y3, &t2);          //y3.add(t2); //y3.norm();//38
        FP_F25519_mul(&x3, &x3, &t3);          //x3.mul(t3);//39
        FP_F25519_sub(&(P->x), &x3, &t1);          //x3.sub(t1);//40
        FP_F25519_mul(&z3, &z3, &t4);          //z3.mul(t4);//41

        FP_F25519_mul(&t1, &t3, &t0);          //t1.mul(t0);//42
        FP_F25519_add(&(P->z), &z3, &t1);          //z3.add(t1);
        FP_F25519_norm(&(P->x));               //x.norm();

        FP_F25519_norm(&(P->y));               //y.norm();
        FP_F25519_norm(&(P->z));               //z.norm();
#endif

#else
    FP_F25519 A, B, C, D, E, F, G, b;

    FP_F25519_mul(&A, &(P->z), &(Q->z));   //A.mul(Q.z);
    FP_F25519_sqr(&B, &A);                 //B.sqr();
    FP_F25519_mul(&C, &(P->x), &(Q->x));   //C.mul(Q.x);
    FP_F25519_mul(&D, &(P->y), &(Q->y));   //D.mul(Q.y);
    FP_F25519_mul(&E, &C, &D);             //E.mul(D);

    if (CURVE_B_I_ED25519 == 0)         //if (ROM.CURVE_B_I==0)
    {
        FP_F25519_rcopy(&b, CURVE_B_ED25519);  //FP b=new FP(new BIG(ROM.CURVE_B));
        FP_F25519_mul(&E, &E, &b);         //E.mul(b);
    }
    else
        FP_F25519_imul(&E, &E, CURVE_B_I_ED25519); //E.imul(ROM.CURVE_B_I);

    FP_F25519_sub(&F, &B, &E);         //F.sub(E);
    FP_F25519_add(&G, &B, &E);         //G.add(E);

#if CURVE_A_ED25519 == 1
        FP_F25519_sub(&E, &D, &C);     //E.sub(C);
#endif
    FP_F25519_add(&C, &C, &D);         //C.add(D);
    FP_F25519_add(&B, &(P->x), &(P->y));   //B.add(y);

    FP_F25519_add(&D, &(Q->x), &(Q->y));   //D.add(Q.y);
    FP_F25519_norm(&B);                //B.norm();
    FP_F25519_norm(&D);                //D.norm();
    FP_F25519_mul(&B, &B, &D);         //B.mul(D);
    FP_F25519_sub(&B, &B, &C);         //B.sub(C);
    FP_F25519_norm(&B);                //B.norm();
    FP_F25519_norm(&F);                //F.norm();
    FP_F25519_mul(&B, &B, &F);         //B.mul(F);
    FP_F25519_mul(&(P->x), &A, &B); //x.mul(B);
    FP_F25519_norm(&G);                //G.norm();

#if CURVE_A_ED25519 == 1
        FP_F25519_norm(&E);            //E.norm();
        FP_F25519_mul(&C, &E, &G);     //C.mul(G);
#endif
#if CURVE_A_ED25519 == -1
        FP_F25519_norm(&C);            //C.norm();
        FP_F25519_mul(&C, &C, &G);     //C.mul(G);
#endif
    FP_F25519_mul(&(P->y), &A, &C); //y.mul(C);
    FP_F25519_mul(&(P->z), &F, &G); //z.mul(G);

#endif
}

/* Set P-=Q */
/* SU=16 */
void  ECP_ED25519_sub(ECP_ED25519 *P, ECP_ED25519 *Q)
{
    ECP_ED25519 NQ;
    ECP_ED25519_copy(&NQ, Q);
    ECP_ED25519_neg(&NQ);
    ECP_ED25519_add(P, &NQ);
}

#endif

#if CURVETYPE_ED25519!=MONTGOMERY
/* constant time multiply by small integer of length bts - use ladder */
void ECP_ED25519_pinmul(ECP_ED25519 *P, int e, int bts)
{
    int i, b;
    ECP_ED25519 R0, R1;

    ECP_ED25519_affine(P);
    ECP_ED25519_inf(&R0);
    ECP_ED25519_copy(&R1, P);

    for (i = bts - 1; i >= 0; i--)
    {
        b = (e >> i) & 1;
        ECP_ED25519_copy(P, &R1);
        ECP_ED25519_add(P, &R0);
        ECP_ED25519_cswap(&R0, &R1, b);
        ECP_ED25519_copy(&R1, P);
        ECP_ED25519_dbl(&R0);
        ECP_ED25519_cswap(&R0, &R1, b);
    }
    ECP_ED25519_copy(P, &R0);
}
#endif

// Point multiplication, multiplies a point P by a scalar e
// This code has no inherent awareness of the order of the curve, or the order of the point.
// The order of the curve will be h.r, where h is a cofactor, and r is a large prime
// Typically P will be of order r (but not always), and typically e will be less than r (but not always)
// A problem can arise if a secret e is a few bits less than r, as the leading zeros in e will leak via a timing attack
// The secret e may however be greater than r (see RFC7748 which combines elimination of a small cofactor h with the point multiplication, using an e>r)
// Our solution is to use as a multiplier an e, whose length in bits is that of the logical OR of e and r, hence allowing e>r while forcing inclusion of leading zeros if e<r. 
// The point multiplication methods used will process leading zeros correctly.

// So this function leaks information about the length of e...
void ECP_ED25519_mul(ECP_ED25519 *P,BIG_256_56 e)
{
    ECP_ED25519_clmul(P,e,e);
}

// .. but this one does not (typically set maxe=r)
// Set P=e*P 
void ECP_ED25519_clmul(ECP_ED25519 *P, BIG_256_56 e, BIG_256_56 maxe)
{
    int max;
    BIG_256_56 cm;
#if CURVETYPE_ED25519==MONTGOMERY
    /* Montgomery ladder */
    int nb, i, b;
    ECP_ED25519 R0, R1, D;
    BIG_256_56_or(cm,e,maxe);
    max=BIG_256_56_nbits(cm);

    if (ECP_ED25519_isinf(P)) return;
    if (BIG_256_56_iszilch(e))
    {
        ECP_ED25519_inf(P);
        return;
    }

    ECP_ED25519_copy(&R0, P);
    ECP_ED25519_copy(&R1, P);
    ECP_ED25519_dbl(&R1);

    ECP_ED25519_copy(&D, P);
    ECP_ED25519_affine(&D);

    nb = max;
    for (i = nb - 2; i >= 0; i--)
    {
        b = BIG_256_56_bit(e, i);
        ECP_ED25519_copy(P, &R1);
        ECP_ED25519_add(P, &R0, &D);
        ECP_ED25519_cswap(&R0, &R1, b);
        ECP_ED25519_copy(&R1, P);
        ECP_ED25519_dbl(&R0);

        ECP_ED25519_cswap(&R0, &R1, b);
    }

    ECP_ED25519_copy(P, &R0);

#else
    /* fixed size windows */
    int i, nb, s, ns;
    BIG_256_56 mt, t;
    ECP_ED25519 Q, W[8], C;
    sign8 w[1 + (NLEN_256_56 * BASEBITS_256_56 + 3) / 4];
    BIG_256_56_or(cm,e,maxe);
    max=BIG_256_56_nbits(cm);

    if (ECP_ED25519_isinf(P)) return;
    if (BIG_256_56_iszilch(e))
    {
        ECP_ED25519_inf(P);
        return;
    }

    /* precompute table */

    ECP_ED25519_copy(&Q, P);
    ECP_ED25519_dbl(&Q);

    ECP_ED25519_copy(&W[0], P);

    for (i = 1; i < 8; i++)
    {
        ECP_ED25519_copy(&W[i], &W[i - 1]);
        ECP_ED25519_add(&W[i], &Q);
    }

//printf("W[1]= ");ECP_output(&W[1]); printf("\n");

    /* make exponent odd - add 2P if even, P if odd */
    BIG_256_56_copy(t, e);
    s = BIG_256_56_parity(t);
    BIG_256_56_inc(t, 1);
    BIG_256_56_norm(t);
    ns = BIG_256_56_parity(t);
    BIG_256_56_copy(mt, t);
    BIG_256_56_inc(mt, 1);
    BIG_256_56_norm(mt);
    BIG_256_56_cmove(t, mt, s);
    ECP_ED25519_cmove(&Q, P, ns);
    ECP_ED25519_copy(&C, &Q);

    nb = 1 + (max + 3) / 4;

    /* convert exponent to signed 4-bit window */
    for (i = 0; i < nb; i++)
    {
        w[i] = BIG_256_56_lastbits(t, 5) - 16;
        BIG_256_56_dec(t, w[i]);
        BIG_256_56_norm(t);
        BIG_256_56_fshr(t, 4);
    }
    w[nb] = BIG_256_56_lastbits(t, 5);

    //ECP_ED25519_copy(P, &W[(w[nb] - 1) / 2]);
    ECP_ED25519_select(P, W, w[i]);
    for (i = nb - 1; i >= 0; i--)
    {
        ECP_ED25519_select(&Q, W, w[i]);
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        ECP_ED25519_add(P, &Q);
    }
    ECP_ED25519_sub(P, &C); /* apply correction */
#endif
}

#if CURVETYPE_ED25519!=MONTGOMERY

// Generic multi-multiplication, fixed 4-bit window, P=Sigma e_i*X_i
void ECP_ED25519_muln(ECP_ED25519 *P,int n,ECP_ED25519 X[],BIG_256_56 e[])
{
    int i,j,k,nb;
    BIG_256_56 t,mt;
    ECP_ED25519 S,R,B[16];
    ECP_ED25519_inf(P);

    BIG_256_56_copy(mt,e[0]);  BIG_256_56_norm(mt);
    for (i=1;i<n;i++)
    { // find biggest
        BIG_256_56_copy(t,e[i]); BIG_256_56_norm(t);
        k=BIG_256_56_comp(t,mt);
        BIG_256_56_cmove(mt,t,(k+1)/2);
    }
    nb=(BIG_256_56_nbits(mt)+3)/4;
    for (int i=nb-1;i>=0;i--)
    { // Pippenger's algorithm
        for (j=0;j<16;j++)
            ECP_ED25519_inf(&B[j]);
        for (j=0;j<n;j++)
        {
            BIG_256_56_copy(mt,e[j]); BIG_256_56_norm(mt);
            BIG_256_56_shr(mt,4*i);
            k=BIG_256_56_lastbits(mt,4);
            ECP_ED25519_add(&B[k],&X[j]);
        }
        ECP_ED25519_inf(&R); ECP_ED25519_inf(&S);
        for (j=15;j>=1;j--)
        {
            ECP_ED25519_add(&R,&B[j]);
            ECP_ED25519_add(&S,&R);
        }
        for (j=0;j<4;j++)
            ECP_ED25519_dbl(P);
        ECP_ED25519_add(P,&S);
    }
}

/* Set P=eP+fQ double multiplication */
/* constant time - as useful for GLV method in pairings */
/* SU=456 */

void ECP_ED25519_mul2(ECP_ED25519 *P, ECP_ED25519 *Q, BIG_256_56 e, BIG_256_56 f)
{
    BIG_256_56 te, tf, mt;
    ECP_ED25519 S, T, W[8], C;
    sign8 w[1 + (NLEN_256_56 * BASEBITS_256_56 + 1) / 2];
    int i, a, b, s, ns, nb;

    BIG_256_56_copy(te, e);
    BIG_256_56_copy(tf, f);

    /* precompute table */
    ECP_ED25519_copy(&W[1], P);
    ECP_ED25519_sub(&W[1], Q); /* P+Q */
    ECP_ED25519_copy(&W[2], P);
    ECP_ED25519_add(&W[2], Q); /* P-Q */
    ECP_ED25519_copy(&S, Q);
    ECP_ED25519_dbl(&S);  /* S=2Q */
    ECP_ED25519_copy(&W[0], &W[1]);
    ECP_ED25519_sub(&W[0], &S);
    ECP_ED25519_copy(&W[3], &W[2]);
    ECP_ED25519_add(&W[3], &S);
    ECP_ED25519_copy(&T, P);
    ECP_ED25519_dbl(&T); /* T=2P */
    ECP_ED25519_copy(&W[5], &W[1]);
    ECP_ED25519_add(&W[5], &T);
    ECP_ED25519_copy(&W[6], &W[2]);
    ECP_ED25519_add(&W[6], &T);
    ECP_ED25519_copy(&W[4], &W[5]);
    ECP_ED25519_sub(&W[4], &S);
    ECP_ED25519_copy(&W[7], &W[6]);
    ECP_ED25519_add(&W[7], &S);

    /* if multiplier is odd, add 2, else add 1 to multiplier, and add 2P or P to correction */

    s = BIG_256_56_parity(te);
    BIG_256_56_inc(te, 1);
    BIG_256_56_norm(te);
    ns = BIG_256_56_parity(te);
    BIG_256_56_copy(mt, te);
    BIG_256_56_inc(mt, 1);
    BIG_256_56_norm(mt);
    BIG_256_56_cmove(te, mt, s);
    ECP_ED25519_cmove(&T, P, ns);
    ECP_ED25519_copy(&C, &T);

    s = BIG_256_56_parity(tf);
    BIG_256_56_inc(tf, 1);
    BIG_256_56_norm(tf);
    ns = BIG_256_56_parity(tf);
    BIG_256_56_copy(mt, tf);
    BIG_256_56_inc(mt, 1);
    BIG_256_56_norm(mt);
    BIG_256_56_cmove(tf, mt, s);
    ECP_ED25519_cmove(&S, Q, ns);
    ECP_ED25519_add(&C, &S);

    BIG_256_56_add(mt, te, tf);
    BIG_256_56_norm(mt);
    nb = 1 + (BIG_256_56_nbits(mt) + 1) / 2;

    /* convert exponent to signed 2-bit window */
    for (i = 0; i < nb; i++)
    {
        a = BIG_256_56_lastbits(te, 3) - 4;
        BIG_256_56_dec(te, a);
        BIG_256_56_norm(te);
        BIG_256_56_fshr(te, 2);
        b = BIG_256_56_lastbits(tf, 3) - 4;
        BIG_256_56_dec(tf, b);
        BIG_256_56_norm(tf);
        BIG_256_56_fshr(tf, 2);
        w[i] = 4 * a + b;
    }
    w[nb] = (4 * BIG_256_56_lastbits(te, 3) + BIG_256_56_lastbits(tf, 3));

    //ECP_ED25519_copy(P, &W[(w[nb] - 1) / 2]);
    ECP_ED25519_select(P, W, w[i]);
    for (i = nb - 1; i >= 0; i--)
    {
        ECP_ED25519_select(&T, W, w[i]);
        ECP_ED25519_dbl(P);
        ECP_ED25519_dbl(P);
        ECP_ED25519_add(P, &T);
    }
    ECP_ED25519_sub(P, &C); /* apply correction */
}

#endif

int ECP_ED25519_generator(ECP_ED25519 *G)
{
    BIG_256_56 x, y;
    BIG_256_56_rcopy(x, CURVE_Gx_ED25519);
#if CURVETYPE_ED25519!=MONTGOMERY
    BIG_256_56_rcopy(y, CURVE_Gy_ED25519);
    return ECP_ED25519_set(G, x, y);
#else
    return ECP_ED25519_set(G, x);
#endif
}

#ifdef HAS_MAIN

int main()
{
    int i;
    ECP_ED25519 G, P;
    csprng RNG;
    BIG_256_56 r, s, x, y, b, m, w, q;
    BIG_256_56_rcopy(x, CURVE_Gx);
#if CURVETYPE_ED25519!=MONTGOMERY
    BIG_256_56_rcopy(y, CURVE_Gy);
#endif
    BIG_256_56_rcopy(m, Modulus_F25519);

    printf("x= ");
    BIG_256_56_output(x);
    printf("\n");
#if CURVETYPE_ED25519!=MONTGOMERY
    printf("y= ");
    BIG_256_56_output(y);
    printf("\n");
#endif
    RNG_seed(&RNG, 3, "abc");

#if CURVETYPE_ED25519!=MONTGOMERY
    ECP_ED25519_set(&G, x, y);
#else
    ECP_ED25519_set(&G, x);
#endif
    if (ECP_ED25519_isinf(&G)) printf("Failed to set - point not on curve\n");
    else printf("set success\n");

    ECP_ED25519_output(&G);

    BIG_256_56_rcopy(r, CURVE_Order); //BIG_dec(r,7);
    printf("r= ");
    BIG_256_56_output(r);
    printf("\n");

    ECP_ED25519_copy(&P, &G);

    ECP_ED25519_mul(&P, r);

    ECP_ED25519_output(&P);
//exit(0);
    BIG_256_56_randomnum(w, &RNG);
    BIG_256_56_mod(w, r);

    ECP_ED25519_copy(&P, &G);
    ECP_ED25519_mul(&P, w);

    ECP_ED25519_output(&P);

    return 0;
}

#endif
