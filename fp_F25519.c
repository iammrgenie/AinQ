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

/* CORE mod p functions */
/* Small Finite Field arithmetic */
/* SU=m, SU is Stack Usage (NOT_SPECIAL Modulus) */

#include "fp_F25519.h"

/* Fast Modular Reduction Methods */

/* r=d mod m */
/* d MUST be normalised */
/* Products must be less than pR in all cases !!! */
/* So when multiplying two numbers, their product *must* be less than MODBITS+BASEBITS*NLEN */
/* Results *may* be one bit bigger than MODBITS */

#if MODTYPE_F25519 == PSEUDO_MERSENNE
/* r=d mod m */

/* Converts from BIG integer to residue form mod Modulus */
void FP_F25519_nres(FP_F25519 *y, BIG_256_56 x)
{
    BIG_256_56 mdls;
    BIG_256_56_rcopy(mdls, Modulus_F25519);
    BIG_256_56_copy(y->g, x);
    BIG_256_56_mod(y->g,mdls);
    y->XES = 1;
}

/* Converts from residue form back to BIG integer form */
void FP_F25519_redc(BIG_256_56 x, FP_F25519 *y)
{
    BIG_256_56_copy(x, y->g);
}

/* reduce a DBIG to a BIG exploiting the special form of the modulus */
void FP_F25519_mod(BIG_256_56 r, DBIG_256_56 d)
{
    BIG_256_56 t, b;
    chunk v, tw;
    BIG_256_56_split(t, b, d, MODBITS_F25519);

    /* Note that all of the excess gets pushed into t. So if squaring a value with a 4-bit excess, this results in
       t getting all 8 bits of the excess product! So products must be less than pR which is Montgomery compatible */

    if (MConst_F25519 < NEXCESS_256_56)
    {
        BIG_256_56_imul(t, t, MConst_F25519);
        BIG_256_56_norm(t);
        BIG_256_56_add(r, t, b);
        BIG_256_56_norm(r);
        tw = r[NLEN_256_56 - 1];
        r[NLEN_256_56 - 1] &= TMASK_F25519;
        r[0] += MConst_F25519 * ((tw >> TBITS_F25519));
    }
    else
    {
        v = BIG_256_56_pmul(t, t, MConst_F25519);
        BIG_256_56_add(r, t, b);
        BIG_256_56_norm(r);
        tw = r[NLEN_256_56 - 1];
        r[NLEN_256_56 - 1] &= TMASK_F25519;
#if CHUNK == 16
        r[1] += muladd_256_56(MConst_F25519, ((tw >> TBITS_F25519) + (v << (BASEBITS_256_56 - TBITS_F25519))), 0, &r[0]);
#else
        r[0] += MConst_F25519 * ((tw >> TBITS_F25519) + (v << (BASEBITS_256_56 - TBITS_F25519)));
#endif
    }
    BIG_256_56_norm(r);
}
#endif

/* This only applies to Curve C448, so specialised (for now) */
#if MODTYPE_F25519 == GENERALISED_MERSENNE

void FP_F25519_nres(FP_F25519 *y, BIG_256_56 x)
{
    BIG_256_56 mdls;
    BIG_256_56_rcopy(mdls, Modulus_F25519);
    BIG_256_56_copy(y->g, x);
    BIG_256_56_mod(y->g,mdls);
    y->XES = 1;
}

/* Converts from residue form back to BIG integer form */
void FP_F25519_redc(BIG_256_56 x, FP_F25519 *y)
{
    BIG_256_56_copy(x, y->g);
}

/* reduce a DBIG to a BIG exploiting the special form of the modulus */
void FP_F25519_mod(BIG_256_56 r, DBIG_256_56 d)
{
    BIG_256_56 t, b;
    chunk carry;
    BIG_256_56_split(t, b, d, MBITS_F25519);

    BIG_256_56_add(r, t, b);

    BIG_256_56_dscopy(d, t);
    BIG_256_56_dshl(d, MBITS_F25519 / 2);

    BIG_256_56_split(t, b, d, MBITS_F25519);

    BIG_256_56_add(r, r, t);
    BIG_256_56_add(r, r, b);
    BIG_256_56_norm(r);
    BIG_256_56_shl(t, MBITS_F25519 / 2);

    BIG_256_56_add(r, r, t);

    carry = r[NLEN_256_56 - 1] >> TBITS_F25519;

    r[NLEN_256_56 - 1] &= TMASK_F25519;
    r[0] += carry;

    r[224 / BASEBITS_256_56] += carry << (224 % BASEBITS_256_56); /* need to check that this falls mid-word */
    BIG_256_56_norm(r);
}

#endif

#if MODTYPE_F25519 == MONTGOMERY_FRIENDLY

/* convert to Montgomery n-residue form */
void FP_F25519_nres(FP_F25519 *y, BIG_256_56 x)
{
    DBIG_256_56 d;
    BIG_256_56 r;
    BIG_256_56_rcopy(r, R2modp_F25519);
    BIG_256_56_mul(d, x, r);
    FP_F25519_mod(y->g, d);
    y->XES = 2;
}

/* convert back to regular form */
void FP_F25519_redc(BIG_256_56 x, FP_F25519 *y)
{
    DBIG_256_56 d;
    BIG_256_56_dzero(d);
    BIG_256_56_dscopy(d, y->g);
    FP_F25519_mod(x, d);
}

/* fast modular reduction from DBIG to BIG exploiting special form of the modulus */
void FP_F25519_mod(BIG_256_56 a, DBIG_256_56 d)
{
    int i;

    for (i = 0; i < NLEN_256_56; i++)
        d[NLEN_256_56 + i] += muladd_256_56(d[i], MConst_F25519 - 1, d[i], &d[NLEN_256_56 + i - 1]);

    BIG_256_56_sducopy(a, d);
    BIG_256_56_norm(a);
}

#endif

#if MODTYPE_F25519 == NOT_SPECIAL

/* convert to Montgomery n-residue form */
void FP_F25519_nres(FP_F25519 *y, BIG_256_56 x)
{
    DBIG_256_56 d;
    BIG_256_56 r;
    BIG_256_56_rcopy(r, R2modp_F25519);
    BIG_256_56_mul(d, x, r);
    FP_F25519_mod(y->g, d);
    y->XES = 2;
}

/* convert back to regular form */
void FP_F25519_redc(BIG_256_56 x, FP_F25519 *y)
{
    DBIG_256_56 d;
    BIG_256_56_dzero(d);
    BIG_256_56_dscopy(d, y->g);
    FP_F25519_mod(x, d);
}


/* reduce a DBIG to a BIG using Montgomery's no trial division method */
/* d is expected to be dnormed before entry */
/* SU= 112 */
void FP_F25519_mod(BIG_256_56 a, DBIG_256_56 d)
{
    BIG_256_56 mdls;
    BIG_256_56_rcopy(mdls, Modulus_F25519);
    BIG_256_56_monty(a, mdls, MConst_F25519, d);
}

#endif

void FP_F25519_from_int(FP_F25519 *x,int a)
{
    BIG_256_56 w;
    if (a<0) BIG_256_56_rcopy(w, Modulus_F25519);
    else BIG_256_56_zero(w); 
    BIG_256_56_inc(w,a); BIG_256_56_norm(w); 
    FP_F25519_nres(x,w);
}

/* test x==0 ? */
/* SU= 48 */
int FP_F25519_iszilch(FP_F25519 *x)
{
    BIG_256_56 m;
    FP_F25519 y;
    FP_F25519_copy(&y,x);
    FP_F25519_reduce(&y);
    FP_F25519_redc(m,&y);
    return BIG_256_56_iszilch(m);
}

int FP_F25519_isunity(FP_F25519 *x)
{
    BIG_256_56 m;
    FP_F25519 y;
    FP_F25519_copy(&y,x);
    FP_F25519_reduce(&y);
    FP_F25519_redc(m,&y);
    return BIG_256_56_isunity(m);
}


void FP_F25519_copy(FP_F25519 *y, FP_F25519 *x)
{
    BIG_256_56_copy(y->g, x->g);
    y->XES = x->XES;
}

void FP_F25519_rcopy(FP_F25519 *y, const BIG_256_56 c)
{
    BIG_256_56 b;
    BIG_256_56_rcopy(b, c);
    FP_F25519_nres(y, b);
}

/* Swap a and b if d=1 */
void FP_F25519_cswap(FP_F25519 *a, FP_F25519 *b, int d)
{
    sign32 t, c = d;
    BIG_256_56_cswap(a->g, b->g, d);

    c = ~(c - 1);
    t = c & ((a->XES) ^ (b->XES));
    a->XES ^= t;
    b->XES ^= t;

}

/* Move b to a if d=1 */
void FP_F25519_cmove(FP_F25519 *a, FP_F25519 *b, int d)
{
    sign32 c = -d;

    BIG_256_56_cmove(a->g, b->g, d);
    a->XES ^= (a->XES ^ b->XES)&c;
}

void FP_F25519_zero(FP_F25519 *x)
{
    BIG_256_56_zero(x->g);
    x->XES = 1;
}

int FP_F25519_equals(FP_F25519 *x, FP_F25519 *y)
{
    FP_F25519 xg, yg;
    FP_F25519_copy(&xg, x);
    FP_F25519_copy(&yg, y);
    FP_F25519_reduce(&xg);
    FP_F25519_reduce(&yg);
    if (BIG_256_56_comp(xg.g, yg.g) == 0) return 1;
    return 0;
}

// Is x lexically larger than p-x?
// return -1 for no, 0 if x=0, 1 for yes
int FP_F25519_islarger(FP_F25519 *x)
{
    BIG_256_56 p,fx,sx;
    if (FP_F25519_iszilch(x)) return 0;
    BIG_256_56_rcopy(p,Modulus_F25519);
    FP_F25519_redc(fx,x);
    BIG_256_56_sub(sx,p,fx);  BIG_256_56_norm(sx); 
    return BIG_256_56_comp(fx,sx);
}

void FP_F25519_toBytes(char *b,FP_F25519 *x)
{
    BIG_256_56 t;
    FP_F25519_redc(t, x);
    BIG_256_56_toBytes(b, t);
}

void FP_F25519_fromBytes(FP_F25519 *x,char *b)
{
    BIG_256_56 t;
    BIG_256_56_fromBytes(t, b);
    FP_F25519_nres(x, t);
}

/* output FP */
/* SU= 48 */
void FP_F25519_output(FP_F25519 *r)
{
    BIG_256_56 c;
    FP_F25519_reduce(r);
    FP_F25519_redc(c, r);
    BIG_256_56_output(c);
}

void FP_F25519_rawoutput(FP_F25519 *r)
{
    BIG_256_56_rawoutput(r->g);
}

#ifdef GET_STATS
int tsqr = 0, rsqr = 0, tmul = 0, rmul = 0;
int tadd = 0, radd = 0, tneg = 0, rneg = 0;
int tdadd = 0, rdadd = 0, tdneg = 0, rdneg = 0;
#endif

#ifdef FUSED_MODMUL

/* Insert fastest code here */

#endif

/* r=a*b mod Modulus */
/* product must be less that p.R - and we need to know this in advance! */
/* SU= 88 */
void FP_F25519_mul(FP_F25519 *r, FP_F25519 *a, FP_F25519 *b)
{
    DBIG_256_56 d;

    if ((sign64)a->XES * b->XES > (sign64)FEXCESS_F25519)
    {
#ifdef DEBUG_REDUCE
        printf("Product too large - reducing it\n");
#endif
        FP_F25519_reduce(a);  /* it is sufficient to fully reduce just one of them < p */
    }

#ifdef FUSED_MODMUL
    FP_F25519_modmul(r->g, a->g, b->g);
#else
    BIG_256_56_mul(d, a->g, b->g);
    FP_F25519_mod(r->g, d);
#endif
    r->XES = 2;
}


/* multiplication by an integer, r=a*c */
/* SU= 136 */
void FP_F25519_imul(FP_F25519 *r, FP_F25519 *a, int c)
{
    int s = 0;

    if (c < 0)
    {
        c = -c;
        s = 1;
    }

#if MODTYPE_F25519==PSEUDO_MERSENNE || MODTYPE_F25519==GENERALISED_MERSENNE
    DBIG_256_56 d;
    BIG_256_56_pxmul(d, a->g, c);
    FP_F25519_mod(r->g, d);
    r->XES = 2;

#else
    //Montgomery
    BIG_256_56 k;
    FP_F25519 f;
    if (a->XES * c <= FEXCESS_F25519)
    {
        BIG_256_56_pmul(r->g, a->g, c);
        r->XES = a->XES * c; // careful here - XES jumps!
    }
    else
    {
        // don't want to do this - only a problem for Montgomery modulus and larger constants
        BIG_256_56_zero(k);
        BIG_256_56_inc(k, c);
        BIG_256_56_norm(k);
        FP_F25519_nres(&f, k);
        FP_F25519_mul(r, a, &f);
    }
#endif

    if (s)
    {
        FP_F25519_neg(r, r);
        FP_F25519_norm(r);
    }
}

/* Set r=a^2 mod m */
/* SU= 88 */
void FP_F25519_sqr(FP_F25519 *r, FP_F25519 *a)
{
    DBIG_256_56 d;

    if ((sign64)a->XES * a->XES > (sign64)FEXCESS_F25519)
    {
#ifdef DEBUG_REDUCE
        printf("Product too large - reducing it\n");
#endif
        FP_F25519_reduce(a);
    }

    BIG_256_56_sqr(d, a->g);
    FP_F25519_mod(r->g, d);
    r->XES = 2;
}

/* SU= 16 */
/* Set r=a+b */
void FP_F25519_add(FP_F25519 *r, FP_F25519 *a, FP_F25519 *b)
{
    BIG_256_56_add(r->g, a->g, b->g);
    r->XES = a->XES + b->XES;
    if (r->XES > FEXCESS_F25519)
    {
#ifdef DEBUG_REDUCE
        printf("Sum too large - reducing it \n");
#endif
        FP_F25519_reduce(r);
    }
}

/* Set r=a-b mod m */
/* SU= 56 */
void FP_F25519_sub(FP_F25519 *r, FP_F25519 *a, FP_F25519 *b)
{
    FP_F25519 n;
    FP_F25519_neg(&n, b);
    FP_F25519_add(r, a, &n);
}

// https://graphics.stanford.edu/~seander/bithacks.html
// constant time log to base 2 (or number of bits in)

static int logb2(unsign32 v)
{
    int r;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    v = v - ((v >> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
    r = (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
    return r;
}

// find appoximation to quotient of a/m
// Out by at most 2.
// Note that MAXXES is bounded to be 2-bits less than half a word
static int quo(BIG_256_56 n, BIG_256_56 m)
{
    int sh;
    chunk num, den;
    int hb = CHUNK / 2;
    if (TBITS_F25519 < hb)
    {
        sh = hb - TBITS_F25519;
        num = (n[NLEN_256_56 - 1] << sh) | (n[NLEN_256_56 - 2] >> (BASEBITS_256_56 - sh));
        den = (m[NLEN_256_56 - 1] << sh) | (m[NLEN_256_56 - 2] >> (BASEBITS_256_56 - sh));
    }
    else
    {
        num = n[NLEN_256_56 - 1];
        den = m[NLEN_256_56 - 1];
    }
    return (int)(num / (den + 1));
}

/* SU= 48 */
/* Fully reduce a mod Modulus */
void FP_F25519_reduce(FP_F25519 *a)
{
    BIG_256_56 m, r;
    int sr, sb, q;
    chunk carry;

    BIG_256_56_rcopy(m, Modulus_F25519);

    BIG_256_56_norm(a->g);

    if (a->XES > 16)
    {
        q = quo(a->g, m);
        carry = BIG_256_56_pmul(r, m, q);
        r[NLEN_256_56 - 1] += (carry << BASEBITS_256_56); // correction - put any carry out back in again
        BIG_256_56_sub(a->g, a->g, r);
        BIG_256_56_norm(a->g);
        sb = 2;
    }
    else sb = logb2(a->XES - 1); // sb does not depend on the actual data

    BIG_256_56_fshl(m, sb);

    while (sb > 0)
    {
// constant time...
        sr = BIG_256_56_ssn(r, a->g, m); // optimized combined shift, subtract and norm
        BIG_256_56_cmove(a->g, r, 1 - sr);
        sb--;
    }

    //BIG_256_56_mod(a->g,m);
    a->XES = 1;
}

void FP_F25519_norm(FP_F25519 *x)
{
    BIG_256_56_norm(x->g);
}

/* Set r=-a mod Modulus */
/* SU= 64 */
void FP_F25519_neg(FP_F25519 *r, FP_F25519 *a)
{
    int sb;
    BIG_256_56 m;

    BIG_256_56_rcopy(m, Modulus_F25519);

    sb = logb2(a->XES - 1);
    BIG_256_56_fshl(m, sb);
    BIG_256_56_sub(r->g, m, a->g);
    r->XES = ((sign32)1 << sb) + 1;

    if (r->XES > FEXCESS_F25519)
    {
#ifdef DEBUG_REDUCE
        printf("Negation too large -  reducing it \n");
#endif
        FP_F25519_reduce(r);
    }

}

/* Set r=a/2. */
/* SU= 56 */
void FP_F25519_div2(FP_F25519 *r, FP_F25519 *a)
{
    BIG_256_56 m;
    BIG_256_56 w;
    BIG_256_56_rcopy(m, Modulus_F25519);
    int pr=BIG_256_56_parity(a->g);

    FP_F25519_copy(r, a);
    BIG_256_56_copy(w,r->g);
    BIG_256_56_fshr(r->g,1);
    BIG_256_56_add(w, w, m);
    BIG_256_56_norm(w);
    BIG_256_56_fshr(w, 1);   
    
    BIG_256_56_cmove(r->g,w,pr);

}

// Could leak size of b
// but not used here with secret exponent b
void FP_F25519_pow(FP_F25519 *r, FP_F25519 *a, BIG_256_56 b)
{
    sign8 w[1 + (NLEN_256_56 * BASEBITS_256_56 + 3) / 4];
    FP_F25519 tb[16];
    BIG_256_56 t;
    int i, nb;

    FP_F25519_copy(r,a);
    FP_F25519_norm(r);
    BIG_256_56_copy(t, b);
    BIG_256_56_norm(t);
    nb = 1 + (BIG_256_56_nbits(t) + 3) / 4;
    /* convert exponent to 4-bit window */
    for (i = 0; i < nb; i++)
    {
        w[i] = BIG_256_56_lastbits(t, 4);
        BIG_256_56_dec(t, w[i]);
        BIG_256_56_norm(t);
        BIG_256_56_fshr(t, 4);
    }

    FP_F25519_one(&tb[0]);
    FP_F25519_copy(&tb[1], r);
    for (i = 2; i < 16; i++)
        FP_F25519_mul(&tb[i], &tb[i - 1], r);

    FP_F25519_copy(r, &tb[w[nb - 1]]);
    for (i = nb - 2; i >= 0; i--)
    {
        FP_F25519_sqr(r, r);
        FP_F25519_sqr(r, r);
        FP_F25519_sqr(r, r);
        FP_F25519_sqr(r, r);
        FP_F25519_mul(r, r, &tb[w[i]]);
    }
    FP_F25519_reduce(r);
}


#if MODTYPE_F25519 == PSEUDO_MERSENNE || MODTYPE_F25519==GENERALISED_MERSENNE

// See eprint paper https://eprint.iacr.org/2018/1038
// If p=3 mod 4 r= x^{(p-3)/4}, if p=5 mod 8 r=x^{(p-5)/8}

static void FP_F25519_fpow(FP_F25519 *r, FP_F25519 *x)
{
    int i, j, k, bw, w, c, nw, lo, m, n, nd, e=PM1D2_F25519;
    FP_F25519 xp[11], t, key;
    const int ac[] = {1, 2, 3, 6, 12, 15, 30, 60, 120, 240, 255};
// phase 1
    FP_F25519_copy(&xp[0], x); // 1
    FP_F25519_sqr(&xp[1], x); // 2
    FP_F25519_mul(&xp[2], &xp[1], x); //3
    FP_F25519_sqr(&xp[3], &xp[2]); // 6
    FP_F25519_sqr(&xp[4], &xp[3]); // 12
    FP_F25519_mul(&xp[5], &xp[4], &xp[2]); // 15
    FP_F25519_sqr(&xp[6], &xp[5]); // 30
    FP_F25519_sqr(&xp[7], &xp[6]); // 60
    FP_F25519_sqr(&xp[8], &xp[7]); // 120
    FP_F25519_sqr(&xp[9], &xp[8]); // 240
    FP_F25519_mul(&xp[10], &xp[9], &xp[5]); // 255

#if MODTYPE_F25519==PSEUDO_MERSENNE
    n = MODBITS_F25519;
#endif
#if MODTYPE_F25519==GENERALISED_MERSENNE  // Goldilocks ONLY
    n = MODBITS_F25519 / 2;
#endif

    n-=(e+1);
    c=(MConst_F25519+(1<<e)+1)/(1<<(e+1));

// need c to be odd
    nd=0;
    while (c%2==0)
    {
        c/=2;
        n-=1;
        nd++;
    }

    bw = 0; w = 1; while (w < c) {w *= 2; bw += 1;}
    k = w - c;

    if (k != 0)
    {
        i = 10; while (ac[i] > k) i--;
        FP_F25519_copy(&key, &xp[i]);
        k -= ac[i];
    }
    while (k != 0)
    {
        i--;
        if (ac[i] > k) continue;
        FP_F25519_mul(&key, &key, &xp[i]);
        k -= ac[i];
    }

// phase 2
    FP_F25519_copy(&xp[1], &xp[2]);
    FP_F25519_copy(&xp[2], &xp[5]);
    FP_F25519_copy(&xp[3], &xp[10]);

    j = 3; m = 8;
    nw = n - bw;
    while (2 * m < nw)
    {
        FP_F25519_copy(&t, &xp[j++]);
        for (i = 0; i < m; i++)
            FP_F25519_sqr(&t, &t);
        FP_F25519_mul(&xp[j], &xp[j - 1], &t);
        m *= 2;
    }

    lo = nw - m;
    FP_F25519_copy(r, &xp[j]);

    while (lo != 0)
    {
        m /= 2; j--;
        if (lo < m) continue;
        lo -= m;
        FP_F25519_copy(&t, r);
        for (i = 0; i < m; i++)
            FP_F25519_sqr(&t, &t);
        FP_F25519_mul(r, &t, &xp[j]);
    }
// phase 3

    if (bw != 0)
    {
        for (i = 0; i < bw; i++ )
            FP_F25519_sqr(r, r);
        FP_F25519_mul(r, r, &key);
    }
#if MODTYPE_F25519==GENERALISED_MERSENNE  // Goldilocks ONLY
    FP_F25519_copy(&key, r);
    FP_F25519_sqr(&t, &key);
    FP_F25519_mul(r, &t, &xp[0]);
    for (i = 0; i < n + 1; i++)
        FP_F25519_sqr(r, r);
    FP_F25519_mul(r, r, &key);
#endif

    for (i=0;i<nd;i++)
        FP_F25519_sqr(r,r);
}

#endif


// calculates r=x^(p-1-2^e)/2^{e+1) where 2^e|p-1
void FP_F25519_progen(FP_F25519 *r,FP_F25519 *x)
{
#if MODTYPE_F25519==PSEUDO_MERSENNE  || MODTYPE_F25519==GENERALISED_MERSENNE
    FP_F25519_fpow(r, x);  
#else
    int e=PM1D2_F25519;
    BIG_256_56 m;
    BIG_256_56_rcopy(m, Modulus_F25519);
    BIG_256_56_dec(m,1);
    BIG_256_56_shr(m,e);
    BIG_256_56_dec(m,1);
    BIG_256_56_fshr(m,1);
    FP_F25519_pow(r,x,m);
#endif
}

/* Is x a QR? return optional hint for fast follow-up square root */
int FP_F25519_qr(FP_F25519 *x,FP_F25519 *h)
{
    FP_F25519 r;
    int i,e=PM1D2_F25519;
    FP_F25519_progen(&r,x);
    if (h!=NULL)
        FP_F25519_copy(h,&r);

    FP_F25519_sqr(&r,&r);
    FP_F25519_mul(&r,x,&r);
    for (i=0;i<e-1;i++ )
        FP_F25519_sqr(&r,&r);


//    for (i=0;i<e;i++)
//        FP_F25519_sqr(&r,&r);
//    FP_F25519_copy(&s,x);
//    for (i=0;i<e-1;i++ )
//        FP_F25519_sqr(&s,&s);
//    FP_F25519_mul(&r,&r,&s);
    
    return FP_F25519_isunity(&r);
}

/* Modular inversion */
void FP_F25519_inv(FP_F25519 *r,FP_F25519 *x,FP_F25519 *h)
{
    int i,e=PM1D2_F25519;
    FP_F25519 s,t;
    FP_F25519_norm(x);
    FP_F25519_copy(&s,x);

    if (h==NULL)
        FP_F25519_progen(&t,x);
    else
        FP_F25519_copy(&t,h);

    for (i=0;i<e-1;i++)
    {  
        FP_F25519_sqr(&s,&s);
        FP_F25519_mul(&s,&s,x);
    }
  
    for (i=0;i<=e;i++)
        FP_F25519_sqr(&t,&t);
    
    FP_F25519_mul(r,&t,&s);
    FP_F25519_reduce(r);
}

// Tonelli-Shanks in constant time
void FP_F25519_sqrt(FP_F25519 *r, FP_F25519 *a, FP_F25519 *h)
{
    int i,j,k,u,e=PM1D2_F25519;
    FP_F25519 v,g,t,b;
    BIG_256_56 m;

    if (h==NULL)
        FP_F25519_progen(&g,a);
    else
        FP_F25519_copy(&g,h);

    BIG_256_56_rcopy(m,ROI_F25519);
    FP_F25519_nres(&v,m);

    FP_F25519_sqr(&t,&g);
    FP_F25519_mul(&t,&t,a);
   
    FP_F25519_mul(r,&g,a);
    FP_F25519_copy(&b,&t);
    for (k=e;k>1;k--)
    {
        for (j=1;j<k-1;j++)
            FP_F25519_sqr(&b,&b);
        u=1-FP_F25519_isunity(&b);
        FP_F25519_mul(&g,r,&v);
        FP_F25519_cmove(r,&g,u);
        FP_F25519_sqr(&v,&v);
        FP_F25519_mul(&g,&t,&v);
        FP_F25519_cmove(&t,&g,u);
        FP_F25519_copy(&b,&t);
    }
// always return +ve square root
    k=FP_F25519_sign(r);
    FP_F25519_neg(&v,r); FP_F25519_norm(&v);
    FP_F25519_cmove(r,&v,k);
}

// Calculate both inverse and square root of x, return QR
int FP_F25519_invsqrt(FP_F25519 *i, FP_F25519 *s, FP_F25519 *x)
{
    FP_F25519 h;
    int qr=FP_F25519_qr(x,&h);
    FP_F25519_sqrt(s,x,&h);
    FP_F25519_inv(i,x,&h);
    return qr;
}

// Two for Price of One - See Hamburg https://eprint.iacr.org/2012/309.pdf
// Calculate inverse of i and square root of s, return QR
int FP_F25519_tpo(FP_F25519 *i, FP_F25519 *s)
{
    int qr;
    FP_F25519 w,t;
    FP_F25519_mul(&w,s,i);
    FP_F25519_mul(&t,&w,i);
    qr=FP_F25519_invsqrt(i,s,&t);
    FP_F25519_mul(i,i,&w);
    FP_F25519_mul(s,s,i);
    return qr;
}

/* SU=8 */
/* set n=1 */
void FP_F25519_one(FP_F25519 *n)
{
    BIG_256_56 b;
    BIG_256_56_one(b);
    FP_F25519_nres(n, b);
}

int FP_F25519_sign(FP_F25519 *x)
{
#ifdef BIG_ENDIAN_SIGN_F25519
    int cp;
    BIG_256_56 m,pm1d2;
    FP_F25519 y;
    BIG_256_56_rcopy(pm1d2, Modulus_F25519);
    BIG_256_56_dec(pm1d2,1);
    BIG_256_56_fshr(pm1d2,1); //(p-1)/2
     
    FP_F25519_copy(&y,x);
    FP_F25519_reduce(&y);
    FP_F25519_redc(m,&y);
    cp=BIG_256_56_comp(m,pm1d2);
    return ((cp+1)&2)>>1;

#else
    BIG_256_56 m;
    FP_F25519 y;
    FP_F25519_copy(&y,x);
    FP_F25519_reduce(&y);
    FP_F25519_redc(m,&y);
    return BIG_256_56_parity(m);
#endif
}

void FP_F25519_rand(FP_F25519 *x,csprng *rng)
{
    BIG_256_56 w,m;
    BIG_256_56_rcopy(m,Modulus_F25519);
    BIG_256_56_randomnum(w,m,rng);
    FP_F25519_nres(x,w);
}


