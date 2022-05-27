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

/* CORE basic functions for BIG type */
/* SU=m, SU is Stack Usage */

#include "big_256_56.h"

/* test a=0? */
int BIG_256_56_iszilch(BIG_256_56 a)
{
    int i;
    chunk d=0;
    for (i = 0; i < NLEN_256_56; i++)
        d|=a[i];
    return (1 & ((d-1)>>BASEBITS_256_56));
}

/* test a=1? */
int BIG_256_56_isunity(BIG_256_56 a)
{
    int i;
    chunk d=0;
    for (i = 1; i < NLEN_256_56; i++)
        d|=a[i];
    return (1 & ((d-1)>>BASEBITS_256_56) & (((a[0]^1)-1)>>BASEBITS_256_56));
}

/* test a=0? */
int BIG_256_56_diszilch(DBIG_256_56 a)
{
    int i;
    chunk d=0;
    for (i = 0; i < DNLEN_256_56; i++)
        d|=a[i];
    return (1 & ((d-1)>>BASEBITS_256_56));
}

/* SU= 56 */
/* output a */
void BIG_256_56_output(BIG_256_56 a)
{
    BIG_256_56 b;
    int i, len;
    len = BIG_256_56_nbits(a);
    if (len % 4 == 0) len /= 4;
    else
    {
        len /= 4;
        len++;
    }
    if (len < MODBYTES_256_56 * 2) len = MODBYTES_256_56 * 2;

    for (i = len - 1; i >= 0; i--)
    {
        BIG_256_56_copy(b, a);
        BIG_256_56_shr(b, i * 4);
        printf("%01x", (unsigned int) b[0] & 15);
    }
}

/* SU= 16 */
void BIG_256_56_rawoutput(BIG_256_56 a)
{
    int i;
    printf("(");
    for (i = 0; i < NLEN_256_56 - 1; i++)
#if CHUNK==64
        printf("%"PRIxMAX",", (uintmax_t) a[i]);
    printf("%"PRIxMAX")", (uintmax_t) a[NLEN_256_56 - 1]);
#else
        printf("%x,", (unsigned int) a[i]);
    printf("%x)", (unsigned int) a[NLEN_256_56 - 1]);
#endif
}

/* Swap a and b if d=1 */
void BIG_256_56_cswap(BIG_256_56 a, BIG_256_56 b, int d)
{
    int i;
    chunk t, c = d;
    c = ~(c - 1);
#ifdef DEBUG_NORM
    for (i = 0; i < NLEN_256_56 + 2; i++)
#else
    for (i = 0; i < NLEN_256_56; i++)
#endif
    {
        t = c & (a[i] ^ b[i]);
        a[i] ^= t;
        b[i] ^= t;
    }
}

/* Move b to a if d=1 */
void BIG_256_56_cmove(BIG_256_56 f, BIG_256_56 g, int d)
{
    int i;
    chunk b = (chunk) - d;
#ifdef DEBUG_NORM
    for (i = 0; i < NLEN_256_56 + 2; i++)
#else
    for (i = 0; i < NLEN_256_56; i++)
#endif
    {
        f[i] ^= (f[i] ^ g[i])&b;
    }
}

/* Move g to f if d=1 */
void BIG_256_56_dcmove(DBIG_256_56 f, DBIG_256_56 g, int d)
{
    int i;
    chunk b = (chunk) - d;
#ifdef DEBUG_NORM
    for (i = 0; i < DNLEN_256_56 + 2; i++)
#else
    for (i = 0; i < DNLEN_256_56; i++)
#endif
    {
        f[i] ^= (f[i] ^ g[i])&b;
    }
}

/* convert BIG to/from bytes */
/* SU= 64 */
void BIG_256_56_toBytes(char *b, BIG_256_56 a)
{
    int i;
    BIG_256_56 c;
    BIG_256_56_copy(c, a);
    BIG_256_56_norm(c);
    for (i = MODBYTES_256_56 - 1; i >= 0; i--)
    {
        b[i] = c[0] & 0xff;
        BIG_256_56_fshr(c, 8);
    }
}

/* SU= 16 */
void BIG_256_56_fromBytes(BIG_256_56 a, char *b)
{
    int i;
    BIG_256_56_zero(a);
    for (i = 0; i < MODBYTES_256_56; i++)
    {
        BIG_256_56_fshl(a, 8);
        a[0] += (int)(unsigned char)b[i];
    }
#ifdef DEBUG_NORM
    a[MPV_256_56] = 1;
    a[MNV_256_56] = 0;
#endif
}

void BIG_256_56_fromBytesLen(BIG_256_56 a, char *b, int s)
{
    int i, len = s;
    BIG_256_56_zero(a);

    if (len > MODBYTES_256_56) len = MODBYTES_256_56;
    for (i = 0; i < len; i++)
    {
        BIG_256_56_fshl(a, 8);
        a[0] += (int)(unsigned char)b[i];
    }
#ifdef DEBUG_NORM
    a[MPV_256_56] = 1;
    a[MNV_256_56] = 0;
#endif
}



/* SU= 88 */
void BIG_256_56_doutput(DBIG_256_56 a)
{
    DBIG_256_56 b;
    int i, len;
    BIG_256_56_dnorm(a);
    len = BIG_256_56_dnbits(a);
    if (len % 4 == 0) len /= 4;
    else
    {
        len /= 4;
        len++;
    }

    for (i = len - 1; i >= 0; i--)
    {
        BIG_256_56_dcopy(b, a);
        BIG_256_56_dshr(b, i * 4);
        printf("%01x", (unsigned int) b[0] & 15);
    }
}


void BIG_256_56_drawoutput(DBIG_256_56 a)
{
    int i;
    printf("(");
    for (i = 0; i < DNLEN_256_56 - 1; i++)
#if CHUNK==64
        printf("%"PRIxMAX",", (uintmax_t) a[i]);
    printf("%"PRIxMAX")", (uintmax_t) a[DNLEN_256_56 - 1]);
#else
        printf("%x,", (unsigned int) a[i]);
    printf("%x)", (unsigned int) a[DNLEN_256_56 - 1]);
#endif
}

/* Copy b=a */
void BIG_256_56_copy(BIG_256_56 b, BIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        b[i] = a[i];
#ifdef DEBUG_NORM
    b[MPV_256_56] = a[MPV_256_56];
    b[MNV_256_56] = a[MNV_256_56];
#endif
}

/* Copy from ROM b=a */
void BIG_256_56_rcopy(BIG_256_56 b, const BIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        b[i] = a[i];
#ifdef DEBUG_NORM
    b[MPV_256_56] = 1;
    b[MNV_256_56] = 0;
#endif
}

/* double length DBIG copy b=a */
void BIG_256_56_dcopy(DBIG_256_56 b, DBIG_256_56 a)
{
    int i;
    for (i = 0; i < DNLEN_256_56; i++)
        b[i] = a[i];
#ifdef DEBUG_NORM
    b[DMPV_256_56] = a[DMPV_256_56];
    b[DMNV_256_56] = a[DMNV_256_56];
#endif
}

/* Copy BIG to bottom half of DBIG */
void BIG_256_56_dscopy(DBIG_256_56 b, BIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56 - 1; i++)
        b[i] = a[i];

    b[NLEN_256_56 - 1] = a[NLEN_256_56 - 1] & BMASK_256_56; /* top word normalized */
    b[NLEN_256_56] = a[NLEN_256_56 - 1] >> BASEBITS_256_56;

    for (i = NLEN_256_56 + 1; i < DNLEN_256_56; i++) b[i] = 0;
#ifdef DEBUG_NORM
    b[DMPV_256_56] = a[MPV_256_56];
    b[DMNV_256_56] = a[MNV_256_56];
#endif
}

/* Copy BIG to top half of DBIG */
void BIG_256_56_dsucopy(DBIG_256_56 b, BIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        b[i] = 0;
    for (i = NLEN_256_56; i < DNLEN_256_56; i++)
        b[i] = a[i - NLEN_256_56];
#ifdef DEBUG_NORM
    b[DMPV_256_56] = a[MPV_256_56];
    b[DMNV_256_56] = a[MNV_256_56];
#endif
}

/* Copy bottom half of DBIG to BIG */
void BIG_256_56_sdcopy(BIG_256_56 b, DBIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        b[i] = a[i];
#ifdef DEBUG_NORM
    b[MPV_256_56] = a[DMPV_256_56];
    b[MNV_256_56] = a[DMNV_256_56];
#endif
}

/* Copy top half of DBIG to BIG */
void BIG_256_56_sducopy(BIG_256_56 b, DBIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        b[i] = a[i + NLEN_256_56];
#ifdef DEBUG_NORM
    b[MPV_256_56] = a[DMPV_256_56];
    b[MNV_256_56] = a[DMNV_256_56];

#endif
}

/* Set a=0 */
void BIG_256_56_zero(BIG_256_56 a)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        a[i] = 0;
#ifdef DEBUG_NORM
    a[MPV_256_56] = a[MNV_256_56] = 0;
#endif
}

void BIG_256_56_dzero(DBIG_256_56 a)
{
    int i;
    for (i = 0; i < DNLEN_256_56; i++)
        a[i] = 0;
#ifdef DEBUG_NORM
    a[DMPV_256_56] = a[DMNV_256_56] = 0;
#endif
}

/* set a=1 */
void BIG_256_56_one(BIG_256_56 a)
{
    int i;
    a[0] = 1;
    for (i = 1; i < NLEN_256_56; i++)
        a[i] = 0;
#ifdef DEBUG_NORM
    a[MPV_256_56] = 1;
    a[MNV_256_56] = 0;
#endif
}



/* Set c=a+b */
/* SU= 8 */
void BIG_256_56_add(BIG_256_56 c, BIG_256_56 a, BIG_256_56 b)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        c[i] = a[i] + b[i];
#ifdef DEBUG_NORM
    c[MPV_256_56] = a[MPV_256_56] + b[MPV_256_56];
    c[MNV_256_56] = a[MNV_256_56] + b[MNV_256_56];
    if (c[MPV_256_56] > NEXCESS_256_56)  printf("add problem - positive digit overflow %d\n", (int)c[MPV_256_56]);
    if (c[MNV_256_56] > NEXCESS_256_56)  printf("add problem - negative digit overflow %d\n", (int)c[MNV_256_56]);

#endif
}

/* Set c=a or b */
void BIG_256_56_or(BIG_256_56 c, BIG_256_56 a, BIG_256_56 b)
{
    int i;
    BIG_256_56_norm(a);
    BIG_256_56_norm(b);
    for (i = 0; i < NLEN_256_56; i++)
        c[i] = a[i] | b[i];
#ifdef DEBUG_NORM
    c[MPV_256_56] = 1;
    c[MNV_256_56] = 0;
#endif
}


/* Set c=c+d */
void BIG_256_56_inc(BIG_256_56 c, int d)
{
    BIG_256_56_norm(c);
    c[0] += (chunk)d;
#ifdef DEBUG_NORM
    c[MPV_256_56] += 1;
#endif
}

/* Set c=a-b */
/* SU= 8 */
void BIG_256_56_sub(BIG_256_56 c, BIG_256_56 a, BIG_256_56 b)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++)
        c[i] = a[i] - b[i];
#ifdef DEBUG_NORM
    c[MPV_256_56] = a[MPV_256_56] + b[MNV_256_56];
    c[MNV_256_56] = a[MNV_256_56] + b[MPV_256_56];
    if (c[MPV_256_56] > NEXCESS_256_56)  printf("sub problem - positive digit overflow %d\n", (int)c[MPV_256_56]);
    if (c[MNV_256_56] > NEXCESS_256_56)  printf("sub problem - negative digit overflow %d\n", (int)c[MNV_256_56]);

#endif
}

/* SU= 8 */

void BIG_256_56_dsub(DBIG_256_56 c, DBIG_256_56 a, DBIG_256_56 b)
{
    int i;
    for (i = 0; i < DNLEN_256_56; i++)
        c[i] = a[i] - b[i];
#ifdef DEBUG_NORM
    c[DMPV_256_56] = a[DMPV_256_56] + b[DMNV_256_56];
    c[DMNV_256_56] = a[DMNV_256_56] + b[DMPV_256_56];
    if (c[DMPV_256_56] > NEXCESS_256_56)  printf("double sub problem - positive digit overflow %d\n", (int)c[DMPV_256_56]);
    if (c[DMNV_256_56] > NEXCESS_256_56)  printf("double sub problem - negative digit overflow %d\n", (int)c[DMNV_256_56]);
#endif
}

void BIG_256_56_dadd(DBIG_256_56 c, DBIG_256_56 a, DBIG_256_56 b)
{
    int i;
    for (i = 0; i < DNLEN_256_56; i++)
        c[i] = a[i] + b[i];
#ifdef DEBUG_NORM
    c[DMPV_256_56] = a[DMPV_256_56] + b[DMNV_256_56];
    c[DMNV_256_56] = a[DMNV_256_56] + b[DMPV_256_56];
    if (c[DMPV_256_56] > NEXCESS_256_56)  printf("double add problem - positive digit overflow %d\n", (int)c[DMPV_256_56]);
    if (c[DMNV_256_56] > NEXCESS_256_56)  printf("double add problem - negative digit overflow %d\n", (int)c[DMNV_256_56]);
#endif
}

/* Set c=c-1 */
void BIG_256_56_dec(BIG_256_56 c, int d)
{
    BIG_256_56_norm(c);
    c[0] -= (chunk)d;
#ifdef DEBUG_NORM
    c[MNV_256_56] += 1;
#endif
}

/* multiplication r=a*c by c<=NEXCESS_256_56 */
void BIG_256_56_imul(BIG_256_56 r, BIG_256_56 a, int c)
{
    int i;
    for (i = 0; i < NLEN_256_56; i++) r[i] = a[i] * c;
#ifdef DEBUG_NORM
    r[MPV_256_56] = a[MPV_256_56] * c;
    r[MNV_256_56] = a[MNV_256_56] * c;
    if (r[MPV_256_56] > NEXCESS_256_56)  printf("int mul problem - positive digit overflow %d\n", (int)r[MPV_256_56]);
    if (r[MNV_256_56] > NEXCESS_256_56)  printf("int mul problem - negative digit overflow %d\n", (int)r[MNV_256_56]);

#endif
}

/* multiplication r=a*c by larger integer - c<=FEXCESS */
/* SU= 24 */
chunk BIG_256_56_pmul(BIG_256_56 r, BIG_256_56 a, int c)
{
    int i;
    chunk ak, carry = 0;
    for (i = 0; i < NLEN_256_56; i++)
    {
        ak = a[i];
        r[i] = 0;
        carry = muladd_256_56(ak, (chunk)c, carry, &r[i]);
    }
#ifdef DEBUG_NORM
    r[MPV_256_56] = 1;
    r[MNV_256_56] = 0;
#endif
    return carry;
}

/* r/=3 */
/* SU= 16 */
int BIG_256_56_div3(BIG_256_56 r)
{
    int i;
    chunk ak, base, carry = 0;
    BIG_256_56_norm(r);
    base = ((chunk)1 << BASEBITS_256_56);
    for (i = NLEN_256_56 - 1; i >= 0; i--)
    {
        ak = (carry * base + r[i]);
        r[i] = ak / 3;
        carry = ak % 3;
    }
    return (int)carry;
}

/* multiplication c=a*b by even larger integer b>FEXCESS, resulting in DBIG */
/* SU= 24 */
void BIG_256_56_pxmul(DBIG_256_56 c, BIG_256_56 a, int b)
{
    int j;
    chunk carry;
    BIG_256_56_dzero(c);
    carry = 0;
    for (j = 0; j < NLEN_256_56; j++)
        carry = muladd_256_56(a[j], (chunk)b, carry, &c[j]);
    c[NLEN_256_56] = carry;
#ifdef DEBUG_NORM
    c[DMPV_256_56] = 1;
    c[DMNV_256_56] = 0;
#endif
}

/* .. if you know the result will fit in a BIG, c must be distinct from a and b */
/* SU= 40 */
void BIG_256_56_smul(BIG_256_56 c, BIG_256_56 a, BIG_256_56 b)
{
    int i, j;
    chunk carry;

    BIG_256_56_zero(c);
    for (i = 0; i < NLEN_256_56; i++)
    {
        carry = 0;
        for (j = 0; j < NLEN_256_56; j++)
        {
            if (i + j < NLEN_256_56)
                carry = muladd_256_56(a[i], b[j], carry, &c[i + j]);
        }
    }
#ifdef DEBUG_NORM
    c[MPV_256_56] = 1;
    c[MNV_256_56] = 0;
#endif

}

/* Set c=a*b */
/* SU= 72 */
//void BIG_256_56_mul(chunk c[restrict DNLEN_256_56],chunk a[restrict NLEN_256_56],chunk b[restrict NLEN_256_56])
void BIG_256_56_mul(DBIG_256_56 c, BIG_256_56 a, BIG_256_56 b)
{
    int i;
#ifdef dchunk
    dchunk t, co;
    dchunk s;
    dchunk d[NLEN_256_56];
    int k;
#endif

#ifdef DEBUG_NORM
    if ((a[MPV_256_56] != 1 && a[MPV_256_56] != 0) || a[MNV_256_56] != 0) printf("First input to mul not normed\n");
    if ((b[MPV_256_56] != 1 && b[MPV_256_56] != 0) || b[MNV_256_56] != 0) printf("Second input to mul not normed\n");
#endif

    /* Faster to Combafy it.. Let the compiler unroll the loops! */

#ifdef COMBA

    /* faster psuedo-Karatsuba method */
#ifdef UNWOUND

#ifdef USE_KARATSUBA

    	d[0]=(dchunk)a[0]*b[0];
	d[1]=(dchunk)a[1]*b[1];
	d[2]=(dchunk)a[2]*b[2];
	d[3]=(dchunk)a[3]*b[3];
	d[4]=(dchunk)a[4]*b[4];

	s=d[0];
	t = s; c[0]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56;
	s+=d[1]; t=co+s +(dchunk)(a[1]-a[0])*(b[0]-b[1]); c[1]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s+=d[2]; t=co+s +(dchunk)(a[2]-a[0])*(b[0]-b[2]); c[2]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s+=d[3]; t=co+s +(dchunk)(a[3]-a[0])*(b[0]-b[3])+(dchunk)(a[2]-a[1])*(b[1]-b[2]); c[3]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s+=d[4]; t=co+s +(dchunk)(a[4]-a[0])*(b[0]-b[4])+(dchunk)(a[3]-a[1])*(b[1]-b[3]); c[4]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 

	s-=d[0]; t=co+s +(dchunk)(a[4]-a[1])*(b[1]-b[4])+(dchunk)(a[3]-a[2])*(b[2]-b[3]); c[5]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s-=d[1]; t=co+s +(dchunk)(a[4]-a[2])*(b[2]-b[4]); c[6]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s-=d[2]; t=co+s +(dchunk)(a[4]-a[3])*(b[3]-b[4]); c[7]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	s-=d[3]; t=co+s ; c[8]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	c[9]=(chunk)co;


#else

    	t=(dchunk)a[0]*b[0]; c[0]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[0]*b[1]+(dchunk)a[1]*b[0]; c[1]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[0]*b[2]+(dchunk)a[1]*b[1]+(dchunk)a[2]*b[0]; c[2]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[0]*b[3]+(dchunk)a[1]*b[2]+(dchunk)a[2]*b[1]+(dchunk)a[3]*b[0]; c[3]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[0]*b[4]+(dchunk)a[1]*b[3]+(dchunk)a[2]*b[2]+(dchunk)a[3]*b[1]+(dchunk)a[4]*b[0]; c[4]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[1]*b[4]+(dchunk)a[2]*b[3]+(dchunk)a[3]*b[2]+(dchunk)a[4]*b[1]; c[5]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[2]*b[4]+(dchunk)a[3]*b[3]+(dchunk)a[4]*b[2]; c[6]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[3]*b[4]+(dchunk)a[4]*b[3]; c[7]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	t=t+(dchunk)a[4]*b[4]; c[8]=(chunk)t & BMASK_256_56; t=t>>BASEBITS_256_56;
	c[9]=(chunk)t;


#endif


#else

#ifndef USE_KARATSUBA

    t=(dchunk)a[0]*b[0];
    c[0]=(chunk)t & BMASK_256_56;
    t = t >> BASEBITS_256_56;
    for (i=1;i<NLEN_256_56;i++)
    {
        k=0; 
        while (k<=i) {t+=(dchunk)a[k]*b[i-k]; k++;}
        c[i]=(chunk)t & BMASK_256_56;
        t = t >> BASEBITS_256_56;
    }

    for (i=NLEN_256_56;i<2*NLEN_256_56-1;i++)
    {
        k=i-(NLEN_256_56-1);
        while (k<=NLEN_256_56-1) {t+=(dchunk)a[k]*b[i-k]; k++;}
        c[i]=(chunk)t & BMASK_256_56;
        t = t >> BASEBITS_256_56;
    }

    c[2 * NLEN_256_56 - 1] = (chunk)t;
#else

    for (i = 0; i < NLEN_256_56; i++)
        d[i] = (dchunk)a[i] * b[i];

    s = d[0];
    t = s;
    c[0] = (chunk)t & BMASK_256_56;
    co = t >> BASEBITS_256_56;

    for (k = 1; k < NLEN_256_56; k++)
    {
        s += d[k];
        t = co + s;
        
        /*for (i = k; i >= 1 + k / 2; i--) This causes a huge slow down! gcc/g++ optimizer problem (I think) */
        for (i=1+k/2;i<=k;i++) t += (dchunk)(a[i] - a[k - i]) * (b[k - i] - b[i]);
        c[k] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
    }
    for (k = NLEN_256_56; k < 2 * NLEN_256_56 - 1; k++)
    {
        s -= d[k - NLEN_256_56];
        t = co + s;
        for (i=1+k/2;i<NLEN_256_56;i++) t += (dchunk)(a[i] - a[k - i]) * (b[k - i] - b[i]);
        c[k] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
    }
    c[2 * NLEN_256_56 - 1] = (chunk)co;
#endif
#endif

#else
    int j;
    chunk carry;
    BIG_256_56_dzero(c);
    for (i = 0; i < NLEN_256_56; i++)
    {
        carry = 0;
        for (j = 0; j < NLEN_256_56; j++)
            carry = muladd_256_56(a[i], b[j], carry, &c[i + j]);

        c[NLEN_256_56 + i] = carry;
    }

#endif

#ifdef DEBUG_NORM
    c[DMPV_256_56] = 1;
    c[DMNV_256_56] = 0;
#endif
}

/* Set c=a*a */
/* SU= 80 */

//void BIG_256_56_sqr(chunk c[restrict DNLEN_256_56],chunk a[restrict NLEN_256_56])
void BIG_256_56_sqr(DBIG_256_56 c, BIG_256_56 a)
{
    int i, j;
#ifdef dchunk
    dchunk t, co;
#endif

#ifdef DEBUG_NORM
    if ((a[MPV_256_56] != 1 && a[MPV_256_56] != 0) || a[MNV_256_56] != 0) printf("Input to sqr not normed\n");
#endif
    /* Note 2*a[i] in loop below and extra addition */

#ifdef COMBA

#ifdef UNWOUND

    
	t=(dchunk)a[0]*a[0]; c[0]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56;
	t= +(dchunk)a[1]*a[0]; t+=t; t+=co; c[1]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t= +(dchunk)a[2]*a[0]; t+=t; t+=co; t+=(dchunk)a[1]*a[1]; c[2]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t= +(dchunk)a[3]*a[0]+(dchunk)a[2]*a[1]; t+=t; t+=co; c[3]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t= +(dchunk)a[4]*a[0]+(dchunk)a[3]*a[1]; t+=t; t+=co; t+=(dchunk)a[2]*a[2]; c[4]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 

	t= +(dchunk)a[4]*a[1]+(dchunk)a[3]*a[2]; t+=t; t+=co; c[5]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t= +(dchunk)a[4]*a[2]; t+=t; t+=co; t+=(dchunk)a[3]*a[3]; c[6]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t= +(dchunk)a[4]*a[3]; t+=t; t+=co; c[7]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
	t=co; t+=(dchunk)a[4]*a[4]; c[8]=(chunk)t&BMASK_256_56; co=t>>BASEBITS_256_56; 
 	c[9]=(chunk)co;


#else


    t = (dchunk)a[0] * a[0];
    c[0] = (chunk)t & BMASK_256_56;
    co = t >> BASEBITS_256_56;

    for (j = 1; j < NLEN_256_56 - 1; )
    {
        t = (dchunk)a[j] * a[0];
        for (i = 1; i < (j + 1) / 2; i++)
        {
            t += (dchunk)a[j - i] * a[i];
        }
        t += t;
        t += co;
        c[j] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
        j++;
        t = (dchunk)a[j] * a[0];
        for (i = 1; i < (j + 1) / 2; i++)
        {
            t += (dchunk)a[j - i] * a[i];
        }
        t += t;
        t += co;
        t += (dchunk)a[j / 2] * a[j / 2];
        c[j] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
        j++;
    }

    for (j = NLEN_256_56 - 1 + NLEN_256_56 % 2; j < DNLEN_256_56 - 3; )
    {
        t = (dchunk)a[NLEN_256_56 - 1] * a[j - NLEN_256_56 + 1];
        for (i = j - NLEN_256_56 + 2; i < (j + 1) / 2; i++)
        {
            t += (dchunk)a[j - i] * a[i];
        }
        t += t;
        t += co;
        c[j] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
        j++;
        t = (dchunk)a[NLEN_256_56 - 1] * a[j - NLEN_256_56 + 1];
        for (i = j - NLEN_256_56 + 2; i < (j + 1) / 2; i++)
        {
            t += (dchunk)a[j - i] * a[i];
        }
        t += t;
        t += co;
        t += (dchunk)a[j / 2] * a[j / 2];
        c[j] = (chunk)t & BMASK_256_56;
        co = t >> BASEBITS_256_56;
        j++;
    }

    t = (dchunk)a[NLEN_256_56 - 2] * a[NLEN_256_56 - 1];
    t += t;
    t += co;
    c[DNLEN_256_56 - 3] = (chunk)t & BMASK_256_56;
    co = t >> BASEBITS_256_56;

    t = (dchunk)a[NLEN_256_56 - 1] * a[NLEN_256_56 - 1] + co;
    c[DNLEN_256_56 - 2] = (chunk)t & BMASK_256_56;
    co = t >> BASEBITS_256_56;
    c[DNLEN_256_56 - 1] = (chunk)co;


#endif

#else
    chunk carry;
    BIG_256_56_dzero(c);
    for (i = 0; i < NLEN_256_56; i++)
    {
        carry = 0;
        for (j = i + 1; j < NLEN_256_56; j++)
            carry = muladd_256_56(a[i], a[j], carry, &c[i + j]);
        c[NLEN_256_56 + i] = carry;
    }

    for (i = 0; i < DNLEN_256_56; i++) c[i] *= 2;

    for (i = 0; i < NLEN_256_56; i++)
        c[2 * i + 1] += muladd_256_56(a[i], a[i], 0, &c[2 * i]);

    BIG_256_56_dnorm(c);
#endif


#ifdef DEBUG_NORM
    c[DMPV_256_56] = 1;
    c[DMNV_256_56] = 0;
#endif

}

/* Montgomery reduction */
//void BIG_256_56_monty(chunk a[restrict NLEN_256_56], chunk md[restrict NLEN_256_56], chunk MC, chunk d[restrict DNLEN_256_56])
void BIG_256_56_monty(BIG_256_56 a, BIG_256_56 md, chunk MC, DBIG_256_56 d)
{
    int i, k;

#ifdef dchunk
    dchunk t, c, s;
    dchunk dd[NLEN_256_56];
    chunk v[NLEN_256_56];
#endif

#ifdef DEBUG_NORM
    if ((d[DMPV_256_56] != 1 && d[DMPV_256_56] != 0) || d[DMNV_256_56] != 0) printf("Input to redc not normed\n");
#endif

#ifdef COMBA

#ifdef UNWOUND

#ifdef USE_KARATSUBA

    	t=d[0]; v[0]=((chunk)t*MC)&BMASK_256_56; t+=(dchunk)v[0]*md[0];  s=0; c=(t>>BASEBITS_256_56);

	t=d[1]+c+s+(dchunk)v[0]*md[1]; v[1]=((chunk)t*MC)&BMASK_256_56; t+=(dchunk)v[1]*md[0];  dd[1]=(dchunk)v[1]*md[1]; s+=dd[1]; c=(t>>BASEBITS_256_56); 
	t=d[2]+c+s+(dchunk)v[0]*md[2]; v[2]=((chunk)t*MC)&BMASK_256_56; t+=(dchunk)v[2]*md[0];  dd[2]=(dchunk)v[2]*md[2]; s+=dd[2]; c=(t>>BASEBITS_256_56); 
	t=d[3]+c+s+(dchunk)v[0]*md[3]+(dchunk)(v[1]-v[2])*(md[2]-md[1]); v[3]=((chunk)t*MC)&BMASK_256_56; t+=(dchunk)v[3]*md[0];  dd[3]=(dchunk)v[3]*md[3]; s+=dd[3]; c=(t>>BASEBITS_256_56); 
	t=d[4]+c+s+(dchunk)v[0]*md[4]+(dchunk)(v[1]-v[3])*(md[3]-md[1]); v[4]=((chunk)t*MC)&BMASK_256_56; t+=(dchunk)v[4]*md[0];  dd[4]=(dchunk)v[4]*md[4]; s+=dd[4]; c=(t>>BASEBITS_256_56); 

	t=d[5]+c+s+(dchunk)(v[1]-v[4])*(md[4]-md[1])+(dchunk)(v[2]-v[3])*(md[3]-md[2]); a[0]=(chunk)t&BMASK_256_56;  s-=dd[1]; c=(t>>BASEBITS_256_56); 
	t=d[6]+c+s+(dchunk)(v[2]-v[4])*(md[4]-md[2]); a[1]=(chunk)t&BMASK_256_56;  s-=dd[2]; c=(t>>BASEBITS_256_56); 
	t=d[7]+c+s+(dchunk)(v[3]-v[4])*(md[4]-md[3]); a[2]=(chunk)t&BMASK_256_56;  s-=dd[3]; c=(t>>BASEBITS_256_56); 
	t=d[8]+c+s; a[3]=(chunk)t&BMASK_256_56;  s-=dd[4]; c=(t>>BASEBITS_256_56); 
	a[4]=d[9]+((chunk)c&BMASK_256_56);


#else

    	t = d[0];
	v[0] = ((chunk)t * MC)&BMASK_256_56;
	t += (dchunk)v[0] * md[0];
	t = (t >> BASEBITS_256_56) + d[1];
	t += (dchunk)v[0] * md[1] ; v[1] = ((chunk)t * MC)&BMASK_256_56; t += (dchunk)v[1] * md[0]; t = (t >> BASEBITS_256_56) + d[2];
	t += (dchunk)v[0] * md[2] + (dchunk)v[1]*md[1]; v[2] = ((chunk)t * MC)&BMASK_256_56; t += (dchunk)v[2] * md[0]; t = (t >> BASEBITS_256_56) + d[3];
	t += (dchunk)v[0] * md[3] + (dchunk)v[1]*md[2]+ (dchunk)v[2]*md[1]; v[3] = ((chunk)t * MC)&BMASK_256_56; t += (dchunk)v[3] * md[0]; t = (t >> BASEBITS_256_56) + d[4];
	t += (dchunk)v[0] * md[4] + (dchunk)v[1]*md[3]+ (dchunk)v[2]*md[2]+ (dchunk)v[3]*md[1]; v[4] = ((chunk)t * MC)&BMASK_256_56; t += (dchunk)v[4] * md[0]; t = (t >> BASEBITS_256_56) + d[5];
	t=t + (dchunk)v[1]*md[4] + (dchunk)v[2]*md[3] + (dchunk)v[3]*md[2] + (dchunk)v[4]*md[1] ; a[0] = (chunk)t & BMASK_256_56; t = (t >> BASEBITS_256_56) + d[6];
	t=t + (dchunk)v[2]*md[4] + (dchunk)v[3]*md[3] + (dchunk)v[4]*md[2] ; a[1] = (chunk)t & BMASK_256_56; t = (t >> BASEBITS_256_56) + d[7];
	t=t + (dchunk)v[3]*md[4] + (dchunk)v[4]*md[3] ; a[2] = (chunk)t & BMASK_256_56; t = (t >> BASEBITS_256_56) + d[8];
	t=t + (dchunk)v[4]*md[4] ; a[3] = (chunk)t & BMASK_256_56; t = (t >> BASEBITS_256_56) + d[9];
	a[4] = (chunk)t & BMASK_256_56;


#endif


#else
#ifndef USE_KARATSUBA 
    t = d[0];
    v[0] = ((chunk)t * MC)&BMASK_256_56;
    t += (dchunk)v[0] * md[0];
    t = (t >> BASEBITS_256_56) + d[1];
   
    for (i = 1; i < NLEN_256_56; i++)
    {
        k=1;
        t += (dchunk)v[0] * md[i];
        while (k<i) {t += (dchunk)v[k]*md[i-k]; k++;}
        v[i] = ((chunk)t * MC)&BMASK_256_56;
        t += (dchunk)v[i] * md[0];
        t = (t >> BASEBITS_256_56) + d[i + 1];
    }
    for (i = NLEN_256_56; i < 2 * NLEN_256_56 - 1; i++)
    {
        k=i-(NLEN_256_56-1);
        while (k<=NLEN_256_56-1) {t += (dchunk)v[k]*md[i-k]; k++;}
        a[i - NLEN_256_56] = (chunk)t & BMASK_256_56;
        t = (t >> BASEBITS_256_56) + d[i + 1];
    }
    a[NLEN_256_56 - 1] = (chunk)t & BMASK_256_56;
#else
    t = d[0];
    v[0] = ((chunk)t * MC)&BMASK_256_56;
    t += (dchunk)v[0] * md[0];
    c = (t >> BASEBITS_256_56) + d[1];
    s = 0;

    for (k = 1; k < NLEN_256_56; k++)
    {
        t = c + s + (dchunk)v[0] * md[k];
        for (i=1+k/2;i<k;i++) t += (dchunk)(v[k - i] - v[i]) * (md[i] - md[k - i]);
        v[k] = ((chunk)t * MC)&BMASK_256_56;
        t += (dchunk)v[k] * md[0];
        c = (t >> BASEBITS_256_56) + d[k + 1];
        dd[k] = (dchunk)v[k] * md[k];
        s += dd[k];
    }
    for (k = NLEN_256_56; k < 2 * NLEN_256_56 - 1; k++)
    {
        t = c + s;
        for (i=1+k/2;i<NLEN_256_56;i++) t += (dchunk)(v[k - i] - v[i]) * (md[i] - md[k - i]);
        a[k - NLEN_256_56] = (chunk)t & BMASK_256_56;
        c = (t >> BASEBITS_256_56) + d[k + 1];
        s -= dd[k - NLEN_256_56 + 1];
    }
    a[NLEN_256_56 - 1] = (chunk)c & BMASK_256_56;

#endif
#endif


#else
    int j;
    chunk m, carry;
    for (i = 0; i < NLEN_256_56; i++)
    {
        if (MC == -1) m = (-d[i])&BMASK_256_56;
        else
        {
            if (MC == 1) m = d[i];
            else m = (MC * d[i])&BMASK_256_56;
        }
        carry = 0;
        for (j = 0; j < NLEN_256_56; j++)
            carry = muladd_256_56(m, md[j], carry, &d[i + j]);
        d[NLEN_256_56 + i] += carry;
    }
    BIG_256_56_sducopy(a, d);
    BIG_256_56_norm(a);

#endif

#ifdef DEBUG_NORM
    a[MPV_256_56] = 1;
    a[MNV_256_56] = 0;
#endif
}

/* General shift left of a by n bits */
/* a MUST be normalised */
/* SU= 32 */
void BIG_256_56_shl(BIG_256_56 a, int k)
{
    int i;
    int n = k % BASEBITS_256_56;
    int m = k / BASEBITS_256_56;

    a[NLEN_256_56 - 1] = ((a[NLEN_256_56 - 1 - m] << n));
    if (NLEN_256_56 >= m + 2) a[NLEN_256_56 - 1] |= (a[NLEN_256_56 - m - 2] >> (BASEBITS_256_56 - n));

    for (i = NLEN_256_56 - 2; i > m; i--)
        a[i] = ((a[i - m] << n)&BMASK_256_56) | (a[i - m - 1] >> (BASEBITS_256_56 - n));
    a[m] = (a[0] << n)&BMASK_256_56;
    for (i = 0; i < m; i++) a[i] = 0;

}

/* Fast shift left of a by n bits, where n less than a word, Return excess (but store it as well) */
/* a MUST be normalised */
/* SU= 16 */
int BIG_256_56_fshl(BIG_256_56 a, int n)
{
    int i;

    a[NLEN_256_56 - 1] = ((a[NLEN_256_56 - 1] << n)) | (a[NLEN_256_56 - 2] >> (BASEBITS_256_56 - n)); /* top word not masked */
    for (i = NLEN_256_56 - 2; i > 0; i--)
        a[i] = ((a[i] << n)&BMASK_256_56) | (a[i - 1] >> (BASEBITS_256_56 - n));
    a[0] = (a[0] << n)&BMASK_256_56;

    return (int)(a[NLEN_256_56 - 1] >> ((8 * MODBYTES_256_56) % BASEBITS_256_56)); /* return excess - only used in ff.c */
}

/* double length left shift of a by k bits - k can be > BASEBITS , a MUST be normalised */
/* SU= 32 */
void BIG_256_56_dshl(DBIG_256_56 a, int k)
{
    int i;
    int n = k % BASEBITS_256_56;
    int m = k / BASEBITS_256_56;

    a[DNLEN_256_56 - 1] = ((a[DNLEN_256_56 - 1 - m] << n)) | (a[DNLEN_256_56 - m - 2] >> (BASEBITS_256_56 - n));

    for (i = DNLEN_256_56 - 2; i > m; i--)
        a[i] = ((a[i - m] << n)&BMASK_256_56) | (a[i - m - 1] >> (BASEBITS_256_56 - n));
    a[m] = (a[0] << n)&BMASK_256_56;
    for (i = 0; i < m; i++) a[i] = 0;

}

/* General shift right of a by k bits */
/* a MUST be normalised */
/* SU= 32 */
void BIG_256_56_shr(BIG_256_56 a, int k)
{
    int i;
    int n = k % BASEBITS_256_56;
    int m = k / BASEBITS_256_56;
    for (i = 0; i < NLEN_256_56 - m - 1; i++)
        a[i] = (a[m + i] >> n) | ((a[m + i + 1] << (BASEBITS_256_56 - n))&BMASK_256_56);
    if (NLEN_256_56 > m)  a[NLEN_256_56 - m - 1] = a[NLEN_256_56 - 1] >> n;
    for (i = NLEN_256_56 - m; i < NLEN_256_56; i++) a[i] = 0;

}

/* Fast combined shift, subtract and norm. Return sign of result */
int BIG_256_56_ssn(BIG_256_56 r, BIG_256_56 a, BIG_256_56 m)
{
    int i, n = NLEN_256_56 - 1;
    chunk carry;
    m[0] = (m[0] >> 1) | ((m[1] << (BASEBITS_256_56 - 1))&BMASK_256_56);
    r[0] = a[0] - m[0];
    carry = r[0] >> BASEBITS_256_56;
    r[0] &= BMASK_256_56;

    for (i = 1; i < n; i++)
    {
        m[i] = (m[i] >> 1) | ((m[i + 1] << (BASEBITS_256_56 - 1))&BMASK_256_56);
        r[i] = a[i] - m[i] + carry;
        carry = r[i] >> BASEBITS_256_56;
        r[i] &= BMASK_256_56;
    }

    m[n] >>= 1;
    r[n] = a[n] - m[n] + carry;
#ifdef DEBUG_NORM
    r[MPV_256_56] = 1;
    r[MNV_256_56] = 0;
#endif
    return ((r[n] >> (CHUNK - 1)) & 1);
}

/* Faster shift right of a by k bits. Return shifted out part */
/* a MUST be normalised */
/* SU= 16 */
int BIG_256_56_fshr(BIG_256_56 a, int k)
{
    int i;
    chunk r = a[0] & (((chunk)1 << k) - 1); /* shifted out part */
    for (i = 0; i < NLEN_256_56 - 1; i++)
        a[i] = (a[i] >> k) | ((a[i + 1] << (BASEBITS_256_56 - k))&BMASK_256_56);
    a[NLEN_256_56 - 1] = a[NLEN_256_56 - 1] >> k;
    return (int)r;
}

/* double length right shift of a by k bits - can be > BASEBITS */
/* SU= 32 */
void BIG_256_56_dshr(DBIG_256_56 a, int k)
{
    int i;
    int n = k % BASEBITS_256_56;
    int m = k / BASEBITS_256_56;
    for (i = 0; i < DNLEN_256_56 - m - 1; i++)
        a[i] = (a[m + i] >> n) | ((a[m + i + 1] << (BASEBITS_256_56 - n))&BMASK_256_56);
    a[DNLEN_256_56 - m - 1] = a[DNLEN_256_56 - 1] >> n;
    for (i = DNLEN_256_56 - m; i < DNLEN_256_56; i++ ) a[i] = 0;
}

/* Split DBIG d into two BIGs t|b. Split happens at n bits, where n falls into NLEN word */
/* d MUST be normalised */
/* SU= 24 */
chunk BIG_256_56_split(BIG_256_56 t, BIG_256_56 b, DBIG_256_56 d, int n)
{
    int i;
    chunk nw, carry = 0;
    int m = n % BASEBITS_256_56;

    if (m == 0)
    {
        for (i = 0; i < NLEN_256_56; i++) b[i] = d[i];
        if (t != b)
        {
            for (i = NLEN_256_56; i < 2 * NLEN_256_56; i++) t[i - NLEN_256_56] = d[i];
            carry = t[NLEN_256_56 - 1] >> BASEBITS_256_56;
            t[NLEN_256_56 - 1] = t[NLEN_256_56 - 1] & BMASK_256_56; /* top word normalized */
        }
        return carry;
    }

    for (i = 0; i < NLEN_256_56 - 1; i++) b[i] = d[i];

    b[NLEN_256_56 - 1] = d[NLEN_256_56 - 1] & (((chunk)1 << m) - 1);

    if (t != b)
    {
        carry = (d[DNLEN_256_56 - 1] << (BASEBITS_256_56 - m));
        for (i = DNLEN_256_56 - 2; i >= NLEN_256_56 - 1; i--)
        {
            nw = (d[i] >> m) | carry;
            carry = (d[i] << (BASEBITS_256_56 - m))&BMASK_256_56;
            t[i - NLEN_256_56 + 1] = nw;
        }
    }
#ifdef DEBUG_NORM
    t[MPV_256_56] = 1;
    t[MNV_256_56] = 0;
    b[MPV_256_56] = 1;
    b[MNV_256_56] = 0;
#endif
    return carry;
}

/* you gotta keep the sign of carry! Look - no branching! */
/* Note that sign bit is needed to disambiguate between +ve and -ve values */
/* normalise BIG - force all digits < 2^BASEBITS */
chunk BIG_256_56_norm(BIG_256_56 a)
{
    int i;
    chunk d, carry;

    carry=a[0]>>BASEBITS_256_56;
    a[0]&=BMASK_256_56;

    for (i = 1; i < NLEN_256_56 - 1; i++)
    {
        d = a[i] + carry;
        a[i] = d & BMASK_256_56;
        carry = d >> BASEBITS_256_56;
    }
    a[NLEN_256_56 - 1] = (a[NLEN_256_56 - 1] + carry);

#ifdef DEBUG_NORM
    a[MPV_256_56] = 1;
    a[MNV_256_56] = 0;
#endif
    return (a[NLEN_256_56 - 1] >> ((8 * MODBYTES_256_56) % BASEBITS_256_56)); /* only used in ff.c */
}

void BIG_256_56_dnorm(DBIG_256_56 a)
{
    int i;
    chunk d, carry;

    carry=a[0]>>BASEBITS_256_56;
    a[0]&=BMASK_256_56;

    for (i = 1; i < DNLEN_256_56 - 1; i++)
    {
        d = a[i] + carry;
        a[i] = d & BMASK_256_56;
        carry = d >> BASEBITS_256_56;
    }
    a[DNLEN_256_56 - 1] = (a[DNLEN_256_56 - 1] + carry);
#ifdef DEBUG_NORM
    a[DMPV_256_56] = 1;
    a[DMNV_256_56] = 0;
#endif
}

/* Compare a and b. Return 1 for a>b, -1 for a<b, 0 for a==b */
/* a and b MUST be normalised before call */
/* sodium constant time implementation */

int BIG_256_56_comp(BIG_256_56 a, BIG_256_56 b)
{
    int i;
    chunk gt=0; chunk eq=1;
    for (i = NLEN_256_56-1; i>=0; i--)
    {
        gt |= ((b[i]-a[i]) >> BASEBITS_256_56) & eq;
        eq &= ((b[i]^a[i])-1) >> BASEBITS_256_56;
    }
    return (int)(gt+gt+eq-1);
}

int BIG_256_56_dcomp(DBIG_256_56 a, DBIG_256_56 b)
{
    int i;
    chunk gt=0; chunk eq=1;
    for (i = DNLEN_256_56-1; i>=0; i--)
    {
        gt |= ((b[i]-a[i]) >> BASEBITS_256_56) & eq;
        eq &= ((b[i]^a[i])-1) >> BASEBITS_256_56;
    }
    return (int)(gt+gt+eq-1);
}

/* return number of bits in a */
/* SU= 8 */
int BIG_256_56_nbits(BIG_256_56 a)
{
    int bts, k = NLEN_256_56 - 1;
    BIG_256_56 t;
    chunk c;
    BIG_256_56_copy(t, a);
    BIG_256_56_norm(t);
    while (k >= 0 && t[k] == 0) k--;
    if (k < 0) return 0;
    bts = BASEBITS_256_56 * k;
    c = t[k];
    while (c != 0)
    {
        c /= 2;
        bts++;
    }
    return bts;
}

/* SU= 8, Calculate number of bits in a DBIG - output normalised */
int BIG_256_56_dnbits(DBIG_256_56 a)
{
    int bts, k = DNLEN_256_56 - 1;
    DBIG_256_56 t;
    chunk c;
    BIG_256_56_dcopy(t, a);
    BIG_256_56_dnorm(t);
    while (k >= 0 && t[k] == 0) k--;
    if (k < 0) return 0;
    bts = BASEBITS_256_56 * k;
    c = t[k];
    while (c != 0)
    {
        c /= 2;
        bts++;
    }
    return bts;
}

void BIG_256_56_ctmod(BIG_256_56 b,BIG_256_56 m,int bd)
{
    int k=bd;
    BIG_256_56 r,c;
    BIG_256_56_copy(c,m);
    BIG_256_56_norm(b);

    BIG_256_56_shl(c,k);
    while (k>=0)
    {
        BIG_256_56_sub(r, b, c);
        BIG_256_56_norm(r);
        BIG_256_56_cmove(b, r, 1 - ((r[NLEN_256_56 - 1] >> (CHUNK - 1)) & 1));
        BIG_256_56_fshr(c, 1);
        k--;
    }
}

/* Set b=b mod c */
/* SU= 16 */
void BIG_256_56_mod(BIG_256_56 b, BIG_256_56 m)
{
    int k=BIG_256_56_nbits(b)-BIG_256_56_nbits(m);
    if (k<0) k=0;
    BIG_256_56_ctmod(b,m,k);
}

// Set a=b mod m in constant time (if bd is known at compile time)
// bd is Max number of bits in b - Actual number of bits in m
void BIG_256_56_ctdmod(BIG_256_56 a, DBIG_256_56 b, BIG_256_56 m, int bd)
{
    int k=bd;
    DBIG_256_56 c,r;
    BIG_256_56_dscopy(c,m);
    BIG_256_56_dnorm(b);

    BIG_256_56_dshl(c,k);
    while (k>=0)
    {
        BIG_256_56_dsub(r, b, c);
        BIG_256_56_dnorm(r);
        BIG_256_56_dcmove(b, r, 1 - ((r[DNLEN_256_56 - 1] >> (CHUNK - 1)) & 1));
        BIG_256_56_dshr(c, 1);
        k--;
    }
    BIG_256_56_sdcopy(a,b);
}


/* Set a=b mod m, b is destroyed. Slow but rarely used. */
void BIG_256_56_dmod(BIG_256_56 a, DBIG_256_56 b, BIG_256_56 m)
{
    int k=BIG_256_56_dnbits(b)-BIG_256_56_nbits(m);
    if (k<0) k=0;
    BIG_256_56_ctdmod(a,b,m,k);
}

// a=b/m  in constant time (if bd is known at compile time)
// bd is Max number of bits in b - Actual number of bits in m
void BIG_256_56_ctddiv(BIG_256_56 a,DBIG_256_56 b,BIG_256_56 m,int bd)
{
    int d,k=bd;
    DBIG_256_56 c,dr;
    BIG_256_56 e,r;
    BIG_256_56_dscopy(c,m);
    BIG_256_56_dnorm(b);

    BIG_256_56_zero(a);
    BIG_256_56_zero(e);
    BIG_256_56_inc(e, 1);

    BIG_256_56_shl(e,k);
    BIG_256_56_dshl(c,k);

    while (k >= 0)
    {
        BIG_256_56_dsub(dr, b, c);
        BIG_256_56_dnorm(dr);
        d = 1 - ((dr[DNLEN_256_56 - 1] >> (CHUNK - 1)) & 1);
        BIG_256_56_dcmove(b, dr, d);

        BIG_256_56_add(r, a, e);
        BIG_256_56_norm(r);
        BIG_256_56_cmove(a, r, d);

        BIG_256_56_dshr(c, 1);
        BIG_256_56_fshr(e, 1);
        k--;
    }
}

/* Set a=b/m,  b is destroyed. Slow but rarely used. */
void BIG_256_56_ddiv(BIG_256_56 a, DBIG_256_56 b, BIG_256_56 m)
{
    int k=BIG_256_56_dnbits(b)-BIG_256_56_nbits(m);
    if (k<0) k=0;
    BIG_256_56_ctddiv(a,b,m,k);
}

// a=a/m  in constant time (if bd is known at compile time)
// bd is Max number of bits in b - Actual number of bits in m
void BIG_256_56_ctsdiv(BIG_256_56 b,BIG_256_56 m,int bd)
{
    int d, k=bd;
    BIG_256_56 e,a,r,c;
    BIG_256_56_norm(b);
    BIG_256_56_copy(a,b);
    BIG_256_56_copy(c,m);
    BIG_256_56_zero(b);
    BIG_256_56_zero(e);
    BIG_256_56_inc(e, 1);

    BIG_256_56_shl(c,k);
    BIG_256_56_shl(e,k);

    while (k >= 0)
    {
        BIG_256_56_sub(r, a, c);
        BIG_256_56_norm(r);
        d = 1 - ((r[NLEN_256_56 - 1] >> (CHUNK - 1)) & 1);
        BIG_256_56_cmove(a, r, d);

        BIG_256_56_add(r, b, e);
        BIG_256_56_norm(r);
        BIG_256_56_cmove(b, r, d);

        BIG_256_56_fshr(c, 1);
        BIG_256_56_fshr(e, 1);

        k--;
    }
}

void BIG_256_56_sdiv(BIG_256_56 b, BIG_256_56 m)
{
    int k=BIG_256_56_nbits(b)-BIG_256_56_nbits(m);
    if (k<0) k=0;
    BIG_256_56_ctsdiv(b,m,k);
}

/* return LSB of a */
int BIG_256_56_parity(BIG_256_56 a)
{
    return a[0] % 2;
}

/* return n-th bit of a */
/* SU= 16 */
int BIG_256_56_bit(BIG_256_56 a, int n)
{
    return (int)((a[n / BASEBITS_256_56] & ((chunk)1 << (n % BASEBITS_256_56))) >> (n%BASEBITS_256_56));
//    if (a[n / BASEBITS_256_56] & ((chunk)1 << (n % BASEBITS_256_56))) return 1;
//    else return 0;
}

/* return last n bits of a, where n is small < BASEBITS */
/* SU= 16 */
int BIG_256_56_lastbits(BIG_256_56 a, int n)
{
    int msk = (1 << n) - 1;
    BIG_256_56_norm(a);
    return ((int)a[0])&msk;
}

/* get 8*MODBYTES size random number */
void BIG_256_56_random(BIG_256_56 m, csprng *rng)
{
    int i, b, j = 0, r = 0;
    int len = 8 * MODBYTES_256_56;

    BIG_256_56_zero(m);
    /* generate random BIG */
    for (i = 0; i < len; i++)
    {
        if (j == 0) r = RAND_byte(rng);
        else r >>= 1;
        b = r & 1;
        BIG_256_56_shl(m, 1);
        m[0] += b;
        j++;
        j &= 7;
    }

#ifdef DEBUG_NORM
    m[MPV_256_56] = 1;
    m[MNV_256_56] = 0;
#endif
}

/* get random BIG from rng, modulo q. Done one bit at a time, so its portable */

void BIG_256_56_randomnum(BIG_256_56 m, BIG_256_56 q, csprng *rng)
{
    int i, b, j = 0, r = 0;
    int n=2 * BIG_256_56_nbits(q);
    DBIG_256_56 d;
    BIG_256_56_dzero(d);
    /* generate random DBIG */
    for (i = 0; i < n; i++)
    {
        if (j == 0) r = RAND_byte(rng);
        else r >>= 1;
        b = r & 1;
        BIG_256_56_dshl(d, 1);
        d[0] += b;
        j++;
        j &= 7;
    }
    /* reduce modulo a BIG. Removes bias */
    BIG_256_56_dmod(m, d, q);
#ifdef DEBUG_NORM
    m[MPV_256_56] = 1;
    m[MNV_256_56] = 0;
#endif
}

/* create randum BIG less than r and less than trunc bits */
void BIG_256_56_randtrunc(BIG_256_56 s, BIG_256_56 r, int trunc, csprng *rng)
{
    BIG_256_56_randomnum(s, r, rng);
    if (BIG_256_56_nbits(r) > trunc)
        BIG_256_56_mod2m(s, trunc);
}


/* Set r=a*b mod m */
/* SU= 96 */
void BIG_256_56_modmul(BIG_256_56 r, BIG_256_56 a1, BIG_256_56 b1, BIG_256_56 m)
{
    DBIG_256_56 d;
    BIG_256_56 a, b;
    BIG_256_56_copy(a, a1);
    BIG_256_56_copy(b, b1);
    BIG_256_56_mod(a, m);
    BIG_256_56_mod(b, m);

    BIG_256_56_mul(d, a, b);
    BIG_256_56_ctdmod(r, d, m, BIG_256_56_nbits(m));
}

/* Set a=a*a mod m */
/* SU= 88 */
void BIG_256_56_modsqr(BIG_256_56 r, BIG_256_56 a1, BIG_256_56 m)
{
    DBIG_256_56 d;
    BIG_256_56 a;
    BIG_256_56_copy(a, a1);
    BIG_256_56_mod(a, m);
    BIG_256_56_sqr(d, a);
    BIG_256_56_ctdmod(r, d, m, BIG_256_56_nbits(m));
}

/* Set r=-a mod m */
/* SU= 16 */
void BIG_256_56_modneg(BIG_256_56 r, BIG_256_56 a1, BIG_256_56 m)
{
    BIG_256_56 a;
    BIG_256_56_copy(a, a1);
    BIG_256_56_mod(a, m);
    BIG_256_56_sub(r, m, a); BIG_256_56_norm(r);
}

/* Set r=a+b mod m */
void BIG_256_56_modadd(BIG_256_56 r, BIG_256_56 a1, BIG_256_56 b1, BIG_256_56 m)
{
    BIG_256_56 a, b;
    BIG_256_56_copy(a, a1);
    BIG_256_56_copy(b, b1);
    BIG_256_56_mod(a, m);
    BIG_256_56_mod(b, m);
    BIG_256_56_add(r,a,b); BIG_256_56_norm(r);
    BIG_256_56_ctmod(r,m,1);
}

/* Set a=a/b mod m */
/* SU= 136 */
void BIG_256_56_moddiv(BIG_256_56 r, BIG_256_56 a1, BIG_256_56 b1, BIG_256_56 m)
{
    DBIG_256_56 d;
    BIG_256_56 z;
    BIG_256_56 a, b;
    BIG_256_56_copy(a, a1);
    BIG_256_56_copy(b, b1);

    BIG_256_56_mod(a, m);
    BIG_256_56_invmodp(z, b, m);

    BIG_256_56_mul(d, a, z);
    BIG_256_56_ctdmod(r, d, m, BIG_256_56_nbits(m));
}

/* Get jacobi Symbol (a/p). Returns 0, 1 or -1 */
/* SU= 216 */
int BIG_256_56_jacobi(BIG_256_56 a, BIG_256_56 p)
{
    int n8, k, m = 0;
    BIG_256_56 t, x, n, zilch, one;
    BIG_256_56_one(one);
    BIG_256_56_zero(zilch);
    if (BIG_256_56_parity(p) == 0 || BIG_256_56_comp(a, zilch) == 0 || BIG_256_56_comp(p, one) <= 0) return 0;
    BIG_256_56_norm(a);
    BIG_256_56_copy(x, a);
    BIG_256_56_copy(n, p);
    BIG_256_56_mod(x, p);

    while (BIG_256_56_comp(n, one) > 0)
    {
        if (BIG_256_56_comp(x, zilch) == 0) return 0;
        n8 = BIG_256_56_lastbits(n, 3);
        k = 0;
        while (BIG_256_56_parity(x) == 0)
        {
            k++;
            BIG_256_56_shr(x, 1);
        }
        if (k % 2 == 1) m += (n8 * n8 - 1) / 8;
        m += (n8 - 1) * (BIG_256_56_lastbits(x, 2) - 1) / 4;
        BIG_256_56_copy(t, n);

        BIG_256_56_mod(t, x);
        BIG_256_56_copy(n, x);
        BIG_256_56_copy(x, t);
        m %= 2;

    }
    if (m == 0) return 1;
    else return -1;
}

/* Arazi and Qi inversion mod 256 */
static int invmod256(int a)
{
    int U, t1, t2, b, c;
    t1 = 0;
    c = (a >> 1) & 1;
    t1 += c;
    t1 &= 1;
    t1 = 2 - t1;
    t1 <<= 1;
    U = t1 + 1;

// i=2
    b = a & 3;
    t1 = U * b;
    t1 >>= 2;
    c = (a >> 2) & 3;
    t2 = (U * c) & 3;
    t1 += t2;
    t1 *= U;
    t1 &= 3;
    t1 = 4 - t1;
    t1 <<= 2;
    U += t1;

// i=4
    b = a & 15;
    t1 = U * b;
    t1 >>= 4;
    c = (a >> 4) & 15;
    t2 = (U * c) & 15;
    t1 += t2;
    t1 *= U;
    t1 &= 15;
    t1 = 16 - t1;
    t1 <<= 4;
    U += t1;

    return U;
}

/* a=1/a mod 2^BIGBITS. This is very fast! */
void BIG_256_56_invmod2m(BIG_256_56 a)
{
    int i;
    BIG_256_56 U, t1, b, c;
    BIG_256_56_zero(U);
    BIG_256_56_inc(U, invmod256(BIG_256_56_lastbits(a, 8)));
    for (i = 8; i < BIGBITS_256_56; i <<= 1)
    {
        BIG_256_56_norm(U);
        BIG_256_56_copy(b, a);
        BIG_256_56_mod2m(b, i);  // bottom i bits of a

        BIG_256_56_smul(t1, U, b);
        BIG_256_56_shr(t1, i); // top i bits of U*b

        BIG_256_56_copy(c, a);
        BIG_256_56_shr(c, i);
        BIG_256_56_mod2m(c, i); // top i bits of a

        BIG_256_56_smul(b, U, c);
        BIG_256_56_mod2m(b, i); // bottom i bits of U*c

        BIG_256_56_add(t1, t1, b);
        BIG_256_56_norm(t1);
        BIG_256_56_smul(b, t1, U);
        BIG_256_56_copy(t1, b); // (t1+b)*U
        BIG_256_56_mod2m(t1, i);               // bottom i bits of (t1+b)*U

        BIG_256_56_one(b);
        BIG_256_56_shl(b, i);
        BIG_256_56_sub(t1, b, t1);
        BIG_256_56_norm(t1);

        BIG_256_56_shl(t1, i);

        BIG_256_56_add(U, U, t1);
    }
    BIG_256_56_copy(a, U);
    BIG_256_56_norm(a);
    BIG_256_56_mod2m(a, BIGBITS_256_56);
}

/* Set r=1/a mod p. Binary method */
/* SU= 240 */
void BIG_256_56_invmodp(BIG_256_56 r, BIG_256_56 a, BIG_256_56 p)
{
    BIG_256_56 u, v, x1, x2, t, one;

    BIG_256_56_mod(a, p);
    if (BIG_256_56_iszilch(a))
    {
        BIG_256_56_zero(r);
        return;
    }

    BIG_256_56_copy(u, a);
    BIG_256_56_copy(v, p);
    BIG_256_56_one(one);
    BIG_256_56_copy(x1, one);
    BIG_256_56_zero(x2);

    while (BIG_256_56_comp(u, one) != 0 && BIG_256_56_comp(v, one) != 0)
    {
        while (BIG_256_56_parity(u) == 0)
        {
            BIG_256_56_fshr(u, 1);
            BIG_256_56_add(t,x1,p);
            BIG_256_56_cmove(x1,t,BIG_256_56_parity(x1));
            BIG_256_56_norm(x1);
            BIG_256_56_fshr(x1,1);
        }
        while (BIG_256_56_parity(v) == 0)
        {
            BIG_256_56_fshr(v, 1);
            BIG_256_56_add(t,x2,p);
            BIG_256_56_cmove(x2,t,BIG_256_56_parity(x2));
            BIG_256_56_norm(x2);
            BIG_256_56_fshr(x2,1);
        }
        if (BIG_256_56_comp(u, v) >= 0)
        {
            BIG_256_56_sub(u, u, v);
            BIG_256_56_norm(u);
            BIG_256_56_add(t,x1,p);
            BIG_256_56_cmove(x1,t,(BIG_256_56_comp(x1,x2)>>1)&1); /* move if x1<x2 */
            BIG_256_56_sub(x1,x1,x2);
            BIG_256_56_norm(x1);
        }
        else
        {
            BIG_256_56_sub(v, v, u);
            BIG_256_56_norm(v);
            BIG_256_56_add(t,x2,p);
            BIG_256_56_cmove(x2,t,(BIG_256_56_comp(x2,x1)>>1)&1); /* move if x2<x1 */
            BIG_256_56_sub(x2,x2,x1);
            BIG_256_56_norm(x2);
        }
    }
    BIG_256_56_copy(r,x1);
    BIG_256_56_cmove(r,x2,BIG_256_56_comp(u,one)&1);
}

/* set x = x mod 2^m */
void BIG_256_56_mod2m(BIG_256_56 x, int m)
{
    int i, wd, bt;
    chunk msk;
    BIG_256_56_norm(x);

    wd = m / BASEBITS_256_56;
    bt = m % BASEBITS_256_56;
    msk = ((chunk)1 << bt) - 1;
    x[wd] &= msk;
    for (i = wd + 1; i < NLEN_256_56; i++) x[i] = 0;
}

// new
/* Convert to DBIG number from byte array of given length */
void BIG_256_56_dfromBytesLen(DBIG_256_56 a, char *b, int s)
{
    int i, len = s;
    BIG_256_56_dzero(a);

    for (i = 0; i < len; i++)
    {
        BIG_256_56_dshl(a, 8);
        a[0] += (int)(unsigned char)b[i];
    }
#ifdef DEBUG_NORM
    a[DMPV_256_56] = 1;
    a[DMNV_256_56] = 0;
#endif
}
