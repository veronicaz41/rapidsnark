#include "fr.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <string>
#include <cstring>
#include <iostream>


static mpz_t q;
static mpz_t zero;
static mpz_t one;
static mpz_t mask;
static size_t nBits;
static bool initialized = false;

#ifndef USE_ASM

static FrElement    Fr_q     = {0, 0x80000000, {0x43e1f593f0000001,0x2833e84879b97091,0xb85045b68181585d,0x30644e72e131a029}};
static uint64_t     Fr_rawq[]  =               {0x43e1f593f0000001,0x2833e84879b97091,0xb85045b68181585d,0x30644e72e131a029, 0};
static FrElement    Fr_R3    = {0, 0x80000000, {0x5e94d8e1b4bf0040,0x2a489cbe1cfbb6b8,0x893cc664a19fcfed,0x0cf8594b7fcc657c}};
static FrRawElement Fr_rawR3 =                 {0x5e94d8e1b4bf0040,0x2a489cbe1cfbb6b8,0x893cc664a19fcfed,0x0cf8594b7fcc657c};
static FrRawElement Fr_rawR2 =                 {0x1bb8e645ae216da7,0x53fe3ab1e35c59e3,0x8c49833d53bb8085,0x0216d0b17f4e44a5};
static uint64_t     Fr_np    = {0xc2e1f593efffffff};


#endif

void Fr_toMpz(mpz_t r, PFrElement pE) {
    FrElement tmp;
    Fr_toNormal(&tmp, pE);
    if (!(tmp.type & Fr_LONG)) {
        mpz_set_si(r, tmp.shortVal);
        if (tmp.shortVal<0) {
            mpz_add(r, r, q);
        }
    } else {
        mpz_import(r, Fr_N64, -1, 8, -1, 0, (const void *)tmp.longVal);
    }
}

void Fr_fromMpz(PFrElement pE, mpz_t v) {
    if (mpz_fits_sint_p(v)) {
        pE->type = Fr_SHORT;
        pE->shortVal = mpz_get_si(v);
    } else {
        pE->type = Fr_LONG;
        for (int i=0; i<Fr_N64; i++) pE->longVal[i] = 0;
        mpz_export((void *)(pE->longVal), NULL, -1, 8, -1, 0, v);
    }
}


bool Fr_init() {
    if (initialized) return false;
    initialized = true;
    mpz_init(q);
    mpz_import(q, Fr_N64, -1, 8, -1, 0, (const void *)Fr_q.longVal);
    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(one, 1);
    nBits = mpz_sizeinbase (q, 2);
    mpz_init(mask);
    mpz_mul_2exp(mask, one, nBits);
    mpz_sub(mask, mask, one);
    return true;
}

void Fr_str2element(PFrElement pE, char const *s) {
    mpz_t mr;
    mpz_init_set_str(mr, s, 10);
    mpz_fdiv_r(mr, mr, q);
    Fr_fromMpz(pE, mr);
    mpz_clear(mr);
}

char *Fr_element2str(PFrElement pE) {
    FrElement tmp;
    mpz_t r;
    if (!(pE->type & Fr_LONG)) {
        if (pE->shortVal>=0) {
            char *r = new char[32];
            sprintf(r, "%d", pE->shortVal);
            return r;
        } else {
            mpz_init_set_si(r, pE->shortVal);
            mpz_add(r, r, q);
        }
    } else {
        Fr_toNormal(&tmp, pE);
        mpz_init(r);
        mpz_import(r, Fr_N64, -1, 8, -1, 0, (const void *)tmp.longVal);
    }
    char *res = mpz_get_str (0, 10, r);
    mpz_clear(r);
    return res;
}

void Fr_idiv(PFrElement r, PFrElement a, PFrElement b) {
    mpz_t ma;
    mpz_t mb;
    mpz_t mr;
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(mr);

    Fr_toMpz(ma, a);
    // char *s1 = mpz_get_str (0, 10, ma);
    // printf("s1 %s\n", s1);
    Fr_toMpz(mb, b);
    // char *s2 = mpz_get_str (0, 10, mb);
    // printf("s2 %s\n", s2);
    mpz_fdiv_q(mr, ma, mb);
    // char *sr = mpz_get_str (0, 10, mr);
    // printf("r %s\n", sr);
    Fr_fromMpz(r, mr);

    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(mr);
}

void Fr_mod(PFrElement r, PFrElement a, PFrElement b) {
    mpz_t ma;
    mpz_t mb;
    mpz_t mr;
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(mr);

    Fr_toMpz(ma, a);
    Fr_toMpz(mb, b);
    mpz_fdiv_r(mr, ma, mb);
    Fr_fromMpz(r, mr);

    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(mr);
}

void Fr_pow(PFrElement r, PFrElement a, PFrElement b) {
    mpz_t ma;
    mpz_t mb;
    mpz_t mr;
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(mr);

    Fr_toMpz(ma, a);
    Fr_toMpz(mb, b);
    mpz_powm(mr, ma, mb, q);
    Fr_fromMpz(r, mr);

    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(mr);
}

void Fr_inv(PFrElement r, PFrElement a) {
    mpz_t ma;
    mpz_t mr;
    mpz_init(ma);
    mpz_init(mr);

    Fr_toMpz(ma, a);
    mpz_invert(mr, ma, q);
    Fr_fromMpz(r, mr);
    mpz_clear(ma);
    mpz_clear(mr);
}

void Fr_div(PFrElement r, PFrElement a, PFrElement b) {
    FrElement tmp;
    Fr_inv(&tmp, b);
    Fr_mul(r, a, &tmp);
}

void Fr_fail() {
    assert(false);
}


RawFr::RawFr() {
    Fr_init();
    fromString(fZero, "0");
    fromString(fOne, "1");
    neg(fNegOne, fOne);
}

RawFr::~RawFr() {
}

void RawFr::fromString(Element &r, std::string s) {
    mpz_t mr;
    mpz_init_set_str(mr, s.c_str(), 10);
    mpz_fdiv_r(mr, mr, q);
    for (int i=0; i<Fr_N64; i++) r.v[i] = 0;
    mpz_export((void *)(r.v), NULL, -1, 8, -1, 0, mr);
    Fr_rawToMontgomery(r.v,r.v);
    mpz_clear(mr);
}

void RawFr::fromUI(Element &r, unsigned long int v) {
    mpz_t mr;
    mpz_init(mr);
    mpz_set_ui(mr, v);
    for (int i=0; i<Fr_N64; i++) r.v[i] = 0;
    mpz_export((void *)(r.v), NULL, -1, 8, -1, 0, mr);
    Fr_rawToMontgomery(r.v,r.v);
    mpz_clear(mr);
}



std::string RawFr::toString(Element &a, uint32_t radix) {
    Element tmp;
    mpz_t r;
    Fr_rawFromMontgomery(tmp.v, a.v);
    mpz_init(r);
    mpz_import(r, Fr_N64, -1, 8, -1, 0, (const void *)(tmp.v));
    char *res = mpz_get_str (0, radix, r);
    mpz_clear(r);
    std::string resS(res);
    free(res);
    return resS;
}

void RawFr::inv(Element &r, Element &a) {
    mpz_t mr;
    mpz_init(mr);
    mpz_import(mr, Fr_N64, -1, 8, -1, 0, (const void *)(a.v));
    mpz_invert(mr, mr, q);


    for (int i=0; i<Fr_N64; i++) r.v[i] = 0;
    mpz_export((void *)(r.v), NULL, -1, 8, -1, 0, mr);

    Fr_rawMMul(r.v, r.v,Fr_rawR3);
    mpz_clear(mr);
}

void RawFr::div(Element &r, Element &a, Element &b) {
    Element tmp;
    inv(tmp, b);
    mul(r, a, tmp);
}

#define BIT_IS_SET(s, p) (s[p>>3] & (1 << (p & 0x7)))
void RawFr::exp(Element &r, Element &base, uint8_t* scalar, unsigned int scalarSize) {
    bool oneFound = false;
    Element copyBase;
    copy(copyBase, base);
    for (int i=scalarSize*8-1; i>=0; i--) {
        if (!oneFound) {
            if ( !BIT_IS_SET(scalar, i) ) continue;
            copy(r, copyBase);
            oneFound = true;
            continue;
        }
        square(r, r);
        if ( BIT_IS_SET(scalar, i) ) {
            mul(r, r, copyBase);
        }
    }
    if (!oneFound) {
        copy(r, fOne);
    }
}

void RawFr::toMpz(mpz_t r, Element &a) {
    Element tmp;
    Fr_rawFromMontgomery(tmp.v, a.v);
    mpz_import(r, Fr_N64, -1, 8, -1, 0, (const void *)tmp.v);
}

void RawFr::fromMpz(Element &r, mpz_t a) {
    for (int i=0; i<Fr_N64; i++) r.v[i] = 0;
    mpz_export((void *)(r.v), NULL, -1, 8, -1, 0, a);
    Fr_rawToMontgomery(r.v, r.v);
}


static bool init = Fr_init();

RawFr RawFr::field;

#ifndef USE_ASM

/*****************************************************************************************
 * ASM Functions to C/C++ using GNU MP Lib Begin
******************************************************************************************/

void Fr_rawAdd(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    uint64_t carry = mpn_add_n(pRawResult, pRawA, pRawB, Fr_N64);

    if(carry || mpn_cmp(pRawResult, Fr_rawq, Fr_N64) >= 0)
    {
        mpn_sub_n(pRawResult, pRawResult, Fr_rawq, Fr_N64);
    }
}

void Fr_rawSub(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    uint64_t carry = mpn_sub_n(pRawResult, pRawA, pRawB, Fr_N64);

    if(carry)
    {
        mpn_add_n(pRawResult, pRawResult, Fr_rawq, Fr_N64);
    }
}

void Fr_rawNeg(FrRawElement pRawResult, FrRawElement pRawA)
{
    const uint64_t zero[Fr_N64] = {0, 0, 0, 0};

    if (mpn_cmp(pRawA, zero, Fr_N64) != 0)
    {
        mpn_sub_n(pRawResult, Fr_rawq, pRawA, Fr_N64);
    }
    else
    {
        mpn_copyi(pRawResult, zero, Fr_N64);
    }
}

void Fr_rawCopy(FrRawElement pRawResult, FrRawElement pRawA)
{
    mpn_copyi(pRawResult, pRawA, Fr_N64);
}

void Fr_copy(PFrElement r, PFrElement a)
{
    r->shortVal = a->shortVal;
    r->type = a->type;
    std::memcpy(r->longVal, a->longVal, sizeof(FrRawElement));
}

int Fr_rawIsEq(FrRawElement pRawA, FrRawElement pRawB)
{
    return mpn_cmp(pRawA, pRawB, Fr_N64) == 0;
}

void Fr_rawMMul(FrRawElement pRawResult, FrRawElement pRawA, FrRawElement pRawB)
{
    const mp_size_t  N = Fr_N64+1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    // FirstLoop 0
    product0[4] = mpn_mul_1(product0, pRawB, Fr_N64, pRawA[0]);

    // Second Loop 0
    np0 = Fr_np * product0[0];
    product1[1] = mpn_addmul_1(product0, mq, N, np0);

    // FirstLoop 1
    product1[4] = mpn_addmul_1(product1,  pRawB, Fr_N64, pRawA[1]);
    mpn_add(product1, product1, N, product0+1, N-1);

    // Second Loop 1
    np0 = Fr_np * product1[0];
    product2[1] = mpn_addmul_1(product1, mq, N, np0);

    // FirstLoop 2
    product2[4] = mpn_addmul_1(product2,  pRawB, Fr_N64, pRawA[2]);
    mpn_add(product2, product2, N, product1+1, N-1);

    // Second Loop 2
    np0 = Fr_np * product2[0];
    product3[1] = mpn_addmul_1(product2, mq, N, np0);

    // FirstLoop 3
    product3[4] = mpn_addmul_1(product3, pRawB, Fr_N64, pRawA[3]);
    mpn_add(product3, product3, N, product2+1, N-1);

    // Second Loop 3
    np0 = Fr_np * product3[0];
    mpn_addmul_1(product3, mq, N, np0);

    mpn_copyi(pRawResult,  product3+1, Fr_N64);

    if (mpn_cmp(pRawResult, mq, Fr_N64) >= 0)
    {
        mpn_sub_n(pRawResult, pRawResult, mq, Fr_N64);
    }
}

void Fr_rawMSquare(FrRawElement pRawResult, FrRawElement pRawA)
{
    Fr_rawMMul(pRawResult, pRawA, pRawA);
}

void Fr_rawMMul1(FrRawElement pRawResult, FrRawElement pRawA, uint64_t pRawB)
{
    const mp_size_t  N = Fr_N64+1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N] = {0};
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    // FirstLoop 0
    product0[4] = mpn_mul_1(product0, pRawA, Fr_N64, pRawB);

    // Second Loop 0
    np0 = Fr_np * product0[0];
    product1[1] = mpn_addmul_1(product0, mq, N, np0);
    mpn_add(product1, product1, N, product0+1, N-1);

    // Second Loop 1
    np0 = Fr_np * product1[0];
    product2[1] = mpn_addmul_1(product1, mq, N, np0);
    mpn_add(product2, product2, N, product1+1, N-1);

    // Second Loop 2
    np0 = Fr_np * product2[0];
    product3[1] = mpn_addmul_1(product2, mq, N, np0);
    mpn_add(product3, product3, N, product2+1, N-1);

    // Second Loop 3
    np0 = Fr_np * product3[0];
    mpn_addmul_1(product3, mq, N, np0);

    mpn_copyi(pRawResult,  product3+1, Fr_N64);

    if (mpn_cmp(pRawResult, mq, Fr_N64) >= 0)
    {
        mpn_sub_n(pRawResult, pRawResult, mq, Fr_N64);
    }
}

void Fr_rawToMontgomery(FrRawElement pRawResult, FrRawElement pRawA)
{
    Fr_rawMMul(pRawResult, pRawA, Fr_rawR2);
}

void Fr_rawFromMontgomery(FrRawElement pRawResult, FrRawElement pRawA)
{
    const mp_size_t  N = Fr_N64+1;
    const uint64_t  *mq = Fr_rawq;

    uint64_t  np0;

    uint64_t  product0[N];
    uint64_t  product1[N] = {0};
    uint64_t  product2[N] = {0};
    uint64_t  product3[N] = {0};

    // FirstLoop 0
    mpn_copyi(product0, pRawA, Fr_N64); product0[4] = 0;

    // Second Loop 0
    np0 = Fr_np * product0[0];
    product1[1] = mpn_addmul_1(product0, mq, N, np0);
    mpn_add(product1, product1, N, product0+1, N-1);

    // Second Loop 1
    np0 = Fr_np * product1[0];
    product2[1] = mpn_addmul_1(product1, mq, N, np0);
    mpn_add(product2, product2, N, product1+1, N-1);

    // Second Loop 2
    np0 = Fr_np * product2[0];
    product3[1] = mpn_addmul_1(product2, mq, N, np0);
    mpn_add(product3, product3, N, product2+1, N-1);

    // Second Loop 3
    np0 = Fr_np * product3[0];
    mpn_addmul_1(product3, mq, N, np0);

    mpn_copyi(pRawResult,  product3+1, Fr_N64);

    if (mpn_cmp(pRawResult, mq, Fr_N64) >= 0)
    {
        mpn_sub_n(pRawResult, pRawResult, mq, Fr_N64);
    }
}

void Fr_toNormal(PFrElement r, PFrElement a)
{
    if (a->type == Fr_LONGMONTGOMERY)
    {
        r->type = Fr_LONG;
        Fr_rawFromMontgomery(r->longVal, a->longVal);
    }
    else
    {
        Fr_copy(r, a);
    }
}

int Fr_rawIsZero(FrRawElement pRawB)
{
    for (int i=0; i<Fr_N64; i++)
    {
        if (pRawB[i] != 0)
        {
            return 0;
        }

    }
    return 1;
}

void Fr_rawSwap(FrRawElement pRawResult, FrRawElement pRawA)
{
    FrRawElement tmp = {0};
    for (int i=0; i<Fr_N64; i++)
    {
        tmp[i] = pRawResult[i];
        pRawResult[i] = pRawA[i];
        pRawA[i] = tmp[i];
    }
}

void rawCopyS2L(PFrElement r, int64_t temp);
void mul_s1s2(PFrElement r, PFrElement a, PFrElement b);
void mul_l1nl2n(PFrElement r, PFrElement a, PFrElement b);
void mul_l1ml2n(PFrElement r,PFrElement a,PFrElement b);
void mul_l1nl2m(PFrElement r, PFrElement a, PFrElement b);
void mul_l1ml2m(PFrElement r,PFrElement a,PFrElement b);

void mul_l1ns2n(PFrElement r,PFrElement a,PFrElement b);
void mul_s1nl2n(PFrElement r,PFrElement a,PFrElement b);
void mul_l1ms2n(PFrElement r,PFrElement a,PFrElement b);
void mul_s1nl2m(PFrElement r,PFrElement a,PFrElement b);

void mul_l1ns2m(PFrElement r,PFrElement a,PFrElement b);
void mul_l1ms2m(PFrElement r,PFrElement a,PFrElement b);
void mul_s1ml2n(PFrElement r,PFrElement a,PFrElement b);
void mul_s1ml2m(PFrElement r,PFrElement a,PFrElement b);


void Fr_mul(PFrElement r, PFrElement a, PFrElement b)
{
    //mpz_import(mr3, Fr_N64, -1, 8, -1, 0, (const void *)Fr_rawR3);

    if (a->type & Fr_LONG) // if (mpz_tstbit (ma, 63)) // 2267 ; Check if is short first operand
    {
        // jc     mul_l1
        if (b->type & Fr_LONG) //if (mpz_tstbit (mb, 63)) // 2293 ; Check if is short second operand
        {
            // mul_l1l2
            if (a->type == Fr_LONGMONTGOMERY) // if (mpz_tstbit (ma, 62)) // 2511 ; check if montgomery first
            {
                // mul_l1ml2
                if (b->type == Fr_LONGMONTGOMERY) //if (mpz_tstbit (mb, 62)) // 2554 ; check if montgomery second
                {
                    mul_l1ml2m(r, a, b);
                }
                else
                {
                    mul_l1ml2n(r, a, b);
                }
            }
            else if (b->type == Fr_LONGMONTGOMERY) //if (mpz_tstbit (mb, 62)) // 2514 ; check if montgomery second
            {
                mul_l1nl2m(r, a, b);
            }
            else
            {
                mul_l1nl2n(r, a, b);
            }
        }
        //mul_l1s2:
        else if (a->type == Fr_LONGMONTGOMERY) //if (mpz_tstbit (ma, 62)) // 2298 ; check if montgomery first
        {
            // mul_l1ms2
            if (b->type == Fr_SHORT) //if (mpz_tstbit (mb, 62)) // 2358 ; check if montgomery second
            {

                mul_l1ms2n(r, a, b);
            }
            else
            {
                mul_l1ms2m(r, a, b);
            }

        }
        // mul_l1ns2
        else
        {
            if (b->type == Fr_SHORTMONTGOMERY) //if (mpz_tstbit (mb, 62)) // 2301 ; check if montgomery second
            {
                mul_l1ns2m(r, a, b);
            }
            else
            {
                mul_l1ns2n(r, a, b);
            }
        }
    }
    else if (b->type & Fr_LONG)//if (mpz_tstbit (mb, 63)) // 2269  ; Check if is short second operand
    {
        // mul_s1l2
        if (a->type == Fr_SHORTMONTGOMERY)//if (mpz_tstbit (ma, 62)) // 2406  ; check if montgomery first
        {
            // mul_s1ml2
            if (b->type == Fr_LONGMONTGOMERY)//if (mpz_tstbit (mb, 62)) // 2479  ; check if montgomery second
            {
                mul_s1ml2m(r, a, b);
            }
            else
            {
                mul_s1ml2n(r,a, b);
            }
        }
        // mul_s1nl2
        else if (b->type == Fr_LONGMONTGOMERY) //if (mpz_tstbit (mb, 62)) // 2409; check if montgomery second
        {
            mul_s1nl2m(r, a, b);
        }
        // mul_s1nl2n
        else
        {
            mul_s1nl2n(r, a, b);
        }
    }
    else
    {
         mul_s1s2(r, a, b);
    }
}

void mul_s1s2(PFrElement r, PFrElement a, PFrElement b)
{
    mpz_t rax;
    mpz_init(rax);

    int64_t temp = (int64_t)a->shortVal * (int64_t)b->shortVal;
    r->longVal[0] = temp;
    mpz_import(rax, 1, -1, 8, -1, 0, (const void *)r);
    // mul_manageOverflow
    if (!mpz_fits_sint_p(rax))
    {
        rawCopyS2L(r, temp);
    }
    else
    {
        r->type = Fr_LONG;
    }
    mpz_clear(rax);
}

void rawCopyS2L(PFrElement pResult, int64_t val)
{
    pResult->type = Fr_LONG;
    pResult->shortVal = 0;

    uint64_t *result = pResult->longVal;

    mpn_zero(result, Fr_N64);
    result[0] = val;

    if (val < 0)
    {
        mpn_com(result+1, result+1, Fr_N64-1);
        mpn_add_n(result, result, Fr_rawq, Fr_N64);
    }
}

void mul_l1nl2n(PFrElement r,PFrElement a,PFrElement b)
{
    FrElement tmp1;
    FrElement tmp2;

    r->type = Fr_LONGMONTGOMERY;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);

    tmp1.type = Fr_LONG;
    tmp2.type = Fr_LONG;
    tmp1.shortVal = 0;
    tmp2.shortVal = 0;
    for (int i=0; i<Fr_N64; i++)
    {
        tmp1.longVal[i] = r->longVal[i];
        tmp2.longVal[i] = Fr_R3.longVal[i];
    }

    Fr_rawMMul(&r->longVal[0], &tmp1.longVal[0], &tmp2.longVal[0]);
}

void mul_l1nl2m(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONG;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_l1ml2m(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONGMONTGOMERY;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_l1ml2n(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONG;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_l1ns2n(PFrElement r,PFrElement a,PFrElement b)
{
    FrElement tmp1;
    FrElement tmp2;
    int32_t   tmp3;

    r->type = Fr_LONGMONTGOMERY;
    if (b->shortVal >= 0)
    {
        // tmp_5:
        Fr_rawMMul1(&r->longVal[0], &a->longVal[0], b->shortVal);
        // tmp_6:
        tmp1.type = Fr_LONG;
        tmp2.type = Fr_LONG;
        tmp1.shortVal = 0;
        tmp2.shortVal = 0;
        for (int i=0; i<Fr_N64; i++)
        {
            tmp1.longVal[i] = r->longVal[i];
            tmp2.longVal[i] = Fr_R3.longVal[i];
        }
        Fr_rawMMul(&r->longVal[0], &tmp1.longVal[0], &tmp2.longVal[0]);
    }
    else
    {
        tmp3 = b->shortVal * (-1);
        Fr_rawMMul1(&r->longVal[0], &a->longVal[0], tmp3);
        Fr_rawNeg(&r->longVal[0], &r->longVal[0]);
        // tmp_6:
        tmp1.type = Fr_LONG;
        tmp2.type = Fr_LONG;
        tmp1.shortVal = 0;
        tmp2.shortVal = 0;
        for (int i=0; i<Fr_N64; i++)
        {
            tmp1.longVal[i] = r->longVal[i];
            tmp2.longVal[i] = Fr_R3.longVal[i];
        }
        Fr_rawMMul(&r->longVal[0], &tmp1.longVal[0], &tmp2.longVal[0]);
    }
}

void mul_s1nl2n(PFrElement r,PFrElement a,PFrElement b)
{
    FrElement tmp1;
    FrElement tmp2;
    int32_t   tmp3;

    r->type = Fr_LONGMONTGOMERY;
    if (a->shortVal >= 0)
    {
        // tmp_9:
        Fr_rawMMul1(&r->longVal[0], &b->longVal[0], a->shortVal);
        // tmp_10:
        tmp1.type = Fr_LONG;
        tmp2.type = Fr_LONG;
        tmp1.shortVal = 0;
        tmp2.shortVal = 0;
        for (int i=0; i<Fr_N64; i++)
        {
            tmp1.longVal[i] = r->longVal[i];
            tmp2.longVal[i] = Fr_R3.longVal[i];
        }
        Fr_rawMMul(&r->longVal[0], &tmp1.longVal[0], &tmp2.longVal[0]);
    }
    else
    {
        tmp3 = a->shortVal * (-1);
        Fr_rawMMul1(&r->longVal[0], &b->longVal[0], tmp3);
        Fr_rawNeg(&r->longVal[0], &r->longVal[0]);
        // tmp_6:
        tmp1.type = Fr_LONG;
        tmp2.type = Fr_LONG;
        tmp1.shortVal = 0;
        tmp2.shortVal = 0;
        for (int i=0; i<Fr_N64; i++)
        {
            tmp1.longVal[i] = r->longVal[i];
            tmp2.longVal[i] = Fr_R3.longVal[i];
        }
        Fr_rawMMul(&r->longVal[0], &tmp1.longVal[0], &tmp2.longVal[0]);
    }
}

void mul_l1ms2n(PFrElement r,PFrElement a,PFrElement b)
{
    FrElement tmp1;
    FrElement tmp2;
    int32_t   tmp3;

    r->type = Fr_LONG;
    if (b->shortVal >= 0)
    {
        // tmp_7:
        Fr_rawMMul1(&r->longVal[0], &a->longVal[0], b->shortVal);
        // tmp_8:
    }
    else
    {
        tmp3 = b->shortVal * (-1);
        Fr_rawMMul1(&r->longVal[0], &a->longVal[0], tmp3);
        Fr_rawNeg(&r->longVal[0], &r->longVal[0]);
        // tmp_8:
    }
}

void mul_s1nl2m(PFrElement r,PFrElement a,PFrElement b)
{
    FrElement tmp1;
    FrElement tmp2;
    int32_t   tmp3;

    r->type = Fr_LONG;
    if (a->shortVal >= 0)
    {
        // tmp_11:
        Fr_rawMMul1(&r->longVal[0], &b->longVal[0], a->shortVal);
        // tmp_12:
    }
    else
    {
        tmp3 = a->shortVal * (-1);
        Fr_rawMMul1(&r->longVal[0], &b->longVal[0], tmp3);
        Fr_rawNeg(&r->longVal[0], &r->longVal[0]);
        // tmp_12:
    }
}

void mul_l1ns2m(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONG;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_l1ms2m(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONGMONTGOMERY;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_s1ml2m(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONGMONTGOMERY;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

void mul_s1ml2n(PFrElement r,PFrElement a,PFrElement b)
{
    r->type = Fr_LONG;
    Fr_rawMMul(&r->longVal[0], &a->longVal[0], &b->longVal[0]);
}

/*****************************************************************************************
 * ASM Functions to C/C++ using GNU MP Lib End
******************************************************************************************/

#endif
