#ifndef SEED_ENCODE_H
#define SEED_ENCODE_H

/* ============================================================= */
/*  Fast 4-bit DNA encoding tables for nShift = 62               */
/*  Only A, C, G, T are accepted. All ambiguous bases → invalid  */
/* ============================================================= */

/* wr_table: contribution to the reverse-complement word (wr) */
static const unsigned long long wr_table[16] = {
    0x0000000000000000ULL,   /* 0x0  invalid */
    0xC000000000000000ULL,   /* 0x1  A */
    0x0000000000000000ULL,   /* 0x2  T */
    0x0000000000000000ULL,  
    0x4000000000000000ULL,   /* 0x4  G */
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x8000000000000000ULL,   /* 0x8  C */
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL,   
    0x0000000000000000ULL    /* 0xF  A|G|C|T → ambiguous */
};

/* w_table: contribution to the forward word (w) */
static const unsigned long w_table[16] = {
    0x0UL,   /* 0x0  invalid */
    0x0UL,   /* 0x1  A */
    0x3UL,   /* 0x2  T */
    0x0UL,   
    0x2UL,   /* 0x4  G */
    0x0UL,   
    0x0UL,   
    0x0UL,   
    0x1UL,   /* 0x8  C */
    0x0UL,   
    0x0UL,  
    0x0UL, 
    0x0UL, 
    0x0UL, 
    0x0UL, 
    0x0UL    /* 0xF  invalid */
};

/* valid_base: 1 = acceptable (only pure A, C, G, T), 0 = invalid */
static const unsigned int valid_base[16] = {
    0, 0xffff, 0xffff, 0,     /* 0x0..0x3 */
    0xffff, 0, 0, 0,     /* 0x4..0x7 */
    0xffff, 0, 0, 0,     /* 0x8..0xB */
    0, 0, 0, 0      /* 0xC..0xF */
};

#endif /* SEED_ENCODE_H */
