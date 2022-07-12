// Author: Dr Konstantinos Sfyrakis

#ifndef CROSSPLATFORMCONVERSION_H_
#define CROSSPLATFORMCONVERSION_H_

#define SMALL_NUM 0.00000001            // anything that avoids division overflow
#define ABS(x) ((x) >= 0 ? (x) : -(x))  // absolute value

#define SWAP_2(x) ((((x)&0xff) << 8) | ((unsigned short)(x) >> 8))
#define SWAP_4(x) \
    (((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((x) >> 24))
#define FIX_SHORT(x) (*(unsigned short*)&(x) = SWAP_2(*(unsigned short*)&(x)))
#define FIX_INT(x) (*(unsigned int*)&(x) = SWAP_4(*(unsigned int*)&(x)))
#define FIX_LONG(x) FIX_INT(x)
#define FIX_FLOAT(x) FIX_INT(x)
#define DOUBLESWAP(x) ByteSwap((unsigned char*)&x, sizeof(x))

inline void ByteSwap(unsigned char* b, int n) {
    int i = 0;
    int j = n - 1;

    while (i < j) {
        std::swap(b[i], b[j]);
        i++, j--;
    }
}

#endif /*CROSSPLATFORMCONVERSION_H_*/
