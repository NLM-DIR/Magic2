#include <stdio.h>
#include "sra_read.h"

int main(void)
{
    SRAObj* sra = SraObjNew("SRR24564586");
    int num_bases = 500;

    while (1) {
        const char* p = SraGetReadBatch(sra, num_bases, 1);
        if (!p) {
            break;
        }
        printf("%s", p);
    }
    SraObjFree(sra);
}
