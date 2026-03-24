#include <stdio.h>
#include "sra_read.h"

int main(void)
{
    SRAObj* sra = SraObjNew("SRR24511885");
    int num_bases = 500;
    const char* p1 = NULL;
    const char* p2 = NULL;

    while (1) {
        SraGetReadBatch(sra, num_bases, 1, 1, &p1, &p2);
        if (!p1) {
            break;
        }
        printf("%s", p1);
        printf("\n---------------------\n\n");
        if (p2) {
            printf("%s", p2);
        }
    }
    SraObjFree(sra);
}
