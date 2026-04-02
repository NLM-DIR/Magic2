#include <stdio.h>
#include "sra_read.h"

int main(void)
{
    SRAReadBatch* sra = SRAReadBatchNew("SRR24511885");
    int num_bases = 500;

    while (1) {
        SraGetReadBatch(sra, num_bases, 1, 1);
        if (!sra->seq) {
            break;
        }
        printf("%s", sra->seq);
        printf("\n---------------------\n\n");
        if (sra->seq2) {
            printf("%s", sra->seq2);
        }

        break;
    }

    printf("Paired: %s\n", sra->is_paired ? "yes" : "no");
    printf("Num bases read: %d\n", sra->num_bases);

    SRAReadBatchFree(sra);
}
