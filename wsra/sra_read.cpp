#include <stdlib.h>
#include <NGS.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>
#include <sstream>
#include <utility>
#include "sra_read.h"

using namespace std;

ngs::ReadIterator SraOpen(const char* accession);
size_t SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss,
               bool& is_paired);
size_t SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
               stringstream& ss2, bool& is_paired);
size_t SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss,
                    bool& is_paired);
size_t SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
                    stringstream& ss2, bool& is_paired);


struct SraObj {
    /* Read iterator used to downlaod batches of reads */
    ngs::ReadIterator it;
    /* Buffers that hold downloaded read sequences in FASTA or FASTQ format.
       There are two buffers for paired reads. */
    pair<string, string> buff;

    SraObj(const char* accession) : it(SraOpen(accession))
    {}
};


SRAReadBatch* SRAReadBatchNew(const char* accession)
{
    SRAReadBatch* retval = new SRAReadBatch();
    if (!retval) {
        return nullptr;
    }
    retval->sra_obj = static_cast<SRAObj*>(new SraObj(accession));
    if (!retval->sra_obj) {
        delete retval;
        return nullptr;
    }
    return retval;
}

SRAReadBatch* SRAReadBatchFree(SRAReadBatch* sra)
{
    if (sra) {
        if (sra->sra_obj) {
            SraObj* obj = static_cast<SraObj*>(sra->sra_obj);
            delete obj;
        }
        delete sra;
    }
    return nullptr;
}



int SraGetReadBatch(SRAReadBatch* sra, long int num_bases, int quality_scores,
                    int split_spot)
{
    if (!sra) {
        return -1;
    }
    size_t num_bases_read = 0;
    bool is_paired;
    SraObj* sra_obj = static_cast<SraObj*>(sra->sra_obj);

    sra->seq = nullptr;
    sra->seq2 = nullptr;
    sra->size = static_cast<long unsigned int>(0);
    sra->size2 = static_cast<long unsigned int>(0);

    if (split_spot) {
        stringstream ss1;
        stringstream ss2;
        if (quality_scores) {
            num_bases_read = SraReadFastq(sra_obj->it, num_bases, ss1, ss2,
                                          is_paired);
        }
        else {
            num_bases_read = SraRead(sra_obj->it, num_bases, ss1, ss2, is_paired);
        }

        sra_obj->buff.first = std::move(ss1.str());
        sra_obj->buff.second = std::move(ss2.str());

        if (!sra_obj->buff.first.empty()) {
            sra->seq = sra_obj->buff.first.c_str();
            sra->size = static_cast<unsigned long int>(sra_obj->buff.first.size());
        }
        if (!sra_obj->buff.second.empty()) {
            sra->seq2 = sra_obj->buff.second.c_str();
            sra->size2 = static_cast<unsigned long int>(sra_obj->buff.second.size());
        }
    }
    else {
        stringstream ss;
        if (quality_scores) {
            num_bases_read = SraReadFastq(sra_obj->it, num_bases, ss, is_paired);
        }
        else {
            num_bases_read = SraRead(sra_obj->it, num_bases, ss, is_paired);
        }

        sra_obj->buff.first = std::move(ss.str());

        if (!sra_obj->buff.first.empty()) {
            sra->seq = sra_obj->buff.first.c_str();
            sra->size = static_cast<unsigned long int>(sra_obj->buff.first.size());
        }
     }

    sra->num_bases = num_bases_read;
    sra->is_paired = is_paired ? 1 : 0;

    return 0;
}

ngs::ReadIterator SraOpen(const char* accession)
{
    ngs::ReadCollection run = ncbi::NGS::openReadCollection(accession);
    ngs::ReadIterator it = run.getReads(ngs::Read::all);

    return it;
}


/* Download SRA reads as FASTA, interleaved for paired reads */
size_t SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss,
                 bool& is_paired)
{
    size_t num_bases = 0;
    is_paired = false;
    while (num_bases < max_bases && it.nextRead()) {
        if (it.nextFragment()) {
            if (it.isPaired()) {
                string read_id = it.getReadId().toString();
                ss << ">" << read_id << "/1" << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();

                ss << ">" << read_id << "/2" << endl;
                if (it.nextFragment()) {
                    string bases(std::move(it.getFragmentBases().toString()));
                    ss << bases << endl;
                    num_bases += bases.length();
                }

                is_paired = true;
            }
            else {
                ss << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();
            }
        }
    }

    return num_bases;
}


/* Download SRA reads as FASTA, split paired reads */
size_t SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
                 stringstream& ss2, bool& is_paired)
{
    size_t num_bases = 0;
    is_paired = false;
    while (num_bases < max_bases && it.nextRead()) {
        if (it.nextFragment()) {
            if (it.isPaired()) {
                string read_id = it.getReadId().toString();
                ss1 << ">" << read_id << "/1" << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss1 << bases << endl;
                num_bases += bases.length();

                ss2 << ">" << read_id << "/2" << endl;
                if (it.nextFragment()) {
                    string bases(std::move(it.getFragmentBases().toString()));
                    ss2 << bases << endl;
                    num_bases += bases.length();
                }

                is_paired = true;
            }
            else {
                ss1 << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss1 << bases << endl;
                num_bases += bases.length();
            }
        }
    }

    return num_bases;
}


/* Download SRA reads as FASTQ, interleaved for paired reads */
size_t SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss,
                    bool& is_paired)
{
    size_t num_bases = 0;
    is_paired = false;
    while (num_bases < max_bases && it.nextRead()) {
        if (it.nextFragment()) {
            if (it.isPaired()) {
                string read_id = it.getReadId().toString();
                ss << "@" << read_id << "/1" << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();

                ss << "+" << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss << qualities << endl;

                ss << "@" << read_id << "/2" << endl;
                if (it.nextFragment()) {
                    string bases(std::move(it.getFragmentBases().toString()));
                    ss << bases << endl;
                    num_bases += bases.length();

                    ss << "+" << endl;
                    string qualities(std::move(it.getFragmentQualities().toString()));
                    ss << qualities << endl;
                }

                is_paired = true;
            }
            else {
                ss << "@" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();

                ss << "+" << it.getReadId().data() << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss << qualities << endl;
            }
        }
    }

    return num_bases;
}


/* Download SRA reads as FASTA, split paired reads */
size_t SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
                    stringstream& ss2, bool& is_paired)
{
    size_t num_bases = 0;
    is_paired = false;
    while (num_bases < max_bases && it.nextRead()) {
        if (it.nextFragment()) {
            if (it.isPaired()) {
                string read_id = it.getReadId().toString();
                ss1 << "@" << read_id << "/1" << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss1 << bases << endl;
                num_bases += bases.length();

                ss1 << "+" << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss1 << qualities << endl;

                ss2 << "@" << read_id << "/2" << endl;
                if (it.nextFragment()) {
                    string bases(std::move(it.getFragmentBases().toString()));
                    ss2 << bases << endl;
                    num_bases += bases.length();

                    ss2 << "+" << endl;
                    string qualities(std::move(it.getFragmentQualities().toString()));
                    ss2 << qualities << endl;
                }

                is_paired = true;
            }
            else {
                ss1 << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss1 << bases << endl;
                num_bases += bases.length();

                ss1 << "+" << it.getReadId().data() << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss1 << qualities << endl;
            }
        }
    }

    return num_bases;
}
