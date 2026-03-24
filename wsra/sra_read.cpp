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
int SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss);
int SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
            stringstream& ss2);
int SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss);
int SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
                 stringstream& ss2);


struct SraObj {
    /* Read iterator used to downlaod batches of reads */
    ngs::ReadIterator it;
    /* Buffers that hold downloaded read sequences in FASTA or FASTQ format.
       There are two buffers for paired reads. */
    pair<string, string> buff;

    SraObj(const char* accession) : it(SraOpen(accession))
    {}
};


SRAObj* SraObjNew(const char* accession)
{
    SraObj* retval = new SraObj(accession);

    return retval;
}


SRAObj* SraObjFree(SRAObj* insra)
{
    SraObj* sra = static_cast<SraObj*>(insra);
    if (sra) {
	delete sra;
    }
    return nullptr;
}


int SraGetReadBatch(SRAObj* insra, long int num_bases, int quality_scores,
                    int split_spot, const char** buff1, const char** buff2)
{
    SraObj* sra = static_cast<SraObj*>(insra);
    if (!buff1) {
        return -1;
    }
    if (split_spot && !buff2) {
        return -1;
    }
    if (!sra) {
        *buff1 = nullptr;
        *buff2 = nullptr;
        return -1;
    }
    if (split_spot) {
        stringstream ss1;
        stringstream ss2;
        if (quality_scores) {
            SraReadFastq(sra->it, num_bases, ss1, ss2);
        }
        else {
            SraRead(sra->it, num_bases, ss1, ss2);
        }

        sra->buff.first = std::move(ss1.str());
        sra->buff.second = std::move(ss2.str());

        *buff1 = !sra->buff.first.empty() ? sra->buff.first.c_str() : nullptr;
        *buff2 = !sra->buff.second.empty() ? sra->buff.second.c_str() : nullptr;
    }
    else {
        stringstream ss;
        if (quality_scores) {
            SraReadFastq(sra->it, num_bases, ss);
        }
        else {
            SraRead(sra->it, num_bases, ss);
        }

        sra->buff.first = std::move(ss.str());
        *buff1 = !sra->buff.first.empty() ? sra->buff.first.c_str() : nullptr;
     }

    return 0;
}

ngs::ReadIterator SraOpen(const char* accession)
{
    ngs::ReadCollection run = ncbi::NGS::openReadCollection(accession);
    ngs::ReadIterator it = run.getReads(ngs::Read::all);

    return it;
}


/* Download SRA reads as FASTA, interleaved for paired reads */
int SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss)
{
    size_t num_bases = 0;
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
            }
            else {
                ss << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();
            }
        }
    }

    return 0;
}


/* Download SRA reads as FASTA, split paired reads */
int SraRead(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
            stringstream& ss2)
{
    size_t num_bases = 0;
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
            }
            else {
                ss1 << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss1 << bases << endl;
                num_bases += bases.length();
            }
        }
    }

    return 0;
}


/* Download SRA reads as FASTQ, interleaved for paired reads */
int SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss)
{
    size_t num_bases = 0;
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
            }
            else {
                ss << ">" << it.getReadId().data() << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();

                ss << "+" << it.getReadId().data() << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss << qualities << endl;
            }
        }
    }

    return 0;
}


/* Download SRA reads as FASTA, split paired reads */
int SraReadFastq(ngs::ReadIterator& it, size_t max_bases, stringstream& ss1,
                 stringstream& ss2)
{
    size_t num_bases = 0;
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

    return 0;
}
