#include <NGS.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>
#include <sstream>
#include "sra_read.h"

using namespace std;

ngs::ReadIterator SraOpen(const char* accession);
int SraRead(ngs::ReadIterator& it, int max_bases, stringstream& ss);
int SraReadFastq(ngs::ReadIterator& it, int max_bases, stringstream& ss);


struct SraObj {
    /* Read iterator used to downlaod batches of reads */
    ngs::ReadIterator it;
    /* Buffer that holds downloaded read sequences in FASTA format */
    string buff;

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


const char* SraGetReadBatch(SRAObj* insra, int num_bases, int quality_scores)
{
    SraObj* sra = static_cast<SraObj*>(insra);
    if (!sra) {
	return nullptr;
    }
    stringstream ss;
    if (quality_scores) {
        SraReadFastq(sra->it, num_bases, ss);
    }
    else {
        SraRead(sra->it, num_bases, ss);
    }
    sra->buff = std::move(ss.str());
    return (!sra->buff.empty() ? sra->buff.c_str() : nullptr);
}

ngs::ReadIterator SraOpen(const char* accession)
{
    ngs::ReadCollection run = ncbi::NGS::openReadCollection(accession);
    ngs::ReadIterator it = run.getReads(ngs::Read::all);

    return it;
}


int SraRead(ngs::ReadIterator& it, int max_bases, stringstream& ss)
{
    size_t num_bases = 0;
    while (num_bases < max_bases && it.nextRead()) {
	if (it.nextFragment()) {
	    if (it.isPaired()) {
		string read_id = it.getReadId().toString();
		ss << ">" << read_id << ".1" << endl;
		string bases(std::move(it.getFragmentBases().toString()));
		ss << bases << endl;
		num_bases += bases.length();

		ss << ">" << read_id << ".2" << endl;
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


int SraReadFastq(ngs::ReadIterator& it, int max_bases, stringstream& ss)
{
    size_t num_bases = 0;
    while (num_bases < max_bases && it.nextRead()) {
        if (it.nextFragment()) {
            if (it.isPaired()) {
                string read_id = it.getReadId().toString();
                ss << "@" << read_id << ".1" << endl;
                string bases(std::move(it.getFragmentBases().toString()));
                ss << bases << endl;
                num_bases += bases.length();

                ss << "+" << read_id << ".1" << endl;
                string qualities(std::move(it.getFragmentQualities().toString()));
                ss << qualities << endl;

                ss << "@" << read_id << ".2" << endl;
                if (it.nextFragment()) {
                    string bases(std::move(it.getFragmentBases().toString()));
                    ss << bases << endl;
                    num_bases += bases.length();

                    ss << "+" << read_id << ".1" << endl;
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
