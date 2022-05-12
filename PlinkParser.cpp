//    Copyright 2019 Lars Wienbrandt, Jan Christian Kässens
//
//    This file is part of SNPInt-GPU.
//
//    SNPInt-GPU is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SNPInt-GPU is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with SNPInt-GPU.  If not, see <https://www.gnu.org/licenses/>.

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include "PlinkParser.h"
#include "SNPDB.h"

using namespace std;

static bool file_exists(const string &filename) {
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}

static void bail_if_not_readable(const string &filename) {
    int ret = access(filename.c_str(), R_OK);
    if(ret == 0)
        return;

    throw system_error(errno, system_category(), filename);
}

PlinkParser::PlinkParser(const string &bim, const string &bed, const string &fam)
    : bimfile(bim), bedfile(bed), famfile(fam)
{}

void PlinkParser::parseFam(vector<SNPDB::SampleInfo> &samples, vector<bool> &ignored, size_t &numcases) {

    if(!file_exists(famfile))
        throw runtime_error("Could not find a FAM file that belongs to the given dataset");

    bail_if_not_readable(famfile);

    ifstream fam(famfile);

    string line, sexpheno;
    size_t line_no = 0;
    size_t missing_pheno_count = 0;
    numcases = 0;
    while(getline(fam, line)) {
        line_no++;

        SNPDB::SampleInfo info;
        istringstream s(line);

        if(! (s >> info.family >> info.within_family >> info.father >> info.mother >> sexpheno) )
            throw runtime_error("Malformed FAM file on line " + to_string(line_no));

        switch(stoi(sexpheno)) {
        case 1:
            info.sex = SNPDB::SexCode::Male;
            break;
        case 2:
            info.sex = SNPDB::SexCode::Female;
            break;
        default:
            info.sex = SNPDB::SexCode::Unknown;
            break;
        }

        if(!(s >> sexpheno))
            throw runtime_error("FAM file did not specify any phenotypes.");

        switch(stoi(sexpheno)) {
        case 1:
            info.pheno = SNPDB::Phenotype::Control;
            break;
        case 2:
            info.pheno = SNPDB::Phenotype::Case;
            numcases++;
            break;
        default:
            info.pheno = SNPDB::Phenotype::Missing;
            break;
        }

        // count but ignore samples with missing phenotypes
        if (info.pheno == SNPDB::Phenotype::Missing) {
            missing_pheno_count++;
            ignored.push_back(true);
        } else {
            samples.push_back(info);
            ignored.push_back(false);
        }
    }
    if (missing_pheno_count)
        cerr << "\nWARNING! Ignored " << missing_pheno_count << " samples with missing phenotype." << endl;
}

void PlinkParser::parseBim(vector<SNPDB::SNPInfo> &snps) {

    if(!file_exists(bimfile))
        throw runtime_error("Could not find BIM file.");
    ifstream bim(bimfile);

    string line;
    size_t line_no = 0;
    while(getline(bim, line)) {
        line_no++;

        SNPDB::SNPInfo info;
        istringstream s(line);

        if(! (s >> info.chromosome >> info.variant_id >> info.pos_cm >> info.pos_bp >> info.alleles[0] >> info.alleles[1]) )
            throw runtime_error("Malformed BIM file on line " + to_string(line_no));

        snps.push_back(info);
    }
}

void PlinkParser::parseBed(SNPDB& db, const vector<bool> &ignored_samples) {

    /*
    Source format: https://www.cog-genomics.org/plink2/formats#bed
    Byte-wise storage, LSBs first, two bits per genotype. Cases and controls possibly unsorted
    */

    if(!file_exists(bedfile)) throw runtime_error("Could not open BED file");

    // set a suitably large read-ahead buffer and open bed file
    ifstream bed;
    bed.open(bedfile, ios_base::binary | ios_base::in);

    // check file magic
    char magic[3];
    bed.read(magic, 3);
    if(!(magic[0] == 0x6C && magic[1] == 0x1B && magic[2] == 0x01))
        throw runtime_error("Unsupported BED format (only SNP-major PLINK BED files can be processed. Consider rebuilding with a recent Plink version (1.9 or later))");

    // load bed file for first phenotype
    // reserve space for exactly one array of samples (i.e. one SNP)
    size_t samplecnt = ignored_samples.size(); // this reflects the size of all samples including ignored ones
    streamsize samples_size = samplecnt / 4 + (samplecnt % 4 ? 1 : 0);
    char *gts = new char[samples_size];

    static array<SNPDB::Genotype, 4> translateGenotypes {
        {SNPDB::Genotype::HomozygousVariant, // Plink 0
        SNPDB::Genotype::Invalid,            // Plink 1
        SNPDB::Genotype::Heterozygous,       // Plink 2
        SNPDB::Genotype::HomozygousWild}     // Plink 3
    };

    vector<bool> casephenos(db.getSampleCount()); // only for not-ignored samples!
    for (unsigned i = 0; i < db.getSampleCount(); i++)
        casephenos[i] = db.getSampleInfo(i).pheno == SNPDB::Phenotype::Case;

    // read bunches of samples (one SNP at a time) and use phenotype #0 for initial database setup
    unsigned unrolled = (samplecnt / 4) * 4;
    while(bed.read(gts, samples_size)) {

        auto curr_is_case = casephenos.begin();
        // iterate over samples
        for(unsigned i = 0; i < unrolled; /* no increment */) { // much faster!!!
            unsigned current_gtblock = gts[i/4];
            unsigned current_gt = current_gtblock & 3;
            if (!ignored_samples[i])
                db.addSample(translateGenotypes[current_gt], *curr_is_case++);
            i++;

            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            if (!ignored_samples[i])
                db.addSample(translateGenotypes[current_gt], *curr_is_case++);
            i++;

            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            if (!ignored_samples[i])
                db.addSample(translateGenotypes[current_gt], *curr_is_case++);
            i++;

            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            if (!ignored_samples[i])
                db.addSample(translateGenotypes[current_gt], *curr_is_case++);
            i++;
        }

        for(unsigned i = unrolled; i < samplecnt; i++) {
            unsigned current_gt = (gts[i / 4] >> (i%4)*2) & 3;
            if (!ignored_samples[i])
                db.addSample(translateGenotypes[current_gt], *curr_is_case++);
        }

        // advance to next SNP
        db.finishSNP();
    }

    delete[] gts;

}

void PlinkParser::parse(SNPDB &db, size_t snp_word_size) {
    size_t sample_count, cases = 0, snp_count;
    vector<bool> ignored_samples;
    {
        vector<SNPDB::SampleInfo> samples;
        parseFam(samples, ignored_samples, cases);
        sample_count = samples.size();
        db.setSampleInfo(move(samples));
    }

    {
        vector<SNPDB::SNPInfo> snps;
        parseBim(snps);
        snp_count = snps.size();
        db.setSNPInfo(move(snps));
    }

    size_t controls = sample_count - cases; // works since we ignore missing phenotypes
    db.initialize(snp_count, cases, controls, snp_word_size);

    parseBed(db, ignored_samples);
}
