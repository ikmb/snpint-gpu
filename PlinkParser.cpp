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

static bool file_exists(const std::string &filename) {
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
}

static void bail_if_not_readable(const std::string &filename) {
    int ret = access(filename.c_str(), R_OK);
    if(ret == 0)
        return;

    throw system_error(errno, system_category(), filename);
}

PlinkParser::PlinkParser(const std::string &bim, const string &bed, const string &fam)
    : bimfile(bim), bedfile(bed), famfile(fam)
{}

void PlinkParser::parseFam(vector<SNPDB::SampleInfo> &samples) {

    if(!file_exists(famfile))
        throw runtime_error("Could not find a FAM file that belongs to the given dataset");

    bail_if_not_readable(famfile);

    ifstream fam(famfile);

    string line, tmp;
    unsigned long line_no = 0;
    int missing_pheno_count = 0;
    while(getline(fam, line)) {
        line_no++;

        SNPDB::SampleInfo info;
        istringstream s(line);

        if(! (s >> info.family >> info.within_family >> info.father >> info.mother >> tmp) )
            throw runtime_error("Malformed FAM file on line " + to_string(line_no));

        switch(std::stoi(tmp)) {
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

        if(!(s >> tmp))
            throw runtime_error("FAM file did not specify any phenotypes.");

        switch(std::stoi(tmp)) {
        case 1:
            info.pheno = SNPDB::Phenotype::Control;
            break;
        case 2:
            info.pheno = SNPDB::Phenotype::Case;
            break;
        default:
            info.pheno = SNPDB::Phenotype::Missing;
            break;
        }

        // count but ignore samples with missing phenotypes
        if (info.pheno == SNPDB::Phenotype::Missing)
            missing_pheno_count++;
        else
            samples.push_back(info);
    }
    if (missing_pheno_count)
        std::cerr << "\nWARNING! Ignored " << missing_pheno_count << " samples with missing phenotype after .fam parsing." << std::endl;
}

void PlinkParser::parseBim(vector<SNPDB::SNPInfo> &snps) {

    if(!file_exists(bimfile))
        throw runtime_error("Could not find BIM file.");
    ifstream bim(bimfile);

    string line;
    unsigned long line_no = 0;
    while(getline(bim, line)) {
        line_no++;

        SNPDB::SNPInfo info;
        istringstream s(line);

        if(! (s >> info.chromosome >> info.variant_id >> info.pos_cm >> info.pos_bp >> info.alleles[0] >> info.alleles[1]) )
            throw runtime_error("Malformed BIM file on line " + to_string(line_no));

        snps.push_back(info);
    }
}

void PlinkParser::parseBed(SNPDB& target) {

    /*
    Source format: https://www.cog-genomics.org/plink2/formats#bed
    Byte-wise storage, LSBs first, two bits per genotype. Cases and controls possibly unsorted
    */

    if(!file_exists(bedfile)) throw runtime_error("Could not open BED file");

    // set a suitably large read-ahead buffer and open bed file
    ifstream bed;
    bed.open(bedfile, std::ios_base::binary | std::ios_base::in);

    // check file magic
    char magic[3];
    bed.read(magic, 3);
    if(!(magic[0] == 0x6C && magic[1] == 0x1B && magic[2] == 0x01))
        throw runtime_error("Unsupported BED format (only SNP-major PLINK BED files can be processed. Consider rebuilding with a recent Plink version (1.9 or later))");

    // load bed file for first phenotype
    // reserve space for exactly one array of samples (i.e. one SNP)
    const streamsize samples_size = target.getSampleCount() / 4 + !!(target.getSampleCount() % 4);
    char *samples = new char[samples_size];

    // read bunches of samples (one SNP at a time) and use phenotype #0 for initial database setup
    unsigned sampleCount = target.getSampleCount();
    unsigned unrolled = (sampleCount / 4) * 4;

    static std::array<SNPDB::Genotype, 4> translateGenotypes {
        {SNPDB::Genotype::HomozygousVariant, // Plink 0
        SNPDB::Genotype::Invalid,            // Plink 1
        SNPDB::Genotype::Heterozygous,       // Plink 2
        SNPDB::Genotype::HomozygousWild}     // Plink 3
    };

    std::vector<bool> casephenos(sampleCount);
    for (unsigned i = 0; i < sampleCount; i++)
        casephenos[i] = target.getSampleInfo(i).pheno == SNPDB::Phenotype::Case;

    while(bed.read(samples, samples_size)) {

        auto curr_is_case = casephenos.begin();
        // iterate over samples
        for(unsigned i = 0; i < unrolled; /* no increment */) {
            unsigned current_gtblock = samples[i/4];
            unsigned current_gt = current_gtblock & 3;
            //SNPDB::Phenotype current_pheno = target.getSampleInfo(i).pheno;
            //target.addSample(translateGenotypes[current_gt], current_pheno);
            target.addSample(translateGenotypes[current_gt], *curr_is_case);

            ++i;
            curr_is_case++;
            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            //current_pheno = target.getSampleInfo(i).pheno;
            //target.addSample(translateGenotypes[current_gt], current_pheno);
            target.addSample(translateGenotypes[current_gt], *curr_is_case);

            ++i;
            curr_is_case++;
            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            //current_pheno = target.getSampleInfo(i).pheno;
            //target.addSample(translateGenotypes[current_gt], current_pheno);
            target.addSample(translateGenotypes[current_gt], *curr_is_case);

            ++i;
            curr_is_case++;
            current_gtblock >>= 2;
            current_gt = current_gtblock & 3;
            //current_pheno = target.getSampleInfo(i).pheno;
            //target.addSample(translateGenotypes[current_gt], current_pheno);
            target.addSample(translateGenotypes[current_gt], *curr_is_case);
            ++i;
            curr_is_case++;
        }

        for(unsigned i = unrolled; i < sampleCount; i++) {
            unsigned current_gt = (samples[i / 4] >> (i%4)*2) & 3;
            //SNPDB::Phenotype current_pheno = target.getSampleInfo(i).pheno;
            //target.addSample(translateGenotypes[current_gt], current_pheno);
            target.addSample(translateGenotypes[current_gt], *curr_is_case);
        }

        // advance to next SNP
        target.finishSNP();
    }

    delete[] samples;

}

void PlinkParser::parse(SNPDB &target, size_t snp_word_size) {
    unsigned long sample_count, cases = 0, snp_count;
    {
        std::vector<SNPDB::SampleInfo> samples;
        parseFam(samples);
        sample_count = samples.size();
        cases = count_if(begin(samples), end(samples), [](const SNPDB::SampleInfo &s) { return s.pheno == SNPDB::Phenotype::Case; });
        target.setSampleInfo(std::move(samples));
    }

    {
        std::vector<SNPDB::SNPInfo> snps;
        parseBim(snps);
        snp_count = snps.size();
        target.setSNPInfo(std::move(snps));
    }

    unsigned long controls = sample_count - cases; // works since we ignore missing phenotypes
    target.initialize(snp_count, cases, controls, snp_word_size);

    parseBed(target);
}
