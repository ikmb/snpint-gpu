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

#ifndef SNPDB_H
#define SNPDB_H

#include <array>
#include <vector>

class SNPDB
{

public:
    using AlleleCounterType = uint16_t;

    // ATTENTION!!!
    // If you change the encoding, you need to change the Boolean condition in addSampleInt
    // to add a new genotype to the buffer correctly!
    enum Genotype {
        HomozygousWild = 0,
        HomozygousVariant = 2,
        Heterozygous = 1,
        Invalid = 3
    };
    const Genotype buffer_init_gt = Invalid;
    const unsigned char buffer_init_value =
            (buffer_init_gt << 6) |
            (buffer_init_gt << 4) |
            (buffer_init_gt << 2) |
             buffer_init_gt;

    enum SexCode {
        Male,
        Female,
        Unknown
    };

    enum Phenotype {
        Control,
        Case,
        Missing
    };

    struct SampleInfo {
        std::string family;
        std::string within_family;
        std::string father;
        std::string mother;
        SexCode sex;
        Phenotype pheno;
    };

    struct SNPInfo {
        std::string chromosome;
        std::string variant_id;
        double pos_cm;
        unsigned long pos_bp;
        std::array<std::string, 2> alleles;
    };

    struct SHMHeader {
        std::uint64_t genotype_size;
        std::uint64_t case_count;
        std::uint64_t control_count;
        std::uint64_t snp_count;
    };

    template<typename T>
    class AlleleInfo {
    public:

        AlleleInfo() : counters{0,0,0,0,0,0} {}
        explicit AlleleInfo(T *counters_) {
            memcpy(this->counters, counters_, sizeof(T)*6);
        }

        inline void inc(bool is_case, SNPDB::Genotype gt, T val = 1) {
        	if (gt != SNPDB::Genotype::Invalid)
        		counters[gt + (is_case? 3 : 0)] += val;
        }

        void import(T *source) {
            memcpy(source, counters, 6 * sizeof(T));
        }

        const T *get() const {
            return counters;
        }


    private:
        T counters[6];
    };

    SNPDB(unsigned minCases = 0, unsigned minControls = 0);

    SNPDB(const SNPDB&) = delete; // copy c'tor
    SNPDB& operator=(const SNPDB&) = delete; // copy assignment

    SNPDB(SNPDB&&) = default; // move c'tor
    SNPDB& operator=(SNPDB&&) & = default; // move assignment

    ~SNPDB();
    //void initialize(unsigned long num_cases, unsigned long num_controls, size_t word_size);
    void initialize(unsigned long num_snps, unsigned long num_cases, unsigned long num_controls, size_t word_size);

    void setSNPInfo(std::vector<SNPInfo> &&snpinfo) {
        snp_info = snpinfo;
    }

    void setSampleInfo(std::vector<SampleInfo> &&sampleinfo_) {
        sample_info = sampleinfo_;
    }

    unsigned char *data() { return buffer; }
    const unsigned char *data() const { return buffer; }

    unsigned long getCaseCount() const { return num_cases; }
    unsigned long getControlCount() const { return num_controls; }
    unsigned long getSampleCount() const { return num_cases + num_controls; }
    unsigned long getCaseCountPadded() const { return num_cases_padded; }
    unsigned long getControlCountPadded() const { return num_controls_padded; }
    unsigned long getSampleCountPadded() const { return num_cases_padded + num_controls_padded; }
    unsigned long getSNPCount() const { return num_snps; }

    const std::vector<SampleInfo> &getSampleInfos() const { return sample_info; }
    const SampleInfo &getSampleInfo(unsigned index) const { return sample_info[index]; }
    const SNPInfo &getSNPInfo(unsigned index) const { return snp_info[index]; }

    size_t getSNPSize() const { return snp_size; }
    size_t getBufferSize() const { return buffer_size; }

    inline void addSample(SNPDB::Genotype genotype, bool is_case) {
        //current_snp_alleles.inc(is_case, genotype);
        addSampleInt(genotype, (is_case ? current_case_offset : current_ctrl_offset));
    }

    void finishSNP();

    const unsigned char * operator[](size_t snp_index) const {
        return buffer + (snp_size * snp_index);
    }

private:
    std::vector<SNPInfo> snp_info;
    std::vector<SampleInfo> sample_info;
    //std::vector<AlleleInfo<AlleleCounterType> > allele_info;

    unsigned long num_snps;
    unsigned long num_cases;
    unsigned long num_controls;
    unsigned long num_cases_padded;
    unsigned long num_controls_padded;
    size_t word_size;
    size_t snp_size;
    size_t case_size;
    typedef struct { unsigned char *bufptr; size_t bit;} genotype_offset;
    size_t current_snp;
    genotype_offset current_case_offset;
    genotype_offset current_ctrl_offset;
    //AlleleInfo<AlleleCounterType> current_snp_alleles;
    bool allow_buffer_resize;
    size_t buffer_size;

    const size_t snp_batch_size = 1024; // allocate this number of SNPs at once
    unsigned char *buffer;
    int buffer_fd = -1;
    unsigned long min_cases;
    unsigned long min_controls;

    inline void addSampleInt(SNPDB::Genotype genotype, genotype_offset& current_offset) {

        // this requires each cell to be pre-initialized with 0xFF (i.e. an invalid genotype)
        // ATTENTION! If you change the encoding of an invalid genotype, this need to be changed as well!
        *(current_offset.bufptr) ^= ((3 ^ genotype) << current_offset.bit);

        // byte advance?
        current_offset.bit += 2;
        if(current_offset.bit == 8) {
            current_offset.bit = 0;
            ++(current_offset.bufptr);
        }
    }

    inline void resetCurrOffsets() {
        // set case offset
        current_case_offset.bufptr = buffer;
        current_case_offset.bit = 0;
        // set ctrl offset
        current_ctrl_offset.bufptr = buffer + case_size;
        current_ctrl_offset.bit = 0;
    }

};

#endif // SNPDB_H
