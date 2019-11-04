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

#include <array>
#include <algorithm>
#include "Method.h"

using namespace std;


// 1. method type
// 2. order
// 3. score fields without decomposition
// 4. additional score fields if decomposition is used
// 5. short name
// 6. description
static const std::array<Method::MethodDescription, Method::TYPE_MAX> methods {{
    { Method::Type::IG2,     2, 1, 3, "ig2",     "Information Gain (2-way)"},
    { Method::Type::IG3,     3, 1, 7, "ig3",     "Information Gain (3-way)"},
    { Method::Type::MI2,     2, 1, 2, "mi2",     "Mutual Information (2-way)"},
    { Method::Type::MI3,     3, 1, 2, "mi3",     "Mutual Information (3-way)"},
    { Method::Type::BOOST,  2, 3, 2, "boost",  "BOOST (2-way)"},
    { Method::Type::LOGLIN, 2, 3, 1, "loglin", "Log-linear regression approximation (2-way)"},
    { Method::Type::LOGREG, 2, 4, 4, "logreg", "Logistic regression (2-way)"},
    { Method::Type::LOGREG3, 3, 4, 4, "logreg3", "Logistic regression (3-way) (EXPERIMENTAL)"},
    { Method::Type::INVALID, 0, 0, 0, "invalid", "(invalid)"},
}};

static const std::string SNPIDHEADER_2WAY = "SNPID A\tSNPID B\tIDX A\tIDX B";
static const std::string SNPIDHEADER_3WAY = "SNPID A\tSNPID B\tSNPID C\tIDX A\tIDX B\tIDX C";

static const std::array<std::string, Method::TYPE_MAX> SCOREHEADERS {
  "\tIG", // IG2
  "\tIG", // IG3
  "\tMI", // MI2
  "\tMI", // MI3
  "\tCHISQ\tERR\tP-VAL", // BOOST
  "\tCHISQ\tERR\tP-VAL", // LOGLIN
  "\tCHISQ\tOR\tP-VAL\tBETA", // LOGREG
  "\tCHISQ\tOR\tP-VAL\tBETA", // LOGREG3
  "" // INVALID
};

static const std::string SCOREHEADER_LD2WAYAPPEND = "\tR2_H\tR2_L";
static const std::string SCOREHEADER_LD3WAYAPPEND = "\tR2_AB\tR2_AC\tR2_BC";

static const std::array<std::string, Method::TYPE_MAX> SCOREHEADER_DECOMPAPPEND {
  "\tMI_AB\tMI_A\tMI_B", // IG2
  "\tMI_ABC\tMI_AB\tMI_AC\tMI_BC\tMI_A\tMI_B\tMI_C", // IG3
  "\tH_ABY\tH_AB", // MI2
  "\tH_ABCY\tH_ABC", // MI3
  "\t#ITER\tKSA\tKSASA", // BOOST
  "\t#ITER", // LOGLIN
  "\t#ITER\tDELTA\tMINDELTA\tRETURN", // LOGREG
  "\t#ITER\tDELTA\tMINDELTA\tRETURN", // LOGREG3
  "" // INVALID
};

const std::array<Method::MethodDescription, Method::Type::TYPE_MAX>& Method::getMethods() {
    return methods;
}

unsigned Method::getOrder() const {
    return methods[type].order;
}

unsigned Method::getScoreFields() const {
    unsigned ret =  methods[type].numScoreFields + (details? methods[type].numScoreComponents : 0);
    if (ld)
        ret += methods[type].order; // 2 LD fields for second order, 3 LD fields for 3rd order
    return ret;
}

const char *Method::getShortName() const {
    return methods[type].shortName;
}

const char *Method::getDescription() const {
    return methods[type].descriptiveName;
}

void Method::printHeader(ostream &out) const {
    if (getOrder() == 2)
        out << SNPIDHEADER_2WAY;
    else
        out << SNPIDHEADER_3WAY;
    out << SCOREHEADERS[type];
    if (ld) {
        if (getOrder() == 2)
            out << SCOREHEADER_LD2WAYAPPEND;
        else
            out << SCOREHEADER_LD3WAYAPPEND;
    }
    if (details)
        out << SCOREHEADER_DECOMPAPPEND[type];
    out << endl;
}

// stream parser for kernel enum
istream &operator>>(istream &in, Method &method) {

    string token;
    in >> token;

    auto it = find_if(begin(methods), end(methods), [&](const Method::MethodDescription& m) { return m.shortName == token; });
    if(it == end(methods))
        in.setstate(ios_base::failbit);
    else
        method = Method(it->type);

    return in;
}

ostream &operator<<(ostream &out, const Method &method) {
    out << method.getShortName();
    return out;
}

bool operator==(const Method& lhs, const Method& rhs) { return lhs.getType() == rhs.getType(); }
bool operator==(const Method& lhs, Method::Type t) { return lhs.getType() == t; }

