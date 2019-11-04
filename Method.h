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

#ifndef METHOD_H
#define METHOD_H

#include <istream>
#include <ostream>
#include <string>
#include <array>
#include <utility>

// 0-based, 1st including, 2nd excluding
typedef std::pair<size_t,size_t> snprange;

class Method {
public:
    enum Type {
        IG2,
        IG3,
        MI2,
        MI3,
        BOOST,
        LOGLIN,
        LOGREG,
        LOGREG3,
        INVALID,
        TYPE_MAX
    };

    struct MethodDescription {
        Method::Type type;              ///< Type ID, see Method::Type
        unsigned order;                 ///< Order of interaction (also: number of ID fields in result set)
        unsigned numScoreFields;        ///< Number of score fields in result set
        unsigned numScoreComponents;    ///< Number of *additional* score fields when decomposed results are requested
        const char *shortName;          ///< A short identifier used for command-line args parsing and file extensions
        const char *descriptiveName;    ///< A description of the method in use, for help output
    };

    static const int maxScoreComponents = 11; // including LD fields!

    static const std::array<MethodDescription, TYPE_MAX>& getMethods();

    Method() : type(INVALID), ld(false), details(false) {}
    explicit Method(Type t) : type(t), ld(false), details(false) {}
    unsigned getOrder() const;
    unsigned getScoreFields() const;
    const char *getShortName() const;
    const char *getDescription() const;
    Type getType() const { return type; }

    void setLDAndDetails(bool ld_, bool details_) { ld = ld_; details = details_; }
    bool isLDIncluded() const { return ld; }
    bool isDetailedComputation() const { return details; }

    void printHeader(std::ostream &out) const;

private:
    Type type;
    bool ld;
    bool details;
};

std::istream &operator>>(std::istream &in, Method &method);
std::ostream &operator<<(std::ostream &out, const Method &method);
bool operator==(const Method& lhs, const Method& rhs);
bool operator==(const Method& lhs, Method::Type t);


#endif // METHOD_H
