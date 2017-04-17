#ifndef L1Trigger_L1TGlobal_GtDefinitions_h
#define L1Trigger_L1TGlobal_GtDefinitions_h

/**
 *
 *
 * Description: enums for the L1 GT.
 *
 * Implementation:
 *    Defines various enums for CondFormats L1 GT. For each enum, define the
 *    lightweight "maps" for enum string label and enum value
 *
 *
 * \author: Vladimir Rekovic,   Brian Winer, OSU   Vasile Mihai Ghete - HEPHY Vienna
 *
 * $Date$
 * $Revision$
 *
 */

// system include files
#include <string>

// user include files

namespace l1t {
/// board types in GT
enum L1GtBoardType {
    MP7,
    BoardNull
};

struct L1GtBoardTypeStringToEnum {
    const char* label;
    L1GtBoardType value;
};

L1GtBoardType l1GtBoardTypeStringToEnum(const std::string&);
std::string l1GtBoardTypeEnumToString(const L1GtBoardType&);


/// condition types
/// TypeNull:  null type - for condition constructor only
/// Type1s :   one particle
/// Type2s :   two particles, same type, no spatial correlations among them
/// Type2wsc : two particles, same type, with spatial correlations among them
/// Type2cor : two particles, different type, with spatial correlations among them
/// Type3s : three particles, same type
/// Type4s : four particles, same type
/// TypeETM, TypeETT, TypeHTT, TypeHTM  : ETM, ETT, HTT, HTM
/// TypeExternal: external conditions (logical result only; definition in L1 GT external systems)
/// Type1sWithOverlapRemoval: two particles - as Type1s plus one particle for overlap removal
/// Type2sWithOverlapRemoval: three particles - as Type2s plus one particle for overlap removal
/// Type3sWithOverlapRemoval: four particles - as Type3s plus one particle for overlap removal
/// Type4sWithOverlapRemoval: five particles - as Type4s plus one particle for overlap removal
/// Type2corWithOverlapRemoval: three partices - as Type2cor plus one particle for overlap removal
enum GtConditionType {
    TypeNull,
    Type1s,
    Type2s,
    Type2wsc,
    Type2cor,
    Type3s,
    Type4s,
    TypeETM,
    TypeETT,
    TypeHTT,
    TypeHTM,
    TypeETMHF,
    TypeTowerCount,
    TypeMinBiasHFP0,
    TypeMinBiasHFM0,
    TypeMinBiasHFP1,
    TypeMinBiasHFM1,
    TypeETTem,
    TypeExternal,
    Type1sWithOverlapRemoval,
    Type2sWithOverlapRemoval,
    Type3sWithOverlapRemoval,
    Type4sWithOverlapRemoval,
    Type2corWithOverlapRemoval

};

struct GtConditionTypeStringToEnum {
    const char* label;
    GtConditionType value;
};

GtConditionType l1GtConditionTypeStringToEnum(const std::string&);
std::string l1GtConditionTypeEnumToString(const GtConditionType&);

/// condition categories
enum GtConditionCategory {
    CondNull,
    CondMuon,
    CondCalo,
    CondEnergySum,
    CondCorrelation,
    CondExternal,
    CondMuonWithOverlapRemoval,
    CondCaloWithOverlapRemoval,
    CondCorrelationWithOverlapRemoval
};

struct GtConditionCategoryStringToEnum {
    const char* label;
    GtConditionCategory value;
};

GtConditionCategory l1GtConditionCategoryStringToEnum(const std::string&);
std::string l1GtConditionCategoryEnumToString(const GtConditionCategory&);

}
#endif
