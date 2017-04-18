/**
 * \class CaloWithOverlapRemovalTemplate
 *
 *
 * Description: L1 Global Trigger calo with overlap removal template.
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *
 * \author: Vladimir Rekovic 
 *
 * $Date$
 * $Revision$
 *
 */

// this class header
#include "L1Trigger/L1TGlobal/interface/CaloWithOverlapRemovalTemplate.h"

// system include files

#include <iostream>
#include <iomanip>

// user include files

//   base class


// forward declarations

// constructors
CaloWithOverlapRemovalTemplate::CaloWithOverlapRemovalTemplate()
        : GlobalCondition()
{

    m_condCategory = l1t::CondCaloWithOverlapRemoval;

    m_cond0Category = l1t::CondNull;
    m_condOvrCategory = l1t::CondNull;
    m_cond0Index = -1;
    m_condOvrIndex = -1;


}

CaloWithOverlapRemovalTemplate::CaloWithOverlapRemovalTemplate(const std::string& cName)
        : GlobalCondition(cName)
{

    m_condCategory = l1t::CondCaloWithOverlapRemoval;

    m_cond0Category = l1t::CondNull;
    m_condOvrCategory = l1t::CondNull;
    m_cond0Index = -1;
    m_condOvrIndex = -1;

}

CaloWithOverlapRemovalTemplate::CaloWithOverlapRemovalTemplate(const std::string& cName, const l1t::GtConditionType& cType)
        : GlobalCondition(cName, l1t::CondCaloWithOverlapRemoval, cType)
{

    int nObjects = nrObjects();

    if (nObjects > 0) {
        m_objectParameter.reserve(nObjects);

        m_objectType.reserve(nObjects);
    }

}

// copy constructor
CaloWithOverlapRemovalTemplate::CaloWithOverlapRemovalTemplate(const CaloWithOverlapRemovalTemplate& cp)
        : GlobalCondition(cp.m_condName)
{
    copy(cp);
}

// destructor
CaloWithOverlapRemovalTemplate::~CaloWithOverlapRemovalTemplate()
{
    // empty now
}

// assign operator
CaloWithOverlapRemovalTemplate& CaloWithOverlapRemovalTemplate::operator= (const CaloWithOverlapRemovalTemplate& cp)
{

    copy(cp);
    return *this;
}

// set the category of the mian condition and overlap-removal condition
void CaloWithOverlapRemovalTemplate::setCond0Category(
        const l1t::GtConditionCategory& condCateg) {

    m_cond0Category = condCateg;
}

void CaloWithOverlapRemovalTemplate::setCondOvrCategory(
        const l1t::GtConditionCategory& condCateg) {

    m_condOvrCategory = condCateg;
}


// set the index of the two sub-conditions in the cor* vector from menu
void CaloWithOverlapRemovalTemplate::setCond0Index(const int& condIndex) {
    m_cond0Index = condIndex;
}

void CaloWithOverlapRemovalTemplate::setCondOvrIndex(const int& condIndex) {
    m_condOvrIndex = condIndex;
}

// setConditionParameter - set the parameters of the condition
void CaloWithOverlapRemovalTemplate::setConditionWithOverlapRemovalParameter(
    const std::vector<ObjectParameter>& objParameter,
    const OverlapRemovalParameter& corrParameter)
{

    m_objectParameter = objParameter;
    m_overlapRemovalParameter = corrParameter;

}

void CaloWithOverlapRemovalTemplate::print(std::ostream& myCout) const
{

    myCout << "\n  CaloWithOverlapRemovalTemplate print..." << std::endl;

    GlobalCondition::print(myCout);

    int nObjects = nrObjects();

    for (int i = 0; i < nObjects; i++) {
        myCout << std::endl;
        myCout << "  Template for object " << i << " [ hex ]" << std::endl;
        myCout << "    etThreshold       = "
        << std::hex << m_objectParameter[i].etLowThreshold << "  " << m_objectParameter[i].etHighThreshold << std::endl;
        myCout << "    indexLow       = "
        << std::hex << m_objectParameter[i].indexLow << std::endl;
        myCout << "    indexHigh      = "
        << std::hex << m_objectParameter[i].indexHigh << std::endl;
        myCout << "    etaRange          = "
        << std::hex << m_objectParameter[i].etaRange << std::endl;
        myCout << "    phiRange          = "
        << std::hex << m_objectParameter[i].phiRange << std::endl;
        myCout << "    isolationLUT      = "
        << std::hex << m_objectParameter[i].isolationLUT << std::endl;
        myCout << "    qualityLUT      = "
        << std::hex << m_objectParameter[i].qualityLUT << std::endl;	
    }

    myCout << "\n  First sub-condition category:  " << m_cond0Category <<  std::endl;
    myCout <<   "  Overlap sub-condition category: " << m_condOvrCategory <<  std::endl;

    myCout << "\n  First sub-condition index:  " << m_cond0Index <<  std::endl;
    myCout <<   "  Overlap sub-condition index: " << m_condOvrIndex <<  std::endl;

    myCout << "\n  OverlapRemoval parameters " << "[ hex ]" <<  std::endl;

    myCout << "    Cut Type:  " << m_overlapRemovalParameter.overlapCutType << std::endl;
    myCout << "    minOverlapRemovalEtaCutValue        = " << std::dec << m_overlapRemovalParameter.minOverlapRemovalEtaCutValue << std::endl;
    myCout << "    maxOverlapRemovalEtaCutValue        = " << std::dec << m_overlapRemovalParameter.maxOverlapRemovalEtaCutValue << std::endl;
    myCout << "    precOverlapRemovalEtaCut            = " << std::dec << m_overlapRemovalParameter.precOverlapRemovalEtaCut     << std::endl;
    myCout << "    minOverlapRemovalPhiCutValue        = " << std::dec << m_overlapRemovalParameter.minOverlapRemovalPhiCutValue << std::endl;
    myCout << "    maxOverlapRemovalPhiCutValue        = " << std::dec << m_overlapRemovalParameter.maxOverlapRemovalPhiCutValue << std::endl;
    myCout << "    precOverlapRemovalPhiCut            = " << std::dec << m_overlapRemovalParameter.precOverlapRemovalPhiCut     << std::endl;
    myCout << "    minOverlapRemovalDRCutValue         = " << std::dec << m_overlapRemovalParameter.minOverlapRemovalDRCutValue  << std::endl;
    myCout << "    maxOverlapRemovalDRCutValue         = " << std::dec << m_overlapRemovalParameter.maxOverlapRemovalDRCutValue  << std::endl;
    myCout << "    precOverlapRemovalDRCut             = " << std::dec << m_overlapRemovalParameter.precOverlapRemovalDRCut      << std::endl;
 

    // reset to decimal output
    myCout << std::dec << std::endl;

}

void CaloWithOverlapRemovalTemplate::copy(const CaloWithOverlapRemovalTemplate& cp)
{

    m_condName     = cp.condName();
    m_condCategory = cp.condCategory();
    m_condType     = cp.condType();
    m_objectType   = cp.objectType();
    m_condGEq      = cp.condGEq();
    m_condChipNr   = cp.condChipNr();
    m_condRelativeBx = cp.condRelativeBx();

    m_cond0Category = cp.cond0Category();
    m_condOvrCategory = cp.condOvrCategory();
    m_cond0Index = cp.cond0Index();
    m_condOvrIndex = cp.condOvrIndex();

    m_objectParameter = *(cp.objectParameter());
    m_overlapRemovalParameter = *(cp.overlapRemovalParameter());

}

// output stream operator
std::ostream& operator<<(std::ostream& os, const CaloWithOverlapRemovalTemplate& result)
{
    result.print(os);
    return os;

}



