#ifndef L1Trigger_L1TGlobal_CaloWithOverlapRemovalTemplate_h
#define L1Trigger_L1TGlobal_CaloWithOverlapRemovalTemplate_h

/**
 * \class CaloWithOverlapRemovalTemplate
 *
 *
 * Description: L1 Global Trigger calo template with overlap removal
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

// system include files
#include <string>
#include <iosfwd>

// user include files

//   base class
#include "L1Trigger/L1TGlobal/interface/GlobalCondition.h"

// forward declarations

// class declaration
class CaloWithOverlapRemovalTemplate : public GlobalCondition
{

public:

    // constructor
    CaloWithOverlapRemovalTemplate();

    // constructor
    CaloWithOverlapRemovalTemplate(const std::string& );

    // constructor
    CaloWithOverlapRemovalTemplate(const std::string&, const l1t::GtConditionType& );

    // copy constructor
    CaloWithOverlapRemovalTemplate( const CaloWithOverlapRemovalTemplate& );

    // destructor
    virtual ~CaloWithOverlapRemovalTemplate();

    // assign operator
    CaloWithOverlapRemovalTemplate& operator= (const CaloWithOverlapRemovalTemplate&);

public:

    /// typedef for a single object template
    struct ObjectParameter
    {
      unsigned int etLowThreshold;
      unsigned int etHighThreshold;
      unsigned int indexLow;
      unsigned int indexHigh;
      unsigned int etaRange;
      unsigned int phiRange;

      unsigned int isolationLUT;
      unsigned int qualityLUT;     

      unsigned int etaWindow1Lower;
      unsigned int etaWindow1Upper;
      unsigned int etaWindow2Lower;
      unsigned int etaWindow2Upper;

      unsigned int phiWindow1Lower;
      unsigned int phiWindow1Upper;
      unsigned int phiWindow2Lower;
      unsigned int phiWindow2Upper;

    };

    /// typedef for correlation parameters
    struct OverlapRemovalParameter
    {
        long long minOverlapRemovalEtaCutValue;
        long long maxOverlapRemovalEtaCutValue; 
        unsigned int precOverlapRemovalEtaCut;
        
        long long minOverlapRemovalPhiCutValue;
        long long maxOverlapRemovalPhiCutValue; 
        unsigned int precOverlapRemovalPhiCut;
        
        long long minOverlapRemovalDRCutValue;
        long long maxOverlapRemovalDRCutValue;
        unsigned int precOverlapRemovalDRCut; 

	int overlapCutType;

    };


public:


    /// get / set the category of the two sub-conditions
    //  one for primary object, and one for overlap removal object
    inline const l1t::GtConditionCategory cond0Category() const {
        return m_cond0Category;
    }

    inline const l1t::GtConditionCategory condOvrCategory() const {
        return m_condOvrCategory;
    }

    void setCond0Category(const l1t::GtConditionCategory&);
    void setCondOvrCategory(const l1t::GtConditionCategory&);

    /// get / set the index of the two sub-conditions (main + overla-removal) in the vector from menu
    inline const int cond0Index() const {
        return m_cond0Index;
    }

    inline const int condOvrIndex() const {
        return m_condOvrIndex;
    }

    void setCond0Index(const int&);
    void setCondOvrIndex(const int&);

    inline const std::vector<ObjectParameter>* objectParameter() const
    {
        return &m_objectParameter;
    }

    inline const OverlapRemovalParameter* overlapRemovalParameter() const
    {
        return &m_overlapRemovalParameter;
    }


    /// set functions
    void setConditionWithOverlapRemovalParameter(const std::vector<ObjectParameter>& objParameter,
        const OverlapRemovalParameter& corrParameter);


    /// print the condition
    virtual void print(std::ostream& myCout) const;

    /// output stream operator
    friend std::ostream& operator<<(std::ostream&, const CaloWithOverlapRemovalTemplate&);

protected:

    /// copy function for copy constructor and operator=
    void copy( const CaloWithOverlapRemovalTemplate& cp);


protected:

    l1t::GtConditionCategory m_cond0Category;
    l1t::GtConditionCategory m_condOvrCategory;
    int m_cond0Index;
    int m_condOvrIndex;

    /// variables containing the parameters
    std::vector<ObjectParameter> m_objectParameter;
    OverlapRemovalParameter m_overlapRemovalParameter;

};

#endif
