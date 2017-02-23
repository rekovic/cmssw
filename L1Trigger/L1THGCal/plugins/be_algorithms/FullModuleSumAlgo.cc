#include "L1Trigger/L1THGCal/interface/HGCalTriggerBackendAlgorithmBase.h"
#include "L1Trigger/L1THGCal/interface/fe_codecs/HGCalTriggerCellBestChoiceCodec.h"
#include "L1Trigger/L1THGCal/interface/fe_codecs/HGCalTriggerCellThresholdCodec.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"

using namespace HGCalTriggerBackend;

template<typename FECODEC, typename DATA>
class FullModuleSumAlgo : public Algorithm<FECODEC> 
{
    public:
        using Algorithm<FECODEC>::name;

    protected:
        using Algorithm<FECODEC>::codec_;

    public:
        FullModuleSumAlgo(const edm::ParameterSet& conf):
            Algorithm<FECODEC>(conf),
            cluster_product_( new l1t::HGCalClusterBxCollection ){}

        virtual void setProduces(edm::EDProducer& prod) const override final 
        {
            prod.produces<l1t::HGCalClusterBxCollection>(name());
        }

        virtual void run(const l1t::HGCFETriggerDigiCollection& coll, const edm::EventSetup& es) override final;

        virtual void putInEvent(edm::Event& evt) override final 
        {
            evt.put(std::move(cluster_product_),name());
        }

        virtual void reset() override final 
        {
            cluster_product_.reset( new l1t::HGCalClusterBxCollection );
        }

    private:
        std::unique_ptr<l1t::HGCalClusterBxCollection> cluster_product_;

};

/*****************************************************************/
void FullModuleSumAlgo::run(const l1t::HGCFETriggerDigiCollection& coll, const edm::EventSetup& es) 
/*****************************************************************/
{
    for( const auto& digi : coll ) 
    {
        HGCalTriggerCellBestChoiceCodec::data_type data;
        data.reset();
        const HGCalDetId& moduleId = digi.getDetId<HGCalDetId>();
        digi.decode(codec_, data);

DEFINE_EDM_PLUGIN(HGCalTriggerBackendAlgorithmFactory, 
        FullModuleSumAlgoBestChoice,
        "FullModuleSumAlgoBestChoice");

DEFINE_EDM_PLUGIN(HGCalTriggerBackendAlgorithmFactory, 
        FullModuleSumAlgoThreshold,
        "FullModuleSumAlgoThreshold");
