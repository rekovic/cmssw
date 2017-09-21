#ifndef L1T_PACKER_STAGE2_REGIONALMUONGMTUNPACKER_H
#define L1T_PACKER_STAGE2_REGIONALMUONGMTUNPACKER_H

#include "EventFilter/L1TRawToDigi/interface/Unpacker.h"
#include "EventFilter/L1TRawToDigi/interface/Block.h"
#include "GMTCollections.h"

namespace l1t {
   namespace stage2 {
      class RegionalMuonGMTUnpacker : public Unpacker {
         public:
            virtual bool unpack(const Block& block, UnpackerCollections *coll) override;

         private:
            static const unsigned int nWords_ = 6; // every link transmits 6 words (3 muons) per bx
            static const unsigned int bxzs_enable_shift_ = 1;
      };
   }
}

#endif
