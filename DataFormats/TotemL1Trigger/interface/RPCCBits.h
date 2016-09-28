/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Leszek Grzanka (braciszek@gmail.com)
 * $$RCSfile: RPCCBits.h,v $: $
 * $Revision: 1.1 $
 * $Date: 2009/01/30 19:05:12 $
 *
 ****************************************************************************/

#ifndef DataFormatsTotemL1TriggerRPCCBits_h
#define DataFormatsTotemL1TriggerRPCCBits_h

#include "DataFormats/TotemL1Trigger/interface/RPCCId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include <iosfwd>
#include <iostream>
#include <cstdlib>
#include <bitset>

class RPCCBits {
 public:
  /// Construct from a packed id. It is required that the Detector part of
  /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
  explicit RPCCBits(RPCCIdRaw id, std::bitset<16> bs) : id_(id){ setBS(bs); };

  /// Construct from fully qualified identifier.
  explicit RPCCBits( RPCCId id, std::bitset<16> bs) { id_ = id.rawId(); setBS(bs); };

  RPCCBits() : id_(0){ 
      reset();
  };
  void reset(){
     std::bitset<16> nullBitset;
     nullBitset.reset();
	 setBS(nullBitset);
  }

  void setId( RPCCIdRaw id ){id_ = id;};
  void setId( RPCCId id ){id_ = id.rawId();};

  inline RPCCIdRaw getId() const {return id_;};

  inline std::bitset<16> getBS() const {
    std::bitset<16> res;
    for( unsigned short i = 0 ; i < 16 ; ++i ){
     res[i] = bs_[i];
   }
    return res;
 };
  
  inline const bool* getRawBS() const
  {
    return bs_; 
  }

  void setBS( std::bitset<16> bs ){
    for( unsigned short i = 0 ; i < 16 ; ++i ){
     bs_[i] = bs[i];
   }
 };

 private:
   RPCCIdRaw id_;
   bool bs_[16];

}; // RPCCBits

// Comparison operators
inline bool operator<( const RPCCBits& one, const RPCCBits& other) {
  if(one.getId() < other.getId())
    return true;
  else if(one.getId() == other.getId())
    return one.getBS().to_ulong() < other.getBS().to_ulong();
  else
    return false;
}

#endif  //DataFormatsTotemL1TriggerRPCCBits_h