/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Leszek Grzanka (braciszek@gmail.com)
 * $$RCSfile: RPCCId.h,v $: $
 * $Revision: 1.1 $
 * $Date: 2009/01/30 19:05:12 $
 *
 ****************************************************************************/

#ifndef DataFormatsTotemL1TriggerRPCCId_h
#define DataFormatsTotemL1TriggerRPCCId_h

#include "DataFormats/DetId/interface/DetId.h"
#include <FWCore/Utilities/interface/Exception.h>

#include <iosfwd>
#include <iostream>

typedef uint32_t RPCCIdRaw;

class RPCCId: public DetId
{
 public:
  RPCCId();

  /// Construct from a packed id. It is required that the Detector part of
  /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
  explicit RPCCId(RPCCIdRaw id);

  /// Construct from fully qualified identifier.
  RPCCId(unsigned int Arm, unsigned int Station,
	 unsigned int RomanPot, unsigned int Direction);

  /// Bit 24 = Arm: 1=z>0 0=z<0
  /// Bits [22:23] Station
  /// Bits [19:21] Roman Pot number
  /// Bits [15:18] CC direction: 0 - V, 1 - U
  inline int Arm() const
    {
      return int((id_>>startArmBit) & 0x1);
    }
  inline int Station() const
    {
      return int((id_>>startStationBit) & 0x3);
    }
  inline int RomanPot() const
    {
      return int((id_>>startRPBit) & 0x7);
    }
  inline int Direction() const
    {
      return int((id_>>startDirBit) & 0xf);
    }

  inline bool IsStripsCoordinateUDirection() const {return Direction()%2;}
  inline bool IsStripsCoordinateVDirection() const {return !IsStripsCoordinateUDirection();}

  static const int startArmBit = 24;
  static const int startStationBit = 22;
  static const int startRPBit = 19;
  static const int startDirBit = 15;
  static const int totem_rp_subdet_id = 3;

 private:
  inline void init(unsigned int Arm, unsigned int Station,
		   unsigned int RomanPot, unsigned int Direction);
}; // RPCCId

std::ostream& operator<<(std::ostream& os, const RPCCId& id);

#endif  //DataFormatsTotemL1TriggerRPCCId_h
