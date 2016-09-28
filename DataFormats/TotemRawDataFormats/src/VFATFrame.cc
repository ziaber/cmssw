/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*  Leszek Grzanka
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"

#include <stdio.h>
#include <cstring>

namespace Totem {

//----------------------------------------------------------------------------------------------------

VFATFrame::VFATFrame(unsigned short *_data)
{
  if (_data)
    setData(_data);
  else
    memset(data, 0, 12 * sizeof(unsigned short));
}

//----------------------------------------------------------------------------------------------------

void VFATFrame::setData(const VFATFrame::word *_data)
{
  memcpy(data, _data, 24);
}

//----------------------------------------------------------------------------------------------------

std::vector<unsigned char> VFATFrame::getActiveChannels() const
{
  std::vector<unsigned char> channels;

  for (int i = 0; i < 8; i++)
  {
    // quick check
    if (!data[1 + i])
      continue;

    // go throug bits
    unsigned short mask;
    char offset;
    for (mask = 1 << 15, offset = 15; mask; mask >>= 1, offset--)
    {
      if (data[1 + i] & mask)
        channels.push_back( i * 16 + offset );
    }
  }

  return channels;
}

//----------------------------------------------------------------------------------------------------

bool VFATFrame::checkCRC() const
{
  unsigned short int crc_fin = 0xffff;

  for (int i = 11; i >= 1; i--)
    crc_fin = crc_calc(crc_fin, data[i]);
  
  return (crc_fin == data[0]);
}

//----------------------------------------------------------------------------------------------------

VFATFrame::word VFATFrame::crc_calc(VFATFrame::word crc_in, VFATFrame::word dato)
{
  word v = 0x0001;
  word mask = 0x0001;    
  bool d=0;
  word crc_temp = crc_in;
  unsigned char datalen = 16;

  for (int i=0; i<datalen; i++)
  {
    if (dato & v)
      d = 1; 
    else
      d = 0;
      
    if ((crc_temp & mask)^d)
      crc_temp = crc_temp>>1 ^ 0x8408;
    else
      crc_temp = crc_temp>>1;
       
    v<<=1;
  }

  return crc_temp;
}

//----------------------------------------------------------------------------------------------------

void VFATFrame::Print(bool binary) const
{
  if (binary)
  {
    for (int i = 0; i < 12; i++)
    {
      const unsigned short &word = data[11 - i];
      unsigned short mask = (1 << 15);
      for (int j = 0; j < 16; j++)
      {
        if (word & mask)
          printf("1");
        else
          printf("0");
        mask = (mask >> 1);
        if ((j + 1) % 4 == 0)
          printf("|");
      }
      printf("\n");
    }
  } else {
    printf("ID = %03x, BC = %04u, EC = %03u, flags = %2u, CRC = %04x ", getChipID(), getBC(), getEC(), getFlags(), getCRC());
    if (checkCRC())
      printf("(  OK), footprint ");
    else
      printf("(FAIL), footprint ");
    if (checkFootprint())
      printf("  OK");
    else
      printf("FAIL");
    printf(", frame = %04x|%04x|%04x|", data[11], data[10], data[9]);
    for (int i = 8; i > 0; i--)
      printf("%04x", data[i]);
    printf("|%04x\n", data[0]);
  }
}

} // namespace
