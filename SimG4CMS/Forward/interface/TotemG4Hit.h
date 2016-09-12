#ifndef Forward_TotemG4Hit_h
#define Forward_TotemG4Hit_h 1
// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemG4Hit
//
/**\class TotemG4Hit TotemG4Hit.h SimG4CMS/Forward/interface/TotemG4Hit.h
 
 Description: Transient Hit class for Totem taken from those for Calorimeters
 
 Usage: One Hit object should be created
   -for each new particle entering the calorimeter
   -for each detector unit (= cristal or fiber or scintillator layer)
   -for each nanosecond of the shower development

   This implies that all hit objects created for a given shower
   have the same value for
   - Entry (= local coordinates of the entrance point of the particle
              in the unit where the shower starts) 
   - the TrackID (= Identification number of the incident particle)
   - the IncidentEnergy (= energy of that particle)
 
*/
//
// Original Author:
//         Created:  Tue May 16 10:14:34 CEST 2006
//
 
// system include files

// user include files

#include "G4VHit.hh"
#include <CLHEP/Vector/ThreeVector.h>
#include <boost/cstdint.hpp>
#include <iostream>

using CLHEP::Hep3Vector;

class TotemG4Hit : public G4VHit {
  
public:

  // ---------- Constructor and destructor -----------------
  TotemG4Hit();
  ~TotemG4Hit();
  TotemG4Hit(const TotemG4Hit &right);

  // ---------- operators ----------------------------------
  const TotemG4Hit& operator=(const TotemG4Hit &right); //to delete?
  int operator==(const TotemG4Hit &){return 0;} //to delete?


  // ---------- member functions ---------------------------
  void         Draw(){}
  void         Print();

  double       getEM() const;
  void         setEM (double e);

  double       getHadr() const;
  void         setHadr (double e);

  public:
  Hep3Vector getEntry() const;
  void setEntry(Hep3Vector xyz);
  Hep3Vector getExit() const;
  void setExit(Hep3Vector xyz);

  void setLocalEntry(const Hep3Vector &theLocalEntryPoint);
  void setLocalExit(const Hep3Vector &theLocalExitPoint);
  Hep3Vector getLocalEntry() const;
  Hep3Vector getLocalExit() const;
  
  double       getIncidentEnergy() const;
  void         setIncidentEnergy (double e);
  
  int          getTrackID() const;
  void         setTrackID (int i);
  
  uint32_t     getUnitID() const;
  void         setUnitID (unsigned int i);
  
  double       getTimeSlice() const;     
  void         setTimeSlice(double d);
  int          getTimeSliceID() const;     
  
  void         addEnergyDeposit(double em, double hd);
  void         addEnergyDeposit(const TotemG4Hit& aHit);
  
  double       getEnergyDeposit() const;
  
  double        getPabs() const;
  double        getTof() const;
  double        getEnergyLoss() const;
  int          getParticleType() const;

  void setPabs(double e);
  void setTof(double e);
  void setEnergyLoss(double e);
  void setParticleType(short i);

  void addEnergyLoss(double e);

  double getThetaAtEntry() const;
  double getPhiAtEntry() const;

  void setThetaAtEntry(double t);
  void setPhiAtEntry(double f);

  double getX() const;
  void setX(double t);
  double getY() const;
  double getZ() const;
  void setY(double t);
  void setZ(double t);

  int getParentId() const;
  double getVx() const;
  double getVy() const;
  double getVz() const;

  void setParentId(int p);
  void setVx(double p);
  void setVy(double p);
  void setVz(double p);

  void set_p_x(double p);
  void set_p_y(double p);
  void set_p_z(double p);

  double get_p_x() const;
  double get_p_y() const;
  double get_p_z() const;

private:
  Hep3Vector entry;             //Entry point
  Hep3Vector exit;    //Exit point
  Hep3Vector local_entry;    //local entry point
  Hep3Vector local_exit;     //local exit point
  double theIncidentEnergy; //Energy of the primary particle
  int theTrackID;        //Identification number of the primary
                                  //particle
  uint32_t theUnitID;         //Totem Unit Number
  double theTimeSlice;      //Time Slice Identification

  double theX;
  double theY;
  double theZ;
  double thePabs;
  double theTof;
  double theEnergyLoss;
  int theParticleType;

  double theThetaAtEntry;
  double thePhiAtEntry;
  Hep3Vector theEntryPoint;
  Hep3Vector theExitPoint;

  int theParentId;
  double theVx;
  double theVy;
  double theVz;

  double p_x, p_y, p_z;
};

std::ostream& operator<<(std::ostream&, const TotemG4Hit&);

#endif

