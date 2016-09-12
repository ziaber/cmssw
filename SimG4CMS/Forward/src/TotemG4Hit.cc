// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemG4Hit
//
// Implementation:
//     <Notes on implementation>
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//

// system include files

// user include files
#include "SimG4CMS/Forward/interface/TotemG4Hit.h"
#include <iostream>

//
// constructors and destructor
//

TotemG4Hit::TotemG4Hit(): MeanPosition(0), theEntryPoint(0), theExitPoint(0) {

  elem              = 0.;
  hadr              = 0.;
  theIncidentEnergy = 0.;
  theTrackID        = -1;
  theUnitID         =  0;
  theTimeSlice      = 0.;

  theX              = 0.;
  theY              = 0.;
  theZ              = 0.;
  thePabs           = 0.;
  theTof            = 0.;
  theEnergyLoss     = 0.;
  theParticleType   = 0;
  theThetaAtEntry   = 0.;
  thePhiAtEntry     = 0.;
  theParentId       = 0;
  theVx             = 0.;
  theVy             = 0.;
  theVz             = 0.;
  thePx             = 0;
  thePy             = 0;
  thePz             = 0;
  theVPx            = 0;
  theVPy            = 0;
  theVPz            = 0;
}

TotemG4Hit::~TotemG4Hit() {}

TotemG4Hit::TotemG4Hit(const TotemG4Hit &right) {

  MeanPosition      = right.MeanPosition;

  elem              = right.elem;
  hadr              = right.hadr;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID        = right.theTrackID;
  theUnitID         = right.theUnitID;
  theTimeSlice      = right.theTimeSlice;

  theX              = right.theX;
  theY              = right.theY;
  theZ              = right.theZ;
  thePabs           = right.thePabs;
  theTof            = right.theTof;
  theEnergyLoss     = right.theEnergyLoss;
  theParticleType   = right.theParticleType;

  theThetaAtEntry   = right.theThetaAtEntry;
  thePhiAtEntry     = right.thePhiAtEntry;
  theEntryPoint     = right.theEntryPoint;
  theExitPoint      = right.theExitPoint;
  theParentId       = right.theParentId;
  theVx             = right.theVx;
  theVy             = right.theVy;
  theVz             = right.theVz;
  thePx             = right.thePx;
  thePy             = right.thePy;
  thePz             = right.thePz;
  theVPx            = right.theVPx;
  theVPy            = right.theVPy;
  theVPz            = right.theVPz;
}


const TotemG4Hit& TotemG4Hit::operator=(const TotemG4Hit &right) {

  MeanPosition      = right.MeanPosition;
  elem              = right.elem;
  hadr              = right.hadr;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID        = right.theTrackID;
  theUnitID         = right.theUnitID;
  theTimeSlice      = right.theTimeSlice;
 
  theX              = right.theX;
  theY              = right.theY;
  theZ              = right.theZ;
  thePabs           = right.thePabs;
  theTof            = right.theTof ;
  theEnergyLoss     = right.theEnergyLoss   ;
  theParticleType   = right.theParticleType ;

  theThetaAtEntry   = right.theThetaAtEntry;
  thePhiAtEntry     = right.thePhiAtEntry;
  theEntryPoint     = right.theEntryPoint;
  theExitPoint      = right.theExitPoint;
  theParentId       = right.theParentId;
  theVx             = right.theVx;
  theVy             = right.theVy;
  theVz             = right.theVz;
  thePx             = right.thePx;
  thePy             = right.thePy;
  thePz             = right.thePz;
  theVPx            = right.theVPx;
  theVPy            = right.theVPy;
  theVPz            = right.theVPz;

  return *this;
}

void TotemG4Hit::addEnergyDeposit(const TotemG4Hit& aHit) {

  elem += aHit.getEM();
  hadr += aHit.getHadr();
}


void TotemG4Hit::Print() {
  std::cout << (*this);
}


Hep3Vector TotemG4Hit::getEntryPoint() const           {return theEntryPoint;}
void       TotemG4Hit::setEntryPoint(Hep3Vector xyz)   { theEntryPoint    = xyz; }
Hep3Vector TotemG4Hit::getExitPoint() const           {return theExitPoint;}
void       TotemG4Hit::setExitPoint(Hep3Vector xyz)   { theExitPoint    = xyz; }

double     TotemG4Hit::getEM() const              {return elem; }
void       TotemG4Hit::setEM (double e)           { elem     = e; }
      
double     TotemG4Hit::getHadr() const            {return hadr; }
void       TotemG4Hit::setHadr (double e)         { hadr     = e; }
      
double     TotemG4Hit::getIncidentEnergy() const  {return theIncidentEnergy; }
void       TotemG4Hit::setIncidentEnergy(double e) {theIncidentEnergy  = e; }

int        TotemG4Hit::getTrackID() const         {return theTrackID; }
void       TotemG4Hit::setTrackID (int i)         { theTrackID = i; }

uint32_t   TotemG4Hit::getUnitID() const          {return theUnitID; }
void       TotemG4Hit::setUnitID (uint32_t i)     { theUnitID = i; }

double     TotemG4Hit::getTimeSlice() const       {return theTimeSlice; }
void       TotemG4Hit::setTimeSlice (double d)    { theTimeSlice = d; }
int        TotemG4Hit::getTimeSliceID() const     {return (int)theTimeSlice;}

void       TotemG4Hit::addEnergyDeposit(double em, double hd) {elem += em;  hadr += hd;}

double     TotemG4Hit::getEnergyDeposit() const   {return elem+hadr;}

float      TotemG4Hit::getPabs() const            {return thePabs;}
float      TotemG4Hit::getTof() const             {return theTof;}
float      TotemG4Hit::getEnergyLoss() const      {return theEnergyLoss;}
int        TotemG4Hit::getParticleType() const    {return theParticleType;}
float      TotemG4Hit::getPx() const {return thePx;}
float      TotemG4Hit::getPy() const {return thePy;}
float      TotemG4Hit::getPz() const {return thePz;}
float      TotemG4Hit::getVPx() const {return theVPx;}
float      TotemG4Hit::getVPy() const {return theVPy;}
float      TotemG4Hit::getVPz() const {return theVPz;}

void       TotemG4Hit::setPabs(float e)           {thePabs = e;}
void       TotemG4Hit::setTof(float e)            {theTof = e;}
void       TotemG4Hit::setEnergyLoss(float e)     {theEnergyLoss = e;}
void       TotemG4Hit::setParticleType(short i)   {theParticleType = i;}
void       TotemG4Hit::setPx(float e) {thePx = e;}
void       TotemG4Hit::setPy(float e) {thePy = e;}
void       TotemG4Hit::setPz(float e) {thePz = e;}
void       TotemG4Hit::setVPx(float e) {theVPx = e;}
void       TotemG4Hit::setVPy(float e) {theVPy = e;}
void       TotemG4Hit::setVPz(float e) {theVPz = e;}

float      TotemG4Hit::getThetaAtEntry() const    {return theThetaAtEntry;}   
float      TotemG4Hit::getPhiAtEntry() const      {return thePhiAtEntry;}

void       TotemG4Hit::setThetaAtEntry(float t)   {theThetaAtEntry = t;}
void       TotemG4Hit::setPhiAtEntry(float f)     {thePhiAtEntry = f ;}

float      TotemG4Hit::getX() const               {return theX;}
void       TotemG4Hit::setX(float t)              {theX = t;}

float      TotemG4Hit::getY() const               {return theY;}
void       TotemG4Hit::setY(float t)              {theY = t;}

float      TotemG4Hit::getZ() const               {return theZ;}
void       TotemG4Hit::setZ(float t)              {theZ = t;}

int        TotemG4Hit::getParentId() const        {return theParentId;}
void       TotemG4Hit::setParentId(int p)         {theParentId = p;}

float      TotemG4Hit::getVx() const              {return theVx;}
void       TotemG4Hit::setVx(float t)             {theVx = t;}

float      TotemG4Hit::getVy() const              {return theVy;}
void       TotemG4Hit::setVy(float t)             {theVy = t;}

float      TotemG4Hit::getVz() const              {return theVz;}
void       TotemG4Hit::setVz(float t)             {theVz = t;}

std::ostream& operator<<(std::ostream& os, const TotemG4Hit& hit) {
  os << " Data of this TotemG4Hit are:\n"
  << " Time slice ID: " << hit.getTimeSliceID() << "\n"
  << " EnergyDeposit = " << hit.getEnergyLoss() << "\n"
  << " Energy of primary particle (ID = " << hit.getTrackID()
  << ") = " << hit.getIncidentEnergy() << " (MeV)" << "\n"
  << " Local entry and exit points in Totem unit number " << hit.getUnitID()
  << " are: " << hit.getEntryPoint() << " (mm)" << hit.getExitPoint() << " (mm)" <<"\n"
  << " Global posizion in Totem unit number " << hit.getUnitID()
  << " are: " << hit.getMeanPosition() << " (mm)" <<"\n"
  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  return os;
}


