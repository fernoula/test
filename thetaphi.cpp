#include <TMath.h> 

#ifndef __THETA_PHI_CPP__
#define __THETA_PHI_CPP__

Float_t theta_from_cart(Float_t dx, Float_t dy, Float_t dz)
{
  Float_t rho_cart = sqrt(dx*dx + dy*dy);
  Float_t theta = 0.0; 

  if(dz==0.0)
    theta = TMath::Pi()/2;
  else if(dz>0.0)
    theta = atan(rho_cart/dz);
  else 
    theta = TMath::Pi() + atan(rho_cart/dz);

  return theta;
    
}

Float_t phi_from_cart(Float_t dx, Float_t dy)
{
  Float_t phi = 0.0;
  
  if(dx > 0.0 && dy >= 0.0)
    phi = atan(dy/dx);
  else if(dx < 0.0 && dy >= 0.0)
    phi = TMath::Pi() + atan(dy/dx);
  else if(dx < 0.0 && dy < 0.0)
    phi = TMath::Pi() + atan(dy/dx);
  else if(dx > 0.0 && dy < 0.0)
    phi = 2*TMath::Pi() + atan(dy/dx);
  else if(dx == 0.0 && dy > 0.0)
    phi = TMath::Pi()*1/2;
  else if (dx == 0.0 && dy < 0.0)
    phi = TMath::Pi()*3/2;

  return phi;
  
}

#endif
