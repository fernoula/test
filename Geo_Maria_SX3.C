
{
  

   gSystem->Load("TwoBodySimulator_kinematics_cpp.so"); 
 

    Simulator* S = new Simulator("18Nesrim/18Ne_in_HeCO2_356Torr_18Nerun_09122018.eloss","18Nesrim/H_in_HeCO2_379Torr_allruns_09122018.eloss","18Nesrim/21Na_in_HeCO2_356Torr_18Nerun_10082018.eloss","18Nesrim/20Ne_in_HeCO2_356Torr_18Nerun_12242018.eloss","H_in_Si.txt","WorldCoordinates2018Dec4","pcr_list_07192018.dat");   //uses the new zr-offsets for the SX3s 2.5 & 14.9 cm 
 
 
  CS *cross = new CS;
  cross->ReadFile("cross.txt");

  Strag *straggling = new Strag;
  straggling->ReadFile("18Nesrim/18Ne.txt");
  
  
   //Mass of the Beam and target
  float mB = 16767.09961; //mass of 18Ne
  float mT = 3727.37892;  //mass of He-4.  
  //Mass of the Products.  
  float mL = 938.27197;  //mass of proton; 
  float mH = 19553.56837; //mass of 21Na
  float mH2 = 18617.72807; //mass of 20Ne
   
  

  Double_t e=72.34;    //Beam energy range max after passing Kapton foil
  Double_t emin=0.55;   //Beam energy range min the E_threshold


  S->SetMasses(mB, mT, mL,mH, mH2);

  S->SetBeamEnergyRange(emin,e);

  
  Int_t Det_ID;
  Float_t Xpc,Ypc,Zpc;
  Float_t X,Y,Z,Xa2,Ya2,Za2, Wx, Wy, Wz;
  Float_t Xa1f,Ya1f,Za1f,Xa2f,Ya2f,Za2f;

  Float_t Xpc2,Ypc2,Zpc2;
  Float_t Xp2,Yp2,Zp2;

  
  TGraph2D* World = new TGraph2D(2);

  World->SetPoint(0,15,15,-5);
  World->SetPoint(1,-15,-15,35);

  TGraph2D *Beam_xyz = new TGraph2D();
  Beam_xyz->SetMarkerColor(2);
  Beam_xyz->SetMarkerSize(1);
  Beam_xyz->SetMarkerStyle(7);

  TGraph *Beam_xyz_profile = new TGraph();
  Beam_xyz_profile->SetMarkerColor(2);
  Beam_xyz_profile->SetMarkerSize(1);
  Beam_xyz_profile->SetMarkerStyle(7);

  TGraph *Beam_xz = new TGraph();
  TGraph *Beam_yz = new TGraph();
  
  TGraph2D* SX3_lines = new TGraph2D();
  TGraph2D* QQQ_lines = new TGraph2D();
  TGraph2D* MISS_lines = new TGraph2D();
  TGraph *qqq_profile = new TGraph();
  SX3_lines->SetMarkerColor(4);
  SX3_lines->SetMarkerSize(1);
  SX3_lines->SetMarkerStyle(7);

  QQQ_lines->SetMarkerColor(6);
  QQQ_lines->SetMarkerSize(1);
  QQQ_lines->SetMarkerStyle(7);
  
  MISS_lines->SetMarkerColor(3);
  MISS_lines->SetMarkerSize(1);
  MISS_lines->SetMarkerStyle(7);


  TGraph2D* SX3_lines_2p = new TGraph2D();
  TGraph2D* QQQ_lines_2p = new TGraph2D(); 
  TGraph *qqq_profile_2p = new TGraph();

  SX3_lines_2p->SetMarkerColor(1);
  SX3_lines_2p->SetMarkerSize(1);
  SX3_lines_2p->SetMarkerStyle(3);

  QQQ_lines_2p->SetMarkerColor(3);
  QQQ_lines_2p->SetMarkerSize(1);
  QQQ_lines_2p->SetMarkerStyle(3);

  //TProfile2D *hprof2d = new TProfile2D("hprof2d","Profile QQQ",100,-20,20,100,-20,20,0,20);
   
  bool GoodHit = 0;
  TPolyLine3D** F;
  F = new TPolyLine3D*[30];
  
  for (Int_t j=6; j<30; j++) {  //j=DetID
    F[j] = new TPolyLine3D();
    for (Int_t m=0; m<5; m++) {
      if (j>5) {
	S-> DrawLines(j, m, Wx, Wy, Wz);
	F[j]->SetPoint(m, Wx, Wy, Wz);
	F[j]->SetLineColor(2);
	// cout<<"j="<<j <<" m="<< m <<" Wx="<<Wx<<" Wy="<<Wy<<" Wz="<<Wz<<endl; 
      }
    }
  } 
  
  //int fwd_point_index = 0;
  int nnn =0, nnn2;
  int mmm =0, mmm2;
  int na =0,ma=0;
  int beam_xyz=0, beam_xy=0, b_xz=0, b_yz=0;
  ///////////////////////////////////////////
 
  int i=0;
  const Int_t k=500;
  
   while(i<k){
     if(!S->GenerateEvent(cross,straggling))
       continue;
     

     if(!S->SimThetaPhiE_in(emin,e))
       continue;    
     
     if(S->ke_L_lab <= 0.5)
       continue;
     
     if(!S->GetQQQ3FwdCoord(X,Y,Z, Xpc, Ypc, Zpc, Det_ID, S->theta_L, S->phi_L)){
       if(!S->GetSX3Coord(X,Y,Z, Xpc, Ypc, Zpc, S->theta_L, S->phi_L)){
	 //MISS_lines->SetPoint(na,X,Y,Z);
	 // na++;
	 continue;
       }
     }


     //---------Beam Profile graphs---------------
     Beam_xyz->SetPoint(beam_xyz,S->xr,S->yr,S->zr);
     beam_xyz++;
     Beam_xyz_profile->SetPoint(beam_xy,S->xr,S->yr);
     beam_xy++;

     Beam_xz->SetPoint(b_xz,S->zr,S->xr);
     b_xz++;
     Beam_yz->SetPoint(b_yz,S->zr,S->yr);
     b_yz++;
     

      if(S->GetSX3Coord(X, Y, Z, Xpc, Ypc, Zpc, S->theta_L, S->phi_L)) {
       SX3_lines->SetPoint(nnn, X, Y, Z);
       nnn++;
          }

      if(S->GetQQQ3FwdCoord(X,Y,Z, Xpc, Ypc, Zpc, Det_ID, S->theta_L, S->phi_L)){
	QQQ_lines->SetPoint(mmm, X, Y, Z);	
	qqq_profile -> SetPoint(mmm,X,Y);
	mmm++;
          }

      /////////////////////////////////////////////////////////////////////////////////////
      //--------------- Second Proton ----------------------------------//

      if(!S->SecondProton(S->ee_H_lab, S->Px_H_lab, S->Py_H_lab, S->Pz_H_lab, S->Ex))
	continue;
           
      if(S->ke_1_lab<=0.5)      
	continue;

    
      if(!S->GetQQQ3FwdCoord(Xp2, Yp2, Zp2, Xpc2, Ypc2, Zpc2, Det_ID, S->theta_h1, S->phi_h1))
	{
	  if(!S->GetSX3Coord(Xp2, Yp2, Zp2, Xpc2, Ypc2, Zpc2, S->theta_h1, S->phi_h1))
	    {
	      continue;
	    }
	}
      
      if(S->GetSX3Coord(Xp2, Yp2, Zp2, Xpc2, Ypc2, Zpc2, S->theta_h1, S->phi_h1)) {
       SX3_lines_2p->SetPoint(nnn2, Xp2, Yp2, Zp2);
       nnn2++;
          }

      if(S->GetQQQ3FwdCoord(Xp2,Yp2,Zp2, Xpc2, Ypc2, Zpc2, Det_ID, S->theta_h1, S->phi_h1)){
	QQQ_lines_2p->SetPoint(mmm2, Xp2, Yp2, Zp2);	
	qqq_profile_2p -> SetPoint(mmm2,Xp2,Yp2);
	mmm2++;
          }
       
     i++;
     cout << "\rDone: " <<setprecision(2) << fixed << showpoint <<  i/double(k)*100 << flush;
     
   }

  

  TCanvas* Display2 = new TCanvas("SX3 Geo","SX3 Geo",0,0,800,800);
  Display2->cd();
  World->Draw("p");
  Beam_xyz->Draw("p same");
  SX3_lines->Draw("p same");
  QQQ_lines->Draw("p same");
  SX3_lines_2p->Draw("p same");
  QQQ_lines_2p->Draw("p same");
  
  //MISS_lines->Draw("p same");

  for (Int_t j=6; j<30; j++)  
    F[j]->Draw("l same");

  Display2->SaveAs("graphs/3d_geo_WITHstraggling_WITHcross_500events.png");

   
  TCanvas *c1 = new TCanvas("c1", "Profile QQQ",200, 10,700,500);
  qqq_profile->SetTitle("XY profile QQQ");
  qqq_profile->GetXaxis()->SetTitle("X coord");
  qqq_profile->GetYaxis()->SetTitle("Y coord");
  qqq_profile->SetMarkerColor(6);
  qqq_profile->SetMarkerSize(1);
  qqq_profile->SetMarkerStyle(7);
  qqq_profile->Draw("AP");
  qqq_profile_2p->SetMarkerColor(1);
  qqq_profile_2p->SetMarkerSize(1);
  qqq_profile_2p->SetMarkerStyle(3);
  qqq_profile_2p->Draw("P same");
  Beam_xyz_profile->Draw("P same");

  c1->SaveAs("graphs/XY_profile_QQQ_WITHstraggling_WITHcross_500events.png");

  TCanvas *c2 = new TCanvas("c2","XZ - YZ distribution",300,100,1500,1500);
  c2->Divide(2,1);
  c2->cd(1);
  Beam_xz->SetTitle("XZ distribution");
  Beam_xz->GetXaxis()->SetTitle("Z coord");
  Beam_xz->GetYaxis()->SetTitle("X coord");
  Beam_xz->SetMarkerColor(6);
  Beam_xz->SetMarkerSize(1);
  Beam_xz->SetMarkerStyle(7);
  Beam_xz->Draw("AP");

  c2->cd(2);
  Beam_yz->SetTitle("YZ distribution");
  Beam_yz->GetXaxis()->SetTitle("Z coord");
  Beam_yz->GetYaxis()->SetTitle("Y coord");
  Beam_yz->SetMarkerColor(6);
  Beam_yz->SetMarkerSize(1);
  Beam_yz->SetMarkerStyle(7);
  Beam_yz->Draw("AP");

  
  c2->SaveAs("graphs/XZ_YZ_distribution_WITHstraggling_WITHcross_500events.png");
  
  
  
}
