#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h

class LeptonEfficiencySF {
public:
  LeptonEfficiencySF(int era=2015) : era_(era) { }
  ~LeptonEfficiencySF() {}

  // pair<SF, pair<stat(SF), syst(SF)> >

  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(id==11) return getLeptonIdEfficiency(pt, eta, id); 

    if(era_==2015) { 
      std::pair<float, std::pair<float,float> > effId  = getLeptonIdEfficiency(pt, eta, id); 
      std::pair<float, std::pair<float,float> > effIso = getLeptonIsoEfficiency(pt, eta, id); 

      eff.first = effId.first * effIso.first; 

      // N.B.: assuming that the errors on ID and ISO scale factors are fully uncorrelated --> certainly not true!! 
      //       more conservative approach would be to treat them as fully correlated (i.e. sum them linearly) 
      eff.second.first  = eff.first * sqrt( (effId.second.first *effId.second.first )/(effId.first*effId.first) + (effIso.second.first *effIso.second.first )/(effIso.first*effIso.first) ); 
      eff.second.second = eff.first * sqrt( (effId.second.second*effId.second.second)/(effId.first*effId.first) + (effIso.second.second*effIso.second.second)/(effIso.first*effIso.first) ); 
    } 

    return eff; 
  }



  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonIdEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(era_==2015) { 

      // Medium ele - ID efficiency scale factors 
      if(id==11) { 
	if( pt<20 ) { 
	  if     ( eta<-2.000 ) { eff.first = 1.25131 ; eff.second.first = 0.0279617 ; eff.second.first = 0.0100; } 
	  else if( eta<-1.566 ) { eff.first = 1.09254 ; eff.second.first = 0.0198698 ; eff.second.first = 0.0100; } 
	  else if( eta<-1.444 ) { eff.first = 1.29412 ; eff.second.first = 0.0666726 ; eff.second.first = 0.0100; } 
	  else if( eta<-0.800 ) { eff.first = 1.32673 ; eff.second.first = 0.0384136 ; eff.second.first = 0.0100; } 
	  else if( eta< 0.000 ) { eff.first = 1.23371 ; eff.second.first = 0.0177938 ; eff.second.first = 0.0100; } 
	  else if( eta< 0.800 ) { eff.first = 1.22889 ; eff.second.first = 0.0236843 ; eff.second.first = 0.0100; } 
	  else if( eta< 1.444 ) { eff.first = 1.30827 ; eff.second.first = 0.0269239 ; eff.second.first = 0.0100; } 
	  else if( eta< 1.566 ) { eff.first = 1.35498 ; eff.second.first = 0.0690439 ; eff.second.first = 0.0100; } 
	  else if( eta< 2.000 ) { eff.first = 1.06283 ; eff.second.first = 0.020136  ; eff.second.first = 0.0100; } 
	  else                  { eff.first = 1.20769 ; eff.second.first = 0.0495958 ; eff.second.first = 0.0100; } 
	} 
	else if( pt<30 ) { 
	  if     ( eta<-2.000 ) { eff.first = 1.0536  ; eff.second.first = 0.00990836; eff.second.first = 0.0100; } 
	  else if( eta<-1.566 ) { eff.first = 0.969453; eff.second.first = 0.0116778 ; eff.second.first = 0.0100; } 
	  else if( eta<-1.444 ) { eff.first = 1.04043 ; eff.second.first = 0.0231465 ; eff.second.first = 0.0100; } 
	  else if( eta<-0.800 ) { eff.first = 1.05939 ; eff.second.first = 0.00871653; eff.second.first = 0.0100; } 
	  else if( eta< 0.000 ) { eff.first = 1.02262 ; eff.second.first = 0.00814801; eff.second.first = 0.0100; } 
	  else if( eta< 0.800 ) { eff.first = 1.01799 ; eff.second.first = 0.00672915; eff.second.first = 0.0100; } 
	  else if( eta< 1.444 ) { eff.first = 1.07742 ; eff.second.first = 0.0276387 ; eff.second.first = 0.0100; } 
	  else if( eta< 1.566 ) { eff.first = 1.18801 ; eff.second.first = 0.026376  ; eff.second.first = 0.0100; } 
	  else if( eta< 2.000 ) { eff.first = 0.980392; eff.second.first = 0.00727898; eff.second.first = 0.0100; } 
	  else                  { eff.first = 1.02649 ; eff.second.first = 0.00972224; eff.second.first = 0.0100; } 
	} 
	else if( pt<40 ) { 
	  if     ( eta<-2.000 ) { eff.first = 1.01235 ; eff.second.first = 0.00496475; eff.second.first = 0.0100; } 
	  else if( eta<-1.566 ) { eff.first = 0.971279; eff.second.first = 0.00466581; eff.second.first = 0.0100; } 
	  else if( eta<-1.444 ) { eff.first = 1.02482 ; eff.second.first = 0.0119536 ; eff.second.first = 0.0100; } 
	  else if( eta<-0.800 ) { eff.first = 0.990897; eff.second.first = 0.00290249; eff.second.first = 0.0100; } 
	  else if( eta< 0.000 ) { eff.first = 0.978616; eff.second.first = 0.00280074; eff.second.first = 0.0100; } 
	  else if( eta< 0.800 ) { eff.first = 0.986164; eff.second.first = 0.00280492; eff.second.first = 0.0100; } 
	  else if( eta< 1.444 ) { eff.first = 1.00654 ; eff.second.first = 0.0029268 ; eff.second.first = 0.0100; } 
	  else if( eta< 1.566 ) { eff.first = 0.990991; eff.second.first = 0.0153776 ; eff.second.first = 0.0100; } 
	  else if( eta< 2.000 ) { eff.first = 0.961892; eff.second.first = 0.0046831 ; eff.second.first = 0.0100; } 
	  else                  { eff.first = 1.0096  ; eff.second.first = 0.00496055; eff.second.first = 0.0100; } 
	} 
	else if( pt<50 ) { 
	  if     ( eta<-2.000 ) { eff.first = 1.00875 ; eff.second.first = 0.00451911; eff.second.first = 0.0100; } 
	  else if( eta<-1.566 ) { eff.first = 0.983313; eff.second.first = 0.00334318; eff.second.first = 0.0100; } 
	  else if( eta<-1.444 ) { eff.first = 0.990155; eff.second.first = 0.0089713 ; eff.second.first = 0.0100; } 
	  else if( eta<-0.800 ) { eff.first = 0.974057; eff.second.first = 0.00164621; eff.second.first = 0.0100; } 
	  else if( eta< 0.000 ) { eff.first = 0.969838; eff.second.first = 0.00161607; eff.second.first = 0.0100; } 
	  else if( eta< 0.800 ) { eff.first = 0.974448; eff.second.first = 0.00162168; eff.second.first = 0.0100; } 
	  else if( eta< 1.444 ) { eff.first = 0.976359; eff.second.first = 0.00165201; eff.second.first = 0.0100; } 
	  else if( eta< 1.566 ) { eff.first = 0.956153; eff.second.first = 0.0399717 ; eff.second.first = 0.0100; } 
	  else if( eta< 2.000 ) { eff.first = 0.978469; eff.second.first = 0.00334706; eff.second.first = 0.0100; } 
	  else                  { eff.first = 1.00125 ; eff.second.first = 0.00451998; eff.second.first = 0.0100; } 
	} 
	else { 
	  if     ( eta<-2.000 ) { eff.first = 1.006   ; eff.second.first = 0.00770486; eff.second.first = 0.0100; } 
	  else if( eta<-1.566 ) { eff.first = 0.985109; eff.second.first = 0.00569682; eff.second.first = 0.0100; } 
	  else if( eta<-1.444 ) { eff.first = 0.948856; eff.second.first = 0.0143331 ; eff.second.first = 0.0100; } 
	  else if( eta<-0.800 ) { eff.first = 0.958097; eff.second.first = 0.00403142; eff.second.first = 0.0100; } 
	  else if( eta< 0.000 ) { eff.first = 0.967489; eff.second.first = 0.00311976; eff.second.first = 0.0100; } 
	  else if( eta< 0.800 ) { eff.first = 0.970819; eff.second.first = 0.00312847; eff.second.first = 0.0100; } 
	  else if( eta< 1.444 ) { eff.first = 0.964773; eff.second.first = 0.00405335; eff.second.first = 0.0100; } 
	  else if( eta< 1.566 ) { eff.first = 0.974394; eff.second.first = 0.0179898 ; eff.second.first = 0.0100; } 
	  else if( eta< 2.000 ) { eff.first = 0.994253; eff.second.first = 0.00573526; eff.second.first = 0.0100; } 
	  else                  { eff.first = 1.00481 ; eff.second.first = 0.00869049; eff.second.first = 0.0100; } 
	} 
      } // end id==11 

      // Medium mu - ID efficiency scale factors 
      if(id==13) { 
	if(fabs(eta)<0.9) { 
	  if     (pt<25) { eff.first = 0.979414; eff.second.first = 0.0027415; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 0.984815; eff.second.first = 0.0012909; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 0.989551; eff.second.first = 0.0004332; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.991315; eff.second.first = 0.0003197; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 0.987545; eff.second.first = 0.0008304; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.991321; eff.second.first = 0.0016589; eff.second.second = 0.0100000; } 
	} 
	else if(fabs(eta)<1.2) { 
	  if     (pt<25) { eff.first = 0.987974; eff.second.first = 0.0040435; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 0.985705; eff.second.first = 0.0022352; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 0.991785; eff.second.first = 0.0007619; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.992362; eff.second.first = 0.0005198; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 0.991781; eff.second.first = 0.0013460; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.992222; eff.second.first = 0.0028402; eff.second.second = 0.0100000; } 
	} 
	else if(fabs(eta)<2.1) { 
	  if     (pt<25) { eff.first = 0.995794; eff.second.first = 0.0021834; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 0.991305; eff.second.first = 0.0012558; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 0.993265; eff.second.first = 0.0004708; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.994487; eff.second.first = 0.0003082; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 0.991223; eff.second.first = 0.0009382; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.995587; eff.second.first = 0.0024502; eff.second.second = 0.0100000; } 
	} 
	else { 
	  if     (pt<25) { eff.first = 0.976729; eff.second.first = 0.0047255; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 0.967756; eff.second.first = 0.0030171; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 0.967950; eff.second.first = 0.0014332; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.965109; eff.second.first = 0.0013090; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 0.957879; eff.second.first = 0.0033303; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.974770; eff.second.first = 0.0106513; eff.second.second = 0.0100000; } 
	} 
      } // end id==13 
    } // end era==2015 

    return eff;
  } 


  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonIsoEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(era_==2015) { 

      // Medium ele - ISO efficiency scale factors 
      if(id==11) { eff.first=1.000; eff.second.first=0.000; eff.second.second=0.000; } // useless, everything done in getLeptonIdEfficiency 


      // Medium mu - ISO efficiency scale factors 
      if(id==13) { 
	if(fabs(eta)<0.9) { 
	  if     (pt<25) { eff.first = 1.004513; eff.second.first = 0.0038308; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 1.000996; eff.second.first = 0.0021793; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 1.000794; eff.second.first = 0.0007745; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.998682; eff.second.first = 0.0004395; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 1.000002; eff.second.first = 0.0007444; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.998554; eff.second.first = 0.0008983; eff.second.second = 0.0100000; } 
	} 
	else if(fabs(eta)<1.2) { 
	  if     (pt<25) { eff.first = 1.004341; eff.second.first = 0.0056790; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 1.003733; eff.second.first = 0.0037925; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 1.004714; eff.second.first = 0.0014141; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 1.001664; eff.second.first = 0.0007240; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 0.997668; eff.second.first = 0.0011800; eff.second.second = 0.0100000; } 
	  else           { eff.first = 1.004428; eff.second.first = 0.0015655; eff.second.second = 0.0100000; } 
	} 
	else if(fabs(eta)<2.1) { 
	  if     (pt<25) { eff.first = 0.996890; eff.second.first = 0.0029755; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 1.000283; eff.second.first = 0.0019757; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 1.003120; eff.second.first = 0.0008158; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 0.999851; eff.second.first = 0.0004164; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 1.000648; eff.second.first = 0.0007262; eff.second.second = 0.0100000; } 
	  else           { eff.first = 1.000566; eff.second.first = 0.0009121; eff.second.second = 0.0100000; } 
	} 
	else { 
	  if     (pt<25) { eff.first = 0.993198; eff.second.first = 0.0058158; eff.second.second = 0.0100000; } 
	  else if(pt<30) { eff.first = 0.993565; eff.second.first = 0.0040015; eff.second.second = 0.0100000; } 
	  else if(pt<40) { eff.first = 0.996759; eff.second.first = 0.0016758; eff.second.second = 0.0100000; } 
	  else if(pt<50) { eff.first = 1.002709; eff.second.first = 0.0010661; eff.second.second = 0.0100000; } 
	  else if(pt<60) { eff.first = 1.004295; eff.second.first = 0.0017420; eff.second.second = 0.0100000; } 
	  else           { eff.first = 0.998136; eff.second.first = 0.0022127; eff.second.second = 0.0100000; } 
	} 
      } // end id==13 
    } // end era==2015 

    return eff;
  } 


private: 
  int era_; 
}; 


#endif
