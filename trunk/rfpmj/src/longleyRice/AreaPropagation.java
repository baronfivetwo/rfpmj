package longleyRice;

public class AreaPropagation {
	//********************************************************
	//* Area Mode Calculations                               *
	//********************************************************

	public void area(long ModVar, double deltaH, double tht_m, double rht_m,
			  double dist_km, int TSiteCriteria, int RSiteCriteria, 
	          double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			  double frq_mhz, int radio_climate, int pol, double pctTime, double pctLoc,
			  double pctConf, double &dbloss, char *strmode, int &errnum) // TODO: fix pointers
	{
		// pol: 0-Horizontal, 1-Vertical
		// TSiteCriteria, RSiteCriteria:
		//		   0 - random, 1 - careful, 2 - very careful
		// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
		//                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
		//                7-Maritime Temperate, Over Sea
		// ModVar: 0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
	    //         1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
	    //         2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
	    //         3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
		// pctTime, pctLoc, pctConf: .01 to .99
		// errnum: 0- No Error.
		//         1- Warning: Some parameters are nearly out of range.
		//                     Results should be used with caution.
		//         2- Note: Default parameters have been substituted for impossible ones.
		//         3- Warning: A combination of parameters is out of range.
		//                     Results are probably invalid.
		//         Other-  Warning: Some parameters are out of range.
		//                          Results are probably invalid.
		// NOTE: strmode is not used at this time.
	  prop_type prop;
	  propv_type propv;
	  propa_type propa;
	  double zt, zl, zc, xlb;
	  double fs;
	  long ivar;
	  double eps, eno, sgm;
	  long ipol;
	  int kst[2];

	  kst[0] = (int) TSiteCriteria;
	  kst[1] = (int) RSiteCriteria;
	  zt = qerfi(pctTime);
	  zl = qerfi(pctLoc);
	  zc = qerfi(pctConf);
	  eps = eps_dielect;
	  sgm = sgm_conductivity;
	  eno = eno_ns_surfref;
	  prop.dh = deltaH;
	  prop.hg[0] = tht_m;  prop.hg[1] = rht_m;
	  propv.klim = (__int32) radio_climate;
	  prop.ens = eno;
	  prop.kwx = 0;
	  ivar = (long) ModVar;
	  ipol = (long) pol;
	  qlrps(frq_mhz, 0.0, eno, ipol, eps, sgm, prop);
	  qlra(kst, propv.klim, ivar, prop, propv);
	  if(propv.lvar<1) propv.lvar = 1;
	  lrprop(dist_km * 1000.0, prop, propa);
	  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
	  xlb = fs + avar(zt, zl, zc, prop, propv);
	  dbloss = xlb;
	  if(prop.kwx==0)
		errnum = 0;
	  else
	    errnum = prop.kwx;
	}


	public double ITMAreadBLoss(long ModVar, double deltaH, double tht_m, double rht_m,
			  double dist_km, int TSiteCriteria, int RSiteCriteria, 
	          double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			  double frq_mhz, int radio_climate, int pol, double pctTime, double pctLoc,
			  double pctConf)
	{
		char strmode[200];
		int errnum;
		double dbloss;
		area(ModVar,deltaH,tht_m,rht_m,dist_km,TSiteCriteria,RSiteCriteria, 
	          eps_dielect,sgm_conductivity,eno_ns_surfref,
			  frq_mhz,radio_climate,pol,pctTime,pctLoc,
			  pctConf,dbloss,strmode,errnum);
		return dbloss;
	}

	public double  ITMVersion()
	{
		return 7.0;
	}
}
