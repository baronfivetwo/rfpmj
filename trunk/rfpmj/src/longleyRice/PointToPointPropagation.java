package longleyRice;

public class PointToPointPropagation {
	public void  point_to_point(double elev[], double tht_m, double rht_m,
	          double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			  double frq_mhz, int radio_climate, int pol, double conf, double rel,
			  double &dbloss, char *strmode, int &errnum) // TODO: Fix Pointers
		// pol: 0-Horizontal, 1-Vertical
		// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
		//                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
		//                7-Maritime Temperate, Over Sea
		// conf, rel: .01 to .99
		// elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
		// errnum: 0- No Error.
		//         1- Warning: Some parameters are nearly out of range.
		//                     Results should be used with caution.
		//         2- Note: Default parameters have been substituted for impossible ones.
		//         3- Warning: A combination of parameters is out of range.
		//                     Results are probably invalid.
		//         Other-  Warning: Some parameters are out of range.
		//                          Results are probably invalid.
	{

	  prop_type   prop;
	  propv_type  propv;
	  propa_type  propa;
	  double zsys=0;
	  double zc, zr;
	  double eno, enso, q;
	  long ja, jb, i, np;
	  double dkm, xkm;
	  double fs;

	  prop.hg[0] = tht_m;   prop.hg[1] = rht_m;
	  propv.klim = radio_climate;
	  prop.kwx = 0;
	  propv.lvar = 5;
	  prop.mdp = -1;
	  zc = qerfi(conf);
	  zr = qerfi(rel);
	  np = (long)elev[0];
	  dkm = (elev[1] * elev[0]) / 1000.0;
	  xkm = elev[1] / 1000.0;
	  eno = eno_ns_surfref;
	  enso = 0.0;
	  q = enso;
	  if(q<=0.0)
	  {
	    ja = 3.0 + 0.1 * elev[0];
		jb = np - ja + 6;
		for(i=ja-1;i<jb;++i)
	      zsys+=elev[i];
	    zsys/=(jb-ja+1);
		q=eno;
	  }
	  propv.mdvar=12;
	  qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,prop);
	  qlrpfl(elev,propv.klim,propv.mdvar,prop,propa,propv);
	  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
	  q = prop.dist - propa.dla;
	  if(int(q)<0.0)
	    strcpy(strmode,"Line-Of-Sight Mode");
	  else
	    { if(int(q)==0.0)
	        strcpy(strmode,"Single Horizon");
	      else if(int(q)>0.0)
	        strcpy(strmode,"Double Horizon");
	      if(prop.dist<=propa.dlsa || prop.dist <= propa.dx)
	        strcat(strmode,", Diffraction Dominant");
	      else if(prop.dist>propa.dx)
	        strcat(strmode, ", Troposcatter Dominant");
	    }
	  dbloss = avar(zr,0.0,zc,prop,propv) + fs;
	  errnum = prop.kwx;
	}


	public void point_to_pointMDH (double elev[], double tht_m, double rht_m,
	          double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			  double frq_mhz, int radio_climate, int pol, double timepct, double locpct, double confpct, 
			  double &dbloss, int &propmode, double &deltaH, int &errnum)
		// pol: 0-Horizontal, 1-Vertical
		// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
		//                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
		//                7-Maritime Temperate, Over Sea
		// timepct, locpct, confpct: .01 to .99
		// elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
		// propmode:  Value   Mode
		//             -1     mode is undefined
		//              0     Line of Sight
		//              5     Single Horizon, Diffraction
		//              6     Single Horizon, Troposcatter
		//              9     Double Horizon, Diffraction
		//             10     Double Horizon, Troposcatter
		// errnum: 0- No Error.
		//         1- Warning: Some parameters are nearly out of range.
		//                     Results should be used with caution.
		//         2- Note: Default parameters have been substituted for impossible ones.
		//         3- Warning: A combination of parameters is out of range.
		//                     Results are probably invalid.
		//         Other-  Warning: Some parameters are out of range.
		//                          Results are probably invalid.
	{

	  prop_type   prop;
	  propv_type  propv;
	  propa_type  propa;
	  double zsys=0;
	  double ztime, zloc, zconf;
	  double eno, enso, q;
	  long ja, jb, i, np;
	  double dkm, xkm;
	  double fs;

	  propmode = -1;  // mode is undefined
	  prop.hg[0] = tht_m;   prop.hg[1] = rht_m;
	  propv.klim = radio_climate;
	  prop.kwx = 0;
	  propv.lvar = 5;
	  prop.mdp = -1;
	  ztime = qerfi(timepct);
	  zloc = qerfi(locpct);
	  zconf = qerfi(confpct);

	  np = (long)elev[0];
	  dkm = (elev[1] * elev[0]) / 1000.0;
	  xkm = elev[1] / 1000.0;
	  eno = eno_ns_surfref;
	  enso = 0.0;
	  q = enso;
	  if(q<=0.0)
	  {
	    ja = 3.0 + 0.1 * elev[0];
		jb = np - ja + 6;
		for(i=ja-1;i<jb;++i)
	      zsys+=elev[i];
	    zsys/=(jb-ja+1);
		q=eno;
	  }
	  propv.mdvar=12;
	  qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,prop);
	  qlrpfl(elev,propv.klim,propv.mdvar,prop,propa,propv);
	  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
	  deltaH = prop.dh;
	  q = prop.dist - propa.dla;
	  if(int(q)<0.0)
	    propmode = 0;  // Line-Of-Sight Mode
	  else
	    { if(int(q)==0.0)
	        propmode = 4;  // Single Horizon
	      else if(int(q)>0.0)
	        propmode = 8;  // Double Horizon
	      if(prop.dist<=propa.dlsa || prop.dist <= propa.dx)
	        propmode += 1; // Diffraction Dominant
	      else if(prop.dist>propa.dx)
	        propmode += 2; // Troposcatter Dominant
	    }
	  dbloss = avar(ztime, zloc, zconf, prop, propv) + fs;      //avar(time,location,confidence)
	  errnum = prop.kwx;
	}

	public void point_to_pointDH (double elev[], double tht_m, double rht_m,
	          double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			  double frq_mhz, int radio_climate, int pol, double conf, double rel,
			  double &dbloss, double &deltaH, int &errnum)
		// pol: 0-Horizontal, 1-Vertical
		// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
		//                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
		//                7-Maritime Temperate, Over Sea
		// conf, rel: .01 to .99
		// elev[]: [num points - 1], [delta dist(meters)], [height(meters) point 1], ..., [height(meters) point n]
		// errnum: 0- No Error.
		//         1- Warning: Some parameters are nearly out of range.
		//                     Results should be used with caution.
		//         2- Note: Default parameters have been substituted for impossible ones.
		//         3- Warning: A combination of parameters is out of range.
		//                     Results are probably invalid.
		//         Other-  Warning: Some parameters are out of range.
		//                          Results are probably invalid.
	{

	  char strmode[100];
	  prop_type   prop;
	  propv_type  propv;
	  propa_type  propa;
	  double zsys=0;
	  double zc, zr;
	  double eno, enso, q;
	  long ja, jb, i, np;
	  double dkm, xkm;
	  double fs;

	  prop.hg[0] = tht_m;   prop.hg[1] = rht_m;
	  propv.klim = radio_climate;
	  prop.kwx = 0;
	  propv.lvar = 5;
	  prop.mdp = -1;
	  zc = qerfi(conf);
	  zr = qerfi(rel);
	  np = (long)elev[0];
	  dkm = (elev[1] * elev[0]) / 1000.0;
	  xkm = elev[1] / 1000.0;
	  eno = eno_ns_surfref;
	  enso = 0.0;
	  q = enso;
	  if(q<=0.0)
	  {
	    ja = 3.0 + 0.1 * elev[0];
		jb = np - ja + 6;
		for(i=ja-1;i<jb;++i)
	      zsys+=elev[i];
	    zsys/=(jb-ja+1);
		q=eno;
	  }
	  propv.mdvar=12;
	  qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,prop);
	  qlrpfl(elev,propv.klim,propv.mdvar,prop,propa,propv);
	  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop.dist / 1000.0);
	  deltaH = prop.dh;
	  q = prop.dist - propa.dla;
	  if(int(q)<0.0)
	    strcpy(strmode,"Line-Of-Sight Mode");
	  else
	    { if(int(q)==0.0)
	        strcpy(strmode,"Single Horizon");
	      else if(int(q)>0.0)
	        strcpy(strmode,"Double Horizon");
	      if(prop.dist<=propa.dlsa || prop.dist <= propa.dx)
	        strcat(strmode,", Diffraction Dominant");
	      else if(prop.dist>propa.dx)
	        strcat(strmode, ", Troposcatter Dominant");
	    }
	  dbloss = avar(zr,0.0,zc,prop,propv) + fs;      //avar(time,location,confidence)
	  errnum = prop.kwx;
	}


}