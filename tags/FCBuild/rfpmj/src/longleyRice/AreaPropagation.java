package longleyRice;


public class AreaPropagation {
	// ********************************************************
	// * Area Mode Calculations *
	// ********************************************************
	// TODO: Fix All Pass By reference
	public static void area(long ModVar, double deltaH, double tht_m,
			double rht_m, double dist_km, int TSiteCriteria, int RSiteCriteria,
			double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			double frq_mhz, int radio_climate, int pol, double pctTime,
			double pctLoc, double pctConf, double dbloss, String strmode,
			int errnum) // TODO: fix pointers
	{
		// pol: 0-Horizontal, 1-Vertical
		// TSiteCriteria, RSiteCriteria:
		// 0 - random, 1 - careful, 2 - very careful
		// radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime
		// Tropical,
		// 4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
		// 7-Maritime Temperate, Over Sea
		// ModVar: 0 - Single: pctConf is "Time/Situation/Location", pctTime,
		// pctLoc not used
		// 1 - Individual: pctTime is "Situation/Location", pctConf is
		// "Confidence", pctLoc not used
		// 2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is
		// "Confidence", pctLoc not used
		// 3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is
		// "Confidence"
		// pctTime, pctLoc, pctConf: .01 to .99
		// errnum: 0- No Error.
		// 1- Warning: Some parameters are nearly out of range.
		// Results should be used with caution.
		// 2- Note: Default parameters have been substituted for impossible
		// ones.
		// 3- Warning: A combination of parameters is out of range.
		// Results are probably invalid.
		// Other- Warning: Some parameters are out of range.
		// Results are probably invalid.
		// NOTE: strmode is not used at this time.
		propType prop = new propType();
		propvType propv = new propvType();
		propaType propa = new propaType();
		double zt, zl, zc, xlb;
		double fs;
		int ivar;
		double eps, eno, sgm;
		int ipol;
		int[] kst = { 0, 0 };

		kst[0] = (int) TSiteCriteria;
		kst[1] = (int) RSiteCriteria;
		LongleyRiceCalculations calc = new LongleyRiceCalculations();
		zt = calc.qerfi(pctTime);
		zl = calc.qerfi(pctLoc);
		zc = calc.qerfi(pctConf);
		eps = eps_dielect;
		sgm = sgm_conductivity;
		eno = eno_ns_surfref;
		prop.setDh(deltaH);
		double[] hg = prop.getHg();
		hg[0] = tht_m;
		hg[1] = rht_m;
		prop.setHg(hg);
		propv.setKlim((int) radio_climate);
		prop.setEns(eno);
		prop.setKwx(0);
		ivar = (int) ModVar;
		ipol = (int) pol;
		calc.qlrps(frq_mhz, 0.0, eno, ipol, eps, sgm, prop);
		calc.qlra(kst, propv.getKlim(), ivar, prop, propv);
		if (propv.getLvar() < 1)
			propv.setLvar(1);
		calc.lrprop(dist_km * 1000.0, prop, propa);
		fs = 32.45 + 20.0 * Math.log10(frq_mhz) + 20.0
				* Math.log10(prop.getDist() / 1000.0);
		xlb = fs + calc.avar(zt, zl, zc, prop, propv);
		dbloss = xlb;
		if (prop.getKwx() == 0)
			errnum = 0;
		else
			errnum = prop.getKwx();
	}

	public double ITMAreadBLoss(long ModVar, double deltaH, double tht_m,
			double rht_m, double dist_km, int TSiteCriteria, int RSiteCriteria,
			double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
			double frq_mhz, int radio_climate, int pol, double pctTime,
			double pctLoc, double pctConf) {
		String strmode = null;
		int errnum = 0;
		double dbloss = 0;
		longleyRice.AreaPropagation.area(ModVar, deltaH, tht_m, rht_m, dist_km,
				TSiteCriteria, RSiteCriteria, eps_dielect, sgm_conductivity,
				eno_ns_surfref, frq_mhz, radio_climate, pol, pctTime, pctLoc,
				pctConf, dbloss, strmode, errnum);
		return dbloss;
	}

	public double ITMVersion() {
		return 7.0;
	}
}
