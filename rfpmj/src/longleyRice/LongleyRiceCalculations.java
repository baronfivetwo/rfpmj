package longleyRice;

import org.apache.commons.math.complex.Complex;

//TODO: Fix All Pass By reference 
public class LongleyRiceCalculations {

	private static final double THIRD = 1.0 / 3.0;

	private int mymin(int i, int j) {
		if (i < j)
			return i;
		else
			return j;
	}

	private int mymax(int i, int j) {
		if (i > j)
			return i;
		else
			return j;
	}

	private double mymin(double a, double b) {
		if (a < b)
			return a;
		else
			return b;
	}

	private double mymax(double a, double b) {
		if (a > b)
			return a;
		else
			return b;
	}

	private double FORTRAN_DIM(double x, double y) { // This performs the
														// FORTRAN DIM function.
														// result is x-y if x is
														// greater than y;
														// otherwise result is
														// 0.0
		if (x > y)
			return x - y;
		else
			return 0.0;
	}

	// TODO: Check base of log in original Source Code for now assume base e
	private double aknfe(double v2) {
		double a;
		if (v2 < 5.76)
			a = 6.02 + 9.11 * Math.sqrt(v2) - 1.27 * v2;
		else
			a = 12.953 + 4.343 * Math.log(v2);
		return a;
	}

	// TODO: Check if pow and exp is the same in java as in cpp
	private double fht(double x, double pk) {
		double w, fhtv;
		if (x < 200.0) {
			w = -Math.log(pk);
			if (pk < 1e-5 || x * Math.pow(w, 3.0) > 5495.0) {
				fhtv = -117.0;
				if (x > 1.0)
					fhtv = 17.372 * Math.log(x) + fhtv;
			} else
				fhtv = 2.5e-5 * x * x / pk - 8.686 * w - 15.0;
		} else {
			fhtv = 0.05751 * x - 4.343 * Math.log(x);
			if (x < 2000.0) {
				w = 0.0134 * x * Math.exp(-0.005 * x);
				fhtv = (1.0 - w) * fhtv + w * (17.372 * Math.log(x) - 117.0);
			}
		}
		return fhtv;
	}

	private double h0f(double r, double et) {
		double[] a = { 25.0, 80.0, 177.0, 395.0, 705.0 };
		double[] b = { 24.0, 45.0, 68.0, 80.0, 105.0 };
		double q, x;
		int it;
		double h0fv;
		it = (int) et;
		if (it <= 0) {
			it = 1;
			q = 0.0;
		} else if (it >= 5) {
			it = 5;
			q = 0.0;
		} else
			q = et - it;
		x = Math.pow(1.0 / r, 2.0);
		h0fv = 4.343 * Math.log((a[it - 1] * x + b[it - 1]) * x + 1.0);
		if (q != 0.0)
			h0fv = (1.0 - q) * h0fv + q * 4.343
					* Math.log((a[it] * x + b[it]) * x + 1.0);
		return h0fv;
	}

	private double ahd(double td) {
		int i;
		double[] a = { 133.4, 104.6, 71.8 };
		double[] b = { 0.332e-3, 0.212e-3, 0.157e-3 };
		double[] c = { -4.343, -1.086, 2.171 };
		if (td <= 10e3)
			i = 0;
		else if (td <= 70e3)
			i = 1;
		else
			i = 2;
		return a[i] + b[i] * td + c[i] * Math.log(td);
	}

	private double adiff(double d, propType prop, propaType propa) {
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag()); // complex<double>
																				// prop_zgnd(prop.zgndreal,prop.zgndimag);
		double wd1 = 0, xd1 = 0, afo = 0, qk = 0, aht = 0, xht = 0;
		double a, q, pk, ds, th, wa, ar, wd, adiffv;
		if (d == 0) {
			double[] hg = prop.getHg();
			double[] he = prop.getHg();
			q = hg[0] * hg[1];
			qk = he[0] * he[1] - q;
			if (prop.getMdp() < 0.0)
				q += 10.0;
			wd1 = Math.sqrt(1.0 + qk / q);
			xd1 = propa.getDla() + propa.getTha() / prop.getGme();
			q = (1.0 - 0.8 * Math.exp(-propa.getDlsa() / 50e3)) * prop.getDh();
			// TODO: Check meaning of *= in cpp to java also check ITM ref
			q *= 0.78 * Math.exp(-Math.pow(q / 16.0, 0.25));
			afo = mymin(
					15.0,
					2.171 * Math.log(1.0 + 4.77e-4 * hg[0] * hg[1]
							* prop.getWn() * q));
			qk = 1.0 / propzgnd.abs();
			aht = 20.0;
			xht = 0.0;
			double[] dl = prop.getDl();
			for (int j = 0; j < 2; ++j) {
				a = 0.5 * Math.pow(dl[j], 2.0) / he[j];
				wa = Math.pow(a * prop.getWn(), THIRD);
				pk = qk / wa;
				q = (1.607 - pk) * 151.0 * wa * dl[j] / a;
				xht += q;
				aht += fht(q, pk);
			}
			adiffv = 0.0;
		} else {
			double[] dl = prop.getDl();
			th = propa.getTha() + d * prop.getGme();
			ds = d - propa.getDla();
			q = 0.0795775 * prop.getWn() * ds * Math.pow(th, 2.0);
			adiffv = aknfe(q * dl[0] / (ds + dl[0]))
					+ aknfe(q * dl[1] / (ds + dl[1]));
			a = ds / th;
			wa = Math.pow(a * prop.getWn(), THIRD);
			pk = qk / wa;
			q = (1.607 - pk) * 151.0 * wa * th + xht;
			ar = 0.05751 * q - 4.343 * Math.log(q) - aht;
			q = (wd1 + xd1 / d)
					* mymin(((1.0 - 0.8 * Math.exp(-d / 50e3)) * prop.getDh() * prop
							.getWn()), 6283.2);
			wd = 25.1 / (25.1 + Math.sqrt(q));
			adiffv = ar * wd + (1.0 - wd) * adiffv + afo;
		}
		return adiffv;
	}

	private double ascat(double d, propType prop, propaType propa) {
		@SuppressWarnings("unused")
		// TODO: propzgnd is needed?
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag()); // complex<double>
																				// prop_zgnd(prop.zgndreal,prop.zgndimag);
		double ad = 0, rr = 0, etq = 0, h0s = 0;
		double h0, r1, r2, z0, ss, et, ett, th, q;
		double ascatv;
		double[] dl = prop.getDl();
		double[] he = prop.getHe();
		double[] the = prop.getThe();
		if (d == 0.0) {
			ad = dl[0] - dl[1];
			rr = he[1] / he[0];
			if (ad < 0.0) {
				ad = -ad;
				rr = 1.0 / rr;
			}
			etq = (5.67e-6 * prop.getEns() - 2.32e-3) * prop.getEns() + 0.031;
			h0s = -15.0;
			ascatv = 0.0;
		} else {
			if (h0s > 15.0)
				h0 = h0s;
			else {
				th = the[0] + the[1] + d * prop.getGme();
				r2 = 2.0 * prop.getWn() * th;
				r1 = r2 * he[0];
				r2 *= he[1];
				if (r1 < 0.2 && r2 < 0.2)
					return 1001.0; // <==== early return
				ss = (d - ad) / (d + ad);
				q = rr / ss;
				ss = mymax(0.1, ss);
				q = mymin(mymax(0.1, q), 10.0);
				z0 = (d - ad) * (d + ad) * th * 0.25 / d;
				et = (etq * Math.exp(-Math.pow(mymin(1.7, z0 / 8.0e3), 6.0)) + 1.0)
						* z0 / 1.7556e3;
				ett = mymax(et, 1.0);
				h0 = (h0f(r1, ett) + h0f(r2, ett)) * 0.5;
				h0 += mymin(h0,
						(1.38 - Math.log(ett)) * Math.log(ss) * Math.log(q)
								* 0.49);
				h0 = FORTRAN_DIM(h0, 0.0);
				if (et < 1.0)
					h0 = et
							* h0
							+ (1.0 - et)
							* 4.343
							* Math.log(Math.pow((1.0 + 1.4142 / r1)
									* (1.0 + 1.4142 / r2), 2.0)
									* (r1 + r2) / (r1 + r2 + 2.8284));
				if (h0 > 15.0 && h0s >= 0.0)
					h0 = h0s;
			}
			h0s = h0;
			th = propa.getTha() + d * prop.getGme();
			ascatv = ahd(th * d) + 4.343
					* Math.log(47.7 * prop.getWn() * Math.pow(th, 4.0)) - 0.1
					* (prop.getEns() - 301.0) * Math.exp(-th * d / 40e3) + h0;
		}
		return ascatv;
	}

	protected double qerfi(double q) {
		double x, t, v;
		double c0 = 2.515516698;
		double c1 = 0.802853;
		double c2 = 0.010328;
		double d1 = 1.432788;
		double d2 = 0.189269;
		double d3 = 0.001308;

		x = 0.5 - q;
		t = mymax(0.5 - Math.abs(x), 0.000001); // TODO: Check if fabs == to
												// Math.abs cast to float
		t = Math.sqrt(-2.0 * Math.log(t));
		v = t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0);
		if (x < 0.0)
			v = -v;
		return v;
	}

	protected void qlrps(double fmhz, double zsys, double en0, int ipol,
			double eps, double sgm, propType prop) {
		double gma = 157e-9;
		prop.setWn(fmhz / 47.7);
		prop.setEns(en0);
		double ens = 0;
		if (zsys != 0.0)
			ens = prop.getEns();
		ens *= Math.exp(-zsys / 9460.0);
		prop.setEns(ens);
		prop.setGme(gma * (1.0 - 0.04665 * Math.exp(prop.getEns() / 179.3)));
		Complex zq;
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag());
		zq = new Complex(eps, 376.62 * sgm / prop.getWn());
		Complex temp = new Complex(1.0, 0); //
		zq.subtract(temp); // Hopefully equivalent to prop_zgnd=sqrt(zq-1.0);
		propzgnd = zq.sqrt(); //
		if (ipol != 0.0)
			propzgnd = propzgnd.divide(zq);

		prop.setZgndreal(propzgnd.getReal());
		prop.setZgndimag(propzgnd.getImaginary());
	}

	private double abq_alos(Complex r) {
		return r.getReal() * r.getReal() + r.getImaginary() * r.getImaginary();
	}

	private double alos(double d, propType prop, propaType propa) {
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag());
		double wls = 0;
		Complex r;
		double s, sps, q;
		double alosv;
		if (d == 0.0) {
			wls = 0.021 / (0.021 + prop.getWn() * prop.getDh()
					/ mymax(10e3, propa.getDlsa()));
			alosv = 0.0;
		} else {
			q = (1.0 - 0.8 * Math.exp(-d / 50e3)) * prop.getDh();
			s = 0.78 * q * Math.exp(-Math.pow(q / 16.0, 0.25));
			double[] he = prop.getHe();
			q = he[0] + he[1];
			sps = q / Math.sqrt(d * d + q * q);
			Complex Csps = new Complex(sps, 0.0);
			r = (Csps.subtract(propzgnd)).divide((Csps.add(propzgnd)))
					.multiply(Math.exp(-mymin(10.0, prop.getWn() * s * sps))); // original
																				// :
																				// r=(sps-prop_zgnd)/(sps+prop_zgnd)*exp(-mymin(10.0,prop.wn*s*sps));
			q = abq_alos(r);
			if (q < 0.25 || q < sps)
				r = r.multiply(Math.sqrt(sps / q));
			alosv = propa.getEmd() * d + propa.getAed();
			q = prop.getWn() * he[0] * he[1] * 2.0 / d;
			if (q > 1.57)
				q = 3.14 - 2.4649 / q;
			alosv = (-4.343
					* Math.log(abq_alos(new Complex(Math.cos(q), -Math.sin(q))
							.add(r)) - alosv) * wls + alosv);
		}
		return alosv;
	}

	protected void qlra(int kst[], int klimx, int mdvarx, propType prop,
			propvType propv) {
		@SuppressWarnings("unused")
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag());
		double q;
		double[] he = prop.getHe();
		double[] hg = prop.getHg();
		for (int j = 0; j < 2; ++j) {
			if (kst[j] <= 0)
				he[j] = hg[j];
			else {
				q = 4.0;
				if (kst[j] != 1)
					q = 9.0;
				if (hg[j] < 5.0)
					q *= Math.sin(0.3141593 * hg[j]);
				he[j] = hg[j]
						+ (1.0 + q)
						* Math.exp(-mymin(20.0,
								2.0 * hg[j] / mymax(1e-3, prop.getDh())));
				prop.setHe(he);
			}
			q = Math.sqrt(2.0 * he[j] / prop.getGme());
			double[] dl = prop.getDl();
			dl[j] = q
					* Math.exp(-0.07
							* Math.sqrt(prop.getDh() / mymax(he[j], 5.0)));
			double[] the = prop.getThe();
			the[j] = (0.65 * prop.getDh() * (q / dl[j] - 1.0) - 2.0 * he[j])
					/ q;
			prop.setDl(dl);
			prop.setThe(the);
		}
		prop.setMdp(1);
		propv.setLvar(mymax(propv.getLvar(), 3));
		if (mdvarx >= 0) {
			propv.setMdvar(mdvarx);
			propv.setLvar(mymax(propv.getLvar(), 4));
		}
		if (klimx > 0) {
			propv.setKlim(klimx);
			propv.setLvar(5);
		}
	}

	protected void lrprop(double d, propType prop, propaType propa) // PaulM_lrprop
	{
		boolean wlos = false, wscat = false;
		double dmin = 0, xae = 0;
		Complex propzgnd = new Complex(prop.getZgndreal(), prop.getZgndimag());
		double a0, a1, a2, a3, a4, a5, a6;
		double d0, d1, d2, d3, d4, d5, d6;
		boolean wq;
		double q;
		int j;

		if (prop.getMdp() != 0) {
			double[] dls = propa.getDls();
			double[] he = prop.getHe();
			for (j = 0; j < 2; j++) {
				dls[j] = Math.sqrt(2.0 * he[j] / prop.getGme());
			}
			propa.setDls(dls);
			double[] hg = prop.getHg();
			double[] dl = prop.getDl();
			double[] the = prop.getThe();

			propa.setDlsa(dls[0] + dls[1]);
			propa.setDla(dl[0] + dl[1]);
			propa.setTha(mymax(the[0] + the[1], -propa.getDla() * prop.getGme()));
			wlos = false;
			wscat = false;
			if (prop.getWn() < 0.838 || prop.getWn() > 210.0) {
				prop.setKwx(mymax(prop.getKwx(), 1));
			}
			for (j = 0; j < 2; j++)
				if (hg[j] < 1.0 || hg[j] > 1000.0) {
					prop.setKwx(mymax(prop.getKwx(), 1));
				}
			for (j = 0; j < 2; j++)
				if (Math.abs(the[j]) > 200e-3 || dl[j] < 0.1 * dls[j]
						|| dl[j] > 3.0 * dls[j]) {
					prop.setKwx(mymax(prop.getKwx(), 3));}
			if (prop.getEns() < 250.0 || prop.getEns() > 400.0
					|| prop.getGme() < 75e-9 || prop.getGme() > 250e-9
					|| propzgnd.getReal() <= Math.abs(propzgnd.getImaginary())
					|| prop.getWn() < 0.419 || prop.getWn() > 420.0) {
				prop.setKwx(4);
			}
			for (j = 0; j < 2; j++)
				// Experiment here see if getter will allow for getHG()[j]
				if (prop.getHg()[j] < 0.5 || prop.getHg()[j] > 3000.0) {
					prop.setKwx(4);
				}
			dmin = Math.abs(prop.getHe()[0] - prop.getHe()[1]) / 200e-3;
			q = adiff(0.0, prop, propa);
			xae = Math.pow(prop.getWn() * Math.pow(prop.getGme(), 2), -THIRD);
			d3 = mymax(propa.getDlsa(), 1.3787 * xae + propa.getDla());
			d4 = d3 + 2.7574 * xae;
			a3 = adiff(d3, prop, propa);
			a4 = adiff(d4, prop, propa);
			propa.setEmd((a4 - a3) / (d4 - d3));
			propa.setAed(a3 - propa.getEmd() * d3);
		}
		if (prop.getMdp() >= 0) {
			prop.setMdp(0);
			prop.setDist(d);
		}
		if (prop.getDist() > 0.0) {
			if (prop.getDist() > 1000e3) {
				prop.setKwx(mymax(prop.getKwx(), 1));
			}
			if (prop.getDist() < dmin) {
				prop.setKwx(mymax(prop.getKwx(), 3));
			}
			if (prop.getDist() < 1e3 || prop.getDist() > 2000e3) {
				prop.setKwx(4);
			}
		}
		if (prop.getDist() < propa.getDlsa()) {
			if (!wlos) {
				q = alos(0.0, prop, propa);
				d2 = propa.getDlsa();
				a2 = propa.getAed() + d2 * propa.getEmd();
				d0 = 1.908 * prop.getWn() * prop.getHe()[0] * prop.getHe()[1];
				if (propa.getAed() >= 0.0) {
					d0 = mymin(d0, 0.5 * propa.getDla());
					d1 = d0 + 0.25 * (propa.getDla() - d0);
				} else
					d1 = mymax(-propa.getAed() / propa.getEmd(),
							0.25 * propa.getDla());
				a1 = alos(d1, prop, propa);
				wq = false;
				if (d0 < d1) {
					a0 = alos(d0, prop, propa);
					q = Math.log(d2 / d0);
					propa.setAk2(mymax(0.0, ((d2 - d0) * (a1 - a0) - (d1 - d0)
							* (a2 - a0))
							/ ((d2 - d0) * Math.log(d1 / d0) - (d1 - d0) * q)));
					wq = propa.getAed() >= 0.0 || propa.getAk2() > 0.0;
					if (wq) {
						propa.setAk1((a2 - a0 - propa.getAk2() * q) / (d2 - d0));
						if (propa.getAk1() < 0.0) {
							propa.setAk1(0.0);
							propa.setAk2(FORTRAN_DIM(a2, a0) / q);
							if (propa.getAk2() == 0.0)
								propa.setAk1(propa.getEmd());
						}
					}
				}
				if (!wq) {
					propa.setAk1(FORTRAN_DIM(a2, a1) / (d2 - d1));
					propa.setAk2(0.0);
					if (propa.getAk1() == 0.0)
						propa.setAk1(propa.getEmd());
				}
				propa.setAel(a2 - propa.getAk1() * d2 - propa.getAk2()
						* Math.log(d2));
				wlos = true;
			}
			if (prop.getDist() > 0.0)
				prop.setAref(propa.getAel() + propa.getAk1() * prop.getDist()
						+ propa.getAk2() * Math.log(prop.getDist()));
		}
		if (prop.getDist() <= 0.0 || prop.getDist() >= propa.getDlsa()) {
			if (!wscat) {
				q = ascat(0.0, prop, propa);
				d5 = propa.getDla() + 200e3;
				d6 = d5 + 200e3;
				a6 = ascat(d6, prop, propa);
				a5 = ascat(d5, prop, propa);
				if (a5 < 1000.0) {
					propa.setEms((a6 - a5) / 200e3);
					propa.setDx(mymax(
							propa.getDlsa(),
							mymax(propa.getDla() + 0.3 * xae
									* Math.log(47.7 * prop.getWn()), (a5
									- propa.getAed() - propa.getEms() * d5)
									/ (propa.getEmd() - propa.getEms()))));
					propa.setAes((propa.getEmd() - propa.getEms())
							* propa.getDx() + propa.getAed());
				} else {
					propa.setEms(propa.getEmd());
					propa.setAes(propa.getAed());
					propa.setDx(10.e6);
				}
				wscat = true;
			}
			if (prop.getDist() > propa.getDx())
				prop.setAref(propa.getAes() + propa.getEms() * prop.getDist());
			else
				prop.setAref(propa.getAed() + propa.getEmd() * prop.getDist());
		}
		prop.setAref(mymax(prop.getAref(), 0.0));
	}


	private double curve(double c1, double c2, double x1, double x2, double x3,
			double de) {
		return (c1 + c2 / (1.0 + Math.pow((de - x2) / x3, 2.0)))
				* Math.pow(de / x1, 2.0) / (1.0 + Math.pow(de / x1, 2.0));
	}

	protected double avar(double zzt, double zzl, double zzc, propType prop,
			propvType propv) {
		int kdv = 0;
		double dexa = 0, de = 0, vmd = 0, vs0 = 0, sgl = 0, sgtm = 0, sgtp = 0, sgtd = 0, tgtd = 0, gm = 0, gp = 0, cv1 = 0, cv2 = 0, yv1 = 0, yv2 = 0, yv3 = 0, csm1 = 0, csm2 = 0, ysm1 = 0, ysm2 = 0, ysm3 = 0, csp1 = 0, csp2 = 0, ysp1 = 0, ysp2 = 0, ysp3 = 0, csd1 = 0, zd = 0, cfm1 = 0, cfm2 = 0, cfm3 = 0, cfp1 = 0, cfp2 = 0, cfp3 = 0;
		double[] bv1 = { -9.67, -0.62, 1.26, -9.21, -0.62, -0.39, 3.15 };
		double[] bv2 = { 12.7, 9.19, 15.5, 9.05, 9.19, 2.86, 857.9 };
		double[] xv1 = { 144.9e3, 228.9e3, 262.6e3, 84.1e3, 228.9e3, 141.7e3,
				2222.e3 };
		double[] xv2 = { 190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3,
				164.8e3 };
		double[] xv3 = { 133.8e3, 143.6e3, 99.8e3, 98.6e3, 143.6e3, 167.4e3,
				116.3e3 };
		double[] bsm1 = { 2.13, 2.66, 6.11, 1.98, 2.68, 6.86, 8.51 };
		double[] bsm2 = { 159.5, 7.67, 6.65, 13.11, 7.16, 10.38, 169.8 };
		double[] xsm1 = { 762.2e3, 100.4e3, 138.2e3, 139.1e3, 93.7e3, 187.8e3,
				609.8e3 };
		double[] xsm2 = { 123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3,
				119.9e3 };
		double[] xsm3 = { 94.5e3, 136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3,
				106.6e3 };
		double[] bsp1 = { 2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43 };
		double[] bsp2 = { 102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19 };
		double[] xsp1 = { 636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3,
				136.2e3 };
		double[] xsp2 = { 134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3,
				188.5e3 };
		double[] xsp3 = { 95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3,
				122.9e3 };
		double[] bsd1 = { 1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518 };
		double[] bzd1 = { 1.282, 2.161, 1.282, 20., 1.282, 1.282, 1.282 };
		double[] bfm1 = { 1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0 };
		double[] bfm2 = { 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0 };
		double[] bfm3 = { 0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0 };
		double[] bfp1 = { 1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0 };
		double[] bfp2 = { 0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0 };
		double[] bfp3 = { 0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0 };
		boolean ws = false, w1 = false;
		double rt = 7.8, rl = 24.0, avarv, q, vs, zt, zl, zc;
		double sgt, yr;
		int temp_klim = propv.getKlim() - 1;

		if (propv.getLvar() > 0) {
			switch (propv.getLvar()) {
			default:
				if (propv.getKlim() <= 0 || propv.getKlim() > 7) {
					propv.setKlim(5);
					temp_klim = 4;
					{
						prop.setKwx(mymax(prop.getKwx(), 2));
					}
				}
				cv1 = bv1[temp_klim];
				cv2 = bv2[temp_klim];
				yv1 = xv1[temp_klim];
				yv2 = xv2[temp_klim];
				yv3 = xv3[temp_klim];
				csm1 = bsm1[temp_klim];
				csm2 = bsm2[temp_klim];
				ysm1 = xsm1[temp_klim];
				ysm2 = xsm2[temp_klim];
				ysm3 = xsm3[temp_klim];
				csp1 = bsp1[temp_klim];
				csp2 = bsp2[temp_klim];
				ysp1 = xsp1[temp_klim];
				ysp2 = xsp2[temp_klim];
				ysp3 = xsp3[temp_klim];
				csd1 = bsd1[temp_klim];
				zd = bzd1[temp_klim];
				cfm1 = bfm1[temp_klim];
				cfm2 = bfm2[temp_klim];
				cfm3 = bfm3[temp_klim];
				cfp1 = bfp1[temp_klim];
				cfp2 = bfp2[temp_klim];
				cfp3 = bfp3[temp_klim];
			case 4:
				kdv = propv.getMdvar();
				ws = kdv >= 20;
				if (ws)
					kdv -= 20;
				w1 = kdv >= 10;
				if (w1)
					kdv -= 10;
				if (kdv < 0 || kdv > 3) {
					kdv = 0;
					prop.setKwx(mymax(prop.getKwx(), 2));
				}
			case 3:
				q = Math.log(0.133 * prop.getWn());
				gm = cfm1 + cfm2 / (Math.pow(cfm3 * q, 2.0) + 1.0);
				gp = cfp1 + cfp2 / (Math.pow(cfp3 * q, 2.0) + 1.0);
			case 2:
				dexa = Math.sqrt(18e6 * prop.getHe()[0])
						+ Math.sqrt(18e6 * prop.getHe()[1])
						+ Math.pow((575.7e12 / prop.getWn()), THIRD);
			case 1:
				if (prop.getDist() < dexa)
					de = 130e3 * prop.getDist() / dexa;
				else
					de = 130e3 + prop.getDist() - dexa;
			}
			vmd = curve(cv1, cv2, yv1, yv2, yv3, de);
			sgtm = curve(csm1, csm2, ysm1, ysm2, ysm3, de) * gm;
			sgtp = curve(csp1, csp2, ysp1, ysp2, ysp3, de) * gp;
			sgtd = sgtp * csd1;
			tgtd = (sgtp - sgtd) * zd;
			if (w1)
				sgl = 0.0;
			else {
				q = (1.0 - 0.8 * Math.exp(-prop.getDist() / 50e3))
						* prop.getDh() * prop.getWn();
				sgl = 10.0 * q / (q + 13.0);
			}
			if (ws)
				vs0 = 0.0;
			else
				vs0 = Math.pow(5.0 + 3.0 * Math.exp(-de / 100e3), 2.0);
			propv.setLvar(0);
		}
		zt = zzt;
		zl = zzl;
		zc = zzc;
		switch (kdv) {
		case 0:
			zt = zc;
			zl = zc;
			break;
		case 1:
			zl = zc;
			break;
		case 2:
			zl = zt;
		}

		if (Math.abs(zt) > 3.1 || Math.abs(zl) > 3.1 || Math.abs(zc) > 3.1) // Original
																			// Line
																			// fabs(zt)>3.1
																			// ||
																			// fabs(zl)>3.1
																			// ||
																			// fabs(zc)>3.1)
		{
			prop.setKwx(mymax(prop.getKwx(), 1));
		}
		if (zt < 0.0)
			sgt = sgtm;
		else if (zt <= zd)
			sgt = sgtp;
		else
			sgt = sgtd + tgtd / zt;
		vs = vs0 + Math.pow(sgt * zt, 2.0) / (rt + zc * zc)
				+ Math.pow(sgl * zl, 2.0) / (rl + zc * zc);
		if (kdv == 0) {
			yr = 0.0;
			propv.setSgc(Math.sqrt(sgt * sgt + sgl * sgl + vs));
		} else if (kdv == 1) {
			yr = sgt * zt;
			propv.setSgc(Math.sqrt(sgl * sgl + vs));
		} else if (kdv == 2) {
			yr = Math.sqrt(sgt * sgt + sgl * sgl) * zt;
			propv.setSgc(Math.sqrt(vs));
		} else {
			yr = sgt * zt + sgl * zl;
			propv.setSgc(Math.sqrt(vs));
		}
		avarv = prop.getAref() - vmd - yr - propv.getSgc() * zc;
		if (avarv < 0.0)
			avarv = avarv * (29.0 - avarv) / (29.0 - 10.0 * avarv);
		return avarv;

	}

	private void hzns(double pfl[], propType prop) {
		boolean wq;
		int np;
		double xi, za, zb, qc, q, sb, sa;

		np = (int) pfl[0];
		xi = pfl[1];
		za = pfl[2] + prop.getHg()[0];
		zb = pfl[np + 2] + prop.getHg()[1];
		qc = 0.5 * prop.getGme();
		q = qc * prop.getDist();
		double[] the = prop.getThe();
		the[1] = (zb - za) / prop.getDist();
		the[0] = the[1] - q;
		the[1] = -the[1] - q;
		prop.setThe(the);
		double[] dl = prop.getDl();
		dl[0] = prop.getDist();
		dl[1] = prop.getDist();
		prop.setDl(dl);
		if (np >= 2) {
			sa = 0.0;
			sb = prop.getDist();
			wq = true;
			for (int i = 1; i < np; i++) {
				sa += xi;
				sb -= xi;
				q = pfl[i + 2] - (qc * sa + prop.getThe()[0]) * sa - za;
				if (q > 0.0) {
					double[] the_temp = prop.getThe();
					double[] dl_temp = prop.getDl();
					the_temp[0] += q / sa;
					dl_temp[0] = sa;
					prop.setThe(the_temp);
					prop.setDl(dl_temp);
					wq = false;
				}
				if (!wq) {
					q = pfl[i + 2] - (qc * sb + prop.getThe()[1]) * sb - zb;
					if (q > 0.0) {
						double[] the_temp = prop.getThe();
						double[] dl_temp = prop.getDl();
						the_temp[1] += q / sb;
						dl_temp[1] = sb;
						prop.setThe(the_temp);
						prop.setDl(dl_temp);
					}
				}
			}
		}
	}

	private void z1sq1(double[] z, double x1, double x2, double z0, double zn) {
		double xn, xa, xb, x, a, b;
		int n, ja, jb;
		xn = z[0];
		xa = (int) FORTRAN_DIM(x1 / z[1], 0.0);
		xb = xn - (int) FORTRAN_DIM(xn, x2 / z[1]);
		if (xb <= xa) {
			xa = FORTRAN_DIM(xa, 1.0);
			xb = xn - FORTRAN_DIM(xn, xb + 1.0);
		}
		ja = (int) xa;
		jb = (int) xb;
		n = jb - ja;
		xa = xb - xa;
		x = -0.5 * xa;
		xb += x;
		a = 0.5 * (z[ja + 2] + z[jb + 2]);
		b = 0.5 * (z[ja + 2] - z[jb + 2]) * x;
		for (int i = 2; i <= n; ++i) {
			++ja;
			x += 1.0;
			a += z[ja + 2];
			b += z[ja + 2] * x;
		}
		a /= xa;
		b = b * 12.0 / ((xa * xa + 2.0) * xa);
		z0 = a - b * xb;
		zn = a + b * (xn - xb);
	}

	private double qtile(int nn, double[] a, int ir) {
		double q = 0, r;
		int m, n, i, j, j1 = 0, i0 = 0, k;
		boolean done = false;
		boolean goto10 = true;

		m = 0;
		n = nn;
		k = mymin(mymax(0, ir), n);
		while (!done) {
			if (goto10) {
				q = a[k];
				i0 = m;
				j1 = n;
			}
			i = i0;
			while (i <= n && a[i] >= q)
				i++;
			if (i > n)
				i = n;
			j = j1;
			while (j >= m && a[j] <= q)
				j--;
			if (j < m)
				j = m;
			if (i < j) {
				r = a[i];
				a[i] = a[j];
				a[j] = r;
				i0 = i + 1;
				j1 = j - 1;
				goto10 = false;
			} else if (i < k) {
				a[k] = a[i];
				a[i] = q;
				m = i + 1;
				goto10 = true;
			} else if (j > k) {
				a[k] = a[j];
				a[j] = q;
				n = j - 1;
				goto10 = true;
			} else
				done = true;
		}
		return q;
	}

	@SuppressWarnings("unused")
	private double qerf(double z) {
		double b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937;
		double b4 = -1.821255987, b5 = 1.330274429;
		double rp = 4.317008, rrt2pi = 0.398942280;
		double t, x, qerfv;
		x = z;
		t = Math.abs(x);
		if (t >= 10.0)
			qerfv = 0.0;
		else {
			t = rp / (t + rp);
			qerfv = Math.exp(-0.5 * x * x) * rrt2pi
					* ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
		}
		if (x < 0.0)
			qerfv = 1.0 - qerfv;
		return qerfv;
	}

	private double d1thx(double pfl[], double x1, double x2) {
		int np, ka, kb, n, k, j;
		double d1thxv, sn, xa, xb;
		double[] s = { 0, 0, 0 };

		np = (int) pfl[0];
		xa = x1 / pfl[1];
		xb = x2 / pfl[1];
		d1thxv = 0.0;
		if (xb - xa < 2.0) // exit out
			return d1thxv;
		ka = (int) (0.1 * (xb - xa + 8.0));
		ka = mymin(mymax(4, ka), 25);
		n = 10 * ka - 5;
		kb = n - ka + 1;
		sn = n - 1;
		s[0] = sn;
		s[1] = 1.0;
		xb = (xb - xa) / sn;
		k = (int) (xa + 1.0);
		xa -= k;
		for (j = 0; j < n; j++) {
			while (xa > 0.0 && k < np) {
				xa -= 1.0;
				++k;
			}
			s[j + 2] = pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xa;
			xa = xa + xb;
		}
		z1sq1(s, 0.0, sn, xa, xb);
		xb = (xb - xa) / sn;
		for (j = 0; j < n; j++) {
			s[j + 2] -= xa;
			xa = xa + xb;
		}
		d1thxv = qtile(n - 1, s, ka - 1) - qtile(n - 1, s, kb - 1);
		d1thxv /= 1.0 - 0.8 * Math.exp(-(x2 - x1) / 50.0e3);
		return d1thxv;
	}

	@SuppressWarnings("null")
	// Don't know if this is wise
	protected void qlrpfl(double[] pfl, int klimx, int mdvarx, propType prop,
			propaType propa, propvType propv) {
		int np, j;
		double[] xl = null;
		double q = 0, za = 0, zb = 0;

		prop.setDist(pfl[0] * pfl[1]);
		np = (int) pfl[0];
		hzns(pfl, prop);
		for (j = 0; j < 2; j++)
			xl[j] = mymin(15.0 * prop.getHg()[j], 0.1 * prop.getDl()[j]);
		xl[1] = prop.getDist() - xl[1];
		prop.setDh(d1thx(pfl, xl[0], xl[1]));
		if (prop.getDl()[0] + prop.getDl()[1] > 1.5 * prop.getDist()) {
			z1sq1(pfl, xl[0], xl[1], za, zb);
			double[] he_temp = prop.getHe();
			he_temp[0] = prop.getHg()[0] + FORTRAN_DIM(pfl[2], za);
			he_temp[1] = prop.getHg()[1] + FORTRAN_DIM(pfl[np + 2], zb);
			prop.setHe(he_temp);
			for (j = 0; j < 2; j++) {
				double[] dl_temp = prop.getDl();
				dl_temp[j] = Math.sqrt(2.0 * prop.getHe()[j] / prop.getGme())
						* Math.exp(-0.07
								* Math.sqrt(prop.getDh()
										/ mymax(prop.getHe()[j], 5.0)));
				prop.setDl(dl_temp);
			}
			q = prop.getDl()[0] + prop.getDl()[1];

			if (q <= prop.getDist()) {
				q = Math.pow(prop.getDist() / q, 2.0);
				for (j = 0; j < 2; j++) {
					double[] dl_temp = prop.getDl();
					double[] he_temp1 = prop.getHe();// messy but needed scope
														// conflict
					he_temp1[j] *= q;
					prop.setHe(he_temp1);
					dl_temp[j] = Math.sqrt(2.0 * prop.getHe()[j]
							/ prop.getGme())
							* Math.exp(-0.07
									* Math.sqrt(prop.getDh()
											/ mymax(prop.getHe()[j], 5.0)));
					prop.setDl(dl_temp);
				}
			}
			for (j = 0; j < 2; j++) {
				q = Math.sqrt(2.0 * prop.getHe()[j] / prop.getGme());
				double[] the_temp = prop.getThe();
				the_temp[j] = (0.65 * prop.getDh()
						* (q / prop.getDl()[j] - 1.0) - 2.0 * prop.getHe()[j])
						/ q;
				prop.setThe(the_temp);
			}
		} else {
			z1sq1(pfl, xl[0], 0.9 * prop.getDl()[0], za, q);
			z1sq1(pfl, prop.getDist() - 0.9 * prop.getDl()[1], xl[1], q, zb);
			double[] he_temp = prop.getHe();
			he_temp[0] = prop.getHg()[0] + FORTRAN_DIM(pfl[2], za);
			he_temp[1] = prop.getHg()[1] + FORTRAN_DIM(pfl[np + 2], zb);
		}
		prop.setMdp(-1);
		propv.setLvar(mymax(propv.getLvar(), 3));
		if (mdvarx >= 0) {
			propv.setMdvar(mdvarx);
			propv.setLvar(mymax(propv.getLvar(), 4));
		}
		if (klimx > 0) {
			propv.setKlim(klimx);
			propv.setLvar(5);
		}
		lrprop(0.0, prop, propa);
	}

	@SuppressWarnings("unused")
	private double deg2rad(double d) {
		return d * 3.1415926535897 / 180.0;
	}

}
