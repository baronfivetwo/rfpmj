package longleyRice;

public class propzgndType {
	private double zgndreal;
	private double zgndimag;

	protected propzgndType(double zgndreal, double zgndimag) {
		this.zgndreal = zgndreal;
		this.zgndimag = zgndimag;

	}

	protected double getZgndreal() {
		return zgndreal;
	}

	protected void setZgndreal(double zgndreal) {
		this.zgndreal = zgndreal;
	}

	protected double getZgndimag() {
		return zgndimag;
	}

	protected void setZgndimag(double zgndimag) {
		this.zgndimag = zgndimag;
	}
}
