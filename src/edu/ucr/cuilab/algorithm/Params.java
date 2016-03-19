package edu.ucr.cuilab.algorithm;

public class Params implements Cloneable {

	private int particles;
	private int neighbor;
	private int seqs;
	private int transOrder;
	private double majority;
	private double threshold;
	private double alpha;
	private double alphaLow;
	private double alphaHigh;
	private double zero;
	private int forTest;
	private double coverage;
	public Params() {
		super();
	}
	
	public Params(Params param) {
		super();
		this.particles = param.getParticles();
		this.neighbor = param.getNeighbor();
		this.seqs = param.getSeqs();
		this.transOrder = param.getTransOrder();
		this.majority = param.getMajority();
		this.threshold = param.getThreshold();
		this.alpha = param.getAlpha();
		this.alphaLow = param.getAlphaLow();
		this.alphaHigh = param.getAlphaHigh();
		this.zero = param.getZero();
		this.forTest = param.getForTest();
		this.coverage = param.getCoverage();
	}

	public Params(int particles, int neighbor, int seqs, int transOrder,
			double majority, double threshold, double alpha, double alphaLow,
			double alphaHigh, double zero, int forTest, double coverage) {
		super();
		this.particles = particles;
		this.neighbor = neighbor;
		this.seqs = seqs;
		this.transOrder = transOrder;
		this.majority = majority;
		this.threshold = threshold;
		this.alpha = alpha;
		this.alphaLow = alphaLow;
		this.alphaHigh = alphaHigh;
		this.zero = zero;
		this.forTest = forTest;
		this.coverage = coverage;
	}

	public double getCoverage() {
		return coverage;
	}

	public void setCoverage(double coverage) {
		this.coverage = coverage;
	}

	public int getForTest() {
		return forTest;
	}

	public void setForTest(int forTest) {
		this.forTest = forTest;
	}

	public int getParticles() {
		return particles;
	}

	public void setParticles(int particles) {
		this.particles = particles;
	}

	public int getNeighbor() {
		return neighbor;
	}

	public void setNeighbor(int neighbor) {
		this.neighbor = neighbor;
	}

	public int getSeqs() {
		return seqs;
	}

	public void setSeqs(int seqs) {
		this.seqs = seqs;
	}

	public int getTransOrder() {
		return transOrder;
	}

	public void setTransOrder(int transOrder) {
		this.transOrder = transOrder;
	}

	public double getMajority() {
		return majority;
	}

	public void setMajority(double majority) {
		this.majority = majority;
	}

	public double getThreshold() {
		return threshold;
	}

	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getAlphaLow() {
		return alphaLow;
	}

	public void setAlphaLow(double alphaLow) {
		this.alphaLow = alphaLow;
	}

	public double getAlphaHigh() {
		return alphaHigh;
	}

	public void setAlphaHigh(double alphaHigh) {
		this.alphaHigh = alphaHigh;
	}

	public double getZero() {
		return zero;
	}

	public void setZero(double zero) {
		this.zero = zero;
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
}
