package edu.ucr.cuilab.algorithm;
import java.util.Comparator;

public class DoubleReadCompareNP implements Comparator<DoubleRead> {

	@Override
	public int compare(DoubleRead arg0, DoubleRead arg1) {
		if (arg0.getIdentification() > arg1.getIdentification()) {
			return -1;
		} else if (arg0.getIdentification() < arg1.getIdentification()){
			return 1;
		} else {
			return 0;
		}
	}

}
