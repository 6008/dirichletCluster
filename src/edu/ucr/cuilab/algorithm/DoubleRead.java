package edu.ucr.cuilab.algorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class DoubleRead {

	private int id;
	private String readA;
	private String readB;
	private int transOrder;
	private List<Integer> countList;
	private List<Integer> accumulateCountList;
	private double GC;
	private int[] startWhere;
	private double[] newGroupStartp;
	private double[] newGroupTransp;

	private double identification;

	public DoubleRead(String readA, String readB, int transOrder,
			List<String> permutationList) {
		this.readA = readA;
		this.readB = readB;
		this.transOrder = transOrder;
		fillCountList(permutationList);
		calculateGC();
		// calculateGCTest();
		calculateStartWhere();
		this.identification = -1.0;
	}

	public DoubleRead(int id, String readA, String readB, int transOrder,
			List<String> permutationList) {
		this.id = id;
		this.readA = readA;
		this.readB = readB;
		this.transOrder = transOrder;
		fillCountList(permutationList);
		calculateGC();
		// calculateGCTest();
		calculateStartWhere();
		this.identification = -1.0;
	}

	public void printStartWhere() {
		System.out.println(Arrays.toString(startWhere));
	}

	public void printStartp() {
		System.out.println(Arrays.toString(newGroupStartp));
	}

	public void printTransp() {
		System.out.println(Arrays.toString(newGroupTransp));
	}

	public void printCountList() {
		System.out.println(countList.toString());
	}

	public void printAccumList() {
		System.out.println(accumulateCountList.toString());
	}

	public void outputAccum() {
		if (accumulateCountList.size() > 0) {
			for (Integer i : accumulateCountList) {
				System.out.print(i);
				System.out.print('\t');
			}
			System.out.println();
		} else {
			System.out.println("No elements");
		}
	}

	public void outputStartp() {
		if (newGroupStartp.length > 0) {
			for (double d : newGroupStartp) {
				System.out.print(d);
				System.out.print('\t');
			}
			System.out.println();
		}
	}

	public void outputCountList() {
		if (countList.size() > 0) {
			for (Integer i : countList) {
				System.out.print(i);
				System.out.print('\t');
			}
			System.out.println();
		} else {
			System.out.println("No elements");
		}
	}

	private int getStartPoint(char c, boolean origin) {
		int startWhere = 0;
		if (origin) {
			switch (c) {
			case 'C':
				startWhere = 1;
				break;
			case 'G':
				startWhere = 2;
				break;
			case 'T':
				startWhere = 3;
				break;
			default:
				startWhere = 0;
				break;
			}
		} else {
			switch (c) {
			case 'A':
				startWhere = 3;
				break;
			case 'C':
				startWhere = 2;
				break;
			case 'G':
				startWhere = 1;
				break;
			default:
				startWhere = 0;
				break;
			}
		}
		return startWhere;
	}

	private void calculateStartWhere() {
		this.startWhere = new int[4];
		this.startWhere[0] = getStartPoint(this.readA.charAt(0), true);
		this.startWhere[1] = getStartPoint(this.readB.charAt(0), true);
		this.startWhere[2] = getStartPoint(
				this.readA.charAt(this.readA.length() - 1), false);
		this.startWhere[3] = getStartPoint(
				this.readB.charAt(this.readB.length() - 1), false);
	}

	private void calculateGC() {
		int gcCount = 0;
		char[] charList = readA.toUpperCase().toCharArray();
		for (char c : charList) {
			if (('C' == c) || ('G' == c)) {
				gcCount++;
			}
		}
		charList = readB.toUpperCase().toCharArray();
		for (char c : charList) {
			if (('C' == c) || ('G' == c)) {
				gcCount++;
			}
		}
		GC = (gcCount + 0.0) / (readA.length() + readB.length());
	}

	private void calculateGCTest() {
		double gcCount = 0.0;
		char[] charList = readA.toUpperCase().toCharArray();
		for (int i = 1; i < charList.length - 1; i++) {
			if (('C' == charList[i]) || ('G' == charList[i])) {
				gcCount += 1.0;
			}
		}

		if (('C' == charList[0]) || ('G' == charList[0])) {
			gcCount += 0.5;
		}
		if (('C' == charList[charList.length - 1])
				|| ('G' == charList[charList.length - 1])) {
			gcCount += 0.5;
		}

		charList = readB.toUpperCase().toCharArray();
		for (int i = 1; i < charList.length - 1; i++) {
			if (('C' == charList[i]) || ('G' == charList[i])) {
				gcCount += 1.0;
			}
		}

		if (('C' == charList[0]) || ('G' == charList[0])) {
			gcCount += 0.5;
		}
		if (('C' == charList[charList.length - 1])
				|| ('G' == charList[charList.length - 1])) {
			gcCount += 0.5;
		}
		GC = gcCount / (readA.length() + readB.length() - 2);
	}

	public String getTransReverse(String origin) {
		char[] charList = origin.toUpperCase().toCharArray();
		for (int i = 0; i < charList.length; i++) {
			switch (charList[i]) {
			case 'A':
				charList[i] = 'T';
				break;
			case 'C':
				charList[i] = 'G';
				break;
			case 'G':
				charList[i] = 'C';
				break;
			case 'T':
				charList[i] = 'A';
				break;
			default:
				break;
			}
		}
		return (new StringBuilder(String.valueOf(charList))).reverse()
				.toString();
	}

	private void countMerForElement(String data, Map<String, Integer> countMap,
			List<String> permutationList) {
		for (int i = 0; i < data.length() - transOrder + 1; i++) {
			String key = data.substring(i, i + transOrder);
			countMap.put(key, countMap.get(key) + 1);
		}
	}

	private List<Integer> countMerForList(String[] dataList,
			List<String> permutationList) {
		List<Integer> integerList = new ArrayList<Integer>(
				permutationList.size());
		Map<String, Integer> countMap = new TreeMap<String, Integer>();
		for (String permutation : permutationList) {
			countMap.put(permutation, 0);
		}
		for (int i = 0; i < dataList.length; i++) {
			countMerForElement(dataList[i], countMap, permutationList);
		}
		for (int i = 0; i < permutationList.size(); i++) {
			integerList.add(countMap.get(permutationList.get(i)));
		}
		return integerList;
	}

	private void fillCountList(List<String> permutationList) {
		String[] dataList = { readA, readB, getTransReverse(readA),
				getTransReverse(readB) };
		countList = countMerForList(dataList, permutationList);
	}
	
	public double calcIdentification() {
		double result = 1.0;
		for (int i = 0; i < 4; i++) {
			result *= this.getNewGroupStartp()[this.getStartWhere()[i]];
		}
		double tempSum = 0.0;
		for (int i = 0; i < this.getCountList().size(); i++) {
			if (this.getNewGroupTransp()[i] < Double.MIN_NORMAL) {
				tempSum += countList.get(i) * DefaultConstants.ZERO;
			} else {
				tempSum += countList.get(i)
						* Math.log(this.getNewGroupTransp()[i]);
			}
		}
		result *= Math.exp(tempSum);
		this.identification = result;
		return result;
	}

	public double getLogIdentification() {
		double result = 0.0;
		for (int i = 0; i < 4; i++) {
			result += Math.log10(this.getNewGroupStartp()[this.getStartWhere()[i]]);
		}
		double tempSum = 0.0;
		for (int i = 0; i < this.getCountList().size(); i++) {
			if (this.getNewGroupTransp()[i] < DefaultConstants.ZERO) {
				tempSum += countList.get(i) * Math.log10(DefaultConstants.ZERO);
			} else {
				tempSum += countList.get(i)
						* Math.log10(this.getNewGroupTransp()[i]);
			}
		}
		return result + tempSum;
	}
	
	public double getIdentification() {
		if (this.identification < 0) {
			double result = 1.0;
			for (int i = 0; i < 4; i++) {
				result *= this.getNewGroupStartp()[this.getStartWhere()[i]];
			}
			double tempSum = 0.0;
			for (int i = 0; i < this.getCountList().size(); i++) {
				if (this.getNewGroupTransp()[i] < DefaultConstants.ZERO) {
					tempSum += countList.get(i) * Math.log(DefaultConstants.ZERO);
				} else {
					tempSum += countList.get(i)
							* Math.log(this.getNewGroupTransp()[i]);
				}
			}
			result *= Math.exp(tempSum);
			this.identification = result;
		}
		return this.identification;
	}

	public void setIdentification(double identification) {
		this.identification = identification;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public String getReadA() {
		return readA;
	}

	public void setReadA(String readA) {
		this.readA = readA;
	}

	public String getReadB() {
		return readB;
	}

	public void setReadB(String readB) {
		this.readB = readB;
	}

	public int getTransOrder() {
		return transOrder;
	}

	public void setTransOrder(int transOrder) {
		this.transOrder = transOrder;
	}

	public List<Integer> getCountList() {
		return countList;
	}

	public void setCountList(List<Integer> countList) {
		this.countList = countList;
	}

	public double getGC() {
		return GC;
	}

	public void setGC(double gC) {
		GC = gC;
	}

	public List<Integer> getAccumulateCountList() {
		return accumulateCountList;
	}

	public void setAccumulateCountList(List<Integer> accumulateCountList) {
		this.accumulateCountList = accumulateCountList;
	}

	public int[] getStartWhere() {
		return startWhere;
	}

	public void setStartWhere(int[] startWhere) {
		this.startWhere = startWhere;
	}

	public double[] getNewGroupStartp() {
		return newGroupStartp;
	}

	public void setNewGroupStartp(double[] newGroupStartp) {
		this.newGroupStartp = newGroupStartp;
	}

	public double[] getNewGroupTransp() {
		return newGroupTransp;
	}

	public void setNewGroupTransp(double[] newGroupTransp) {
		this.newGroupTransp = newGroupTransp;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((readA == null) ? 0 : readA.hashCode());
		result = prime * result + ((readB == null) ? 0 : readB.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		DoubleRead other = (DoubleRead) obj;
		if (readA == null) {
			if (other.readA != null)
				return false;
		} else if (!readA.equals(other.readA))
			return false;
		if (readB == null) {
			if (other.readB != null)
				return false;
		} else if (!readB.equals(other.readB))
			return false;
		return true;
	}

}
