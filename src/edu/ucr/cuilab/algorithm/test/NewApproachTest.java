package edu.ucr.cuilab.algorithm.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.ucr.cuilab.algorithm.DefaultConstants;
import edu.ucr.cuilab.algorithm.DirichletClusterSingle;
import edu.ucr.cuilab.algorithm.DoubleRead;
import edu.ucr.cuilab.algorithm.DoubleReadCompareGC;
import edu.ucr.cuilab.algorithm.DoubleReadCompareID;
import edu.ucr.cuilab.algorithm.Params;

public class NewApproachTest {
	// Get all possible pair
	private static List<String> getPermutations(char[] charList, int depth) {
		List<String> stringList = new ArrayList<String>();
		if (1 == depth) {
			for (int i = 0; i < charList.length; i++) {
				stringList.add(String.valueOf(charList[i]));
			}
		} else {
			List<String> subStringList = getPermutations(charList, depth - 1);
			for (int i = 0; i < charList.length; i++) {
				for (int j = 0; j < subStringList.size(); j++) {
					stringList.add(String.valueOf(charList[i])
							+ subStringList.get(j));
				}
			}
		}
		return stringList;
	}

	private static double parse(double defaultValue, String input,
			String descript) {
		double result = defaultValue;
		try {
			result = Double.parseDouble(input);
		} catch (NullPointerException npe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		} catch (NumberFormatException nfe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		}
		return result;
	}

	private static int parse(int defaultValue, String input, String descript) {
		int result = defaultValue;
		try {
			result = Integer.parseInt(input);
		} catch (NullPointerException npe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		} catch (NumberFormatException nfe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		}
		return result;
	}

	public static void frontEnd(String[] args) throws Exception {

		int transOrder = DefaultConstants.TRANSORDER;
		int particles = DefaultConstants.PARTICLES;
		int neighbor = DefaultConstants.NEIGHBOR;
		int seqs = 0;
		double alpha = DefaultConstants.ALPHA;
		double alphaLow = DefaultConstants.ALPHALOW;
		double alphaHigh = DefaultConstants.ALPHAHIGH;
		double threshold = DefaultConstants.THRESHOLD;
		double zero = DefaultConstants.ZERO;
		double majority = DefaultConstants.MAJORITY;
		int clusterNumForTest = DefaultConstants.TESTNUM;
		double coverage = DefaultConstants.COVERAGE;

		List<Double> alphaList = new ArrayList<Double>();
//		 String input = "/home/xinping/Desktop/6008/5_2345Borr_5W.txt";
//		 String input = "/home/xinping/Desktop/6008/Data/5_3910AcaryEhrl.txt";
//		 String input = "/home/xinping/Desktop/6008/Data/5_unrelated_test.txt";
//		 String input = "/home/xinping/Desktop/6008/Data/AcaryEhrlNew.txt";
		 String input = "/home/xinping/Desktop/6008/Data/AcaryEhrlNew_5W.txt";
//		String input = "/home/xinping/Desktop/6008/5_2345Borr.txt";
		String output = null;

		Date date = new Date();

		Params params = new Params(particles, neighbor, seqs, transOrder,
				majority, threshold, alpha, alphaLow, alphaHigh, zero,
				clusterNumForTest, coverage);

		if (0 == args.length) {
			printHelp();
			alphaList.add(0.00000001);
			alphaList.add(0.000001);
			// alphaList.add(0.000075);
			// alphaList.add(0.1);
			// return;
		}

		int pos = 0;
		while (pos < args.length - 1) {
			if (args[pos].charAt(0) == '-') {
				switch (args[pos].charAt(1)) {
				case 'a':
				case 'A':
					params.setAlpha(parse(DefaultConstants.ALPHA,
							args[pos + 1], "Parameter alpha"));
					pos += 2;
				case 'c':
				case 'C':
					params.setCoverage(parse(DefaultConstants.COVERAGE,
							args[pos + 1], "Parameter coverage"));
					pos += 2;
					break;
				case 'i':
				case 'I':
					input = args[pos + 1];
					pos += 2;
					break;
				case 'l':
				case 'L':
					params.setAlphaLow(parse(DefaultConstants.ALPHALOW,
							args[pos + 1], "Parameter alpha_low"));
					pos += 2;
					break;
				case 'm':
				case 'M':
					params.setMajority(parse(DefaultConstants.MAJORITY,
							args[pos + 1], "Parameter majority"));
					pos += 2;
					break;
				case 'n':
				case 'N':
					params.setNeighbor(parse(DefaultConstants.NEIGHBOR,
							args[pos + 1], "Parameter neighbor_number"));
					pos += 2;
					break;
				case 'o':
				case 'O':
					params.setTransOrder(parse(DefaultConstants.TRANSORDER,
							args[pos + 1], "Parameter transorder"));
					pos += 2;
					break;
				case 'p':
				case 'P':
					params.setParticles(parse(DefaultConstants.PARTICLES,
							args[pos + 1], "Parameter particle_amount"));
					pos += 2;
					break;
				case 't':
				case 'T':
					params.setThreshold(parse(DefaultConstants.THRESHOLD,
							args[pos + 1], "Parameter threshold"));
					pos += 2;
					break;
				case 'u':
				case 'U':
					output = args[pos + 1];
					pos += 2;
					break;
				case 'z':
				case 'Z':
					params.setZero(parse(DefaultConstants.ZERO, args[pos + 1],
							"Parameter zero"));
					pos += 2;
					break;
				case '5':
					params.setForTest(parse(DefaultConstants.TESTNUM,
							args[pos + 1], "Parameter forTest"));
					pos += 2;
				default:
					pos++;
					break;
				}
			} else {
				pos++;
			}
		}

		if (null == output) {
			output = input + "." + date.toString() + ".output";
		}

		// mainJob(input, output, params, alphaList);
		mainJobTest(input, output, params, alphaList);
	}

	public static void printHelp() {
		System.out.println("Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " -i abundance_species_equal.txt");
		System.out.println("The file contains pair-end sequence\n");
		System.out.println("More Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " [OPTION] -i abundance_species_equal.txt\n");
		System.out.println("Options:");
		// System.out
		// .println("\t-a:\tThe value followed will be the alpha parameter");
		// System.out.println("\t   \tDefault value " + DefaultConstants.ALPHA);
		// System.out.println("\t\tExample 1: -a3 0.0000001 0.00001 0.001");
		// System.out.println("\t\tExample 2: -a 0.0000001");
		// System.out
		// .println("\t-h:\tThe value followed will be the alpha_high parameter");
		// System.out
		// .println("\t   \tDefault value " + DefaultConstants.ALPHAHIGH);
		System.out
				.println("\t-a:\tThe value followed will be the alpha parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ALPHA);
		System.out
				.println("\t-c:\tThe value followed will be the coverage parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.COVERAGE);
		System.out
				.println("\t-l:\tThe value followed will be the alpha_low parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ALPHALOW);
		System.out
				.println("\t-m:\tThe value followed will be the majority parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.MAJORITY);
		System.out
				.println("\t-n:\tThe value followed will be the neighbor_num parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.NEIGHBOR);
		System.out
				.println("\t-o:\tThe value followed will be the transform_order parameter");
		System.out.println("\t   \tDefault value "
				+ DefaultConstants.TRANSORDER);
		System.out
				.println("\t-p:\tThe value followed will be particle_number parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.PARTICLES);
		System.out
				.println("\t-t:\tThe value followed will be the threshold parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.THRESHOLD);
		System.out
				.println("\t-u:\tThe string followed will be the output file");
		System.out
				.println("\t   \tDefault value is the input file name with additional time");
		System.out.println("\t-z:\tThe value followed will be zero parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ZERO);
	}

	private static List<String> readFile(String inputFile) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(
				inputFile)));
		List<String> readList = new ArrayList<String>();
		String line = null;
		while (null != (line = br.readLine())) {
			readList.add(line);
		}
		br.close();
		return readList;
	}

//	private static double[] startPCounts(List<Integer> counts) {
//		double[] probs = new double[counts.size() / 4];
//		double sum = 0.0;
//		for (int i = 0; i < probs.length; i++) {
//			probs[i] = counts.get(i * 4) + counts.get(i * 4 + 1)
//					+ counts.get(i * 4 + 2) + counts.get(i * 4 + 3);
//			sum += probs[i];
//		}
//		for (int i = 0; i < probs.length; i++) {
//			probs[i] /= sum;
//		}
//		return probs;
//	}
//
//	private static double[] transPCounts(List<Integer> counts) {
//		double[] probs = new double[counts.size()];
//		for (int i = 0; i < probs.length; i += 4) {
//			double sum = counts.get(i) + counts.get(i + 1) + counts.get(i + 2)
//					+ counts.get(i + 3);
//			probs[i] = counts.get(i) / sum;
//			probs[i + 1] = counts.get(i + 1) / sum;
//			probs[i + 2] = counts.get(i + 2) / sum;
//			probs[i + 3] = counts.get(i + 3) / sum;
//		}
//		return probs;
//	}
	
	private static double[] startPCounts(List<Integer> counts) {
		double[] probs = new double[DefaultConstants.CHARSETSIZE];
		int length = counts.size() / DefaultConstants.CHARSETSIZE;
		double sum = 0.0;
		for (int i = 0; i < probs.length; i++) {
			probs[i] = 0.0;
			for (int j = 0; j < length; j++) {
				probs[i] += counts.get(i * length + j);
			}
			sum += probs[i];
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		return probs;
	}

	private static double[] transPCounts(List<Integer> counts) {
		double[] probs = new double[counts.size()];
		int length = counts.size() / DefaultConstants.CHARSETSIZE;
		for (int i = 0; i < probs.length; i += length) {
			double sum = 0.0;
			for (int j = 0; j < length; j++) {
				sum += counts.get(i + j);
			}
			for (int j = 0; j < length; j++) {
				probs[i + j] = counts.get(i + j) / sum;
			}
		}
		return probs;
	}

	private static void newGroup(List<DoubleRead> doubleReadList, int neighbor) {
		Collections.sort(doubleReadList, new DoubleReadCompareGC());
		for (int i = 0; i < doubleReadList.size(); i++) {
			int nflag = 1;
			int pflag = 1;
			List<Integer> accumCountList = new ArrayList<Integer>(
					doubleReadList.get(i).getCountList());

			while ((nflag <= neighbor / 2) && (i - nflag >= 0)) {
				List<Integer> temp = doubleReadList.get(i - nflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				nflag++;
			}

			while ((pflag <= neighbor / 2)
					&& (i + pflag < doubleReadList.size())) {
				List<Integer> temp = doubleReadList.get(i + pflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				pflag++;
			}

			doubleReadList.get(i).setAccumulateCountList(accumCountList);
			doubleReadList.get(i).setNewGroupStartp(
					startPCounts(accumCountList));
			doubleReadList.get(i).setNewGroupTransp(
					transPCounts(accumCountList));
			doubleReadList.get(i).getIdentification();

		}
//		Collections.sort(doubleReadList, new DoubleReadCompareID());
//		Collections.sort(doubleReadList, new DoubleReadCompareNP());
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File(
					"/home/xinping/Desktop/6008/metric_output"
							+ (new Date()).toString())));
			for (DoubleRead dr : doubleReadList) {
				double[] transp = transPCounts(dr.getCountList());
				double[] startp = startPCounts(dr.getCountList());
				double result = 0.0;
				for (int i = 0; i < 4; i++) {
					result += Math.log10(startp[dr.getStartWhere()[i]]);
				}
				double tempSum = 0.0;
				for (int i = 0; i < dr.getCountList().size(); i++) {
					if (transp[i] < DefaultConstants.ZERO) {
						tempSum += dr.getCountList().get(i)
								* Math.log10(DefaultConstants.ZERO);
					} else {
						tempSum += dr.getCountList().get(i)
								* Math.log10(transp[i]);
					}
				}
				result += tempSum;
				pw.print(dr.getId());
				pw.print('\t');
				pw.print(dr.getGC());
				pw.print('\t');
				pw.print(result);
				pw.print('\t');
				pw.print(dr.getLogIdentification());
				pw.print('\t');
				pw.print(dr.getCountList().get(0));
				pw.print('\t');
				pw.print(dr.getCountList().get(3));
				pw.print('\t');
				pw.print(dr.getCountList().get(12));
				pw.print('\t');
				pw.print(dr.getCountList().get(5));
				pw.print('\t');
				pw.print(dr.getCountList().get(6));
				pw.print('\t');
				pw.println(dr.getCountList().get(9));
//				pw.println(dr.getCountList().toString());
			}
			pw.close();
		} catch (Exception e) {

		}
	}

	private static List<Map<Integer, Integer>> updateGroupMapList(
			List<Map<Integer, Integer>> groupMapList, int[][] z, int j) {
		for (int i = 0; i < groupMapList.size(); i++) {
			Map<Integer, Integer> temp = groupMapList.get(i);
			if (temp.containsKey(z[j][i])) {
				temp.put(z[j][i], temp.get(z[j][i]) + 1);
			} else {
				temp.put(z[j][i], 1);
			}
		}
		return groupMapList;
	}

	private static List<Map<Integer, List<Integer>>> updateCountLists(
			List<DoubleRead> doubleReadList, int[][] z, int elem,
			List<Map<Integer, List<Integer>>> originCountLists) {

		int j = doubleReadList.get(elem).getId();
		for (int i = 0; i < originCountLists.size(); i++) {
			Map<Integer, List<Integer>> listMap = originCountLists.get(i);
			if (listMap.containsKey(z[j][i])) {
				List<Integer> temp = new ArrayList<Integer>(
						listMap.get(z[j][i]));
				for (int k = 0; k < temp.size(); k++) {
					temp.set(k, temp.get(k)
							+ doubleReadList.get(elem).getCountList().get(k));
				}
				listMap.put(z[j][i], temp);
			} else {
				listMap.put(z[j][i],
						new ArrayList<Integer>(doubleReadList.get(elem)
								.getCountList()));
			}
			originCountLists.set(i, listMap);
		}

		return originCountLists;
	}

	private static int[] dirichletClusterSingle(
			List<DoubleRead> doubleReadList, List<Set<Integer>> overlapList,
			Params params, int zModeLower) {
		System.out.println((new Date()).toString()
				+ "\tDirichletClusterSingle \nBegin");
		int[][] z = new int[params.getSeqs()][params.getParticles()];
		double[] w = new double[params.getParticles()];
		for (int i = 0; i < params.getSeqs(); i++) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[i][j] = 0;
			}
		}
		for (int i = 0; i < params.getParticles(); i++) {
			w[i] = 1.0 / (0.0 + params.getParticles());
		}
		List<Map<Integer, List<Integer>>> tempAccumCountList = new ArrayList<Map<Integer, List<Integer>>>(
				params.getParticles());
		List<Map<Integer, Integer>> mapCountList = new ArrayList<Map<Integer, Integer>>(
				params.getParticles());
		for (int i = 0; i < params.getParticles(); i++) {
			tempAccumCountList.add(new TreeMap<Integer, List<Integer>>());
			mapCountList.add(new TreeMap<Integer, Integer>());
		}
		newGroup(doubleReadList, params.getNeighbor());

		System.out.println((new Date()).toString()
				+ "\tNew Group calculation complete");

		for (Integer elem : overlapList.get(doubleReadList.get(0).getId())) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[elem][j] = 1;
			}
		}
		int group = 1;
		tempAccumCountList = updateCountLists(doubleReadList, z, 0,
				tempAccumCountList);
		mapCountList = updateGroupMapList(mapCountList, z, doubleReadList
				.get(0).getId());
		Map<Integer, Set<Integer>> overallResult = new HashMap<Integer, Set<Integer>>();

		long[] timeArr = new long[3];
		// try {
		// PrintWriter pw = new PrintWriter(new FileWriter(new File(
		// "/home/xinping/Desktop/6008/metric_output"
		// + (new Date()).toString())));
		// for (DoubleRead dr : doubleReadList) {
		// pw.print(dr.calcIdentification());
		// pw.print("\t");
		// pw.println(dr.getId());
		// }
		// pw.close();
		// } catch (Exception e) {
		//
		// }

		for (int i = 1; i < doubleReadList.size(); i++) {
			// System.out.println(i);
			int elem = doubleReadList.get(i).getId();
			Date d0 = new Date();
//			System.out.println(i);
			z[elem] = DirichletClusterSingle.clusterOneSeq(mapCountList,
					tempAccumCountList, doubleReadList, i, z, w,
					params.getParticles(), params.getAlpha(), i, group);

			w = DirichletClusterSingle.updateWeights(tempAccumCountList,
					doubleReadList, elem, z, w, params.getParticles(), group);
			// System.out.println(Arrays.toString(w));
			double eff = 0.0;
			for (double d : w) {
				eff += d * d;
			}
			if (eff > 1 / (params.getThreshold() * params.getParticles())) {
				int[] resampleIndex = DirichletClusterSingle.sampleInt(w,
						params.getParticles(), 0);
				int[] temp = new int[params.getParticles()];
				for (int j = 0; j < params.getParticles(); j++) {
					w[j] = 1.0 / (0.0 + params.getParticles());
					temp[j] = z[elem][resampleIndex[j]];
				}
				z[elem] = temp;
			}

			for (int j = 0; j < z[elem].length; j++) {
				if (z[elem][j] > group) {
					group = z[elem][j];
					break;
				}
			}

			tempAccumCountList = updateCountLists(doubleReadList, z, i,
					tempAccumCountList);
			mapCountList = updateGroupMapList(mapCountList, z, elem);

		}

		int[] result = DirichletClusterSingle
				.compressResult(DirichletClusterSingle.getMode(z));
		for (int i = 0; i < result.length; i++) {
			if (overallResult.containsKey(result[i])) {
				overallResult.get(result[i]).add(i);
			} else {
				Set<Integer> temp = new TreeSet<Integer>();
				temp.add(i);
				overallResult.put(result[i], temp);
			}
		}
		System.out.println("Overall");
		for (Integer key : overallResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			for (Integer temp : overallResult.get(key)) {
				int i = doubleReadList.get(temp).getId() / params.getForTest()
						+ 1;
				if (tempMap.containsKey(i)) {
					tempMap.put(i, tempMap.get(i) + 1);
				} else {
					tempMap.put(i, 1);
				}
			}
			for (Integer tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey)
						/ (overallResult.get(key).size() + 0.0));
			}
			System.out.println();
			System.out.println();
		}
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] + zModeLower;
		}
		return result;
	}

	private static void dirichletClusterTest(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, int[] zMode, Params params) {
		params.setAlpha(0.000001);

		// dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		// params.setAlphaHigh(1);
		// dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		// params.setAlphaHigh(0.1);
		// dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		// params.setAlphaHigh(0.01);
		// dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		// params.setAlphaHigh(0.001);
		dirichletClusterSingle(doubleReadList, overlapList, params, 0);
	}

	public static void mainJobTest(String inputFile, String outputFile,
			Params params, List<Double> alphaList) throws Exception {
		System.out.println((new Date()).toString()
				+ "\tProgram start... \nReading dataset...");

		List<String> readList = readFile(inputFile);
		params.setSeqs(readList.size() / 2);

		char[] charList = { 'A', 'C', 'G', 'T' };
		List<String> permutationList = getPermutations(charList,
				params.getTransOrder() + 1);

		System.out.println((new Date()).toString()
				+ "\tBegin to calculate overlapping infos");

		// List<Set<Integer>> overlapList =
		// EncodeBinTask.encodeMainJob(readList,
		// 32);
		List<Set<Integer>> overlapList = new ArrayList<Set<Integer>>();
		for (int i = 0; i < params.getSeqs(); i++) {
			Set<Integer> temp = new TreeSet<Integer>();
			temp.add(i);
			overlapList.add(temp);
		}

		System.out.println((new Date()).toString() + "\tComplete");

		List<DoubleRead> doubleReadList = new ArrayList<DoubleRead>();
		for (int i = 0; i < params.getSeqs(); i++) {
			DoubleRead dr = new DoubleRead(i, readList.get(i * 2),
					readList.get(i * 2 + 1), params.getTransOrder() + 1,
					permutationList);
			doubleReadList.add(dr);
			// dr.printCountList();
		}

		System.out
				.println((new Date()).toString() + "\tMain processor started");

		int[] zMode = new int[params.getSeqs()];
		for (int i = 0; i < params.getSeqs(); i++) {
			zMode[i] = 1;
		}
		params.setForTest(params.getSeqs() / params.getForTest());
		dirichletClusterTest(doubleReadList, overlapList, zMode, params);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			frontEnd(args);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
