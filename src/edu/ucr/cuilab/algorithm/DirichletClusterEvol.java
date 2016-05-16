package edu.ucr.cuilab.algorithm;

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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class DirichletClusterEvol {

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
		// String input = "/home/xinping/Desktop/6008/5_2345Borr.txt.sample";
		String input = "/home/xinping/Desktop/6008/5_test.txt";
//		 String input = "/home/xinping/Desktop/6008/5_2345Borr.txt";
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

	private static List<Set<Integer>> updateOverlap(
			List<DoubleRead> doubleReadList, List<Set<Integer>> overlapList) {
		Map<Integer, Integer> reverseMap = new HashMap<Integer, Integer>();
		for (int i = 0; i < doubleReadList.size(); i++) {
			reverseMap.put(doubleReadList.get(i).getId(), i);
		}
		List<Set<Integer>> resultList = new ArrayList<Set<Integer>>();
		for (int i = 0; i < overlapList.size(); i++) {
			Set<Integer> temp = new HashSet<Integer>();
			for (Integer id : overlapList.get(i)) {
				if (reverseMap.containsKey(id)) {
					temp.add(reverseMap.get(id));
				}
			}
			if (temp.size() > 0) {
				resultList.add(temp);
			}
		}
		Collections.sort(resultList, new ListCompareLength());
		return resultList;
	}

	public static int[] dirichletClusterSingleV1(
			List<DoubleRead> doubleReadList, List<Set<Integer>> overlapList,
			Params params, int zModeLower) {
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

		Map<Integer, String> idTagMap = new HashMap<Integer, String>();
		for (int i = 0; i < params.getSeqs(); i++) {
			idTagMap.put(i, "");
		}

		List<Map<Integer, List<Integer>>> tempAccumCountList = new ArrayList<Map<Integer, List<Integer>>>(
				params.getParticles());
		List<Map<Integer, Integer>> mapCountList = new ArrayList<Map<Integer, Integer>>(
				params.getParticles());

		for (int i = 0; i < params.getParticles(); i++) {
			tempAccumCountList.add(new TreeMap<Integer, List<Integer>>());
			mapCountList.add(new TreeMap<Integer, Integer>());
		}
		for (Integer elem : overlapList.get(0)) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[elem][j] = 1;
			}
			idTagMap.put(elem, idTagMap.get(elem) + "Init");
		}
		tempAccumCountList = DirichletClusterSingle.updateCountLists(
				doubleReadList, z, overlapList.get(0), tempAccumCountList);
		mapCountList = DirichletClusterSingle.updateGroupMapList(mapCountList,
				z, overlapList.get(0));
		DirichletClusterSingle.newGroup(doubleReadList, params.getNeighbor());
		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();
		int accumSeqCount = overlapList.get(0).size();
		int overlapInTotal = 0;
		int selectedOverlap = 0;
		Map<Integer, Set<Integer>> overlapResult = new HashMap<Integer, Set<Integer>>();
		Map<Integer, Set<Integer>> overallResult = new HashMap<Integer, Set<Integer>>();

		List<Integer> checkAgainList = new ArrayList<Integer>();

		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() > 1) {
				overlapInTotal++;
				int[][] tempZLower = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaLow(), accumSeqCount);
				int[][] tempZUpper = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaHigh(), accumSeqCount);

				boolean isMarjorityLower = DirichletClusterSingle
						.checkMajorityVote(tempZLower, params.getMajority());
				boolean isMarjorityUpper = DirichletClusterSingle
						.checkMajorityVote(tempZUpper, params.getMajority());

				int[] voteLower = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZLower));
				int[] voteUpper = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZUpper));

				int majorityVoteLower = voteLower[0];
				int majorityVoteUpper = voteUpper[0];
				if (isMarjorityLower && isMarjorityUpper
						&& (majorityVoteLower == majorityVoteUpper)) {

					for (Integer j : overlapList.get(i)) {
						for (int k = 0; k < params.getParticles(); k++) {
							z[j][k] = majorityVoteLower;
						}
						idTagMap.put(j, idTagMap.get(j) + "Majority");
					}

					tempAccumCountList = DirichletClusterSingle
							.updateCountLists(doubleReadList, z,
									overlapList.get(i), tempAccumCountList);
					mapCountList = DirichletClusterSingle.updateGroupMapList(
							mapCountList, z, overlapList.get(i));

					accumSeqCount += overlapList.get(i).size();
					if (overlapResult.containsKey(majorityVoteLower)) {
						overlapResult.get(majorityVoteLower).add(i);
					} else {
						Set<Integer> temp = new HashSet<Integer>();
						temp.add(i);
						overlapResult.put(majorityVoteLower, temp);
					}

					selectedOverlap++;
				} else {
					checkAgainList.add(i);
				}
			}
		}

		System.out.println();
		System.out.print(params.getAlphaHigh());
		System.out.print("\t");
		System.out.print(doubleReadList.size());
		System.out.print("\t");
		System.out.print(overlapInTotal);
		System.out.print("\t");
		System.out.print(selectedOverlap);
		System.out.print("\t");
		System.out.print((selectedOverlap + 0.0) / overlapInTotal);
		System.out.print("\t");
		System.out.println(String
				.valueOf(((selectedOverlap + 0.0) / overlapInTotal) > params
						.getCoverage()));
		System.out.println("Overlap");
		for (Integer key : overlapResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			int size = 0;
			for (Integer groupIndex : overlapResult.get(key)) {
				size += overlapList.get(groupIndex).size();
				for (Integer temp : overlapList.get(groupIndex)) {
					int i = doubleReadList.get(temp).getId()
							/ params.getForTest() + 1;
					if (tempMap.containsKey(i)) {
						tempMap.put(i, tempMap.get(i) + 1);
					} else {
						tempMap.put(i, 1);
					}
				}
			}
			for (Integer tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey) / (size + 0.0));
			}
			System.out.println();
		}

		int goodKey = 0;
		List<Integer> removeKeyList = new ArrayList<Integer>();
		for (Integer key : overlapResult.keySet()) {
			int sumGroup = overlapResult.get(key).size();
			int sumSeqs = 0;
			for (Integer id : overlapResult.get(key)) {
				sumSeqs += overlapList.get(id).size();
			}
			if (sumSeqs >= 50) {
				goodKey++;
			} else {
				for (Integer id : overlapResult.get(key)) {
					for (Integer seq : overlapList.get(id)) {
						for (int i = 0; i < params.getParticles(); i++) {
							z[seq][i] = 0;
						}
					}
					rmList.add(overlapList.get(id));
					selectedOverlap--;
					accumSeqCount -= overlapList.get(id).size();
				}
				removeKeyList.add(key);
			}
		}
		for (Integer key : removeKeyList) {
			overlapResult.remove(key);
		}

		mapCountList = DirichletClusterSingle.rebuildGroupMapList(z,
				params.getParticles());
		tempAccumCountList = DirichletClusterSingle.rebuildCountLists(
				doubleReadList, z, params.getParticles());

		System.out.println(Arrays.toString(DirichletClusterSingle.getMode(z)));
		
		for (Integer i : checkAgainList) {
			int[][] tempZLower = DirichletClusterSingle.clusterOverlapSeqs(
					mapCountList, tempAccumCountList, doubleReadList,
					overlapList.get(i), z, w, params.getParticles(),
					params.getAlphaLow() * 100, accumSeqCount);
			int[][] tempZUpper = DirichletClusterSingle.clusterOverlapSeqs(
					mapCountList, tempAccumCountList, doubleReadList,
					overlapList.get(i), z, w, params.getParticles(),
					params.getAlphaHigh() / 10, accumSeqCount);

			boolean isMarjorityLower = DirichletClusterSingle
					.checkMajorityVote(tempZLower, params.getMajority());
			boolean isMarjorityUpper = DirichletClusterSingle
					.checkMajorityVote(tempZUpper, params.getMajority());

			int[] voteLower = DirichletClusterSingle
					.getMode(DirichletClusterSingle.getMode(tempZLower));
			int[] voteUpper = DirichletClusterSingle
					.getMode(DirichletClusterSingle.getMode(tempZUpper));

			int majorityVoteLower = voteLower[0];
			int majorityVoteUpper = voteUpper[0];
			if (isMarjorityLower && isMarjorityUpper
					&& (majorityVoteLower == majorityVoteUpper)) {

				for (Integer j : overlapList.get(i)) {
					for (int k = 0; k < params.getParticles(); k++) {
						z[j][k] = majorityVoteLower;
					}
					idTagMap.put(j, idTagMap.get(j) + "Majority");
				}

				tempAccumCountList = DirichletClusterSingle.updateCountLists(
						doubleReadList, z, overlapList.get(i),
						tempAccumCountList);
				mapCountList = DirichletClusterSingle.updateGroupMapList(
						mapCountList, z, overlapList.get(i));

				accumSeqCount += overlapList.get(i).size();
				if (overlapResult.containsKey(majorityVoteLower)) {
					overlapResult.get(majorityVoteLower).add(i);
				} else {
					Set<Integer> temp = new HashSet<Integer>();
					temp.add(i);
					overlapResult.put(majorityVoteLower, temp);
				}

				selectedOverlap++;
			} else {
				rmList.add(overlapList.get(i));
			}
		}

		System.out.println();
		System.out.print(params.getAlphaHigh());
		System.out.print("\t");
		System.out.print(doubleReadList.size());
		System.out.print("\t");
		System.out.print(overlapInTotal);
		System.out.print("\t");
		System.out.print(selectedOverlap);
		System.out.print("\t");
		System.out.print((selectedOverlap + 0.0) / overlapInTotal);
		System.out.print("\t");
		System.out.println(String
				.valueOf(((selectedOverlap + 0.0) / overlapInTotal) > params
						.getCoverage()));
		System.out.println("Overlap");
		for (Integer key : overlapResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			int size = 0;
			for (Integer groupIndex : overlapResult.get(key)) {
				size += overlapList.get(groupIndex).size();
				for (Integer temp : overlapList.get(groupIndex)) {
					int i = doubleReadList.get(temp).getId()
							/ params.getForTest() + 1;
					if (tempMap.containsKey(i)) {
						tempMap.put(i, tempMap.get(i) + 1);
					} else {
						tempMap.put(i, 1);
					}
				}
			}
			for (Integer tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey) / (size + 0.0));
			}
			System.out.println();
		}

		goodKey = 0;
		removeKeyList = new ArrayList<Integer>();
		for (Integer key : overlapResult.keySet()) {
			int sumGroup = overlapResult.get(key).size();
			int sumSeqs = 0;
			for (Integer id : overlapResult.get(key)) {
				sumSeqs += overlapList.get(id).size();
			}
			if ((sumSeqs >= 50)) {
				goodKey++;
			} else {
				for (Integer id : overlapResult.get(key)) {
					for (Integer seq : overlapList.get(id)) {
						for (int i = 0; i < params.getParticles(); i++) {
							z[seq][i] = 0;
						}
					}
					rmList.add(overlapList.get(id));
					selectedOverlap--;
					accumSeqCount -= overlapList.get(id).size();
				}
				removeKeyList.add(key);
				System.out.println(key);
			}
		}
		for (Integer key : removeKeyList) {
			overlapResult.remove(key);
		}

		mapCountList = DirichletClusterSingle.rebuildGroupMapList(z,
				params.getParticles());
		tempAccumCountList = DirichletClusterSingle.rebuildCountLists(
				doubleReadList, z, params.getParticles());
//		System.out.println(Arrays.toString(DirichletClusterSingle.getMode(z)));

		System.out.print("Removed:");
		Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
		int size = 0;
		for (Set<Integer> overlap : rmList) {
			size += overlap.size();
			for (Integer elem : overlap) {
				int i = doubleReadList.get(elem).getId() / params.getForTest()
						+ 1;
				if (tempMap.containsKey(i)) {
					tempMap.put(i, tempMap.get(i) + 1);
				} else {
					tempMap.put(i, 1);
				}
			}
		}
		for (Integer tempKey : tempMap.keySet()) {
			System.out.print("\t");
			System.out.print(tempKey);
			System.out.print("\t");
			System.out.print(tempMap.get(tempKey));
			System.out.print("\t");
			System.out.print(tempMap.get(tempKey) / (size + 0.0));
		}
		System.out.println();

		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() == 1) {
				for (Integer elem : overlapList.get(i)) {
					idTagMap.put(elem, idTagMap.get(elem) + "Single");
					// update
					int[] temp = DirichletClusterSingle.clusterOverlapSeqs(
							mapCountList, tempAccumCountList, doubleReadList,
							overlapList.get(i), z, w, params.getParticles(),
							0.1, accumSeqCount)[0];
					int majorityVoteLower = DirichletClusterSingle
							.getMode(temp)[0];
					if (overlapResult.containsKey(majorityVoteLower)) {
						for (int k = 0; k < params.getParticles(); k++) {
							z[elem][k] = majorityVoteLower;
						}

						// w = DirichletClusterSingle.updateWeights(
						// tempAccumCountList, doubleReadList, elem, z, w,
						// params.getParticles());
						//
						// double eff = 0.0;
						// for (double d : w) {
						// eff += d * d;
						// }
						// if (eff > 1 / (params.getThreshold() * params
						// .getParticles())) {
						// int[] resampleIndex =
						// DirichletClusterSingle.sampleInt(
						// w, params.getParticles(), 0);
						// int[] temp = new int[params.getParticles()];
						// for (int j = 0; j < params.getParticles(); j++) {
						// w[j] = 1.0 / (0.0 + params.getParticles());
						// temp[j] = z[elem][resampleIndex[j]];
						// }
						// z[elem] = temp;
						// }

						tempAccumCountList = DirichletClusterSingle
								.updateCountLists(doubleReadList, z,
										overlapList.get(i), tempAccumCountList);
						mapCountList = DirichletClusterSingle
								.updateGroupMapList(mapCountList, z,
										overlapList.get(i));
					} else {
//						System.out.println(overlapList.get(i).toString());
//						System.out.println(Arrays.toString(temp));
					}
				}
			}
		}

		System.out.println(accumSeqCount);
		System.out.println("Overall");
//		int[] result = DirichletClusterSingle
//				.compressResult(DirichletClusterSingle.getMode(z));
		int[] result = DirichletClusterSingle.getMode(z);
		for (int i = 0; i < result.length; i++) {
			if (overallResult.containsKey(result[i])) {
				overallResult.get(result[i]).add(i);
			} else {
				Set<Integer> temp = new TreeSet<Integer>();
				temp.add(i);
				overallResult.put(result[i], temp);
			}
		}
		for (Integer key : overallResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			tempMap = new HashMap<Integer, Integer>();
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
//		System.out.println(Arrays.toString(DirichletClusterSingle.getMode(z)));

		for (Set<Integer> overlap : rmList) {
			int[][] tempZLower = DirichletClusterSingle.clusterOverlapSeqs(
					mapCountList, tempAccumCountList, doubleReadList, overlap,
					z, w, params.getParticles(), 0.001,
					accumSeqCount);
			boolean isMarjorityLower = DirichletClusterSingle
					.checkMajorityVote(tempZLower, params.getMajority());

			int[] voteLower = DirichletClusterSingle
					.getMode(DirichletClusterSingle.getMode(tempZLower));

			int majorityVoteLower = voteLower[0];
			if (isMarjorityLower
					&& overlapResult.containsKey(majorityVoteLower)) {

				for (Integer j : overlap) {
					for (int k = 0; k < params.getParticles(); k++) {
						z[j][k] = majorityVoteLower;
					}
					idTagMap.put(j, idTagMap.get(j) + "Majority");
				}

				tempAccumCountList = DirichletClusterSingle.updateCountLists(
						doubleReadList, z, overlap, tempAccumCountList);
				mapCountList = DirichletClusterSingle.updateGroupMapList(
						mapCountList, z, overlap);

				accumSeqCount += overlap.size();

				selectedOverlap++;
			} else {
//				System.out.println(overlap.toString());
//				for (int[] temp : tempZLower) {
//					System.out.println(Arrays.toString(temp));
//				}
			}
		}
		


		System.out.println(accumSeqCount);
		System.out.println("Overall");
//		result = DirichletClusterSingle
//				.compressResult(DirichletClusterSingle.getMode(z));
		result = DirichletClusterSingle.getMode(z);
		overallResult.clear();
		for (int i = 0; i < result.length; i++) {
			if (overallResult.containsKey(result[i])) {
				overallResult.get(result[i]).add(i);
			} else {
				Set<Integer> temp = new TreeSet<Integer>();
				temp.add(i);
				overallResult.put(result[i], temp);
			}
		}
		for (Integer key : overallResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			tempMap = new HashMap<Integer, Integer>();
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
		//System.out.println(Arrays.toString(DirichletClusterSingle.getMode(z)));
		return null;
	}

	public static int[] dirichletClusterSingle(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, Params params, int zModeLower) {
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

		Map<Integer, String> idTagMap = new HashMap<Integer, String>();
		for (int i = 0; i < params.getSeqs(); i++) {
			idTagMap.put(i, "");
		}

		List<Map<Integer, List<Integer>>> tempAccumCountList = new ArrayList<Map<Integer, List<Integer>>>(
				params.getParticles());
		List<Map<Integer, Integer>> mapCountList = new ArrayList<Map<Integer, Integer>>(
				params.getParticles());

		for (int i = 0; i < params.getParticles(); i++) {
			tempAccumCountList.add(new TreeMap<Integer, List<Integer>>());
			mapCountList.add(new TreeMap<Integer, Integer>());
		}
		for (Integer elem : overlapList.get(0)) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[elem][j] = 1;
			}
			idTagMap.put(elem, idTagMap.get(elem) + "Init");
		}
		tempAccumCountList = DirichletClusterSingle.updateCountLists(
				doubleReadList, z, overlapList.get(0), tempAccumCountList);
		mapCountList = DirichletClusterSingle.updateGroupMapList(mapCountList,
				z, overlapList.get(0));
		DirichletClusterSingle.newGroup(doubleReadList, params.getNeighbor());
		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();
		int accumSeqCount = overlapList.get(0).size();
		int overlapInTotal = 0;
		int selectedOverlap = 0;
		Map<Integer, Set<Integer>> overlapResult = new HashMap<Integer, Set<Integer>>();
		Map<Integer, Set<Integer>> overallResult = new HashMap<Integer, Set<Integer>>();

		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() > 1) {
				overlapInTotal++;
				int[][] tempZLower = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaLow(), accumSeqCount);
				int[][] tempZUpper = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaHigh(), accumSeqCount);

				boolean isMarjorityLower = DirichletClusterSingle
						.checkMajorityVote(tempZLower, params.getMajority());
				boolean isMarjorityUpper = DirichletClusterSingle
						.checkMajorityVote(tempZUpper, params.getMajority());

				int[] voteLower = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZLower));
				int[] voteUpper = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZUpper));

				int majorityVoteLower = voteLower[0];
				int majorityVoteUpper = voteUpper[0];
				if (isMarjorityLower && isMarjorityUpper
						&& (majorityVoteLower == majorityVoteUpper)) {

					for (Integer j : overlapList.get(i)) {
						for (int k = 0; k < params.getParticles(); k++) {
							z[j][k] = majorityVoteLower;
						}
						idTagMap.put(j, idTagMap.get(j) + "Majority");
					}

					tempAccumCountList = DirichletClusterSingle
							.updateCountLists(doubleReadList, z,
									overlapList.get(i), tempAccumCountList);
					mapCountList = DirichletClusterSingle.updateGroupMapList(
							mapCountList, z, overlapList.get(i));

					accumSeqCount += overlapList.get(i).size();
					if (overlapResult.containsKey(majorityVoteLower)) {
						overlapResult.get(majorityVoteLower).add(i);
					} else {
						Set<Integer> temp = new HashSet<Integer>();
						temp.add(i);
						overlapResult.put(majorityVoteLower, temp);
					}

					selectedOverlap++;
				} else {
					rmList.add(overlapList.get(i));
				}
			}
		}

		System.out.println();
		System.out.print(params.getAlphaHigh());
		System.out.print("\t");
		System.out.print(doubleReadList.size());
		System.out.print("\t");
		System.out.print(overlapInTotal);
		System.out.print("\t");
		System.out.print(selectedOverlap);
		System.out.print("\t");
		System.out.print((selectedOverlap + 0.0) / overlapInTotal);
		System.out.print("\t");
		System.out.println(String
				.valueOf(((selectedOverlap + 0.0) / overlapInTotal) > params
						.getCoverage()));
		System.out.println("Overlap");
		for (Integer key : overlapResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			for (Integer groupIndex : overlapResult.get(key)) {
				for (Integer temp : overlapList.get(groupIndex)) {
					int i = doubleReadList.get(temp).getId()
							/ params.getForTest() + 1;
					if (tempMap.containsKey(i)) {
						tempMap.put(i, tempMap.get(i) + 1);
					} else {
						tempMap.put(i, 1);
					}
				}
			}
			for (Integer tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey)
						/ (overlapResult.get(key).size() + 0.0));
			}
			System.out.println();
		}

		int goodKey = 0;
		for (Integer key : overlapResult.keySet()) {
			int sumGroup = overlapResult.get(key).size();
			int sumSeqs = 0;
			for (Integer id : overlapResult.get(key)) {
				sumSeqs += overlapList.get(id).size();
			}
			if ((sumGroup * 20 >= overlapInTotal) && (sumSeqs >= 100)) {
				goodKey++;
			} else {
				for (Integer id : overlapResult.get(key)) {
					for (Integer seq : overlapList.get(id)) {
						for (int i = 0; i < params.getParticles(); i++) {
							z[seq][i] = 0;
						}
					}
					rmList.add(overlapList.get(id));
					selectedOverlap--;
					accumSeqCount -= overlapList.get(id).size();
				}
				// overlapResult.remove(key);
			}
		}

		if (goodKey < 2) {
			return null;
		}

		double noNewGroupAlpha = 0.00000001;

		mapCountList = DirichletClusterSingle.rebuildGroupMapList(z,
				params.getParticles());
		tempAccumCountList = DirichletClusterSingle.rebuildCountLists(
				doubleReadList, z, params.getParticles());

		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() == 1) {
				for (Integer elem : overlapList.get(i)) {
					idTagMap.put(elem, idTagMap.get(elem) + "Single");
					// update
					z[elem] = DirichletClusterSingle.clusterOverlapSeqs(
							mapCountList, tempAccumCountList, doubleReadList,
							overlapList.get(i), z, w, params.getParticles(),
							noNewGroupAlpha, accumSeqCount)[0];
					w = DirichletClusterSingle.updateWeights(
							tempAccumCountList, doubleReadList, elem, z, w,
							params.getParticles());

					double eff = 0.0;
					for (double d : w) {
						eff += d * d;
					}
					if (eff > 1 / (params.getThreshold() * params
							.getParticles())) {
						int[] resampleIndex = DirichletClusterSingle.sampleInt(
								w, params.getParticles(), 0);
						int[] temp = new int[params.getParticles()];
						for (int j = 0; j < params.getParticles(); j++) {
							w[j] = 1.0 / (0.0 + params.getParticles());
							temp[j] = z[elem][resampleIndex[j]];
						}
						z[elem] = temp;
					}

					tempAccumCountList = DirichletClusterSingle
							.updateCountLists(doubleReadList, z,
									overlapList.get(i), tempAccumCountList);
					mapCountList = DirichletClusterSingle.updateGroupMapList(
							mapCountList, z, overlapList.get(i));
				}
			}
		}

		for (Set<Integer> seqList : rmList) {
			// update
			int[][] tempZ = DirichletClusterSingle
					.clusterOverlapSeqs(mapCountList, tempAccumCountList,
							doubleReadList, seqList, z, w,
							params.getParticles(), noNewGroupAlpha,
							accumSeqCount);
			int majorityVote = DirichletClusterSingle
					.getMode(DirichletClusterSingle.getMode(tempZ))[0];
			for (Integer seqId : seqList) {
				idTagMap.put(seqId, idTagMap.get(seqId) + "Remove");
				for (int k = 0; k < params.getParticles(); k++) {
					z[seqId][k] = majorityVote;
				}
			}
			accumSeqCount += seqList.size();

			tempAccumCountList = DirichletClusterSingle.updateCountLists(
					doubleReadList, z, seqList, tempAccumCountList);
			mapCountList = DirichletClusterSingle.updateGroupMapList(
					mapCountList, z, seqList);
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
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] + zModeLower;
		}

		System.out.println();
		System.out.print(params.getAlphaHigh());
		System.out.print("\t");
		System.out.print(doubleReadList.size());
		System.out.print("\t");
		System.out.print(overlapInTotal);
		System.out.print("\t");
		System.out.print(selectedOverlap);
		System.out.print("\t");
		System.out.print((selectedOverlap + 0.0) / overlapInTotal);
		System.out.print("\t");
		System.out.println(String
				.valueOf(((selectedOverlap + 0.0) / overlapInTotal) > params
						.getCoverage()));
		System.out.println("Overlap");
		for (Integer key : overlapResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			for (Integer groupIndex : overlapResult.get(key)) {
				for (Integer temp : overlapList.get(groupIndex)) {
					int i = doubleReadList.get(temp).getId()
							/ params.getForTest() + 1;
					if (tempMap.containsKey(i)) {
						tempMap.put(i, tempMap.get(i) + 1);
					} else {
						tempMap.put(i, 1);
					}
				}
			}
			for (Integer tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey)
						/ (overlapResult.get(key).size() + 0.0));
			}
			System.out.println();
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
		return result;
		// if (((selectedOverlap + 0.0) / overlapInTotal) >
		// params.getCoverage()) {
		// return result;
		// } else {
		// return null;
		// }
	}

	private static void dirichletClusterTest(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, int[] zMode, Params params)
			throws CloneNotSupportedException {
		params.setAlphaHigh(10);
		params.setAlpha(0.000001);
		params.setAlphaLow(0.000001);
		dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		params.setAlphaHigh(1);
		dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		params.setAlphaHigh(0.1);
		dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		params.setAlphaHigh(0.01);
		dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
		params.setAlphaHigh(0.001);
		dirichletClusterSingleV1(doubleReadList, overlapList, params, 0);
	}

	private static int[] dirichletCluster(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, int[] zMode, Params params)
			throws CloneNotSupportedException {
		int[] updateZMode = zMode.clone();
		int tempZModeLower = 0;
		Set<Integer> zModeSet = new TreeSet<Integer>();
		for (int i = 0; i < updateZMode.length; i++) {
			zModeSet.add(updateZMode[i]);
		}
		double[] alphaHighArr = new double[] { 10, 5, 2, 1, 0.5, 0.1, 0.01,
				0.001 };
		System.out
				.println((new Date()).toString()
						+ "\tDirichlet Cluster start with alpha = "
						+ params.getAlpha());

		for (Integer elem : zModeSet) {

			System.out.println((new Date()).toString() + "\tzMode = " + elem);

			List<DoubleRead> subDoubleReadList = new ArrayList<DoubleRead>();
			for (int i = 0; i < zMode.length; i++) {
				if (elem.equals(zMode[i])) {
					subDoubleReadList.add(doubleReadList.get(i));
				}
			}

			System.out.println((new Date()).toString()
					+ "\tUpdating overlap infomation");

			List<Set<Integer>> subOverlapList = updateOverlap(
					subDoubleReadList, overlapList);
			Params tempParams = (Params) params.clone();
			tempParams.setSeqs(subDoubleReadList.size());
			tempParams.setAlpha(params.getAlpha());
			tempParams.setAlphaLow(params.getAlphaLow());

			System.out.println((new Date()).toString()
					+ "\tDirichlet single round begin");

			int[] subResult = new int[subDoubleReadList.size()];
			for (int i = 0; i < subResult.length; i++) {
				subResult[i] = tempZModeLower + 1;
			}

			boolean flag = false;

			for (double alphaHigh : alphaHighArr) {
				tempParams.setAlphaHigh(alphaHigh);
				int[] tempSubResult = dirichletClusterSingle(subDoubleReadList,
						subOverlapList, tempParams, tempZModeLower);
				if (!flag && (null != tempSubResult)) {
					subResult = tempSubResult;
					flag = true;
					break;
				}
			}

			// System.out.print("subResult");
			for (int i = 0; i < subResult.length; i++) {
				updateZMode[subDoubleReadList.get(i).getId()] = subResult[i];
				// System.out.print(" ");
				// System.out.print(subDoubleReadList.get(i).getId());
				if (subResult[i] > tempZModeLower) {
					tempZModeLower = subResult[i];
				}
			}
			// System.out.println();
			System.out.println((new Date()).toString()
					+ "\tDirichlet single round complete");

		}

		System.out.println((new Date()).toString()
				+ "\tDirichlet Cluster complete");

		return updateZMode;
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

		List<Set<Integer>> overlapList = EncodeBinTask.encodeMainJob(readList,
				32);

		System.out.println((new Date()).toString() + "\tComplete");

		// for (Set<Integer> overlapSet:overlapList) {
		// System.out.println(overlapSet.toString());
		// }

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

	public static void mainJob(String inputFile, String outputFile,
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

		List<Set<Integer>> overlapList = EncodeBinTask.encodeMainJob(readList,
				32);

		System.out.println((new Date()).toString() + "\tComplete");

		// for (Set<Integer> overlapSet:overlapList) {
		// System.out.println(overlapSet.toString());
		// }

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
		int oldZModeMax = 1;

		// TEST
		params.setAlpha(DefaultConstants.ALPHA);
		params.setForTest(params.getSeqs() / params.getForTest());
		while (true) {
			zMode = dirichletCluster(doubleReadList, overlapList, zMode, params);
			int zModeMax = oldZModeMax;
			for (int j = 0; j < zMode.length; j++) {
				if (zModeMax < zMode[j]) {
					zModeMax = zMode[j];
				}
			}

			if (oldZModeMax == zModeMax) {
				break;
			}
			oldZModeMax = zModeMax;
			// break;
		}

		System.out.println((new Date()).toString()
				+ "\tComplete\nWriting result...");

		PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFile)));
		for (int i = 0; i < zMode.length; i++) {
			pw.println(zMode[i]);
		}
		pw.close();

	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		frontEnd(args);
	}

}
