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
import java.util.TreeSet;

public class DirichletCluster {

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

		List<Double> alphaList = new ArrayList<Double>();
		String input = "/home/xinping/Desktop/6008/5_2345Borr.txt.sample";
		// String input = "/home/xinping/Desktop/6008/test.txt";
		// String input = "/home/xinping/Desktop/6008/5_2345Borr.txt";
		String output = null;

		Date date = new Date();

		Params params = new Params(particles, neighbor, seqs, transOrder,
				majority, threshold, alpha, alphaLow, alphaHigh, zero,
				DefaultConstants.TESTNUM, DefaultConstants.COVERAGE);

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
					if (args[pos].length() == 2) {
						params.setAlpha(parse(DefaultConstants.ALPHA,
								args[pos + 1], "Parameter alpha"));
						pos += 2;
					} else {
						int alphaListSize = Integer.parseInt(args[pos]
								.substring(2));

						for (int i = 1; i <= alphaListSize; i++) {
							alphaList.add(Double.valueOf(args[pos + i]));
						}

						pos += (alphaListSize + 1);
					}
					break;
				case 'h':
				case 'H':
					params.setAlphaHigh(parse(DefaultConstants.ALPHAHIGH,
							args[pos + 1], "Parameter alpha_high"));
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

		mainJob(input, output, params, alphaList);
	}

	public static void printHelp() {
		System.out.println("Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " -i abundance_species_equal.txt");
		System.out.println("The file contains pair-end sequence\n");
		System.out.println("More Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " [OPTION] abundance_species_equal.txt\n");
		System.out.println("Options:");
		System.out
				.println("\t-a:\tThe value followed will be the alpha parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ALPHA);
		System.out.println("\t\tExample 1: -a3 0.0000001 0.00001 0.001");
		System.out.println("\t\tExample 2: -a 0.0000001");
		System.out
				.println("\t-h:\tThe value followed will be the alpha_high parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.ALPHAHIGH);
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

	private static int[] dirichletCluster(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, int[] zMode, Params params)
			throws CloneNotSupportedException {
		int[] updateZMode = zMode.clone();
		int tempZModeLower = 0;
		Set<Integer> zModeSet = new TreeSet<Integer>();
		for (int i = 0; i < updateZMode.length; i++) {
			zModeSet.add(updateZMode[i]);
		}

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
			// System.out.println("subLength = " + subDoubleReadList.size());

			System.out.println((new Date()).toString()
					+ "\tUpdating overlap infomation");

			List<Set<Integer>> subOverlapList = updateOverlap(
					subDoubleReadList, overlapList);
			Params tempParams = (Params) params.clone();
			tempParams.setSeqs(subDoubleReadList.size());

			System.out.println((new Date()).toString()
					+ "\tDirichlet single round begin");

			// for (Set<Integer> set:subOverlapList) {
			// System.out.println(set.toString());
			// }

			int[] subResult = DirichletClusterSingle.dirichletClusterSingle(
					subDoubleReadList, subOverlapList, tempParams,
					tempZModeLower);

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
			System.out.println("subResult");
			System.out.println(tempParams.getAlpha());
			System.out.println(elem);
			System.out.println(tempZModeLower);
			System.out.println(Arrays.toString(subResult));
		}

		System.out.println((new Date()).toString()
				+ "\tDirichlet Cluster complete");

		return updateZMode;
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

		if (alphaList.size() > 0) {

			for (int i = 0; i < alphaList.size(); i++) {
				params.setAlpha(alphaList.get(i));
				// System.out.println(alphaList.get(i));
				zMode = dirichletCluster(doubleReadList, overlapList, zMode,
						params);
				int zModeMax = oldZModeMax;
				for (int j = 0; j < zMode.length; j++) {
					if (zModeMax < zMode[j]) {
						zModeMax = zMode[j];
					}
				}

				// if (oldZModeMax == zModeMax) {
				// break;
				// }
				oldZModeMax = zModeMax;
			}
		} else {
			zMode = dirichletCluster(doubleReadList, overlapList, zMode, params);
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
