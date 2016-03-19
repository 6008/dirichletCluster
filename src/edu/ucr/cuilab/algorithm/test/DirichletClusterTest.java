package edu.ucr.cuilab.algorithm.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.ucr.cuilab.algorithm.DefaultConstants;
import edu.ucr.cuilab.algorithm.DoubleRead;
import edu.ucr.cuilab.algorithm.DoubleReadCompareGC;
import edu.ucr.cuilab.algorithm.DoubleReadCompareID;
import edu.ucr.cuilab.algorithm.EncodeBinTask;
import edu.ucr.cuilab.algorithm.ListCompareLength;
import edu.ucr.cuilab.algorithm.Params;

public class DirichletClusterTest {

	private static final int NUM = 1000;

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

		}
		Collections.sort(doubleReadList, new DoubleReadCompareID());
	}

	private static int[] compressResult(int[] datas) {
		Map<Integer, Integer> keyMap = new TreeMap<Integer, Integer>();
		int[] result = new int[datas.length];
		int key = 1;
		for (int i : datas) {
			if (!keyMap.containsKey(i)) {
				keyMap.put(i, key);
				key++;
			}
		}
		for (int i = 0; i < datas.length; i++) {
			result[i] = keyMap.get(datas[i]);
		}
		return result;
	}

	// return {mode, length of mode}
	private static int[] getMode(int[] datas) {
		Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
		for (int i : datas) {
			if (tempMap.containsKey(i)) {
				tempMap.put(i, tempMap.get(i) + 1);
			} else {
				tempMap.put(i, 1);
			}
		}

		int[] result = new int[] { 0, 0 };

		for (Integer key : tempMap.keySet()) {
			if (tempMap.get(key) >= result[1]) {
				result[0] = key;
				result[1] = tempMap.get(key);
			}
		}

		return result;
	}

	private static int[] getMode(int[][] datas) {
		int[] result = new int[datas.length];

		for (int i = 0; i < datas.length; i++) {
			result[i] = getMode(datas[i])[0];
		}

		return result;
	}

	private static boolean checkMajorityVote(int[][] datas, double majority) {
		return (majority * datas.length < getMode(getMode(datas))[1]);
	}

	private static int max(int[][] data) {
		int result = 0;
		for (int i = 0; i < data.length; i++) {
			for (int d : data[i]) {
				if (d > result) {
					result = d;
				}
			}
		}
		return result;
	}

	private static List<Map<Integer, Integer>> updateGroupMapList(
			List<Map<Integer, Integer>> groupMapList, int[][] z,
			Set<Integer> seqList) {
		for (int i = 0; i < groupMapList.size(); i++) {
			Map<Integer, Integer> temp = groupMapList.get(i);
			for (Integer j : seqList) {
				if (temp.containsKey(z[j][i])) {
					temp.put(z[j][i], temp.get(z[j][i]) + 1);
				} else {
					temp.put(z[j][i], 1);
				}
			}
		}
		return groupMapList;
	}

	private static double[] priorP(List<Map<Integer, Integer>> groupMapList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha,
			int completeSeqCount) {
		int group = max(z);
		double[][] pList = new double[particles][group + 1];
		double[] colSum = new double[group + 1];
		for (int j = 0; j < group + 1; j++) {
			colSum[j] = 0.0;
		}
		for (int i = 0; i < particles; i++) {
			Map<Integer, Integer> groupMap = groupMapList.get(i);
			int[] groupCol = new int[groupMap.keySet().size()];
			int temp = 0;
			for (Integer id : groupMap.keySet()) {
				groupCol[temp] = id;
				temp++;
			}
			for (int j = 0; j < group + 1; j++) {
				pList[i][j] = w[i] * alpha / (alpha + completeSeqCount)
						/ (group + 1 - groupCol.length);
			}
			for (int j : groupCol) {
				pList[i][j - 1] = w[i] * groupMap.get(j)
						/ (completeSeqCount + alpha);
			}
			for (int j = 0; j < group + 1; j++) {
				colSum[j] += pList[i][j];
			}
		}
		return colSum;
	}

	private static double[] startPCounts(List<Integer> counts) {
		double[] probs = new double[counts.size() / 4];
		double sum = 0.0;
		for (int i = 0; i < probs.length; i++) {
			probs[i] = counts.get(i * 4) + counts.get(i * 4 + 1)
					+ counts.get(i * 4 + 2) + counts.get(i * 4 + 3);
			sum += probs[i];
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		return probs;
	}

	private static double[] transPCounts(List<Integer> counts) {
		double[] probs = new double[counts.size()];
		for (int i = 0; i < probs.length; i += 4) {
			double sum = counts.get(i) + counts.get(i + 1) + counts.get(i + 2)
					+ counts.get(i + 3);
			probs[i] = counts.get(i) / sum;
			probs[i + 1] = counts.get(i + 1) / sum;
			probs[i + 2] = counts.get(i + 2) / sum;
			probs[i + 3] = counts.get(i + 3) / sum;
		}
		return probs;
	}

	private static List<Map<Integer, List<Integer>>> updateCountLists(
			List<DoubleRead> doubleReadList, int[][] z, Set<Integer> seqList,
			List<Map<Integer, List<Integer>>> originCountLists) {

		for (int i = 0; i < originCountLists.size(); i++) {
			Map<Integer, List<Integer>> listMap = originCountLists.get(i);
			for (Integer j : seqList) {
				if (listMap.containsKey(z[j][i])) {
					List<Integer> temp = new ArrayList<Integer>(
							listMap.get(z[j][i]));
					for (int k = 0; k < temp.size(); k++) {
						temp.set(k, temp.get(k)
								+ doubleReadList.get(j).getCountList().get(k));
					}
					listMap.put(z[j][i], temp);
				} else {
					listMap.put(z[j][i], doubleReadList.get(j).getCountList());
				}
			}
			originCountLists.set(i, listMap);
		}

		return originCountLists;
	}

	private static double[][][] startP(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			int particles) {
		int group = max(z);
		double[][][] p = new double[particles][group + 1][];
		double[] newStartP = doubleReadList.get(seqIndex).getNewGroupStartp();
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				p[k][j] = new double[newStartP.length];
			}
			for (Integer j : countLists.get(k).keySet()) {
				p[k][j - 1] = startPCounts(countLists.get(k).get(j));
			}
			p[k][group] = newStartP;
		}
		return p;
	}

	private static double[][][] transP(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			int particles) {
		int group = max(z);
		double[] newTransP = doubleReadList.get(seqIndex).getNewGroupTransp();
		double[][][] p = new double[particles][group + 1][];
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				p[k][j] = new double[newTransP.length];
			}
			for (Integer j : countLists.get(k).keySet()) {
				p[k][j - 1] = transPCounts(countLists.get(k).get(j));
			}
			p[k][group] = newTransP;
		}
		return p;
	}

	private static double[] postPTemp(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles) {

		int group = max(z);

		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles);
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles);

		double[][] p = new double[particles][group + 1];
		double[] result = new double[group + 1];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();
		List<Integer> countList = doubleReadList.get(seqIndex).getCountList();
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group + 1; j++) {
				double sumTempStartP = 0.0;
				double sumTempTransP = 0.0;
				for (int i = 0; i < tempTransP[k][j].length; i++) {
					sumTempTransP += tempTransP[k][j][i];
				}
				for (int i = 0; i < tempStartP[k][j].length; i++) {
					sumTempStartP += tempStartP[k][j][i];
				}
				if ((sumTempTransP > Double.MIN_NORMAL)
						&& (sumTempStartP > Double.MIN_NORMAL)) {
					double[] temp = tempTransP[k][j];
					double tempSum = 0.0;
					for (int i = 0; i < temp.length; i++) {
						if (temp[i] < Double.MIN_NORMAL) {
							tempSum += countList.get(i) * DefaultConstants.ZERO;
						} else {
							tempSum += countList.get(i) * Math.log(temp[i]);
						}
					}
					p[k][j] = w[k] * Math.exp(tempSum);

					for (int i = 0; i < startWhere.length; i++) {
						p[k][j] *= tempStartP[k][j][startWhere[i]];
					}
				} else {
					p[k][j] = 0.0;
				}
			}
		}

		for (int i = 0; i < particles; i++) {
			for (int j = 0; j < group + 1; j++) {
				result[j] += p[i][j];
			}
		}
		return result;
	}

	private static double[] postP(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount) {
		int group = max(z);
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(countLists, doubleReadList, seqIndex, z,
				w, particles);

		double[] tempPriorP = priorP(groupCount, seqIndex, z, w, particles,
				alpha, completeSeqCount);
		double sum = 0.0;
		for (int i = 0; i < p.length; i++) {
			p[i] = tempPostP[i] * tempPriorP[i];
			sum += p[i];
		}
		if (sum < Double.MIN_NORMAL) {
			sum = Double.MIN_NORMAL;
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}

	private static int[] sampleInt(double[] probList, int particles, int start) {
		int[] samples = new int[particles];
		double[] accumList = new double[probList.length];
		accumList[0] = probList[0];
		for (int i = 1; i < probList.length; i++) {
			accumList[i] = probList[i] + accumList[i - 1];
		}
		for (int i = 0; i < probList.length; i++) {
			accumList[i] = accumList[i] / accumList[probList.length - 1];
		}
		Random rand = new Random();
		for (int i = 0; i < particles; i++) {
			double r = rand.nextDouble();
			int j = 0;
			while ((j < accumList.length) && (accumList[j] < r)) {
				j++;
			}
			samples[i] = j + start;
		}
		return samples;
	}

	private static int[] clusterOneSeq(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount) {
		double[] posterior = postP(groupCount, countLists, doubleReadList,
				seqIndex, z, w, particles, alpha, completeSeqCount);

		int[] sample = sampleInt(posterior, particles, 1);
		return sample;
	}

	private static int[][] clusterOverlapSeqs(
			List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, Set<Integer> seqList, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount) {
		int[][] tempZ = new int[seqList.size()][particles];
		int flag = 0;
		for (Integer elem : seqList) {
			tempZ[flag] = clusterOneSeq(groupCount, countLists, doubleReadList,
					elem, z, w, particles, alpha, completeSeqCount);
			flag++;
		}
		return tempZ;
	}

	private static double[] updateWeights(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles) {
		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles);
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles);

		double[] p = new double[particles];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		List<Integer> countList = doubleReadList.get(seqIndex).getCountList();
		for (int k = 0; k < particles; k++) {
			double[] temp = tempTransP[k][z[seqIndex][k] - 1];
			double tempSum = 0.0;
			for (int i = 0; i < temp.length; i++) {
				if (temp[i] < Double.MIN_NORMAL) {
					tempSum += countList.get(i) * DefaultConstants.ZERO;
				} else {
					tempSum += countList.get(i) * Math.log(temp[i]);
				}
			}
			p[k] = w[k] * Math.exp(tempSum);
			for (int i = 0; i < startWhere.length; i++) {
				p[k] *= tempStartP[k][z[seqIndex][k] - 1][startWhere[i]];
			}
		}

		double sum = 0.0;
		for (double elem : p) {
			sum += elem;
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}

	public static int[] dirichletClusterSingle(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, Params params, int zModeLower,
			PrintWriter pw) {
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
		tempAccumCountList = updateCountLists(doubleReadList, z,
				overlapList.get(0), tempAccumCountList);
		mapCountList = updateGroupMapList(mapCountList, z, overlapList.get(0));

		newGroup(doubleReadList, params.getNeighbor());

		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();

		int accumSeqCount = overlapList.get(0).size();

		int overlapInTotal = 0;
		Map<Integer, Set<Integer>> overlapResult = new HashMap<Integer, Set<Integer>>();

		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() > 1) {
				overlapInTotal++;
				if (tempAccumCountList.get(0).size() > 1) {
				}
				int[][] tempZLower = clusterOverlapSeqs(mapCountList,
						tempAccumCountList, doubleReadList, overlapList.get(i),
						z, w, params.getParticles(), params.getAlphaLow(),
						accumSeqCount);
				int[][] tempZUpper = clusterOverlapSeqs(mapCountList,
						tempAccumCountList, doubleReadList, overlapList.get(i),
						z, w, params.getParticles(), params.getAlphaHigh(),
						accumSeqCount);

				boolean isMarjorityLower = checkMajorityVote(tempZLower,
						params.getMajority());
				boolean isMarjorityUpper = checkMajorityVote(tempZUpper,
						params.getMajority());

				int[] voteLower = getMode(getMode(tempZLower));
				int[] voteUpper = getMode(getMode(tempZUpper));

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

					tempAccumCountList = updateCountLists(doubleReadList, z,
							overlapList.get(i), tempAccumCountList);
					mapCountList = updateGroupMapList(mapCountList, z,
							overlapList.get(i));

					accumSeqCount += overlapList.get(i).size();
					if (overlapResult.containsKey(majorityVoteLower)) {
						overlapResult.get(majorityVoteLower).addAll(
								overlapList.get(i));
					} else {
						overlapResult
								.put(majorityVoteLower, overlapList.get(i));
					}

				} else {
					rmList.add(overlapList.get(i));
				}
			} else {
				for (Integer elem : overlapList.get(i)) {
					idTagMap.put(elem, idTagMap.get(elem) + "Single");
					z[elem] = clusterOneSeq(mapCountList, tempAccumCountList,
							doubleReadList, elem, z, w, params.getParticles(),
							params.getAlpha(), accumSeqCount);
					w = updateWeights(tempAccumCountList, doubleReadList, elem,
							z, w, params.getParticles());

					double eff = 0.0;
					for (double d : w) {
						eff += d * d;
					}

					if (eff > 1 / (params.getThreshold() * params
							.getParticles())) {
						int[] resampleIndex = sampleInt(w,
								params.getParticles(), 0);
						int[] temp = new int[params.getParticles()];
						for (int j = 0; j < params.getParticles(); j++) {
							w[j] = 1.0 / (0.0 + params.getParticles());
							temp[j] = z[elem][resampleIndex[j]];
						}
						z[elem] = temp;
					}

					tempAccumCountList = updateCountLists(doubleReadList, z,
							overlapList.get(i), tempAccumCountList);
					mapCountList = updateGroupMapList(mapCountList, z,
							overlapList.get(i));
				}

				accumSeqCount += overlapList.get(i).size();
			}

		}

		// pw.println(String.valueOf(doubleReadList.size()) + "\t"
		// + String.valueOf(overlapInTotal) + "\t"
		// + String.valueOf(rmList.size()));
		pw.print(params.getAlphaHigh());
		pw.print("\t");
		pw.print(doubleReadList.size());
		pw.print("\t");
		if (overlapInTotal > 0) {
			pw.print((overlapInTotal - rmList.size() + 0.0) / overlapInTotal);
		} else {
			pw.print(0);
		}
		double powerEntropy = 0.0;
		int seqNum = 1;
		for (Integer key : overlapResult.keySet()) {
			// pw.println(key);
			// pw.println(overlapResult.get(key).toString());
			Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
			for (Integer temp : overlapResult.get(key)) {
				int i = doubleReadList.get(temp).getId() / NUM + 1;
				if (tempMap.containsKey(i)) {
					tempMap.put(i, tempMap.get(i) + 1);
				} else {
					tempMap.put(i, 1);
				}
			}
			seqNum += overlapResult.get(key).size();
			double entropy = 0.0;
			for (Integer temp : tempMap.keySet()) {
				System.out.println(String.valueOf(temp) + "\t"
						+ String.valueOf(tempMap.get(temp)));
				entropy -= tempMap.get(temp)
						* Math.log((tempMap.get(temp) + 0.0)
								/ overlapResult.get(key).size());
			}
			powerEntropy += entropy;
		}
		pw.print("\t");
		pw.println(powerEntropy / seqNum);
		for (Set<Integer> seqList : rmList) {

			int[][] tempZ = clusterOverlapSeqs(mapCountList,
					tempAccumCountList, doubleReadList, seqList, z, w,
					params.getParticles(), params.getAlpha(), accumSeqCount);
			int majorityVote = getMode(getMode(tempZ))[0];
			for (Integer seqId : seqList) {
				idTagMap.put(seqId, idTagMap.get(seqId) + "Remove");
				for (int k = 0; k < params.getParticles(); k++) {
					z[seqId][k] = majorityVote;
				}
			}
			accumSeqCount += seqList.size();
			tempAccumCountList = updateCountLists(doubleReadList, z, seqList,
					tempAccumCountList);
			mapCountList = updateGroupMapList(mapCountList, z, seqList);
		}

		int[] result = compressResult(getMode(z));
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] + zModeLower;
		}

		return result;
	}

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
			List<Set<Integer>> overlapList, int[] zMode, Params params,
			PrintWriter pw) throws CloneNotSupportedException {
		int[] updateZMode = zMode.clone();
		int tempZModeLower = 0;
		Set<Integer> zModeSet = new TreeSet<Integer>();
		for (int i = 0; i < updateZMode.length; i++) {
			zModeSet.add(updateZMode[i]);
		}

		for (Integer elem : zModeSet) {
			List<DoubleRead> subDoubleReadList = new ArrayList<DoubleRead>();
			for (int i = 0; i < zMode.length; i++) {
				if (elem.equals(zMode[i])) {
					subDoubleReadList.add(doubleReadList.get(i));
				}
			}

			List<Set<Integer>> subOverlapList = updateOverlap(
					subDoubleReadList, overlapList);
			Params tempParams = new Params(params);
			tempParams.setSeqs(subDoubleReadList.size());
			int[] subResult = dirichletClusterSingle(subDoubleReadList,
					subOverlapList, tempParams, tempZModeLower, pw);
			for (int i = 0; i < subResult.length; i++) {
				updateZMode[subDoubleReadList.get(i).getId()] = subResult[i];
				if (subResult[i] > tempZModeLower) {
					tempZModeLower = subResult[i];
				}
			}
		}
		return updateZMode;
	}

	public static void frontEnd() throws Exception {
		String input_4 = "/home/xinping/Desktop/6008/5_2345Borr.txt.sample";
		String input_3 = "/home/xinping/Desktop/6008/5_2345Borr.txt";
		String input_1 = "/home/xinping/Desktop/6008/5_test.txt";
		String input = "/home/xinping/Desktop/6008/5_2345Borr_5W.txt";
		List<String> readList = readFile(input);
		int seqs = readList.size() / 2;
		int transOrder = DefaultConstants.TRANSORDER;
		int particles = DefaultConstants.PARTICLES;
		int neighbor = DefaultConstants.NEIGHBOR;

		double alpha = DefaultConstants.ALPHA;
		double alphaLow = DefaultConstants.ALPHALOW;
		// double alphaHigh = DefaultConstants.ALPHAHIGH;
		double threshold = DefaultConstants.THRESHOLD;
		double zero = DefaultConstants.ZERO;
		double majority = DefaultConstants.MAJORITY;
		char[] charList = { 'A', 'C', 'G', 'T' };
		List<String> permutationList = getPermutations(charList, transOrder + 1);
		List<Double> alphaList = new ArrayList<Double>();
		alphaList.add(0.00000001);
		alphaList.add(0.000001);
		alphaList.add(0.0001);
		List<Set<Integer>> overlapList = EncodeBinTask.encodeMainJob(readList,
				32);
		List<DoubleRead> doubleReadList = new ArrayList<DoubleRead>();
		for (int i = 0; i < seqs; i++) {
			DoubleRead dr = new DoubleRead(i, readList.get(i * 2),
					readList.get(i * 2 + 1), transOrder + 1, permutationList);
			doubleReadList.add(dr);
		}

		double[] alphaHighArr = new double[] { 10, 5, 2, 1, 0.5, 0.1, 0.01,
				0.001 };
		String output = input + "alphaHigh.output";
		PrintWriter pw = new PrintWriter(new FileWriter(new File(output)));
		for (double alphaHigh : alphaHighArr) {

			Params params = new Params(particles, neighbor, seqs, transOrder,
					majority, threshold, alpha, alphaLow, alphaHigh, zero,
					DefaultConstants.TESTNUM, DefaultConstants.COVERAGE);

			int[] zMode = new int[params.getSeqs()];
			for (int i = 0; i < params.getSeqs(); i++) {
				zMode[i] = 1;
			}
			int oldZModeMax = 1;

			if (alphaList.size() > 0) {

				for (int i = 0; i < alphaList.size(); i++) {
					params.setAlpha(alphaList.get(i));
					zMode = dirichletCluster(doubleReadList, overlapList,
							zMode, params, pw);
					int zModeMax = oldZModeMax;
					for (int j = 0; j < zMode.length; j++) {
						if (zModeMax < zMode[j]) {
							zModeMax = zMode[j];
						}
					}
					oldZModeMax = zModeMax;
				}
			} else {
				zMode = dirichletCluster(doubleReadList, overlapList, zMode,
						params, pw);
			}
		}
		pw.close();
	}

	public static void main(String[] args) throws Exception {
		frontEnd();
	}

}
