package edu.ucr.cuilab.algorithm;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class DirichletClusterSingle {

	public static long[] timeArr = new long[2];

	public static void newGroup(List<DoubleRead> doubleReadList, int neighbor) {
		Collections.sort(doubleReadList, new DoubleReadCompareGC());
		for (int i = 0; i < doubleReadList.size(); i++) {
			// System.out.print(doubleReadList.get(i).getId());
			// System.out.print('\t');
			// System.out.println(doubleReadList.get(i).getGC());
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
			// System.out.println(doubleReadList.get(i).getId() + 1);
			// System.out.println(Arrays.toString(doubleReadList.get(i).getNewGroupStartp()));

		}
		Collections.sort(doubleReadList, new DoubleReadCompareID());
		
	}

	public static int[] compressResult(int[] datas) {
		Map<Integer, Integer> keyMap = new TreeMap<Integer, Integer>();
		keyMap.put(0, 0);
		int[] result = new int[datas.length];
		int key = 1;
		for (int i : datas) {
			if (i != 0) {
				if (!keyMap.containsKey(i)) {
					keyMap.put(i, key);
					key++;
				}
			}
		}
		for (int i = 0; i < datas.length; i++) {
			result[i] = keyMap.get(datas[i]);
		}
		return result;
	}

	// return {mode, length of mode}
	public static int[] getMode(int[] datas) {
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

	public static int[] getMode(int[][] datas) {
		int[] result = new int[datas.length];

		for (int i = 0; i < datas.length; i++) {
			result[i] = getMode(datas[i])[0];
		}

		return result;
	}

	public static boolean checkMajorityVote(int[][] datas, double majority) {
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

	private static Map<Integer, Integer> groupCol(int[][] data, int k) {
		Map<Integer, Integer> dataSet = new HashMap<Integer, Integer>();
		for (int i = 0; i < data.length; i++) {
			if (data[i][k] != 0) {
				if (dataSet.containsKey(data[i][k])) {
					dataSet.put(data[i][k], dataSet.get(data[i][k]) + 1);
				} else {
					dataSet.put(data[i][k], 1);
				}
			}
		}
		return dataSet;
	}

	private static double[] priorP(int seqIndex, int[][] z, double[] w,
			int particles, double alpha, int completeSeqCount) {
		int group = max(z);
		double[][] pList = new double[particles][group + 1];
		double[] colSum = new double[group + 1];
		for (int j = 0; j < group + 1; j++) {
			colSum[j] = 0.0;
		}
		for (int i = 0; i < particles; i++) {
			Map<Integer, Integer> groupMap = groupCol(z, i);
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

	public static List<Map<Integer, Integer>> updateGroupMapList(
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

	public static List<Map<Integer, Integer>> rebuildGroupMapList(int[][] z,
			int particles) {
		List<Map<Integer, Integer>> result = new ArrayList<Map<Integer, Integer>>();
		for (int i = 0; i < particles; i++) {
			result.add(groupCol(z, i));
		}
		return result;
	}

	// private static double[] priorP(List<Map<Integer, Integer>> groupMapList,
	// int seqIndex, int[][] z, double[] w, int particles, double alpha,
	// int completeSeqCount) {
	// int group = max(z);
	// double[][] pList = new double[particles][group + 1];
	// double[] colSum = new double[group + 1];
	// for (int j = 0; j < group + 1; j++) {
	// colSum[j] = 0.0;
	// }
	// for (int i = 0; i < particles; i++) {
	// Map<Integer, Integer> groupMap = groupCol(z, i);
	// int[] groupCol = new int[groupMap.keySet().size()];
	// int temp = 0;
	// for (Integer id : groupMap.keySet()) {
	// groupCol[temp] = id;
	// temp++;
	// }
	// for (int j = 0; j < group + 1; j++) {
	// pList[i][j] = w[i] * alpha / (alpha + completeSeqCount)
	// / (group + 1 - groupCol.length);
	// }
	// for (int j : groupCol) {
	// pList[i][j - 1] = w[i] * groupMap.get(j)
	// / (completeSeqCount + alpha);
	// }
	// for (int j = 0; j < group + 1; j++) {
	// colSum[j] += pList[i][j];
	// }
	// }
	// return colSum;
	// }

	private static double[] priorP(List<Map<Integer, Integer>> groupMapList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha,
			int completeSeqCount, int group) {
		// int group = max(z);
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

	private static double[] startPCounts(int[] counts) {
		double[] probs = new double[counts.length / 4];
		double sum = 0.0;
		for (int i = 0; i < probs.length; i++) {
			probs[i] = counts[i * 4] + counts[i * 4 + 1] + counts[i * 4 + 2]
					+ counts[i * 4 + 3];
			sum += probs[i];
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		return probs;
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

	private static double[] transPCounts(int[] counts) {
		double[] probs = new double[counts.length];
		for (int i = 0; i < probs.length; i += 4) {
			double sum = counts[i] + counts[i + 1] + counts[i + 2]
					+ counts[i + 3];
			probs[i] = counts[i] / sum;
			probs[i + 1] = counts[i + 1] / sum;
			probs[i + 2] = counts[i + 2] / sum;
			probs[i + 3] = counts[i + 3] / sum;
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

	public static List<Map<Integer, List<Integer>>> rebuildCountLists(
			List<DoubleRead> doubleReadList, int[][] z, int particles) {

		List<Map<Integer, List<Integer>>> originCountLists = new ArrayList<Map<Integer, List<Integer>>>();
		for (int j = 0; j < particles; j++) {
			Map<Integer, List<Integer>> countMap = new HashMap<Integer, List<Integer>>();
			for (int i = 0; i < doubleReadList.size(); i++) {
				if (z[i][j] != 0) {
					if (countMap.containsKey(z[i][j])) {
						List<Integer> temp = new ArrayList<Integer>(
								countMap.get(z[i][j]));
						for (int k = 0; k < temp.size(); k++) {
							temp.set(k,
									temp.get(k)
											+ doubleReadList.get(i)
													.getCountList().get(k));
						}
						countMap.put(z[i][j], temp);
					} else {
						countMap.put(z[i][j], new ArrayList<Integer>(
								doubleReadList.get(i).getCountList()));
					}
				}
			}
			originCountLists.add(countMap);
		}
		return originCountLists;
	}

	public static List<Map<Integer, List<Integer>>> updateCountLists(
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
					listMap.put(z[j][i], new ArrayList<Integer>(doubleReadList
							.get(j).getCountList()));
				}
			}
			originCountLists.set(i, listMap);
		}

		return originCountLists;
	}

	private static double[][][] startP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, int particles) {
		int group = max(z);
		double[][][] p = new double[particles][group + 1][];
		double[] newStartP = doubleReadList.get(seqIndex).getNewGroupStartp();
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				List<Integer> clusterIndex = new ArrayList<Integer>();
				for (int i = 0; i < z.length; i++) {
					if (z[i][k] == j + 1) {
						clusterIndex.add(i);
					}
				}
				if (clusterIndex.size() > 0) {
					List<Integer> temp = new ArrayList<Integer>(doubleReadList
							.get(clusterIndex.get(0)).getCountList());
					for (int i = 1; i < clusterIndex.size(); i++) {
						for (int t = 0; t < temp.size(); t++) {
							temp.set(
									t,
									temp.get(t)
											+ doubleReadList
													.get(clusterIndex.get(i))
													.getCountList().get(t));
						}
					}
					p[k][j] = startPCounts(temp);
				} else {
					p[k][j] = new double[newStartP.length];
				}
			}
			p[k][group] = newStartP;
		}
		return p;
	}

	private static double[][][] transP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, int particles) {

		int group = max(z);
		double[] newTransP = doubleReadList.get(seqIndex).getNewGroupTransp();
		double[][][] p = new double[particles][group + 1][];
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				List<Integer> clusterIndex = new ArrayList<Integer>();
				for (int i = 0; i < z.length; i++) {
					if (z[i][k] == j + 1) {
						clusterIndex.add(i);
					}
				}
				if (clusterIndex.size() > 0) {
					List<Integer> temp = new ArrayList<Integer>(doubleReadList
							.get(clusterIndex.get(0)).getCountList());
					for (int i = 1; i < clusterIndex.size(); i++) {
						for (int t = 0; t < temp.size(); t++) {
							temp.set(
									t,
									temp.get(t)
											+ doubleReadList
													.get(clusterIndex.get(i))
													.getCountList().get(t));
						}
					}
					p[k][j] = transPCounts(temp);
				} else {
					p[k][j] = new double[newTransP.length];
				}

			}
			p[k][group] = newTransP;
		}

		return p;
	}

	private static double[] postPTemp(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles) {

		int group = max(z);
		double[][][] tempTransP = transP(doubleReadList, seqIndex, z, particles);

		double[][][] tempStartP = startP(doubleReadList, seqIndex, z, particles);

		// System.out.println(Arrays.toString(tempStartP[0][1]));
		// System.out.println(Arrays.toString(tempTransP[0][0]));

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

	private static double[][][] startP(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			int particles, int group) {
		// int group = max(z);
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
			int particles, int group) {
		// int group = max(z);
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
			double[] w, int particles, int group) {

		// int group = max(z);

		Date date1 = new Date();
		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles, group);
		Date date2 = new Date();
		// System.out.println("transP " + (date2.getTime() - date1.getTime()));
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles, group);
		Date date3 = new Date();
		// System.out.println("tempStartP " + (date3.getTime() -
		// date2.getTime()));
		timeArr[0] += date2.getTime() - date1.getTime();
		timeArr[1] += date3.getTime() - date2.getTime();

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
							tempSum += countList.get(i)
									* Math.log(DefaultConstants.ZERO);
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

		// System.out.println("postPTemp "
		// + ((new Date()).getTime() - date3.getTime()));

		return result;
	}

	private static double[] postP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha,
			int completeSeqCount) {
		int group = max(z);
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(doubleReadList, seqIndex, z, w,
				particles);

		// System.out.println(seqIndex);
		// System.out.println(Arrays.toString(tempPostP));
		double[] tempPriorP = priorP(seqIndex, z, w, particles, alpha,
				completeSeqCount);
		// System.out.println(Arrays.toString(tempPriorP));
		double sum = 0.0;
		for (int i = 0; i < p.length; i++) {
			p[i] = tempPostP[i] * tempPriorP[i] + Double.MIN_NORMAL;
			sum += p[i];
		}
		// if (sum < Double.MIN_NORMAL) {
		// sum = Double.MIN_NORMAL;
		// }
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}

	private static double[] postP(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount) {
		int group = max(z);
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(countLists, doubleReadList, seqIndex, z,
				w, particles, group);

		// Date date = new Date();
		// System.out.println(seqIndex);
		// System.out.println(doubleReadList.get(seqIndex).getCountList()
		// .toString());
		// System.out.println(Arrays.toString(tempPostP));
		double[] tempPriorP = priorP(groupCount, seqIndex, z, w, particles,
				alpha, completeSeqCount, group);
		// System.out.println(Arrays.toString(tempPriorP));
		// System.out.println("priorP "
		// + ((new Date()).getTime() - date.getTime()));
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

	private static double[] postP(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount,
			int group) {
		// int group = max(z);
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(countLists, doubleReadList, seqIndex, z,
				w, particles, group);

		// Date date = new Date();
		// System.out.println(seqIndex);
		// System.out.println(doubleReadList.get(seqIndex).getCountList()
		// .toString());
		// System.out.println(Arrays.toString(tempPostP));
		double[] tempPriorP = priorP(groupCount, seqIndex, z, w, particles,
				alpha, completeSeqCount, group);
		// System.out.println(Arrays.toString(tempPriorP));
		// System.out.println("priorP "
		// + ((new Date()).getTime() - date.getTime()));
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

	public static int[] sampleInt(double[] probList, int particles, int start) {
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

	public static int[] clusterOneSeq(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha,
			int completeSeqCount) {
		double[] posterior = postP(doubleReadList, seqIndex, z, w, particles,
				alpha, completeSeqCount);

		// System.out.println(Arrays.toString(posterior));
		int[] sample = sampleInt(posterior, particles, 1);
		// System.out.println(Arrays.toString(sample));
		return sample;
	}

	public static int[] clusterOneSeq(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount) {
		double[] posterior = postP(groupCount, countLists, doubleReadList,
				seqIndex, z, w, particles, alpha, completeSeqCount);

		// Date date = new Date();
		// System.out.println(Arrays.toString(posterior));
		int[] sample = sampleInt(posterior, particles, 1);
		// System.out.println(Arrays.toString(sample));
		// System.out.println("sampleInt "
		// + ((new Date()).getTime() - date.getTime()));
		return sample;
	}

	public static int[] clusterOneSeq(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount,
			int group) {
		double[] posterior = postP(groupCount, countLists, doubleReadList,
				seqIndex, z, w, particles, alpha, completeSeqCount, group);

		// Date date = new Date();
		// System.out.println(Arrays.toString(posterior));
		int[] sample = sampleInt(posterior, particles, 1);
		// System.out.println(Arrays.toString(sample));
		// System.out.println("sampleInt "
		// + ((new Date()).getTime() - date.getTime()));
		return sample;
	}

	public static int[][] clusterOverlapSeqs(
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

	public static int[][] clusterOverlapSeqs(List<DoubleRead> doubleReadList,
			Set<Integer> seqList, int[][] z, double[] w, int particles,
			double alpha, int completeSeqCount) {
		int[][] tempZ = new int[seqList.size()][particles];
		int flag = 0;
		for (Integer elem : seqList) {
			tempZ[flag] = clusterOneSeq(doubleReadList, elem, z, w, particles,
					alpha, completeSeqCount);
			flag++;
		}
		return tempZ;
	}

	// May Be Something Wrong
	public static double[] updateWeights(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, int group) {

		// int group = max(z);
		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles, group);
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles, group);

		double[] p = new double[particles];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		List<Integer> countList = doubleReadList.get(seqIndex).getCountList();
		for (int k = 0; k < particles; k++) {
			double[] temp = tempTransP[k][z[seqIndex][k] - 1];
			// System.out.println(z[seqIndex][k]);
			// System.out.println(Arrays.toString(temp));
			double tempSum = 0.0;
			for (int i = 0; i < temp.length; i++) {
				if (temp[i] < Double.MIN_NORMAL) {
					tempSum += countList.get(i)
							* Math.log(DefaultConstants.ZERO);
				} else {
					tempSum += countList.get(i) * Math.log(temp[i]);
				}
			}
			p[k] = w[k] * Math.exp(tempSum);
			for (int i = 0; i < startWhere.length; i++) {
				p[k] *= tempStartP[k][z[seqIndex][k] - 1][startWhere[i]];
			}
			p[k] += Double.MIN_NORMAL;
		}
		// System.out.println(Arrays.toString(tempStartP[0][0]));
		// System.out.println(Arrays.toString(tempStartP[0][1]));
		// System.out.println(Arrays.toString(tempTransP[0][0]));
		// System.out.println(Arrays.toString(tempTransP[0][1]));
		// System.out.println(countList.toString());
		// System.out.println(Arrays.toString(p));

		double sum = 0.0;
		for (double elem : p) {
			sum += elem;
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}

	public static double[] updateWeights(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles) {
		double[][][] tempTransP = transP(doubleReadList, seqIndex, z, particles);
		double[][][] tempStartP = startP(doubleReadList, seqIndex, z, particles);

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
		tempAccumCountList = updateCountLists(doubleReadList, z,
				overlapList.get(0), tempAccumCountList);
		mapCountList = updateGroupMapList(mapCountList, z, overlapList.get(0));
		System.out.println(tempAccumCountList.get(0).get(1).toString());
		System.out.println((new Date()).toString() + "\tnewGroup start");

		newGroup(doubleReadList, params.getNeighbor());

		System.out.println((new Date()).toString() + "\tnewGroup complete");

		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();

		int accumSeqCount = overlapList.get(0).size();

		for (int i = 1; i < overlapList.size(); i++) {
			System.out.println((new Date()).toString() + "\tOverlap Group: "
					+ i);
			System.out.println(overlapList.get(i).toString());
			if (overlapList.get(i).size() > 0) {
				// System.out.println(Arrays.toString(transPCounts(tempAccumCountList.get(0).get(1))));
				if (tempAccumCountList.get(0).size() > 1) {
					System.out.println(Arrays
							.toString(startPCounts(tempAccumCountList.get(0)
									.get(2))));
				}
				// update
				// int[][] tempZLower = clusterOverlapSeqs(mapCountList,
				// tempAccumCountList, doubleReadList, overlapList.get(i),
				// z, w, params.getParticles(), params.getAlphaLow(),
				// accumSeqCount);
				// int[][] tempZUpper = clusterOverlapSeqs(mapCountList,
				// tempAccumCountList, doubleReadList, overlapList.get(i),
				// z, w, params.getParticles(), params.getAlphaHigh(),
				// accumSeqCount);
				// origin
				int[][] tempZLower = clusterOverlapSeqs(doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaLow(), accumSeqCount);
				int[][] tempZUpper = clusterOverlapSeqs(doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaHigh(), accumSeqCount);

				boolean isMarjorityLower = checkMajorityVote(tempZLower,
						params.getMajority());
				boolean isMarjorityUpper = checkMajorityVote(tempZUpper,
						params.getMajority());

				int[] voteLower = getMode(getMode(tempZLower));
				int[] voteUpper = getMode(getMode(tempZUpper));

				int majorityVoteLower = voteLower[0];
				int majorityVoteUpper = voteUpper[0];
				// System.out.println(isMarjorityLower && isMarjorityUpper
				// && (majorityVoteLower == majorityVoteUpper));
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
					// System.out.println(Arrays.toString(getMode(z)));

					accumSeqCount += overlapList.get(i).size();
					System.out.println("Group :" + majorityVoteLower);

				} else {
					rmList.add(overlapList.get(i));
					System.out.println("Removed");
				}
				// System.out.println();
			} else {
				for (Integer elem : overlapList.get(i)) {
					idTagMap.put(elem, idTagMap.get(elem) + "Single");
					// update
					// z[elem] = clusterOneSeq(mapCountList, tempAccumCountList,
					// doubleReadList, elem, z, w, params.getParticles(),
					// params.getAlpha(), accumSeqCount);
					// w = updateWeights(tempAccumCountList, doubleReadList,
					// elem,
					// z, w, params.getParticles());
					// origin
					z[elem] = clusterOneSeq(doubleReadList, elem, z, w,
							params.getParticles(), params.getAlpha(),
							accumSeqCount);
					w = updateWeights(doubleReadList, elem, z, w,
							params.getParticles());
					System.out.println("Single");
					System.out.println(Arrays.toString(z[elem]));
					// System.out.println(Arrays.toString(w));
					double eff = 0.0;
					for (double d : w) {
						eff += d * d;
					}
					// System.out.println(eff);
					// System.out.println(Arrays.toString(w));
					// System.out.println(eff);
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
						// System.out.println(Arrays.toString(resampleIndex));
						System.out.println("Resample");
						System.out.println(Arrays.toString(z[elem]));
					}

					tempAccumCountList = updateCountLists(doubleReadList, z,
							overlapList.get(i), tempAccumCountList);
					mapCountList = updateGroupMapList(mapCountList, z,
							overlapList.get(i));
					// System.out.println();
					// System.out.println(Arrays.toString(z[elem]));
					// System.out.println(Arrays.toString(w));
				}

				accumSeqCount += overlapList.get(i).size();
			}

		}

		System.out.println((new Date()).toString()
				+ "\tbegin to process removed groups");

		for (Set<Integer> seqList : rmList) {
			// update
			// int[][] tempZ = clusterOverlapSeqs(mapCountList,
			// tempAccumCountList, doubleReadList, seqList, z, w,
			// params.getParticles(), params.getAlpha(), accumSeqCount);
			// origin
			int[][] tempZ = clusterOverlapSeqs(doubleReadList, seqList, z, w,
					params.getParticles(), params.getAlpha(), accumSeqCount);
			int majorityVote = getMode(getMode(tempZ))[0];
			for (Integer seqId : seqList) {
				idTagMap.put(seqId, idTagMap.get(seqId) + "Remove");
				for (int k = 0; k < params.getParticles(); k++) {
					z[seqId][k] = majorityVote;
				}
			}
			accumSeqCount += seqList.size();
			System.out.println(seqList.toString());
			System.out.println(majorityVote);

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
}
