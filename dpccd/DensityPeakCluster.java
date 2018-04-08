package edu.fzu.bigdatalab.algorithms.communitydiscovery.dpccd;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import edu.fzu.bigdatalab.algorithms.communitydiscovery.exception.GraphException;
import edu.fzu.bigdatalab.algorithms.communitydiscovery.graph.Graph;
import edu.fzu.bigdatalab.algorithms.communitydiscovery.graph.GraphTools;

public class DensityPeakCluster {

	public static Graph graph; // 声明一个图
	/** 图中节点集合 */
	public static Set<Integer> vertexSet;
	/** 节点的邻接链表 */
	public static Map<Integer, Set<Integer>> adjList;
	/** 局部密度Map ：<index,densitycount> */
	private static HashMap<Integer, Integer> densityCountMap;
	/** 由大到小排序的Density list */
	private static ArrayList<Map.Entry<Integer, Integer>> sortedDensityList;
	/** deltaMap:<index, delta> */
	private static HashMap<Integer, Double> deltaMap;
	/** 每个样本的最近邻：<sampleIndex, nearestNeighborIndex> */
	private static HashMap<Integer, Integer> nearestNeighborMap;
	/** 样本对距离：<"index1 index2", distance> */
	private static HashMap<String, Double> pairDistanceMap;
	/** 最大样本距离 */
	private static double maxDistance;
	/** 最小样本距离 */
	private static double minDistance;
	/** 选取的簇中心 */
	private static ArrayList<Integer> centerList;
	/** 划分的聚类结果<sampleIndex, clusterIndex> */
	private static HashMap<Integer, Integer> clusterMap;
	/** 最终的社区划分结果，key=社区编号，value=社区节点集合 */
	private static Map<Integer, Set<Integer>> resultCommunityMap;

	/**
	 * 根据节点的局部密度和比节点密度更大的最小距离进行聚类，二者分别表示决策图中的横纵坐标
	 * 
	 * @param deltaThreshold
	 *            比节点密度大且距离节点最近的距离
	 * @param rhoThreshold
	 *            节点的密度
	 */
	public static void clustering(double deltaThreshold, double rhoThreshold) {
		centerList = new ArrayList<Integer>();
		clusterMap = new HashMap<Integer, Integer>();
		// get centers
		for (Map.Entry<Integer, Double> deltaEntry : deltaMap.entrySet()) {
			if (deltaEntry.getValue() >= deltaThreshold && densityCountMap.get(deltaEntry.getKey()) >= rhoThreshold) {
				centerList.add(deltaEntry.getKey());
				clusterMap.put(deltaEntry.getKey(), deltaEntry.getKey());
			}
		}
		// calculate clusters，注意：一定要按照密度由大到小逐个划分簇（从高局部密度到低局部密度）
		for (Map.Entry<Integer, Integer> candidate : sortedDensityList) {
			if (!centerList.contains(candidate.getKey())) {
				// 将最近邻居的类别索引作为该样本的类别索引
				if (clusterMap.containsKey(nearestNeighborMap.get(candidate.getKey()))) {
					clusterMap.put(candidate.getKey(), clusterMap.get(nearestNeighborMap.get(candidate.getKey())));
				} else {
					clusterMap.put(candidate.getKey(), -1);
				}
			}
		}
	}

	public static void calDelta() {
		// 局部密度由大到小排序
		sortedDensityList = new ArrayList<Map.Entry<Integer, Integer>>(densityCountMap.entrySet());
		Collections.sort(sortedDensityList, new Comparator<Map.Entry<Integer, Integer>>() {

			@Override
			public int compare(Entry<Integer, Integer> o1, Entry<Integer, Integer> o2) {
				if (o1.getValue() > o2.getValue())
					return -1;
				else if (o1.getValue() < o2.getValue()) {
					return 1;
				}
				return 0;
			}
		});
		nearestNeighborMap = new HashMap<Integer, Integer>(vertexSet.size());
		deltaMap = new HashMap<Integer, Double>(vertexSet.size());
		for (int i = 0; i < sortedDensityList.size(); i++) {
			if (i == 0) {
				nearestNeighborMap.put(sortedDensityList.get(i).getKey(), -1);
				deltaMap.put(sortedDensityList.get(i).getKey(), maxDistance);
			} else {
				double minDij = Double.MAX_VALUE;
				int index = 0;
				for (int j = 0; j < i; j++) {
					double dis = getDistanceFromIndex(sortedDensityList.get(i).getKey(),
							sortedDensityList.get(j).getKey());
					if (dis < minDij) {
						index = j;
						minDij = dis;
					}
				}
				nearestNeighborMap.put(sortedDensityList.get(i).getKey(), sortedDensityList.get(index).getKey());
				deltaMap.put(sortedDensityList.get(i).getKey(), minDij);
			}
		}

		// 输出样本索引+样本局部密度+最近邻索引+delta值
		// System.out.println("输出样本索引 样本局部密度 最近邻索引 delta值");
		// for (Map.Entry<Integer, Integer> entry : sortedDensityList) {
		// System.out.println(entry.getKey() + " " + entry.getValue() + " " +
		// nearestNeighborMap.get(entry.getKey())
		// + " " + deltaMap.get(entry.getKey()));
		// }
	}

	/**
	 * 根据索引获得两个样本间距离
	 * 
	 * @param index1
	 * @param index2
	 * @return
	 */
	private static double getDistanceFromIndex(int index1, int index2) {
		if (pairDistanceMap.containsKey(index1 + " " + index2)) {
			return pairDistanceMap.get(index1 + " " + index2);
		} else if (pairDistanceMap.containsKey(index2 + " " + index1)) {
			return pairDistanceMap.get(index2 + " " + index1);
		} else
			return 1;
	}

	/**
	 * 计算局部密度
	 */
	public static void calRho(double dcThreshold) {

		densityCountMap = new HashMap<Integer, Integer>(vertexSet.size());
		// 初始化为0
		// for (int i = 1; i <= vertexSet.size(); i++) {
		// densityCountMap.put(i, 0);
		// }
		for (Integer node : vertexSet) {
			densityCountMap.put(node, 0);
		}
		for (Map.Entry<String, Double> diss : pairDistanceMap.entrySet()) {
			if (diss.getValue() < dcThreshold) {
				String[] segs = diss.getKey().split(" ");
				int[] indexs = new int[2];
				indexs[0] = Integer.parseInt(segs[0]);
				indexs[1] = Integer.parseInt(segs[1]);
				for (int i = 0; i < indexs.length; i++) {
					densityCountMap.put(indexs[i], densityCountMap.get(indexs[i]) + 1);
				}
			}
		}
	}

	/**
	 * 计算截断距离
	 * 
	 * @return
	 */
	public static double findDC() {
		double tmpMax = maxDistance;
		double tmpMin = minDistance;
		double dc = 0.5 * (tmpMax + tmpMin);
		for (int iteration = 0; iteration < 100; iteration++) {
			int neighbourNum = 0;
			for (Map.Entry<String, Double> dis : pairDistanceMap.entrySet()) {
				if (dis.getValue() < dc)
					neighbourNum += 2;
			}
			double neighborPercentage = neighbourNum / Math.pow(vertexSet.size(), 2);
			if (neighborPercentage >= 0.01 && neighborPercentage <= 0.02)
				break;
			if (neighborPercentage > 0.02) {
				tmpMax = dc;
				dc = 0.5 * (tmpMax + tmpMin);
			}
			if (neighborPercentage < 0.01) {
				tmpMin = dc;
				dc = 0.5 * (tmpMax + tmpMin);
			}

		}
		return dc;
	}

	/**
	 * 计算两个样本的高斯距离
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	private double twoSampleDistance(Sample a, Sample b) {
		double[] aData = a.getAttributes();
		double[] bData = b.getAttributes();
		double distance = 0.0;
		for (int i = 0; i < aData.length; i++) {
			distance += Math.pow(aData[i] - bData[i], 2);
		}
		return 1 - Math.exp(distance * (-0.5));
	}

	public static ArrayList<Integer> getCenterList() {
		return centerList;
	}

	private static void calculateDistance(Graph graph) throws GraphException {
		pairDistanceMap = new HashMap<>();
		maxDistance = Double.MIN_VALUE;
		minDistance = Double.MAX_VALUE;
		Map<Integer, Set<Integer>> adjList = graph.getAdjLists();
		double distance = 0;
		for (Integer point : vertexSet) {
			for (Integer adj : adjList.get(point)) {
				if (!pairDistanceMap.containsKey(point + " " + adj)
						&& !pairDistanceMap.containsKey(adj + " " + point)) {
					distance = 1 - StructSimilarity.calculateStructSimilarity(point, adj, graph);
					pairDistanceMap.put(point + " " + adj, distance);
					pairDistanceMap.put(adj + " " + point, distance);
					if (distance > maxDistance)
						maxDistance = distance;
					if (distance < minDistance)
						minDistance = distance;
				}
				for (Integer indrectAdj : adjList.get(adj)) {
					if (!pairDistanceMap.containsKey(point + " " + indrectAdj)
							&& !pairDistanceMap.containsKey(indrectAdj + " " + point)) {
						distance = 1 - StructSimilarity.calculateStructSimilarity(point, indrectAdj, graph);
						pairDistanceMap.put(point + " " + indrectAdj, distance);
						pairDistanceMap.put(indrectAdj + " " + point, distance);
						if (distance > maxDistance)
							maxDistance = distance;
						if (distance < minDistance)
							minDistance = distance;
					}
				}
			}
		}

	}

	/**
	 * 根据节点所属的社区序号把节点划分到对应的社区中
	 */
	private static void partionCommunity() {
		resultCommunityMap = new HashMap<>();
		for (Integer point : clusterMap.keySet()) {
			Integer csNumber = clusterMap.get(point);
			if (!resultCommunityMap.containsKey(csNumber)) {
				Set<Integer> communityNodeSet = new HashSet<>();
				communityNodeSet.add(point);
				resultCommunityMap.put(csNumber, communityNodeSet);
			} else {
				resultCommunityMap.get(csNumber).add(point);
			}
		}

	}

	public static void main(String[] args) throws GraphException, IOException {
		int n = 8;
		int disThreshold = 60;
		int densityThreshold = 9;
		for (disThreshold = 10; disThreshold <= 96; disThreshold += 2) {
			for (densityThreshold = 5; densityThreshold <= 50; densityThreshold += 1) {
				for (double mu = 0.8; mu <= 0.8; mu = mu + 0.1) {
					for (int time = 2; time <= n; time++) {
						graph = GraphTools.loadGraphFromFile("data\\real\\as-733\\as_t" + time + ".txt");
						long start = System.currentTimeMillis();
						/** 人工数据集 */
						// graph =
						// GraphTools.loadGraphFromFile("data\\artificial\\1k_u\\network1k_"
						// + mu + "_" + time + ".txt");
						vertexSet = graph.getVertices();
						adjList = graph.getAdjLists();

						/** 安然数据集 */
						// graph =
						// GraphTools.loadGraphFromFile("data\\real\\enrondata\\enron_t"
						// + time + ".txt");
						calculateDistance(graph);
						double dc = findDC();
						System.out.println("dc=" + dc);
						// System.out.println(dc);
						calRho(dc);
						calDelta();
						clustering(disThreshold * 1.0 / 100, densityThreshold);
						// System.out.println(getCenterList());
						partionCommunity();
						long end = System.currentTimeMillis();
						GetFinalCommunity.writeCommunity(resultCommunityMap, time, mu, disThreshold, densityThreshold,
								(end - start));

					}
				}
			}
		}

	}

}
