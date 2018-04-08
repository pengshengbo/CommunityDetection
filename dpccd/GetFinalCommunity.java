/**
 * 
 */
package edu.fzu.bigdatalab.algorithms.communitydiscovery.dpccd;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * @author psb
 * @create 2018年3月22日
 * 
 */
public class GetFinalCommunity {
	// public static Map<Integer, Set<Integer>> communityPast = new
	// HashMap<Integer, Set<Integer>>(); // 记录t-1时刻社区编号和属于该社区的节点
	// public static Map<Integer, Set<Integer>> communityNow = new
	// HashMap<Integer, Set<Integer>>(); // 记录t时刻社区编号和属于该社区的节点
	/** 记录节点和节点在t-1时刻所属社区 */
	// public static Map<Integer, Set<Integer>> pointInCommunityPast = new
	// HashMap<>();

	public static void writeCommunity(Map<Integer, Set<Integer>> resultCommunity, int time, double mu,
			double disThreshold, int densityThreshold, long duration) {
		// File communityFile = new File("results\\dbcd\\real\\as-733\\t_" +
		// time + "_" + u + "_" + minpt + ".txt");
		// File communityFile = new File("results\\dbcd\\real\\enron\\t_" + time
		// + "_" + u + "_" + minpt + ".txt");
		File communityFile = new File("results\\dpccd\\artificial\\1k\\dpccd_cm1k_" + mu + "_t" + time + "_"
				+ disThreshold + "_" + densityThreshold + ".txt");
		if (communityFile.exists()) {
			communityFile.delete();
		} else {
			try {
				communityFile.createNewFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		try {
			// FileWriter writerCommunityFile = new FileWriter(
			// "results\\dbcd\\real\\as-733\\t_" + time + "_" + u + "_" + minpt
			// + ".txt", false);
			// FileWriter writerCommunityFile = new FileWriter(
			// "results\\dbcd\\real\\enron\\t_" + time + "_" + u + "_" + minpt +
			// ".txt", false);
			FileWriter writerCommunityFile = new FileWriter("results\\dpccd\\artificial\\1k\\dpccd_cm1k_" + mu + "_t"
					+ time + "_" + disThreshold + "_" + densityThreshold + ".txt", false);
			for (Integer csnum : resultCommunity.keySet()) {
				Set<Integer> communityNodes = resultCommunity.get(csnum);
				if (communityNodes.isEmpty()) {
					continue;
				}
				Set<Integer> communitySet = new HashSet<Integer>();
				for (Iterator<Integer> it = communityNodes.iterator(); it.hasNext();) {
					Integer p = it.next();
					communitySet.add(p);
					writerCommunityFile.write(p + " ");// 写入一个数+一个空格
				}
				writerCommunityFile.write("\r\n");// 写入一个数+一个换行
			}
			writerCommunityFile.close();
			System.out.println("社区数目=" + resultCommunity.size() + ",disThreshold =" + disThreshold
					+ ",densityThreshold=" + densityThreshold + ",time=" + duration);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
