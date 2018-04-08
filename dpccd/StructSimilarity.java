/**
 * 
 */
package edu.fzu.bigdatalab.algorithms.communitydiscovery.dpccd;

import java.util.HashSet;
import java.util.Set;

import edu.fzu.bigdatalab.algorithms.communitydiscovery.exception.GraphException;
import edu.fzu.bigdatalab.algorithms.communitydiscovery.graph.Graph;

/**
 * @author psb
 * @create 2018年3月8日
 * 
 */
public class StructSimilarity {

	public static double calculateStructSimilarity(Integer p, Integer q, Graph graph) throws GraphException {

		Set<Integer> setP = new HashSet<>(graph.getNeighbors(p));
		setP.add(p);
		Set<Integer> setQ = new HashSet<>(graph.getNeighbors(q));
		setQ.add(q);
		Set<Integer> intersection = new HashSet<Integer>();
		intersection.addAll(setP);
		/** 求两集合的交集 */
		intersection.retainAll(setQ);
		return (double) (Math.round(intersection.size() * 1.0 / Math.sqrt(setP.size() * setQ.size() * 1.0) * 10000)
				/ 10000.0);
		// return (double) (Math.round((1 - intersection.size() * 1.0 /
		// Math.sqrt(setP.size() * setQ.size() * 1.0)) * 10000) / 10000.0);

	}

}
