package utility;

import java.util.ArrayList;
import java.util.List;

public class LogSum {
	
	public static double addLogSummand(double logA, double logB){
		if(Double.NEGATIVE_INFINITY == logA && Double.NEGATIVE_INFINITY == logB){
			return Double.NEGATIVE_INFINITY;
		}
		double maxVal = Math.max(logA, logB);
		return maxVal + Math.log(Math.exp(logA - maxVal) + Math.exp(logB - maxVal));
	}

	private static List<Double> lfacts = new ArrayList<Double>();
	public static double lfactorial(int k) {
		if(k == 0){
			return 0d;
		}
		if(lfacts.size() >= k){
			return lfacts.get(k-1);
		} else if(k == 1){
			lfacts.add(0d);
		}
		else{
			//recursively build this guy
			lfacts.add(Math.log(k) + lfactorial(k-1));
		}
		return(lfactorial(k));
	}
	
}
