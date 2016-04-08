package dataStructures;

import java.util.Arrays;


//Contains information about relationship
public class Relationship {
	
	private final double meiosisRate; //alpha
	private final double[] marginalProbs; //omega
	private final Path[] paths; // (up,down,numVisit) e.g.(0,0,0) = unrelated

	
	public Relationship(double meiosisRate, double[] transitionRates, Path[] paths){
		assert meiosisRate > 0;
		
		double sum = 0;
		for (double transitionRate : transitionRates){
			assert transitionRate >= 0 && transitionRate <= 1;
			sum += transitionRate;
		}
		assert Math.abs(sum - 1d) < 1e-13;
		this.marginalProbs = Arrays.copyOf(transitionRates, transitionRates.length);
		this.paths = paths;

		
		this.meiosisRate = meiosisRate;
		
	}
	
	
	public Path getPath(){
		return paths[0];
	}
	
	

	public double getMeiosisRate(){
		return meiosisRate;
	}
	
	public double[] getMarginalProbs(){
		return marginalProbs;
	}
	
	public Path[] getAllPaths(){
		return paths;
	}
}
