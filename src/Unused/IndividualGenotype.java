package Unused;

import java.util.Arrays;

//A sequence of unphased SNP data, with the distance between each SNP
public class IndividualGenotype {
	
	private final String[] genotypes;	//List of strings corresponding to genotype at that snp, e.g. "AA", "AT". Each string is clearly ordered here, but everyone that uses it will unorder it
	private final double[] distances;	//List of distances in terms of expected number of recombinations between each SNP, +infty for completely unlinked
	
	public IndividualGenotype(String[] snps, double[] distances){
		//Have some assertions to check that these are nice
		assert(distances.length == snps.length - 1);		//One fewer breakpoint than SNP
		for(String snp : snps){
			assert snp.length() == 2;
			assert snp.charAt(0) == 'A' || snp.charAt(0) == 'C' || snp.charAt(0) == 'G' || snp.charAt(0) == 'T';
			assert snp.charAt(1) == 'A' || snp.charAt(1) == 'C' || snp.charAt(1) == 'G' || snp.charAt(1) == 'T';
		}
		for(double distance : distances){
			assert (distance >= 0);
		}
		
		this.genotypes = Arrays.copyOf(snps, snps.length);
		this.distances = Arrays.copyOf(distances, distances.length);
	}
	
	public int getNumSNPs(){
		return this.genotypes.length;
	}
	
	public String getSNP(int index){
		return genotypes[index];
	}
	
	public double getDistance(int index){
		return distances[index];
	}
	
	public double getDistance (int startPosIdx, int endPosIdx){
		
		double sum = 0d;
		for (int i=startPosIdx; i<endPosIdx; i++){
			sum += distances[i];
		}
		
		return sum;
		
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(distances);
		result = prime * result + Arrays.hashCode(genotypes);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		IndividualGenotype other = (IndividualGenotype) obj;
		if (!Arrays.equals(distances, other.distances))
			return false;
		if (!Arrays.equals(genotypes, other.genotypes))
			return false;
		return true;
	}
	
	
	
}
