package Unused;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.PoissonDistribution;

import dataStructures.Pedigree;
import dataStructures.SimplePair;
import utility.ArrayUtility;
import utility.LogSum;

//TODO: this class is shit and needs to be fixed if I'm going to use it for anything other than preliminary results for Rasmus's grant
public class LanderGreenLikelihood {
	
	private final double[][] stationaryDistributions;
	private final double genotypingErrorRate;
	
	public LanderGreenLikelihood(double[][] stationaryDistributions, double genotypingErrorRate){
		this.stationaryDistributions = new double[stationaryDistributions.length][4];
		for(int loc = 0; loc < stationaryDistributions.length; loc++){
			assert Math.abs(ArrayUtility.sum(stationaryDistributions[loc]) - 1d) < ArrayUtility.EPSILON;
			assert stationaryDistributions[loc].length == 4;
			this.stationaryDistributions[loc] = Arrays.copyOf(stationaryDistributions[loc], 4);
		}
		 
		
		this.genotypingErrorRate = genotypingErrorRate;
	}
	
	//I thought we could just do sum-product over individuals, but you can't because that requires
	//knowing which alleles are on each haplotype in your parents (to determine if a recombination is
	//required in the transition matrix).
	//TODO: incorportae genotyping error
	public double logLikelihood(Pedigree pedigree){
		//get founders and non-founders
		List<Individual> founders = new ArrayList<Individual>(pedigree.getFounders());	//TODO: this is wrong (consider a pedigree of only siblings. their parents,
																						//are founders, but aren't indivuals, for one, and for two aren't even returned
																						//by the current method
		
		Set<PedigreeNode> nonfounderSet = new HashSet<PedigreeNode>();
		for(Individual founder : founders){
			nonfounderSet.addAll(founder.getDescendants());
		}
		for(Individual founder : founders){
			assert nonfounderSet.contains(founder);
		}
		List<PedigreeNode> nonfounders = new ArrayList<PedigreeNode>(nonfounderSet);
		//make possible inheritance vectors
		double tester = 0d;
		for(int loc = 0; loc < founders.get(0).getGenotype().getNumSNPs(); loc++){
			tester += oneIVlogEmission(loc, new ArrayList<List<SimplePair<Boolean,Boolean>>>(), 0, nonfounders, founders);
		}
		if(nonfounders.size() == 0){
			return tester;
		}
		System.out.println("Just the founders: " + tester);
//		}
		List<List<SimplePair<Boolean, Boolean>>> inheritanceVector = getInheritanceVector(nonfounders.size());
		
//		if(nonfounders.size() == 0){
		
		
		int[][] hammingDists = getNumRecombinations(inheritanceVector);
		
		
		double[] forwardProbs = new double[inheritanceVector.get(0).size()];
		//initialize:
		//TODO: NOTE THIS IS SUPER BROKEN ON SEVERAL LEVELS, but should work
		//on pedigrees where the founders are individuals in the pedigree and there are no non-specified descendents of the founders
		for(int iv = 0; iv < inheritanceVector.get(0).size(); iv++){
			double toAdd = -Math.log(inheritanceVector.get(0).size());	//a priori all inheritance vectors are equally likely
			toAdd += oneIVlogEmission(0, inheritanceVector, iv, nonfounders, founders);
			forwardProbs[iv] = toAdd;
		}
		
		//recurse
		for(int loc = 1; loc < founders.get(0).getGenotype().getNumSNPs(); loc++){
			double[] prevForwardProbs = Arrays.copyOf(forwardProbs, forwardProbs.length);
			for(int ivCurr = 0; ivCurr < inheritanceVector.get(0).size(); ivCurr++){
				double toAdd = oneIVlogEmission(loc, inheritanceVector, ivCurr, nonfounders, founders);
				forwardProbs[ivCurr] = Double.NEGATIVE_INFINITY;
				for(int ivPrev = 0; ivPrev < inheritanceVector.get(0).size(); ivPrev++){
					double trans = getLogTransitionProb(hammingDists[ivCurr][ivPrev], inheritanceVector.size() * 2, founders.get(0).getGenotype().getDistance(loc - 1));
					forwardProbs[ivCurr] = LogSum.addLogSummand(forwardProbs[ivCurr], toAdd + trans + prevForwardProbs[ivPrev]);
				}
			}
		}
		
		double toReturn = Double.NEGATIVE_INFINITY;
		for(double toAdd : forwardProbs){
			toReturn = LogSum.addLogSummand(toReturn, toAdd);
		}
		
		return toReturn;
	}
	
	//independent poissons, checking if even.
	private double getLogTransitionProb(int numRecos, int totalN, double geneticDist) {
		double oneP = 1d - Math.exp(-2 * geneticDist);
		oneP = oneP / 2d;
		return numRecos*Math.log(oneP) + (totalN - numRecos) * Math.log(1-oneP);
	}

	private double oneIVlogEmission(int locus, List<List<SimplePair<Boolean, Boolean>>> inheritanceVector, int iv, List<PedigreeNode> nonfounders, List<Individual> founders){
		double toAdd = 0d;
		for(Individual founder : founders){
			String founderString = founder.getGenotype().getSNP(locus);
			toAdd += Math.log(getAlleleProb(locus, founderString.charAt(0))) + Math.log(getAlleleProb(locus, founderString.charAt(1)));
			toAdd += founderString.charAt(0) == founderString.charAt(1) ? 0d : Math.log(2d);
		}
		for(int indIdx = 0; indIdx < inheritanceVector.size(); indIdx++){
			assert nonfounders.get(indIdx) instanceof Individual;
			Individual ind = (Individual) nonfounders.get(indIdx);
			String indString = ind.getGenotype().getSNP(locus);
			String inheritedString = getInheritedString(ind, inheritanceVector.get(indIdx).get(iv), locus);
			if(getIBS(indString, inheritedString) != 2){
				toAdd += Double.NEGATIVE_INFINITY;
				break;
			}
		}
		return toAdd;
	}
	
	private String getInheritedString(Individual ind, SimplePair<Boolean,Boolean> iv, int loc){
		String toReturn = "";
		int firstParIdx = iv.getFirst() ?  1 : 0;
		int secondParIdx = iv.getSecond() ?  1 : 0;
		toReturn += ((Individual) ind.getParents().get(0)).getGenotype().getSNP(loc).charAt(firstParIdx);
		toReturn += ((Individual) ind.getParents().get(1)).getGenotype().getSNP(loc).charAt(secondParIdx);
		return toReturn;
	}
	
	//This is huge, and why Lander Green is sooo slow O(2^{2n})
	//[individual][possible
	private List<List<SimplePair<Boolean, Boolean>>> getInheritanceVector(int numIndividuals){
		if(numIndividuals == 0){
			return new ArrayList<List<SimplePair<Boolean,Boolean>>>();	//empty!
		}
		if(numIndividuals == 1){
			List<List<SimplePair<Boolean, Boolean>>> toReturn = new ArrayList<List<SimplePair<Boolean,Boolean>>>();
			toReturn.add(new ArrayList<SimplePair<Boolean, Boolean>>());
			toReturn.get(0).add(new SimplePair<Boolean,Boolean>(false, false));
			toReturn.get(0).add(new SimplePair<Boolean,Boolean>(false, true));
			toReturn.get(0).add(new SimplePair<Boolean,Boolean>(true, false));
			toReturn.get(0).add(new SimplePair<Boolean,Boolean>(true, true));
			return toReturn;
		}
		List<List<SimplePair<Boolean, Boolean>>> toCopy = getInheritanceVector(numIndividuals - 1);
		List<List<SimplePair<Boolean, Boolean>>> toReturn = new ArrayList<List<SimplePair<Boolean,Boolean>>>();
		for(int i = 0; i < numIndividuals; i++){
			toReturn.add(new ArrayList<SimplePair<Boolean, Boolean>>());
			if(i < (numIndividuals - 1)){
				for(int rep = 0; rep < 4; rep++){
					//need to quadruple this guy
					toReturn.get(i).addAll(toCopy.get(i));
				}
			} else {
				for(int rep = 0; rep < 4; rep++){
					for(int idx = 0; idx < toCopy.get(0).size(); idx++){
						if(rep == 0){
							toReturn.get(i).add(new SimplePair<Boolean,Boolean>(false, false));
						} else if(rep == 1){
							toReturn.get(i).add(new SimplePair<Boolean,Boolean>(false, true));
						} else if(rep == 2){
							toReturn.get(i).add(new SimplePair<Boolean,Boolean>(true, false));
						} else {
							toReturn.get(i).add(new SimplePair<Boolean,Boolean>(true, true));
						}
					}
				}
			}
		}
		return toReturn;
	}
	
	private double getAlleleProb(int loc, char allele){
		if (allele == 'A') return this.stationaryDistributions[loc][0];
		if (allele == 'C') return this.stationaryDistributions[loc][1]	;
		if (allele == 'G') return this.stationaryDistributions[loc][2];
		if (allele == 'T') return this.stationaryDistributions[loc][3];
		throw new RuntimeException("Tried to access illegal allele");
	}
	
	private int getIBS(String genotype1, String genotype2){
		if(genotype1.charAt(0) == genotype2.charAt(0)){
			if(genotype1.charAt(1) == genotype2.charAt(1)){
				return 2;	//strings are the same
			}
			return 1;	//first character are the same, second characters differ
		} 

		if(genotype1.charAt(0) == genotype2.charAt(1)){
			if(genotype1.charAt(1) == genotype2.charAt(0)){
				return 2;	//strings are the same but reversed
			}
			return 1;	//first character in one is last character in the other
		}
		//first character doesn't equal either character in the other
		//if the second character equals either of the characters in the other string, ibs = 1 (since first character doesn't match), else 0
		return genotype1.charAt(1) == genotype2.charAt(0) || genotype1.charAt(1) == genotype2.charAt(1) ? 1 : 0;

	}
	
	//glorified hamming distance calculator (note this step might not be necessary because if we were smart
	//we would know all of these values as we constructed the initial table.)
	private int[][] getNumRecombinations(List<List<SimplePair<Boolean, Boolean>>> inheritanceVector){
		int[][] toReturn = new int[inheritanceVector.get(0).size()][inheritanceVector.get(0).size()];
		for(int i = 0; i < inheritanceVector.get(0).size(); i++){
			for(int j = i + 1; j < inheritanceVector.get(0).size(); j++){
				int numRecos = 0;
				for(List<SimplePair<Boolean,Boolean>> individual : inheritanceVector){
					if(individual.get(i).getFirst() != individual.get(j).getFirst()){
						numRecos += 1;
					}
					if(individual.get(i).getSecond() != individual.get(j).getSecond()){
						numRecos += 1;
					}
				}
				assert numRecos > 0;	//we implicitly have avoided calculating any of these where we know they're zero (the diagonals)
				toReturn[i][j] = numRecos;
				toReturn[j][i] = numRecos;	//symmetry!
			}
		}

		
		return toReturn;
	}
}
