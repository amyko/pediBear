package Unused;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import utility.ArrayUtility;
import utility.DataParser;
import dataStructures.Genotype;
import dataStructures.Node;
import dataStructures.EmissionKey;
import dataStructures.Path;
import dataStructures.Relationship;
import dataStructures.SimplePair;


//This Class uses the Forward Algorithm to compute the likelihood of a pair of genotypes (SNPs and distances in terms of expected number of recombinations)
//given a relationship and a sequencing error rate
//hidden states=ibd between two individuals, observed states=ibs between two individuals at each site
//genotype unordered, individual ordered



public class PairwiseLikelihoodCoreStream {
	
	
	//machine precision
	private static final double TOL = 1e-12;

	//info file columns
	private static final int POS = 0;
	private static final int LDPOS = 1;
	private static final int U = 2;
	private static final int u = 3;
	private static final int V = 4;
	private static final int v = 5;
	private static final int pA = 6;
	private static final int pC = 7;
	private static final int pG = 8;
	private static final int pT = 9;
	private static final int A = 10;
	private static final int B = 11;
	private static final int C = 12;
	private static final int D = 13;
	
	//pairwise likelihood file columns
	private static final int UP = 1;
	private static final int DOWN = 2;
	private static final int NUMVISIT = 3;
	
	// model parameters and data structures
	private int back;
	private final double genotypingErrorRate;
	private final double recombRate;
	private int numIndiv;
	private final Map<Path, double[][]> pairwiseLkhd = new HashMap<Path, double[][]>(); //path --> pairwise likelihood table
	private final List<Double> marginals = new ArrayList<Double>(); // index --> marginal likelihood

	
	public PairwiseLikelihoodCoreStream(double genotypingErrorRate, double recombRate, int back, int numIndiv){
		this.genotypingErrorRate = genotypingErrorRate;
		this.recombRate = recombRate;
		this.back = back;
		this.numIndiv = numIndiv;
	}
	
	////////////////// GETTER AND SETTER FOR LIKELIHOOD /////////////////////////// WORKS
	
	public double getLikelihood(Node ind1, Node ind2, Path path){
		
		if (!pairwiseLkhd.containsKey(path)){
			System.out.println(String.format("%d\t%d\t%s", ind1.getIndex(), ind2.getIndex(), path.toString()));
			throw new RuntimeException("key not available");
		}
		if (pairwiseLkhd.get(path)==null){
			System.out.println(String.format(path.toString()));
			throw new RuntimeException("Likelihood not computed for this key");
		}
		
		return pairwiseLkhd.get(path)[ind1.getIndex()][ind2.getIndex()];
	}
	

	public double getMarginal(Node ind){ 
		
		int key = ind.getIndex();
		
		if (key < 0 || key >= marginals.size()){
			System.out.println(key);
			throw new RuntimeException("key not available");
		}
		
		return marginals.get(ind.getIndex());

	}
	
	
	public void setLikelihoods(String filePath) throws IOException{ //works
	
		//clear old values
		pairwiseLkhd.clear();
		
		//open files
		BufferedReader reader = DataParser.openReader(filePath);
		
		String line;
		Path path = null;
		int i=-1;
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//update relationship and individual index
			if(fields[0].equals(">")){

				i = 0;				
				double[][] lkhdTable = new double[numIndiv][numIndiv];
				path = new Path(Integer.parseInt(fields[UP]), Integer.parseInt(fields[DOWN]), Integer.parseInt(fields[NUMVISIT]));
				pairwiseLkhd.put(path, lkhdTable);
				
				continue;
			}
			
			//read likelihood
			for(int j=i+1; j<numIndiv; j++){ 
				
				pairwiseLkhd.get(path)[i][j] = Double.parseDouble(fields[j]);
				pairwiseLkhd.get(path)[j][i] = Double.parseDouble(fields[j]);
				
			}
			
			//increment row
			i++;
		}
		
		//close file
		reader.close();
		
	}
	
	
	public void setMarginals(String filePath) throws IOException{ //works
		
		//clear old values
		marginals.clear();
		
		//open file
		BufferedReader reader = DataParser.openReader(filePath);
		//skip header
		reader.readLine();
		
		//read marginals
		String[] fields = reader.readLine().split("\t");
		for(int i=0; i<numIndiv; i++) 
			marginals.add(i, Double.parseDouble(fields[i]));
		
		//close file
		reader.close();
		
	}

	
	
	

	//////////////////////////// PRECOMPUTE LIKELIHOOD ///////////////////////////////////	
	public double[] computeMarginal(String genoPath, String infoPath, int[] indCols) throws IOException{ //works
		
		int numIndiv = indCols.length;
		double[] toReturn = new double[numIndiv];
		int numSnp = DataParser.countLines(infoPath)-1;
		
		// data structures
		Map<Integer,String[]> posToGeno = new HashMap<Integer,String[]>(back); //holds genotypes for the last 50 snps
		Map<Integer, String[]> posToInfo = new HashMap<Integer, String[]>(back); //holds info line for the last 50 snps
		
		//open files
		BufferedReader genoFile = DataParser.openReader(genoPath);
		BufferedReader infoFile = DataParser.openReader(infoPath);
		//skip headers
		genoFile.readLine();
		infoFile.readLine();
		
		
		//First snp
		String[] geno = genoFile.readLine().split("\t");
		String[] info = infoFile.readLine().split("\t");
		
		Map<Genotype, Double> oneLocusGenoProbMap = computePossibleGenotypeProbWithError(info);
		
		for (int i=0; i<numIndiv; i++){
			String g = geno[indCols[i]];
			toReturn[i] += Math.log(oneLocusGenoProbMap.get(new Genotype(g)));
		}
		
		//update prev info
		int prevPos = Integer.parseInt(info[POS]);
		posToGeno.put(prevPos, geno);
		posToInfo.put(prevPos, info);
		
		for(int locus = 1; locus < numSnp; locus++){ //for every snp
			
			//read current
			info = infoFile.readLine().split("\t");
			geno = genoFile.readLine().split("\t");
			int ldPos = Integer.parseInt(info[LDPOS]);
			int currPos = Integer.parseInt(info[POS]);
			
			//read ld
			String[] ldGeno = posToGeno.get(ldPos);
			String[] ldInfo = posToInfo.get(ldPos);
			
			//compute all possible genotype probs
			Map<SimplePair<Genotype,Genotype>, Double> conditionalGenoProbMap = computePossibleConditionalGenotypeProbWithError(ldInfo, info);
			
			
			//compute likelihood
			for (int i=0; i<numIndiv; i++){
				int col = indCols[i];
				
				String geno1 = ldGeno[col];
				String geno2 = geno[col];
				
				//compute conditional allele probability
				double lkhd = conditionalGenoProbMap.get(new SimplePair<Genotype,Genotype>(new Genotype(geno1), new Genotype(geno2)));
				
				toReturn[i] += Math.log(lkhd);
					
			}
			
			//update prev info
			prevPos = currPos;
			posToGeno.put(prevPos, geno);
			posToInfo.put(prevPos, info);
			
			
			//remove useless info
			if(posToGeno.size()>back){
				int minKey = Collections.min(posToGeno.keySet());
				posToGeno.remove(minKey);
				posToInfo.remove(minKey);
			}
		}
		
		return toReturn;
	}
	
	
	public double[][] forwardAlgorithm(String genoPath, String infoPath, int[] indCols, Relationship rel) throws IOException{
		
		int numIndiv = indCols.length;
		double[][] toReturn = new double[numIndiv][numIndiv];
		int numSnp = DataParser.countLines(infoPath)-1;
		//int numSnp = 1000;
		
		//open file
		BufferedReader genoFile = DataParser.openReader(genoPath);
		BufferedReader infoFile = DataParser.openReader(infoPath);
		genoFile.readLine(); //skip header
		infoFile.readLine(); //skip header
		
		// data structures
		Map<Integer,String[]> posToGeno = new HashMap<Integer,String[]>(back); //holds genotypes for the last 50 snps
		Map<Integer, String[]> posToInfo = new HashMap<Integer, String[]>(back); //holds info line for the last 50 snps

		//previous forward probabilities (only 3 distinct ibd states)
		double[][][] alpha = new double[numIndiv][numIndiv][3];
		
		//initialization:
		//read file
		String[] geno = genoFile.readLine().split("\t");
		String[] info = infoFile.readLine().split("\t");
	
		Map<EmissionKey, Double> emissionMap = computePossibleOneLocusEmissionWithError(info);
		
		//initial alpha
		for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
			for(int i=0; i<numIndiv; i++){
				int i1 = indCols[i];
				for(int j=i+1; j<numIndiv; j++){
					int i2 = indCols[j];
					EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(geno[i1]), new Genotype(geno[i2])), null, ibd);
					alpha[i][j][ibd] = rel.getMarginalProbs()[ibd] * emissionMap.get(key);
				}
			}
		}

		// scale alphas
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				if(ArrayUtility.sum(alpha[i][j]) == 0d){
					toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
				}
				
				double scalingFactor = 1d / ArrayUtility.sum(alpha[i][j]);

				for(int ibd = 0; ibd < alpha[0][0].length; ibd++) alpha[i][j][ibd] = alpha[i][j][ibd] * scalingFactor;//scale alpha
				
				//add to lkhd
				toReturn[i][j] -= Math.log(scalingFactor);

			}
		}
		
	
		//update previous data
		int prevPos = Integer.parseInt(info[POS]);
		posToGeno.put(prevPos,geno);
		posToInfo.put(prevPos, info);
		

		//recursion:
		for(int locus = 1; locus < numSnp; locus++){
			if (locus%10000==0) System.out.println(locus);
			
			double[][][] tempAlpha = new double[numIndiv][numIndiv][3];
			
			//read info
			info = infoFile.readLine().split("\t");
			int ldPos = Integer.parseInt(info[LDPOS]);
			String[] ldInfo = posToInfo.get(ldPos);
			
			//read geno
			geno = genoFile.readLine().split("\t");
			String[] ldGeno = posToGeno.get(ldPos);
			
			//read relevant data
			int currPos = Integer.parseInt(info[POS]);
			double dist = (currPos - prevPos) * recombRate;

			
			//compute all possible emission probs
			Map<EmissionKey, Double> conditionalEmissionMap = computePossibleConditionalEmissionWithError(ldInfo, info);
			
			
			for (int ibdNew = 0; ibdNew < alpha[0][0].length; ibdNew++){
				//transition
				for(int ibdPrev = 0; ibdPrev < alpha[0][0].length; ibdPrev++){	
					
					for(int i=0; i<numIndiv; i++){
						for(int j=i+1; j<numIndiv; j++){
							double transitionDensity = transitionDensity(ibdNew, ibdPrev, dist, rel);
							tempAlpha[i][j][ibdNew] += alpha[i][j][ibdPrev] * transitionDensity;
						}
					}

				}
				
				
				//emission	
				for(int i=0; i<numIndiv; i++){
					int i1 = indCols[i];
					for(int j=i+1; j<numIndiv; j++){
						int i2 = indCols[j];
						EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(ldGeno[i1]), new Genotype(ldGeno[i2])), new SimplePair<Genotype,Genotype>(new Genotype(geno[i1]), new Genotype(geno[i2])), ibdNew);
						double emissionDensity = conditionalEmissionMap.get(key); //should never be zero
						//double emissionDensity = emissionDensityWithError(info, geno[i1], geno[i2], ibdNew); //should never be zero		
						tempAlpha[i][j][ibdNew] *= emissionDensity;
						
						
					}
				}
				

			}
			
			
			//scale alphas
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					if(ArrayUtility.sum(tempAlpha[i][j]) == 0d){
						toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
					}
					
					double scalingFactor = 1d / ArrayUtility.sum(tempAlpha[i][j]);

					for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
						alpha[i][j][ibd] = tempAlpha[i][j][ibd] * scalingFactor;//scale alpha
					}
					
					//add to lkhd
					toReturn[i][j] -= Math.log(scalingFactor);
					
				}
			}

			
			//add previous info
			prevPos = currPos;
			posToGeno.put(prevPos,geno);
			posToInfo.put(prevPos, info);
			

			//remove useless
			if(posToGeno.size()>back){
				int minKey = Collections.min(posToGeno.keySet());
				posToGeno.remove(minKey);
				posToInfo.remove(minKey);
			}
			
		}
		
		return toReturn;
		
	}	
	
	
	////////// EMISSION PROBABILITY /////////

	//returns map (geno1, geno2) --> emission prob with error for all possible genotype combinations
	private Map<EmissionKey, Double> computePossibleOneLocusEmissionWithError(String[] info){
		
		//read info
		String[] possibleGenotypes = getPossibleGenotypes(info[V], info[v]);
		double[] oneLocusFreq = getOneLocusFreq(info);
		
		
		//possible one-locus emission density WITHOUT error
		Map<EmissionKey, Double> oneLocusEmissionWithoutError = new HashMap<EmissionKey, Double>(27);
		for(String geno1 : possibleGenotypes){
			for(String geno2 : possibleGenotypes){
				for(int ibd=0; ibd<3; ibd++){
					
					EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(geno1), new Genotype(geno2)), null, ibd);
					Double val = emissionDensityWithoutError(oneLocusFreq, geno1, geno2, ibd);
					oneLocusEmissionWithoutError.put(key, val);
					
				}
			}
		}
		
		//precompute error probs
		Map<SimplePair<Genotype,Genotype>, Double> errorMap = computePossibleErrorProbs(info);
		
		//possible one-locus emission WITH error
		Map<EmissionKey, Double> oneLocusEmissionWithError = new HashMap<EmissionKey, Double>(27);
		for(String geno1 : possibleGenotypes){
			for(String geno2 : possibleGenotypes){
				for(int ibd=0; ibd<3; ibd++){
					
					EmissionKey key =  new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(geno1), new Genotype(geno2)), null, ibd);
					double val = emissionDensityWithError(oneLocusEmissionWithoutError, errorMap, possibleGenotypes, geno1, geno2, ibd);
					oneLocusEmissionWithError.put(key, val);
					
				}
			}
		}
		
		
		return oneLocusEmissionWithError;
		
	}
     		
	
	//no marginalizing assume genotype1 and genotype2 are the real deal; individuals are ordered
	private double emissionDensityWithoutError(double[] oneLocusFreq, String genotype1, String genotype2, int ibd){
		int ibs = getIBS(genotype1, genotype2);
		if( ibs < ibd){
			return 0d;
		}
		
		double delta1 = genotype1.charAt(0) == genotype1.charAt(1) ? 1d : 0d;
		double delta2 = genotype2.charAt(0) == genotype2.charAt(1) ? 1d : 0d;	//we need these deltas for the combinatorial coefficients in a lot of these calculations
		
		if(ibd == 2){	//correct
			return (2d - delta1) * getAlleleProb(oneLocusFreq, genotype1.charAt(0)) * getAlleleProb(oneLocusFreq, genotype1.charAt(1));
		}
		if(ibd == 1){
			//need to do some fancy stuff to deal with this case
			if(ibs == 2){	//correct
				if(delta1==1){
					return Math.pow(getAlleleProb(oneLocusFreq, genotype1.charAt(0)),3);
				}
				else{
					return Math.pow(getAlleleProb(oneLocusFreq, genotype1.charAt(0)),2) * getAlleleProb(oneLocusFreq, genotype1.charAt(1)) + Math.pow(getAlleleProb(oneLocusFreq, genotype1.charAt(1)), 2) * getAlleleProb(oneLocusFreq, genotype1.charAt(0));	
				}
			}
			if(ibs == 1){ 	//correct
				char sharedAllele;
				char loneAllele1;
				char loneAllele2;
				if(genotype1.charAt(0) == genotype2.charAt(0)){
					sharedAllele = genotype1.charAt(0);
					loneAllele1 = genotype1.charAt(1);
					loneAllele2 = genotype2.charAt(1);
				} else if(genotype1.charAt(0) == genotype2.charAt(1)){
					sharedAllele = genotype1.charAt(0);
					loneAllele1 = genotype1.charAt(1);
					loneAllele2 = genotype2.charAt(0);
				} else if(genotype1.charAt(1) == genotype2.charAt(1)){
					sharedAllele = genotype1.charAt(1);
					loneAllele1 = genotype1.charAt(0);
					loneAllele2 = genotype2.charAt(0);
				} else if(genotype1.charAt(1) == genotype2.charAt(0)){
					sharedAllele = genotype1.charAt(1);
					loneAllele1 = genotype1.charAt(0);
					loneAllele2 = genotype2.charAt(1);
				} else {
					throw new RuntimeException("Something bad happened while computing an emission prob");
				}
				return getAlleleProb(oneLocusFreq, sharedAllele) * getAlleleProb(oneLocusFreq, loneAllele1) * getAlleleProb(oneLocusFreq, loneAllele2);	//combinatorial factors should exactly cancel
			}
			throw new RuntimeException("ibd didn't match IBS!");
		}
		if(ibd == 0){	//correct
			return (2d - delta1) * (2d - delta2) * getAlleleProb(oneLocusFreq, genotype1.charAt(0)) * getAlleleProb(oneLocusFreq, genotype1.charAt(1)) * getAlleleProb(oneLocusFreq, genotype2.charAt(0)) * getAlleleProb(oneLocusFreq, genotype2.charAt(1));
		}
		throw new RuntimeException("Illegal IBD value encountered");
	}
	
	
	//Returns xi(ibs(snp1 snp2) | ibd), the emission probability of ibs given ibd assumes that snp1 and snp2 may have seq. errors and marginalizes
	private double emissionDensityWithError(Map<EmissionKey, Double> oneLocusEmissionWithoutError, Map<SimplePair<Genotype,Genotype>, Double> errorMap, String[] possibleGenotypes, String genotype1, String genotype2, int ibd){
		
		double toReturn = 0d;
		for(String ind1TrueGenotype : possibleGenotypes){
			SimplePair<Genotype,Genotype> errorKey1 = new SimplePair<Genotype,Genotype>(new Genotype(ind1TrueGenotype), new Genotype(genotype1));
			
			for(String ind2TrueGenotype: possibleGenotypes){
				SimplePair<Genotype,Genotype> errorKey2 = new SimplePair<Genotype,Genotype>(new Genotype(ind2TrueGenotype), new Genotype(genotype2));
				
				toReturn += errorMap.get(errorKey1) * errorMap.get(errorKey2) * oneLocusEmissionWithoutError.get(new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(ind1TrueGenotype), new Genotype(ind2TrueGenotype)), null, ibd));
			}
		}

		return toReturn;
	}
	
	
	
	//compute two-locus emission prob for all possible genotype combinations
	private Map<EmissionKey, Double> computePossibleConditionalEmissionWithError(String[] ldInfo, String[] info){
		
		//read info
		String[] possibleGenotypesPrevPos = getPossibleGenotypes(ldInfo[V], ldInfo[v]);
		String[] possibleGenotypesCurrPos = getPossibleGenotypes(info[V], info[v]);
		double[] twoLocusFreq = getTwoLocusFreq(info);
		String[] possibleTwoLocusAlleles = getPossibleTwoLocusAlleles(info);
		
		//get one-locus emission with error for LD pos
		Map<EmissionKey, Double> oneLocusEmissionWithError = computePossibleOneLocusEmissionWithError(ldInfo);
		
		
		//possible two-locus emission without error
		Map<EmissionKey, Double> twoLocusEmissionWithoutError = new HashMap<EmissionKey, Double>(243);		
		for(String genoPrevPos1 : possibleGenotypesPrevPos){
			for(String genoPrevPos2 : possibleGenotypesPrevPos){
				
				SimplePair<Genotype,Genotype> genoPrevPos = new SimplePair<Genotype,Genotype>(new Genotype(genoPrevPos1), new Genotype(genoPrevPos2));
				
				for(String genoCurrPos1 : possibleGenotypesCurrPos){
					for(String genoCurrPos2 : possibleGenotypesCurrPos){
						
						SimplePair<Genotype,Genotype> genoCurrPos = new SimplePair<Genotype,Genotype>(new Genotype(genoCurrPos1), new Genotype(genoCurrPos2));
						
						for(int ibd=0; ibd<3; ibd++){
							
							EmissionKey key = new EmissionKey(genoPrevPos, genoCurrPos, ibd);
							double val = twoLocusEmissionDensityWithoutError(twoLocusFreq, possibleTwoLocusAlleles, genoCurrPos1, genoCurrPos2, genoPrevPos1, genoPrevPos2, ibd);
							twoLocusEmissionWithoutError.put(key,val);
							
						}
					}	
				}
			}
		}
		
		//precompute error maps
		Map<SimplePair<Genotype,Genotype>, Double> errorMapPrevPos = computePossibleErrorProbs(ldInfo);
		Map<SimplePair<Genotype,Genotype>, Double> errorMapCurrPos = computePossibleErrorProbs(info);
		
		//possible conditional emission with error
		Map<EmissionKey, Double> conditionalEmissionWithError = new HashMap<EmissionKey, Double>(243);	
		for(String genoPrevPos1 : possibleGenotypesPrevPos){
			for(String genoPrevPos2 : possibleGenotypesPrevPos){
				
				SimplePair<Genotype,Genotype> genoPrevPos = new SimplePair<Genotype,Genotype>(new Genotype(genoPrevPos1), new Genotype(genoPrevPos2));
				
				for(String genoCurrPos1 : possibleGenotypesCurrPos){
					for(String genoCurrPos2 : possibleGenotypesCurrPos){
						
						SimplePair<Genotype,Genotype> genoCurrPos = new SimplePair<Genotype,Genotype>(new Genotype(genoCurrPos1), new Genotype(genoCurrPos2));
						
						for(int ibd=0; ibd<3; ibd++){
							
							EmissionKey key = new EmissionKey(genoPrevPos, genoCurrPos, ibd);
							double val = conditionalEmissionDensityWithError(twoLocusEmissionWithoutError, oneLocusEmissionWithError, errorMapPrevPos, errorMapCurrPos, possibleGenotypesPrevPos, possibleGenotypesCurrPos, genoCurrPos1, genoCurrPos2, genoPrevPos1, genoPrevPos2, ibd);
							conditionalEmissionWithError.put(key,val);
							
						}
					}	
				}
			}
		}
		

		
		return conditionalEmissionWithError;
		
	}
	

	//no marginalizing assume genotype1 and genotype2 are the real deal
	public double twoLocusEmissionDensityWithoutError(double[] twoLocusFreq, String[] possibleTwoLocusAlleles, String genotypeInd1CurrPos, String genotypeInd2CurrPos, String genotypeInd1PrevPos, String genotypeInd2PrevPos, int ibd){
		
		int ibsCurrPos = getIBS(genotypeInd1CurrPos, genotypeInd2CurrPos);
		int ibsPrevPos = getIBS(genotypeInd1PrevPos, genotypeInd2PrevPos);

		
		if (ibsCurrPos < ibd || ibsPrevPos < ibd){
			return 0d;
		}

		double firstSampleProb = twoLocusGenotypeProbWithoutError(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos, genotypeInd1CurrPos);
	
		
		if (ibd==2){ //prob of sampling the first guy //correct
			return firstSampleProb;
		}
		
		else if (ibd==0){ //prob of sampling two guys independently //correct
			return firstSampleProb * twoLocusGenotypeProbWithoutError(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd2PrevPos, genotypeInd2CurrPos);
		}
		
		else if (ibd==1){ //fix this
			
			double delta1 = genotypeInd1PrevPos.charAt(0)==genotypeInd1PrevPos.charAt(1)? 1d : 0d;
			double delta2 = genotypeInd1CurrPos.charAt(0)==genotypeInd1CurrPos.charAt(1)? 1d : 0d;
			
			if (ibsPrevPos==2 && ibsCurrPos==2){
				
				if(delta1+delta2==2){
					return Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(0)),3);
				}
				else if(delta1+delta2==1){
					return getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(0)) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(1), genotypeInd1CurrPos.charAt(1)) * (getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd2PrevPos.charAt(0), genotypeInd2CurrPos.charAt(0)) + getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd2PrevPos.charAt(1), genotypeInd2CurrPos.charAt(1)) );
				}
				else{
					return Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(0)),2) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(1), genotypeInd1CurrPos.charAt(1)) +
							Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(1), genotypeInd1CurrPos.charAt(1)),2) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(0)) +
							Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(1)),2) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(1), genotypeInd1CurrPos.charAt(0)) +
							Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(1), genotypeInd1CurrPos.charAt(0)),2) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeInd1PrevPos.charAt(0), genotypeInd1CurrPos.charAt(1));
				
				}

			}
			
			else if (ibsPrevPos==1 || ibsCurrPos==1){
				char sharedPrevPos;
				char lone1PrevPos;
				char lone2PrevPos;
				char sharedCurrPos;
				char lone1CurrPos;
				char lone2CurrPos;
				
				if(genotypeInd1PrevPos.charAt(0)==genotypeInd2PrevPos.charAt(0)){
					sharedPrevPos = genotypeInd1PrevPos.charAt(0);
					lone1PrevPos = genotypeInd1PrevPos.charAt(1);
					lone2PrevPos = genotypeInd2PrevPos.charAt(1);
				}
				else if(genotypeInd1PrevPos.charAt(0)==genotypeInd2PrevPos.charAt(1)){
					sharedPrevPos = genotypeInd1PrevPos.charAt(0);
					lone1PrevPos = genotypeInd1PrevPos.charAt(1);
					lone2PrevPos = genotypeInd2PrevPos.charAt(0);
				}
				else if(genotypeInd1PrevPos.charAt(1)==genotypeInd2PrevPos.charAt(0)){
					sharedPrevPos = genotypeInd1PrevPos.charAt(1);
					lone1PrevPos = genotypeInd1PrevPos.charAt(0);
					lone2PrevPos = genotypeInd2PrevPos.charAt(1);
				}
				else if(genotypeInd1PrevPos.charAt(1)==genotypeInd2PrevPos.charAt(1)){
					sharedPrevPos = genotypeInd1PrevPos.charAt(1);
					lone1PrevPos = genotypeInd1PrevPos.charAt(0);
					lone2PrevPos = genotypeInd2PrevPos.charAt(0);
				}
				else throw new RuntimeException("something wrong with ibd");
				
				
				if(genotypeInd1CurrPos.charAt(0)==genotypeInd2CurrPos.charAt(0)){
					sharedCurrPos = genotypeInd1CurrPos.charAt(0);
					lone1CurrPos = genotypeInd1CurrPos.charAt(1);
					lone2CurrPos = genotypeInd2CurrPos.charAt(1);
				}
				else if(genotypeInd1CurrPos.charAt(0)==genotypeInd2CurrPos.charAt(1)){
					sharedCurrPos = genotypeInd1CurrPos.charAt(0);
					lone1CurrPos = genotypeInd1CurrPos.charAt(1);
					lone2CurrPos = genotypeInd2CurrPos.charAt(0);
				}
				else if(genotypeInd1CurrPos.charAt(1)==genotypeInd2CurrPos.charAt(0)){
					sharedCurrPos = genotypeInd1CurrPos.charAt(1);
					lone1CurrPos = genotypeInd1CurrPos.charAt(0);
					lone2CurrPos = genotypeInd2CurrPos.charAt(1);
				}
				else if(genotypeInd1CurrPos.charAt(1)==genotypeInd2CurrPos.charAt(1)){
					sharedCurrPos = genotypeInd1CurrPos.charAt(1);
					lone1CurrPos = genotypeInd1CurrPos.charAt(0);
					lone2CurrPos = genotypeInd2CurrPos.charAt(0);
				}
				else throw new RuntimeException("something wrong with ibd");
				
				
				if(ibsPrevPos==1 && ibsCurrPos==1){
					return getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone1PrevPos, lone1CurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone2PrevPos, lone2CurrPos);
				}

				else if(ibsPrevPos==1 && ibsCurrPos==2){
					if(delta2==1){
						return getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone1PrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone2PrevPos, sharedCurrPos);
					}
					else{
						double toReturn = 0d;
						for (int i=0; i<2; i++){
							sharedCurrPos = genotypeInd1CurrPos.charAt(i);
							lone1CurrPos = genotypeInd1CurrPos.charAt((i+1)%2);
							lone2CurrPos = sharedCurrPos==genotypeInd2CurrPos.charAt(0)? genotypeInd2CurrPos.charAt(1) : genotypeInd2CurrPos.charAt(0);
							toReturn += getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone1PrevPos, lone1CurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone2PrevPos, lone2CurrPos);
						}		
						return toReturn;
					}
				}
				
				else if (ibsPrevPos==2 && ibsCurrPos==1){
					if(delta1==1){
						return getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, lone1CurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, lone2CurrPos);
					}
					else{
						double toReturn = 0d;
						for (int i=0; i<2; i++){
							sharedPrevPos = genotypeInd1PrevPos.charAt(i);
							lone1PrevPos = genotypeInd1PrevPos.charAt((i+1)%2);
							lone2PrevPos = sharedPrevPos==genotypeInd2PrevPos.charAt(0)? genotypeInd2PrevPos.charAt(1) : genotypeInd2PrevPos.charAt(0);
							toReturn += getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, sharedPrevPos, sharedCurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone1PrevPos, lone1CurrPos) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, lone2PrevPos, lone2CurrPos);
						}		
						return toReturn;
					}
				}
			}
			
		}
		
		throw new RuntimeException("Illegal ibd combinations");

	}
	
	
	private double conditionalEmissionDensityWithError(Map<EmissionKey, Double> twoLocusEmissionWithoutError ,Map<EmissionKey, Double> oneLocusEmissionWithError, Map<SimplePair<Genotype,Genotype>, Double> errorMapPrevPos,Map<SimplePair<Genotype,Genotype>, Double> errorMapCurrPos, String[] possibleGenotypesPrevPos, String[] possibleGenotypesCurrPos, String genotypeInd1CurrPos, String genotypeInd2CurrPos, String genotypeInd1PrevPos, String genotypeInd2PrevPos, int ibd){

		double twoLocusEmissionDensity = 0d;
		for(String trueGenotypeInd1CurrPos : possibleGenotypesCurrPos){
			SimplePair<Genotype,Genotype> errorKey1 = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd1CurrPos), new Genotype(genotypeInd1CurrPos));
			double errorProbInd1CurrPos = errorMapCurrPos.get(errorKey1);
			
			for(String trueGenotypeInd2CurrPos: possibleGenotypesCurrPos){
				SimplePair<Genotype,Genotype> errorKey2 = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd2CurrPos), new Genotype(genotypeInd2CurrPos));
				double errorProbInd2CurrPos = errorMapCurrPos.get(errorKey2);
				SimplePair<Genotype,Genotype> genoCurrPos = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd1CurrPos), new Genotype(trueGenotypeInd2CurrPos));
				
				for(String trueGenotypeInd1PrevPos: possibleGenotypesPrevPos){
					SimplePair<Genotype,Genotype> errorKey3 = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd1PrevPos), new Genotype(genotypeInd1PrevPos));
					double errorProbInd1PrevPos  = errorMapPrevPos.get(errorKey3);
					
					for(String trueGenotypeInd2PrevPos: possibleGenotypesPrevPos){
						SimplePair<Genotype,Genotype> errorKey4 = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd2PrevPos), new Genotype(genotypeInd2PrevPos));
						double errorProbInd2PrevPos  = errorMapPrevPos.get(errorKey4);
						
						SimplePair<Genotype,Genotype> genoPrevPos = new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeInd1PrevPos), new Genotype(trueGenotypeInd2PrevPos));
						
						twoLocusEmissionDensity += errorProbInd1CurrPos * errorProbInd2CurrPos * errorProbInd1PrevPos * errorProbInd2PrevPos * twoLocusEmissionWithoutError.get(new EmissionKey(genoPrevPos, genoCurrPos, ibd));
					
					}
				}
			}
		}
		
		
		double emissionProbPrevPos = oneLocusEmissionWithError.get(new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(genotypeInd1PrevPos), new Genotype(genotypeInd2PrevPos)), null, ibd));
		double toReturn = twoLocusEmissionDensity/emissionProbPrevPos;	

		if(Math.abs(toReturn)<TOL) toReturn = 0d;
		if(Math.abs(toReturn-1)<TOL) toReturn = 1d;
		
		if (toReturn >1d || toReturn<0d){
			System.out.println("IN EMISSION");
			//double twoNoError = twoLocusEmissionDensityWithoutError(info, genotypeInd1CurrPos, genotypeInd2CurrPos, genotypeInd1PrevPos, genotypeInd2PrevPos, ibd, currPos, prevPos);
			//double oneNoError = emissionDensityWithoutError(prevInfo, genotypeInd1PrevPos, genotypeInd2PrevPos, ibd);
			//System.out.println(twoNoError/oneNoError);
		}
		
		assert toReturn <= 1d && toReturn >= 0d;
			
		return toReturn;

	}
	

	
	private Map<SimplePair<Genotype,Genotype>, Double> computePossibleErrorProbs(String[] info){
		
		Map<SimplePair<Genotype,Genotype>, Double> toReturn = new HashMap<SimplePair<Genotype,Genotype>, Double>(9);
		String[] possibleGenotypes = getPossibleGenotypes(info[V], info[v]);
		
		for(String trueGenotype : possibleGenotypes){
			for(String observedGenotype: possibleGenotypes){
				toReturn.put(new SimplePair<Genotype,Genotype>(new Genotype(trueGenotype), new Genotype(observedGenotype)), errorProbability(trueGenotype, observedGenotype));
			}
		}
		
		return toReturn;
	}
	
	//likelihood of the observedGenotype given the true genotype 
	private double errorProbability(String trueGenotype, String observedGenotype){

		double delta = observedGenotype.charAt(0) == observedGenotype.charAt(1) ? 1d : 0d;
		int ibs = getIBS(observedGenotype, trueGenotype);
		
		if (trueGenotype.charAt(0) == trueGenotype.charAt(1)){//correct
			return (2d - delta) * Math.pow(this.genotypingErrorRate / 3d, 2 - ibs) * Math.pow(1 - this.genotypingErrorRate, ibs);
		} else {
			if(ibs == 2){
				return Math.pow(1 - this.genotypingErrorRate, 2) + Math.pow(this.genotypingErrorRate / 3d, 2);
			} else {
				return (2d-delta)* Math.pow(this.genotypingErrorRate / 3d, 2 - ibs) * Math.pow(1 - this.genotypingErrorRate, ibs);
			}
		}
	}
	
	
	
	/////// GENOTYPE PROBABILITY //////////
	private Map<Genotype, Double> computePossibleGenotypeProbWithError(String[] info){
		
		String[] possibleGenotypes = getPossibleGenotypes(info[V], info[v]);
		double[] oneLocusFreq = getOneLocusFreq(info);
		
		Map<Genotype, Double> genotypeProbWithoutError = new HashMap<Genotype, Double>(3);
		for(String geno : possibleGenotypes){
			genotypeProbWithoutError.put(new Genotype(geno), genotypeProbWithoutError(oneLocusFreq, geno));
		}
		
		Map<Genotype, Double> genotypeProbWithError = new HashMap<Genotype, Double>(3);
		for(String geno: possibleGenotypes){
			genotypeProbWithError.put(new Genotype(geno), genotypeProbWithError(genotypeProbWithoutError, possibleGenotypes, geno));
		}
		
		return genotypeProbWithError;
		
	}
	
	
	private double genotypeProbWithoutError(double[] oneLocusFreq, String g){
		
		double delta = g.charAt(0)==g.charAt(1)? 1:0;
		double lkhd = (2-delta) * getAlleleProb(oneLocusFreq, g.charAt(0)) * getAlleleProb(oneLocusFreq, g.charAt(1));

		return lkhd;
	}
	
	
	private double genotypeProbWithError(Map<Genotype, Double> genotypeProbWithoutError, String[] possibleGenotypes, String genotype){
		
		double toReturn = 0d;
		for(String trueGenotype : possibleGenotypes){
				toReturn += errorProbability(trueGenotype, genotype) *  genotypeProbWithoutError.get(new Genotype(trueGenotype));
		}

		return toReturn;
	}
	
	
	private Map<SimplePair<Genotype,Genotype>, Double> computePossibleConditionalGenotypeProbWithError(String[] ldInfo, String[] info){
		
		//read info
		double[] twoLocusFreq = getTwoLocusFreq(info);
		String[] possibleTwoLocusAlleles = getPossibleTwoLocusAlleles(info);
		String[] possibleGenotypesPrevPos = getPossibleGenotypes(ldInfo[V], ldInfo[v]);
		String[] possibleGenotypesCurrPos = getPossibleGenotypes(info[V], info[v]);
		
		// get one locus genotype probs
		Map<Genotype, Double> oneLocusGenotypeProbWithError = computePossibleGenotypeProbWithError(ldInfo);
		
		
		//two locus geno prob WITHOUT error
		Map<SimplePair<Genotype,Genotype>, Double> twoLocusGenotypeProbWithoutError = new HashMap<SimplePair<Genotype,Genotype>, Double>(9);
		for(String prevGeno : possibleGenotypesPrevPos){
			for(String currGeno : possibleGenotypesCurrPos){
				twoLocusGenotypeProbWithoutError.put(new SimplePair<Genotype, Genotype>(new Genotype(prevGeno), new Genotype(currGeno)), twoLocusGenotypeProbWithoutError(twoLocusFreq, possibleTwoLocusAlleles, prevGeno, currGeno));
			}
		}
		
		//WITH error
		Map<SimplePair<Genotype,Genotype>, Double>conditionalGenotypeProbWithError = new HashMap<SimplePair<Genotype,Genotype>, Double>(9);
		for(String prevGeno : possibleGenotypesPrevPos){
			for(String currGeno : possibleGenotypesCurrPos){
				conditionalGenotypeProbWithError.put(new SimplePair<Genotype,Genotype>(new Genotype(prevGeno), new Genotype(currGeno)), twoLocusGenotypeProbWithError(twoLocusGenotypeProbWithoutError, oneLocusGenotypeProbWithError, possibleGenotypesPrevPos, possibleGenotypesCurrPos, prevGeno, currGeno));
			}
		}
		
		return conditionalGenotypeProbWithError;
		
	}
	
	
	private double twoLocusGenotypeProbWithoutError (double[] twoLocusFreq, String[] possibleTwoLocusAlleles, String genotypeLoc1, String genotypeLoc2){ //probability of observing two-locus genotypes in one person
		
		double toReturn;
		
		int deltaLocus1 = genotypeLoc1.charAt(0) == genotypeLoc1.charAt(1) ? 1 : 0;
		int deltaLocus2 = genotypeLoc2.charAt(0) == genotypeLoc2.charAt(1) ? 1 : 0;
		
		if (deltaLocus1==1 && deltaLocus2==1){ //both loci are homozygous
			toReturn =  Math.pow(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(0), genotypeLoc2.charAt(0)),2);
		}
		else if (deltaLocus1 + deltaLocus2 == 1){ //one locus heterozygous
			toReturn =  2 * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(0), genotypeLoc2.charAt(0)) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(1), genotypeLoc2.charAt(1));
		}
		else if (deltaLocus1 + deltaLocus2 == 0){ //both loci heterozygous	
			toReturn = 2*(getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(0), genotypeLoc2.charAt(0)) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(1), genotypeLoc2.charAt(1)) + 
					getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(0), genotypeLoc2.charAt(1)) * getTwoLocusAlleleFreq(twoLocusFreq, possibleTwoLocusAlleles, genotypeLoc1.charAt(1), genotypeLoc2.charAt(0)));
		}
		
		else throw new RuntimeException("somethiong went wrong while computing two locus genotype emissions");

		return toReturn;
		
	}
	
	
	private double twoLocusGenotypeProbWithError (Map<SimplePair<Genotype,Genotype>, Double> twoLocusGenotypeProbWithoutError, Map<Genotype, Double> oneLocusGenotypeProbWithError, String[] possibleGenotypes1, String[] possibleGenotypes2, String genotypeLoc1, String genotypeLoc2){ 
		
		double twoLocusGenotypeProb = 0d;
		for(String trueGenotypeLoc1 : possibleGenotypes1){
			for(String trueGenotypeLoc2: possibleGenotypes2){
				twoLocusGenotypeProb += errorProbability(trueGenotypeLoc1, genotypeLoc1) * errorProbability(trueGenotypeLoc2, genotypeLoc2) * twoLocusGenotypeProbWithoutError.get(new SimplePair<Genotype,Genotype>(new Genotype(trueGenotypeLoc1), new Genotype(trueGenotypeLoc2)));
			}
		}
		
		double toReturn = twoLocusGenotypeProb/oneLocusGenotypeProbWithError.get(new Genotype(genotypeLoc1));
		
		return toReturn;
		
	}
	
	
	
	

	

	
	//////// MISCELLANEOUS HELPER FUNCTIONS ////////	
	private double getAlleleProb(double[] oneLocusFreq, char allele){
		
		if (allele == 'A') return oneLocusFreq[0];
		if (allele == 'C') return oneLocusFreq[1];
		if (allele == 'G') return oneLocusFreq[2];
		if (allele == 'T') return oneLocusFreq[3];
		throw new RuntimeException("Tried to access illegal allele");
	}
	
	
	private double getTwoLocusAlleleFreq(double[] TwoLocusFreq, String[] possibleTwoLocusAlleles, char a1, char a2){

		String observed = ""+a1+a2;
		
		if (observed.equals(possibleTwoLocusAlleles[0])) return TwoLocusFreq[0];
		else if (observed.equals(possibleTwoLocusAlleles[1])) return TwoLocusFreq[1];
		else if (observed.equals(possibleTwoLocusAlleles[2])) return TwoLocusFreq[2];
		else if (observed.equals(possibleTwoLocusAlleles[3])) return TwoLocusFreq[3];
		else return 0d;
	}	
	
	
	
	private double[] getOneLocusFreq(String[] info){
		
		return new double[]{Double.parseDouble(info[pA]), Double.parseDouble(info[pC]), Double.parseDouble(info[pG]), Double.parseDouble(info[pT])};
		
	}
	
	
	private double[] getTwoLocusFreq(String[] info){
		
		return new double[]{Double.parseDouble(info[A]), Double.parseDouble(info[B]), Double.parseDouble(info[C]), Double.parseDouble(info[D])};
	}
	
	
	private String[] getPossibleTwoLocusAlleles(String[] info){
		return new String[]{info[u]+info[v], info[u]+info[V], info[U]+info[v], info[U]+info[V]};
	}

	
	//returns the number of characters the two strings have in common (unorders them, first)
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
	
	
	private String[] getPossibleGenotypes(String A, String a){
		
		assert !A.equals(a);
		
		String AA = A+A;
		String Aa = A+a;
		String aa = a+a;
		
		return new String[]{AA,Aa,aa};
	}
	
	
	/////// TRANSITION PROBABILITY ///////
	//Returns phi(ibdNew | ibdOld, l), the transition probability of moving from ibdOld to ibdNew over a sequence with l expected recombinations
	private double transitionDensity(int ibdNew, int ibdOld, double rhoBetween, Relationship rel){
		
		//special case for parent-offspring
		if(rel.getMarginalProbs()[1] == 1d){
			if(ibdOld == 1d && ibdNew ==1d) return 1d;
			else return 0d;
		}

		double toReturn = -1d;
		//double alpha = -rel.getMeiosisRate() * Math.log(1-rhoBetween); //rhoBetween cancels out
		double alpha = rel.getMeiosisRate() * rhoBetween;
		if(ibdOld == 0){
			double probTo2 = Math.exp(-alpha * rel.getMarginalProbs()[1]) * rel.getMarginalProbs()[2] / (rel.getMarginalProbs()[1] - 1d);
			probTo2 += Math.exp(-alpha) * rel.getMarginalProbs()[1];
			probTo2 += Math.exp(-alpha) * rel.getMarginalProbs()[0] * rel.getMarginalProbs()[1] / (rel.getMarginalProbs()[1] - 1d);
			probTo2 += rel.getMarginalProbs()[2];
			double probTo1 = (1d - Math.exp(-alpha)) * rel.getMarginalProbs()[1];
			if(ibdNew == 2){
				toReturn = probTo2;
			}
			else if(ibdNew == 1){
				toReturn =  probTo1;
			}
			else if(ibdNew == 0){
				toReturn =  1d - probTo2 - probTo1;
			}
		}
		
		else if(ibdOld == 1){
			if(ibdNew == 0){
				toReturn =  (1d - Math.exp(-alpha)) * rel.getMarginalProbs()[0];
			}
			if(ibdNew == 1){
				toReturn = (1d - Math.exp(-alpha)) * rel.getMarginalProbs()[1] + Math.exp(-alpha);
			}
			if(ibdNew == 2){
				toReturn = (1d - Math.exp(-alpha)) * rel.getMarginalProbs()[2];
			}
		}
		
		else if(ibdOld == 2){
			double probTo0 = Math.exp(-alpha * rel.getMarginalProbs()[1]) * rel.getMarginalProbs()[0] / (rel.getMarginalProbs()[1] - 1d);
			probTo0 += Math.exp(-alpha) * rel.getMarginalProbs()[1];
			probTo0 += Math.exp(-alpha) * rel.getMarginalProbs()[2] * rel.getMarginalProbs()[1] / (rel.getMarginalProbs()[1] - 1d);
			probTo0 += rel.getMarginalProbs()[0];
			double probTo1 = (1d - Math.exp(-alpha)) * rel.getMarginalProbs()[1];
			if(ibdNew == 0){
				toReturn = probTo0;
			}
			else if(ibdNew == 1){
				toReturn = probTo1;
			}
			else if(ibdNew == 2){
				toReturn = 1d - probTo1 - probTo0;
			}
		}
		
		else throw new RuntimeException("Invalid transition");
		
		//correct for machine precision
		if (Math.abs(toReturn) < TOL){
			toReturn = 0d;
		}
		if (Math.abs(toReturn-1d) < TOL){
			toReturn = 1d;
		}
		
		assert toReturn >= 0d && toReturn <= 1d;
		
		return toReturn;
	}
	
	


	 
}