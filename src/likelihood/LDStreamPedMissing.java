package likelihood;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;


import utility.Cubic;
import utility.DataParser;

//input
//tped, tfam

//output
// currPos ldPos U u V v pA pC pG pT A B C D
//U = a1 at ldPos
//V = a1 at currPos

public class LDStreamPedMissing {

	//input data columns
	private static int CHROM = 0;
	private static int GENPOS = 2;
	private static int PHYPOS = 3;

	
	public static void writeLdOutfile(String fileName, String outPath, double back, boolean conditional) throws IOException{
		
		int numIndiv = DataParser.countLines(fileName+".tfam"); 

		//open files
		BufferedReader tpedfile = DataParser.openReader(fileName+".tped");
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//write header for outfile
		writer.write("CHROM\tPOS\tLDPOS\tA\ta\tB\tb\tp_A\tp_C\tp_G\tp_T\tp_ab\tp_aB\tp_Ab\tp_AB\n");
		

		
		//cubic solver
		Cubic solver = new Cubic();
		
		//data structures
		List<Integer> prevPosList = new ArrayList<Integer>();
		List<Double> prevDistList = new ArrayList<Double>();
		List<String[]> prevGenotypeList = new ArrayList<String[]>();
		List<double[]> prevStationaryDistList = new ArrayList<double[]>();
		List<char[]> prevAlleleTypeList = new ArrayList<char[]>();

		
		//read data
		String tpedline;
		int prevChrom = 1;
		int currChrom = 1;
		while ((tpedline = tpedfile.readLine())!=null){

			String[] tpedFields = tpedline.split("\\s");

			
			// update current info
			prevChrom = currChrom;
			currChrom = Integer.parseInt(tpedFields[CHROM]);
			int currPos = Integer.parseInt(tpedFields[PHYPOS]);
			double currDist = Double.parseDouble(tpedFields[GENPOS]);
			if(currChrom!=prevChrom){
				prevDistList.clear();
				prevPosList.clear();
				prevGenotypeList.clear();
				prevStationaryDistList.clear();
				prevAlleleTypeList.clear();
			}
					
			if(prevPosList.contains(currPos)) continue; //skip if it was already processed
			
					
			String[] currGenotypes = new String[numIndiv];
			for (int i=0; i<numIndiv; i++){				
				currGenotypes[i] = tpedFields[2*i+4] + tpedFields[2*i+5];
			}
			double[] currStationaryDist = estimateStationaryDist(currGenotypes);
			
			
			char[] currAlleleTypes = getAlleles(currStationaryDist);

			
			//before processing current position, remove any previous candidate markers that are too far away
			//if only one markers is left, use that no matter what
		
			// update previous info			
			while(prevDistList.size()>1){
				
				if(currDist - prevDistList.get(0) > back){
					prevDistList.remove(0);
					prevPosList.remove(0);
					prevGenotypeList.remove(0);
					prevStationaryDistList.remove(0);
					prevAlleleTypeList.remove(0);
				}
				
				else
					break;
				
			}
				
			
			// process this line; write to outfile
			processInfo(solver, writer, currChrom, currPos, prevPosList, currGenotypes, prevGenotypeList, currStationaryDist, prevStationaryDistList, currAlleleTypes, prevAlleleTypeList);
			
			
			//add current info
			if(conditional){
				prevDistList.add(currDist);
				prevPosList.add(currPos);
				prevGenotypeList.add(currGenotypes);
				prevStationaryDistList.add(currStationaryDist);
				prevAlleleTypeList.add(currAlleleTypes);
		
			}
			
		
		}
		
		
		tpedfile.close();
		writer.close();
		
		
	}
	
	
	private static void processInfo (Cubic solver, PrintWriter writer, int currChrom, int currPos, List<Integer> prevPosList, String[] currGenotypes, List<String[]> prevGenotypeList, double[] currStationaryDist, List<double[]> prevStationaryDistList, char[] currAlleleTypes, List<char[]> prevAlleleTypeList){
		
		
		if (prevGenotypeList.size()==0){
			writeResults(writer, currChrom, currPos, -1, new char[]{'X','X'}, currAlleleTypes, currStationaryDist, new double[]{-1,-1,-1,-1}); //return if there is no previous snps to condition on
			return;
		}
		
		
		//initialize
		int bestPos = prevPosList.get(0);
		double bestLd = 0d;
		double[] bestTwoLocusFreq = new double[4];
		char[] bestPrevAlleleTypes = new char[2];
		
		//consider every prevPos
		for (int i=0; i<prevPosList.size(); i++){
			
			int prevPos = prevPosList.get(i);
			String[] prevGenotypes = prevGenotypeList.get(i);
			char[] prevAlleleTypes = prevAlleleTypeList.get(i);
			double[] prevStationaryDist = prevStationaryDistList.get(i);
			
			
			//estimate two locus freq {uv, uV, Uv, UV}
			double[] twoLocusFreq = estimateTwoLocusFreq(solver, prevGenotypes, currGenotypes, prevAlleleTypes[0], prevAlleleTypes[1], currAlleleTypes[0], currAlleleTypes[1]);

			
			//compute LD
			char A = prevAlleleTypes[0];
			char B = currAlleleTypes[0];
			double pA = getAlleleFreq(A, prevStationaryDist);
			double pB = getAlleleFreq(B, currStationaryDist);
			double pAB = twoLocusFreq[3];
			
			double currLd = computeLD(pA, pB, pAB);
			
			if (currLd >= bestLd){
				bestLd = currLd;
				bestPos = prevPos;
				bestTwoLocusFreq = twoLocusFreq;
				bestPrevAlleleTypes = prevAlleleTypes;
			}
		}

		// write results
		writeResults(writer, currChrom, currPos, bestPos, bestPrevAlleleTypes, currAlleleTypes, currStationaryDist, bestTwoLocusFreq);

		
	}
	

	
	
	private static double[] estimateTwoLocusFreq(Cubic solver, String[] genotypesAtLocus1, String[] genotypesAtLocus2, char U, char u, char V, char v){


		int[] genotypeCounts = countGenotypes(genotypesAtLocus1, genotypesAtLocus2, U, u, V, v);
		

		int numGeneCopies = 0;
		for(int cnt : genotypeCounts) numGeneCopies += cnt;
		numGeneCopies  = 2 * numGeneCopies;

		double cA = 2*genotypeCounts[0] + genotypeCounts[1] + genotypeCounts[3];
		double cB = 2*genotypeCounts[2] + genotypeCounts[1] + genotypeCounts[5];	
		double cC = 2*genotypeCounts[6] + genotypeCounts[7] + genotypeCounts[3]; 
		double cD = 2*genotypeCounts[8] + genotypeCounts[7] + genotypeCounts[5];
		double e = genotypeCounts[4];
		double A;
		double B;
		double C;
		double D;
		
		
		if (e==0){ //if there's no doubly heterozygous individuals, no need to compute p
			A = cA/numGeneCopies;
			B = cB/numGeneCopies;
			C = cC/numGeneCopies;
			D = cD/numGeneCopies;
		}
		
		else{
			
			double a = -2*Math.pow(e,2);
			double b = 3*Math.pow(e,2) + (-cA + cB + cC - cD)*e;
			double c = -Math.pow(e,2) + (cA-cB-cC+cD)*e - cA*cD- cB*cC;
			double d = cA*cD; 

			solver.solve(a,b,c,d);
			assert solver.roots.size() > 0;
		
			//handle multiple roots: return the root with the highest likelihood; assert 0<p<1
			double p = pickBestP(solver.roots, genotypeCounts, cA,cB,cC,cD, numGeneCopies);
	
			if (p<0d){
				System.out.println("negative p");
				solver.solve(a,b,c,d);
			}
			
			A = (cA + p*genotypeCounts[4])/numGeneCopies;
			B = (cB + (1-p)*genotypeCounts[4])/numGeneCopies;		
			C = (cC + (1-p)*genotypeCounts[4])/numGeneCopies;
			D = (cD + p*genotypeCounts[4])/numGeneCopies;
		}
		
		assert A>0d && B>0d && C>0d && D>0d && Math.abs(1-(A+B+C+D)) < 1e-10;
		
		return new double[]{A,B,C,D};

		

	}
	
	
	private static void writeResults(PrintWriter writer, int currChrom, int currPos, int ldPos, char[] Uu, char[] Vv, double[] freq, double[] twoLocusFreq){
		
		//chrom
		writer.write(String.format("%d\t",currChrom));
		
		//pos
		writer.write(String.format("%d\t%d\t",currPos,ldPos));
		
		//alleles
		for (char a: Uu){
			writer.write(String.format("%c\t",a));
		}
		for (char a: Vv){
			writer.write(String.format("%c\t",a));
		}
		
		//freq
		for(double f:freq){
			
			writer.write(String.format("%f\t",f));
		}
		for(double f : twoLocusFreq){
			writer.write(String.format("%f\t",f));
		}
		
		writer.write("\n");
		
	}

	
	private static double[] estimateStationaryDist (String[] genotypes){

		
		int numIndiv = genotypes.length;
		int numAlleleCopies = 0;
		double[] toReturn = new double[4];
		
		// count alleles
		for (int ind=0; ind<numIndiv; ind++){
			char g1 = genotypes[ind].charAt(0);
			char g2 = genotypes[ind].charAt(1);
			
			//increment copies
			numAlleleCopies += 2;
				
			if (g1=='A') toReturn[0]++;
			else if (g1=='C') toReturn[1]++;
			else if (g1=='G') toReturn[2]++;
			else if (g1=='T') toReturn[3]++;
			else if(g1=='0') numAlleleCopies--;
			else throw new RuntimeException("Illegal allele in LD");
				
			if (g2=='A') toReturn[0]++;
			else if (g2=='C') toReturn[1]++;
			else if (g2=='G') toReturn[2]++;
			else if (g2=='T') toReturn[3]++;	
			else if(g2=='0') numAlleleCopies--;
			else throw new RuntimeException("Illegal allele in LD");
			}

		// divide by total number of allele copies
		for (int allele=0; allele<4; allele++){
			toReturn[allele] = toReturn[allele]/numAlleleCopies;
		}

		
		return toReturn;
		
	}
	
	
	//counts different genotype combinations
	private static int[] countGenotypes(String[] genotypesAtLocus1, String[] genotypesAtLocus2, char U, char u, char V, char v){
		
		int numIndiv = genotypesAtLocus1.length;
		
		//possible genotypes
		String uu = ""+u+u;
		String vv = ""+v+v;
		String UU = ""+U+U;
		String VV = ""+V+V;

		
		int[] genotypeCounts = new int[9];
		
		for (int ind=0; ind < numIndiv; ind++){
			
			String g1 = genotypesAtLocus1[ind];
			String g2 = genotypesAtLocus2[ind];
			
			//skip if missing data
			if(g1.charAt(0)=='0' || g1.charAt(1)=='0' || g2.charAt(0)=='0' || g2.charAt(1)=='0')
				continue;
			
			
			if (g1.equals(uu)){
				if (g2.equals(vv)) genotypeCounts[0]++;
				else if (g2.equals(VV)) genotypeCounts[2]++;
				else genotypeCounts[1]++;
			}
			
			else if (g1.equals(UU)){
				if (g2.equals(vv)) genotypeCounts[6]++;
				else if (g2.equals(VV)) genotypeCounts[8]++;
				else genotypeCounts[7]++;
			}	
			
			else{
				if (g2.equals(vv)) genotypeCounts[3]++;
				else if (g2.equals(VV)) genotypeCounts[5]++;
				else genotypeCounts[4]++;
			}
			
			
		}
		

		return genotypeCounts;
		
	}
	
	
	private static double pickBestP(List<Double> roots, int[] genotypeCounts, double cA, double cB, double cC, double cD, int numGeneCopies){
		
		double bestP = -1;
		double highestLikelihood = Double.NEGATIVE_INFINITY;
		
		for (double p : roots){
			
			if (p==Double.NaN || p < 0d || p > 1d) continue;
			
			double A = (cA + p*genotypeCounts[4])/numGeneCopies;
			double B = (cB + (1-p)*genotypeCounts[4])/numGeneCopies;		
			double C = (cC + (1-p)*genotypeCounts[4])/numGeneCopies;
			double D = (cD + p*genotypeCounts[4])/numGeneCopies;
			double[] marginals = getMarginals(A,B,C,D);
			
			double candidateLikelihood = computeLogLikelihood(genotypeCounts, marginals);
			
			if (candidateLikelihood > highestLikelihood){
				highestLikelihood = candidateLikelihood;
				bestP = p;
			}
		}
		
		assert bestP!=-1;
		
		return bestP;
		
	}
	
	
	private static double computeLogLikelihood(int[] counts, double[] marginals){
		
		double toReturn = 0d;
		
		for (int idx=0; idx < counts.length; idx++){
			if (marginals[idx] > 0){ 
				toReturn += (counts[idx]*Math.log(marginals[idx]));
			}
		}
		
		return toReturn;

	}
	
	
	private static double[] getMarginals (double A, double B, double C, double D){ //probability of observing a particular genotype
		
		double[] marginals = new double[9];
		
		marginals[0] = Math.pow(A,2);
		marginals[1] = 2*A*B;
		marginals[2] = Math.pow(B,2);
		marginals[3] = 2*A*C;
		marginals[4] = 2*A*D + 2*B*C;
		marginals[5] = 2*B*D;
		marginals[6] = Math.pow(C,2);
		marginals[7] = 2*C*D;
		marginals[8] = Math.pow(D,2);
		
		return marginals;
		
	}
	
	
	//TODO what about for fixed sites???
	private static double computeLD (double pA, double pB, double pAB){
		
		if(pA==0 || pB==0 || pA==1 || pB==1) return 0;
		
		double D = pAB - pA*pB;
		if (D==0d) return 0;
	
		
		double denom = pA * (1-pA) * pB * (1-pB);

		double rSquared = Math.pow(D,2)/denom;
		
		assert rSquared >= 0d && rSquared <= 1d;
		
		return rSquared;
		
		
	}
	
	
	private static double getAlleleFreq(char a, double[] af){
		
		if (a=='A') return af[0];
		else if (a=='C') return af[1];
		else if (a=='G') return af[2];
		else if (a=='T') return af[3];
		else if(a=='X') return 0;
		else throw new RuntimeException("Illegal allele");
	}
	

	//given frequencies, return alleles with nonzero entries
	private static char[] getAlleles(double[] freq){
		
		char[] toReturn = new char[2];
		int i=0;
		
		if(freq[0] > 0) toReturn[i++] = 'A';
		if(freq[1] > 0) toReturn[i++] = 'C';
		if(freq[2] > 0) toReturn[i++] = 'G';
		if(freq[3] > 0) toReturn[i++] = 'T';
		
		//if only one allele is observed, add X for the alternate allele
		if(i==1) toReturn[i] = 'X';
		
		return toReturn;
		
	}


	
}