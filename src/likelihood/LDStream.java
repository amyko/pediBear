package likelihood;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import utility.Cubic;
import utility.DataParser;

//hapmap input data format:
//pos a1 a2 geno1 geno2 ...

//output
// currPos ldPos U u V v pA pC pG pT A B C D
//U = a1 at ldPos
//V = a1 at currPos

public class LDStream {

	//input data columns
	private static int POS = 0;
	private static int A1 = 1;
	private static int A2 = 2;
	private static int INDIV = 3;
	private static int HEADER_LENGTH = 3;
	
	
	
	public static void addError(String inPath, String outPath, double epsilon, Random rGen) throws IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//write header
		String[] header = reader.readLine().split("\t");
		DataParser.writeStringArray(writer, header);

		String line;
		while((line = reader.readLine())!=null){
			
			//read line
			String[] fields = line.split("\\s");
			
			//for each genotype
			for(int i=0; i<fields.length; i++){
				
				
				for(int j=0; j<2; j++){
					
					if(rGen.nextDouble() < epsilon){
						if(fields[i].charAt(j)=='A')
							writer.write("T");
						else
							writer.write("A");
					}
					else{
						writer.write(fields[i].charAt(j));
					}
					
				}
				
				writer.write("\t");
				
				
			}
			
			writer.write("\n");
				
			
		}
		
		
		//close files
		reader.close();
		writer.close();
		
	}
	
	
	public static void thin(String inPath, String outPath, double thresh, Random rGen) throws IOException{
	
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//write header
		String[] header = reader.readLine().split("\t");
		DataParser.writeStringArray(writer, header);
		int numSnp = DataParser.countLines(inPath) - 1;
		
		int numAdded = 0;
		String line;
		while((line = reader.readLine())!=null){
			
			if(rGen.nextDouble() < thresh){
				writer.write(line);
				writer.write("\n");
				numAdded++;
			}
				
			
		}
		
		
		//close files
		reader.close();
		writer.close();
		
		System.out.println(String.format("Proportion kept: %f", numAdded/(double)numSnp));
	}
	
	
	
	//run left sliding window and prune away snps with r_sqr greater than threshold value
	//input: geno file
	//output: pruned geno file
	public static void pruneLD(String inPath, String outPath, int windowSize, double thresh) throws IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//write header
		String[] header = reader.readLine().split("\t");
		DataParser.writeStringArray(writer, header);
		
		//data structures
		Map<Integer, Double> posToFreq = new HashMap<Integer,Double>();
		Map<Integer, String[]> posToFields = new HashMap<Integer,String[]>();
		Cubic solver = new Cubic();
		
		int numAdded = 0;
		int numSnp = DataParser.countLines(inPath) - 1;
		//int numSnp = 10;
		
		//sliding window
		for(int locus=1; locus<numSnp; locus++){
			
			//if(locus%10000==0) System.out.println(locus);
			
			//read data
			String[] currFields = reader.readLine().split("\t");
			String[] currGenotypes = Arrays.copyOfRange(currFields, INDIV, currFields.length);
			
			//estimate freq
			double[] stationaryDist = estimateStationaryDist(currGenotypes);
			double pA = getAlleleFreq(currFields[A1].charAt(0), stationaryDist);
			
			//for every previous position in active list
			boolean add = true;
			for(int prevPos : posToFreq.keySet()){
				
				//get prev info
				String[] prevFields = posToFields.get(prevPos);
				String[] prevGenotypes = Arrays.copyOfRange(prevFields, INDIV, prevFields.length);
				double pB = posToFreq.get(prevPos);
				
				//estimate haplo freq
				double pAB = estimateTwoLocusFreq(solver, prevGenotypes, currGenotypes, prevFields[A1].charAt(0), prevFields[A2].charAt(0), currFields[A1].charAt(0), currFields[A2].charAt(0))[3];
				
				//compute LD
				double rSqr = computeLD(pA, pB, pAB);
				if(rSqr > thresh){
					add = false;
					break;
				}
			}
			
			//skip if ineligible
			if(add==false) continue;
			
			//add this line
			numAdded++;
			int currPos = Integer.parseInt(currFields[POS]);
			posToFreq.put(currPos, pA);
			posToFields.put(currPos, currFields);
			
			//remove from active list snps that are farther than the block size
			while(true){
				int key = Collections.min(posToFreq.keySet());
				
				if (currPos-key <= windowSize) break;
				
				posToFreq.remove(key);
				posToFields.remove(key);
				
			}
			
			//write this line to file
			DataParser.writeStringArray(writer, currFields);

		}
		
		//close files
		reader.close();
		writer.close();
		
		System.out.println(String.format("Proportion kept: %f", numAdded/(double)numSnp));
		
	}
	
	//input: geno file
	//output: info file
	public static void writeLdOutfile(String inPath, String outPath, int back) throws IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//write header for outfile
		writer.write("POS\tLDPOS\tU\tu\tV\tv\tpA\tpC\tpG\tpT\tA\tB\tC\tD\n");
		
		//num indiv
		int numIndiv = reader.readLine().split("\\s").length-HEADER_LENGTH; //skip header
		
		//TODO temporary solution for header error
		numIndiv = 100;
		
		//cubic solver
		Cubic solver = new Cubic();
		
		//data structures
		List<Integer> prevPosList = new ArrayList<Integer>();
		List<String[]> prevGenotypeList = new ArrayList<String[]>();
		List<double[]> prevStationaryDistList = new ArrayList<double[]>();
		List<char[]> prevAlleleTypeList = new ArrayList<char[]>();

		
		//read data
		String line;
		while ((line = reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			// update current info
			int currPos = Integer.parseInt(fields[POS]);
					
			if(prevPosList.contains(currPos)) continue; //skip if it was already processed
			
			char[] currAlleleTypes = new char[]{fields[A1].charAt(0), fields[A2].charAt(0)};
			
			String[] currGenotypes = new String[numIndiv];
			for (int i=0; i<numIndiv; i++){				
				currGenotypes[i] = fields[i+HEADER_LENGTH];
			}
			double[] currStationaryDist = estimateStationaryDist(currGenotypes);
			
			// process this line; write to outfile
			processInfo(solver, writer, back, currPos, prevPosList, currGenotypes, prevGenotypeList, currStationaryDist, prevStationaryDistList, currAlleleTypes, prevAlleleTypeList);
			
			// update previous info
			while(prevPosList.size()>0){
				
				if(currPos - prevPosList.get(0) > back){
					prevPosList.remove(0);
					prevGenotypeList.remove(0);
					prevStationaryDistList.remove(0);
					prevAlleleTypeList.remove(0);
				}
				
				else
					break;
				
			}
			
			/*
			if (prevPosList.size() >= back){ //remove the leftmost snp
				prevPosList.remove(0);
				prevGenotypeList.remove(0);
				prevStationaryDistList.remove(0);
				prevAlleleTypeList.remove(0);
			}
			*/
			
			//add new info
			prevPosList.add(currPos);
			prevGenotypeList.add(currGenotypes);
			prevStationaryDistList.add(currStationaryDist);
			prevAlleleTypeList.add(currAlleleTypes);
		}
		
		reader.close();
		writer.close();
		
		
	}
	
	
	private static void processInfo (Cubic solver, PrintWriter writer, int back, int currPos, List<Integer> prevPosList, String[] currGenotypes, List<String[]> prevGenotypeList, double[] currStationaryDist, List<double[]> prevStationaryDistList, char[] currAlleleTypes, List<char[]> prevAlleleTypeList){
		
		
		assert prevPosList.size() <= back; // make sure we don't go back farther than specified
		if (prevGenotypeList.size()==0){
			writeResults(writer, currPos, -1, new char[]{'X','X'}, currAlleleTypes, currStationaryDist, new double[]{-1,-1,-1,-1}); //return if there is no previous snps to condition on
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
			//double[] twoLocusFreq = estimateTwoLocusFreqPhased(prevGenotypes, currGenotypes,prevAlleleTypes[0], prevAlleleTypes[1], currAlleleTypes[0], currAlleleTypes[1]);

			
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
		writeResults(writer, currPos, bestPos, bestPrevAlleleTypes, currAlleleTypes, currStationaryDist, bestTwoLocusFreq);

		
	}
	
	
	private static double[] estimateTwoLocusFreqPhased(String[] genotypesAtLocus1, String[] genotypesAtLocus2, char U, char u, char V, char v){
	
		double[] toReturn = new double[4];
		int numIndiv = genotypesAtLocus1.length;
		
		String A = ""+u+v;
		String B = ""+u+V;
		String C = ""+U+v;
		String D = ""+U+V;
		
		for(int i=0; i<numIndiv; i++){
			String h1 = ""+genotypesAtLocus1[i].charAt(0) + genotypesAtLocus2[i].charAt(0);
			String h2 = ""+genotypesAtLocus1[i].charAt(1) + genotypesAtLocus2[i].charAt(1);
			
			if(h1.equals(A)) toReturn[0]++;
			else if(h1.equals(B)) toReturn[1]++;
			else if(h1.equals(C)) toReturn[2]++;
			else if(h1.equals(D)) toReturn[3]++;
			else throw new RuntimeException("Illegal haplotype!");
			
			if(h2.equals(A)) toReturn[0]++;
			else if(h2.equals(B)) toReturn[1]++;
			else if(h2.equals(C)) toReturn[2]++;
			else if(h2.equals(D)) toReturn[3]++;
			else throw new RuntimeException("Illegal haplotype!");
			
		}
		
		for(int i=0; i<toReturn.length; i++){
			toReturn[i] = toReturn[i]/(numIndiv*2.0);
		}
		
		if (Math.abs(toReturn[0]+toReturn[1]+toReturn[2]+toReturn[3]-1d)>1e-12){
			System.out.println("TWO LOCUS PHASE WRONG");
		}
		
		return toReturn;
	
	}
	
	
	private static double[] estimateTwoLocusFreq(Cubic solver, String[] genotypesAtLocus1, String[] genotypesAtLocus2, char U, char u, char V, char v){


		int numGeneCopies = 2 * genotypesAtLocus1.length;

		int[] genotypeCounts = countGenotypes(genotypesAtLocus1, genotypesAtLocus2, U, u, V, v);

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
	
	
	private static void writeResults(PrintWriter writer, int currPos, int ldPos, char[] Uu, char[] Vv, double[] freq, double[] twoLocusFreq){
		
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
		for(double f : freq){
			writer.write(String.format("%f\t",f));
		}
		for(double f : twoLocusFreq){
			writer.write(String.format("%f\t",f));
		}
		
		writer.write("\n");
		
	}

	
	private static double[] estimateStationaryDist (String[] genotypes){

		
		int numIndiv = genotypes.length;
		int numAlleleCopies = 2*numIndiv;
		double[] toReturn = new double[4];
		
		// count alleles
		for (int ind=0; ind<numIndiv; ind++){
			char g1 = genotypes[ind].charAt(0);
			char g2 = genotypes[ind].charAt(1);
				
			if (g1=='A') toReturn[0]++;
			else if (g1=='C') toReturn[1]++;
			else if (g1=='G') toReturn[2]++;
			else if (g1=='T') toReturn[3]++;
			else throw new RuntimeException("Illegal allele in LD");
				
			if (g2=='A') toReturn[0]++;
			else if (g2=='C') toReturn[1]++;
			else if (g2=='G') toReturn[2]++;
			else if (g2=='T') toReturn[3]++;	
			else throw new RuntimeException("Illegal allele in LD");
			}

		// divide by total number of allele copies
		for (int allele=0; allele<4; allele++){
			toReturn[allele] = toReturn[allele]/numAlleleCopies;
		}

		
		return toReturn;
		
	}
	
	
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
		
		/*
		//sanity check
		int sum=0;
		for (int c : genotypeCounts){
			sum += c;
		}
		assert sum == numIndiv;
		 */

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
	
	
	private static double[] getMarginals (double A, double B, double C, double D){ //probabiltiy of observing a particular genotype
		
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
	
	
	private static double computeLD (double pA, double pB, double pAB){
		
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
		else throw new RuntimeException("Illegal allele");
	}
	


	
}
