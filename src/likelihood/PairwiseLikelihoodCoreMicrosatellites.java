package likelihood;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import utility.DataParser;
import dataStructures.Node;
import dataStructures.PairwiseKey;
import dataStructures.Path;
import dataStructures.Relationship;


public class PairwiseLikelihoodCoreMicrosatellites {
	

	//pairwise likelihood file columns
	private static final int UP = 1;
	private static final int DOWN = 2;
	private static final int NUMVISIT = 3;
	
	// model parameters and data structures
	private int numIndiv;
	private double epsilon1;
	private double epsilon2;
	
	//output
	private final Map<Path, double[][]> pairwiseLkhd = new HashMap<Path, double[][]>(); //path --> pairwise likelihood table
	private final List<Double> marginals = new ArrayList<Double>(); // index --> marginal likelihood
	
	//helper data structures
	private final Map<Integer, String[]> possibleGenotypes = new HashMap<Integer, String[]>(); //number of alleles --> possible genotypes
	private PairwiseKey pairwiseKey = new PairwiseKey("11","11","11","11"); // (geno1, geno2)

	//TODO different error rates at each locus
	//TODO filter out any cases of 01 or 10
	public PairwiseLikelihoodCoreMicrosatellites(int numIndiv, double epsilon1, double epsilon2){
		

		this.numIndiv = numIndiv;
		this.epsilon1 = epsilon1;
		this.epsilon2 = epsilon2;
		
		
		
		
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

		
		while((line=reader.readLine())!=null){

			
			String[] fields = line.split("\\s+");
			
			//update relationship
			if(fields[0].equals(">")){
				
				int up = Integer.parseInt(fields[UP]);
				int down = Integer.parseInt(fields[DOWN]);
				int nVisit = Integer.parseInt(fields[NUMVISIT]);
				
				double[][] lkhdTable = new double[numIndiv][numIndiv];
				path = new Path(up, down, nVisit);

				pairwiseLkhd.put(path, lkhdTable);

				continue;
				
			}

			//read likelihood
			int i = Integer.parseInt(fields[0]);
			int j = Integer.parseInt(fields[1]);
			double lkhd = Double.parseDouble(fields[2]);

			if(i<numIndiv && j<numIndiv){
				pairwiseLkhd.get(path)[i][j] = lkhd;
				pairwiseLkhd.get(path)[j][i] = lkhd;
			}
	
			

		}
		
		//close file
		reader.close();
		
	}

	
	
	public void setMarginals(String filePath) throws IOException{ //works
		
		//clear old values
		marginals.clear();
		
		//open file
		BufferedReader reader = DataParser.openReader(filePath);

		int i = 0;
		String line;
		while((line = reader.readLine())!=null){
			
			String fields[] = line.split("\n");
			
			marginals.add(i , Double.parseDouble(fields[0]));
			i++;
			
		}
		
		//close file
		reader.close();
		
	}

	
	
	

	//////////////////////////// PRECOMPUTE LIKELIHOOD ///////////////////////////////////	
	//marginal probability of an individual
	//missing allele is encoded as 0; but only one missing allele is not allowed (i.e. 0x, x0)
	public double[] computeMarginalMicrosat(String genoPath, String freqPath, int[] ids) throws IOException{ //works

		//data structures 
		double[] toReturn = new double[numIndiv];
		Map<String, Double> genotypeProbWithError = new HashMap<String, Double>(); 
		Map<PairwiseKey, Double> errorProb = new HashMap<PairwiseKey, Double>();

		//open files
		BufferedReader genoFile = DataParser.openReader(genoPath); //tped file
		BufferedReader freqFile = DataParser.openReader(freqPath); //holds frequencies of each allele
			
		String genoLine;
		String freqLine;
		
		while((genoLine=genoFile.readLine())!=null && (freqLine=freqFile.readLine())!=null){
		
			//read fields
			String[] geno = genoLine.split("\\s");
			String[] freqString = freqLine.split("\\s");
					
			//compute info
			int k = freqString.length;
			String[] possibleGenotypes = getPossibleGenotypes(k);
			double[] freq = getOneLocusFreq(freqString);
			double e1 = epsilon1 / (1+epsilon1);
			double e2 = epsilon2 / (k-1);
			
			//clear map
			genotypeProbWithError.clear();
			errorProb.clear();
			
			//compute the probability of the genotype for all individuals
			for (int i=0; i<numIndiv; i++){
				
				//order g so that the lower number allele comes first
				String[] orderedAlleles = orderAlleles(geno[2*ids[i]+4], geno[2*ids[i]+5]);
				
				//skip missing data
				if(orderedAlleles[0].equals("0") && orderedAlleles[1].equals("0")) continue;
				
				toReturn[i] += Math.log(getGenotypeProbWithError(orderedAlleles[0], orderedAlleles[1], genotypeProbWithError, freq, possibleGenotypes, e1, e2, epsilon2, errorProb));
				
			}
		
			
		}
		
		
		return toReturn;
		
	}
	

	
	
	//pairwise likelihood for independent sites
	public double[][][] computePairwiseMicrosat(String genoPath, String freqPath, int[] ids, List<Relationship> rel) throws IOException{
		
		int numRel = rel.size();
		double[][][] toReturn = new double[numRel][numIndiv][numIndiv];

		//(geno1, geno2)--> probability of [ (geno1, geno2 | ibd=0), (geno1, geno2 | ibd=1), (geno1, geno2 | ibd=2) ] 
		Map<PairwiseKey, double[]> pairwiseGenotypeProbWithError = new HashMap<PairwiseKey, double[]>(); 
		Map<PairwiseKey, Double> errorProb = new HashMap<PairwiseKey, Double>(); 

		//open files
		BufferedReader genoFile = DataParser.openReader(genoPath); //tped file
		BufferedReader freqFile = DataParser.openReader(freqPath); //holds frequencies of each allele
			
		String genoLine;
		String freqLine;
		
		
		while((genoLine=genoFile.readLine())!=null && (freqLine=freqFile.readLine())!=null){
		
			//read fields
			String[] geno = genoLine.split("\\s");
			String[] freqString = freqLine.split("\\s");
			
			
			//compute info
			int k = freqString.length;
			String[] possibleGenotypes = getPossibleGenotypes(k);
			double[] freq = getOneLocusFreq(freqString);
			double e1 = epsilon1 / (1+epsilon1);
			double e2 = epsilon2 / (k-1);
			
			
			//clear map
			pairwiseGenotypeProbWithError.clear();
			errorProb.clear();
			
			
			//for every pair	
			for(int i=0; i<numIndiv; i++){
				
				String obs1_ind1 = geno[2*ids[i]+4];
				String obs2_ind1 = geno[2*ids[i]+5];
				
				//skip missing data
				if(obs1_ind1.equals("0") && obs2_ind1.equals("0")) continue;
				
				for(int j=i+1; j<numIndiv; j++){
					
					String obs1_ind2 = geno[2*ids[j]+4];
					String obs2_ind2 = geno[2*ids[j]+5];
					
					if(obs1_ind2.equals("0") && obs2_ind2.equals("0")) continue;
					
					//get pairwise genotype probability
					double[] pairwiseProb = getPairwiseGenotypeProbWithError(obs1_ind1, obs2_ind1, obs1_ind2, obs2_ind2, pairwiseGenotypeProbWithError, freq, possibleGenotypes, e1, e2, epsilon2, errorProb);
					
					for(int r=0; r<numRel; r++){
						
						toReturn[r][i][j] += Math.log(rel.get(r).getMarginalProbs()[0] * pairwiseProb[0] + rel.get(r).getMarginalProbs()[1] * pairwiseProb[1] + rel.get(r).getMarginalProbs()[2] * pairwiseProb[2]);
					
					}
					
				}
			
			}
		
		}
			

		
		return toReturn;
		

	}	

	
	
	
	
	////////// PAIRWISE GENOTOYPE PROBABILITY /////////
	private double[] getPairwiseGenotypeProbWithError(String obs1_ind1, String obs2_ind1, String obs1_ind2, String obs2_ind2, Map<PairwiseKey, double[]> pairwiseGenotypeProbWithError, double[] freq, String[] possibleGenotypes, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb){	
		
		this.pairwiseKey.setKey(obs1_ind1, obs2_ind1, obs1_ind2, obs2_ind2);
		
		//if already in table, return stored value
		if(pairwiseGenotypeProbWithError.containsKey(pairwiseKey))
			return pairwiseGenotypeProbWithError.get(pairwiseKey);
		
	
		// compute, store, return
		double[] toReturn = new double[3];
		
		for(int ibd=0; ibd<3; ibd++){
			
			toReturn[ibd] = pairwiseGenotypeProbWithError(freq, obs1_ind1, obs2_ind1, obs1_ind2, obs2_ind2, ibd, e1, e2, epsilon2, errorProb);
			
			
		}
		
		pairwiseGenotypeProbWithError.put(new PairwiseKey(obs1_ind1, obs2_ind1, obs1_ind2, obs2_ind2), toReturn);
	
		return toReturn;
		
		
		
	}
     	
	
	
	// P(obs1, obs2 | error) = sum over all true genotypes { P(true1, true2| ibd) P(obs1| true1) P(obs2 | true2) }
	private double pairwiseGenotypeProbWithError(double[] freq, String obs1_ind1, String obs2_ind1, String obs1_ind2, String obs2_ind2, int ibd, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb){
		
		double toReturn = 0d;
		int k = freq.length;
		
		 //for all possible genotypes of ind1
		for(int i=1; i<=k; i++) {
			
			String true1_ind1 = i+""; //first allele
			
			
			for(int j=i; j<=k; j++) {
				
				String true2_ind1 = j+""; //second allele
				
				double error1 = getErrorProbability(true1_ind1, true2_ind1, obs1_ind1, obs2_ind1, e1, e2, epsilon2, errorProb);
				
				
				//for all possible genotyeps of ind2
				for(int i2=1; i2<=k; i2++) {
					
					String true1_ind2 = i2+""; //first allele
					
					
					for(int j2=i2; j2<=k; j2++) {
						
						String true2_ind2 = j2+""; //second allele
						
						double error2 = getErrorProbability(true1_ind2, true2_ind2, obs1_ind2, obs2_ind2, e1, e2, epsilon2, errorProb);
						
						toReturn += error1 * error2 * pairwiseGenotypeProbWithoutError(freq, true1_ind1, true2_ind1, true1_ind2, true2_ind2, ibd);
						
						
					}
					
					
				}
				
				
				
			}
			
			
		}
		
		
	
		return toReturn;
	}
	
	
	
	//no marginalizing assume genotype1 and genotype2 are the real deal; individuals are ordered
	private double pairwiseGenotypeProbWithoutError(double[] oneLocusFreq, String a1, String a2, String b1, String b2, int ibd){//works
		
		int ibs = getIBS(a1, a2, b1, b2);
		
		if( ibs < ibd){
			return 0d;
		}
		
		double delta1 = a1.equals(a2) ? 1d : 0d;
		double delta2 = b1.equals(b2) ? 1d : 0d;	//we need these deltas for the combinatorial coefficients in a lot of these calculations
		
		if(ibd == 2){	//correct
			return (2d - delta1) * getAlleleProb(oneLocusFreq, a1) * getAlleleProb(oneLocusFreq, a2);
		}
		if(ibd == 1){
			//need to do some fancy stuff to deal with this case
			if(ibs == 2){	//correct
				if(delta1==1){
					return Math.pow(getAlleleProb(oneLocusFreq, a1),3);
				}
				else{
					return Math.pow(getAlleleProb(oneLocusFreq, a1),2) * getAlleleProb(oneLocusFreq, a2) + Math.pow(getAlleleProb(oneLocusFreq, a2), 2) * getAlleleProb(oneLocusFreq, a1);	
				}
			}
			if(ibs == 1){ 	//correct
				String sharedAllele;
				String loneAllele1;
				String loneAllele2;
				if(a1.equals(b1)){
					sharedAllele = a1;
					loneAllele1 = a2;
					loneAllele2 = b2;
				} else if(a1.equals(b2)){
					sharedAllele = a1;
					loneAllele1 = a2;
					loneAllele2 = b1;
				} else if(a2.equals(b2)){
					sharedAllele = a2;
					loneAllele1 = a1;
					loneAllele2 = b1;
				} else if(a2.equals(b1)){
					sharedAllele = a2;
					loneAllele1 = a1;
					loneAllele2 = b2;
				} else {
					throw new RuntimeException("Something bad happened while computing an emission prob");
				}
				return getAlleleProb(oneLocusFreq, sharedAllele) * getAlleleProb(oneLocusFreq, loneAllele1) * getAlleleProb(oneLocusFreq, loneAllele2);	//combinatorial factors should exactly cancel
			}
			throw new RuntimeException("ibd didn't match IBS!");
		}
		if(ibd == 0){	//correct
			return (2d - delta1) * (2d - delta2) * getAlleleProb(oneLocusFreq, a1) * getAlleleProb(oneLocusFreq, a2) * getAlleleProb(oneLocusFreq, b1) * getAlleleProb(oneLocusFreq, b2);
		}
		throw new RuntimeException("Illegal IBD value encountered");
	}
	
	
	
	
	
	/////// GENOTYPE PROBABILITY //////////
	//assume uniform error rate for each locus
	private double getGenotypeProbWithError(String a1, String a2, Map<String, Double> genotypeProbWithError, double[] freq, String[] possibleGenotypes, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb){ //works
		
		//if this genotype prob has already been computed, return stored value
		if(genotypeProbWithError.containsKey(a1+a2))
			return genotypeProbWithError.get(a1+a2);
		
		//else compute, and store, and return
		double toReturn = genotypeProbWithError(a1, a2, freq, e1, e2, epsilon2, errorProb);
		
		genotypeProbWithError.put(a1+a2, toReturn);

		
		return toReturn;
		
	}
	
	
	
	// P(Obs | error) = sum over P(Obs | True, error) x P(True | error)
	private double genotypeProbWithError(String obs1, String obs2, double[] freq, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb){ //works
		
		double toReturn = 0d;
		
		//iterate through all possible true genotypes
		for(int i=1; i<=freq.length; i++) { //first true allele
			
			String true1 = ""+i;
			
			for(int j=i; j<=freq.length; j++) { //second true allele
				
				String true2 = ""+j;
				
				toReturn += getErrorProbability(true1, true2, obs1, obs2, e1, e2, epsilon2, errorProb) *  genotypeProbWithoutError(freq, true1, true2);
				
				

				
			}
			
		}

		return toReturn;
	}
	
	
	
	//two alleles given by a1 and a2
	private double genotypeProbWithoutError(double[] freq, String a1, String a2){//works
		
		double delta = a1.equals(a2)? 1:0;
		double lkhd = (2-delta) * getAlleleProb(freq, a1) * getAlleleProb(freq, a2);

		return lkhd;
		
	}
	
	


	///////// ERROR PROBABILITY ///////////////////
	private double getErrorProbability(String true1, String true2, String obs1, String obs2, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb) {
	
		this.pairwiseKey.setKey(true1, true2, obs1, obs2);
		
		//if this was already computed
		if(errorProb.containsKey(pairwiseKey))
			return errorProb.get(this.pairwiseKey);
		
		PairwiseKey key = new PairwiseKey(true1, true2, obs1, obs2);
		double val = errorProbability(true1, true2, obs1, obs2, e1, e2, epsilon2);
		
		
		errorProb.put(key, val);
		
		return val;
		
		
	}
	
	
	
	
	//likelihood of the observedGenotype given the true genotype (model by Wang (2003))
	//when dropout happens, heterozygotes looks like homozygotes. e.g. 12-->11
	//sequencing error happens after dropout in a normal way
	//true genotypes given by (true1, true2). Observed is (obs1, obs2)
	private double errorProbability(String true1, String true2, String obs1, String obs2, double e1, double e2, double epsilon2){ //works

		double obsHomo = obs1.equals(obs2) ? 1d : 0d;
		double trueHomo = true1.equals(true2) ? 1d : 0d;
		int ibs = getIBS(true1, true2, obs1, obs2);
		

		//no gene drop outs
		if(trueHomo==1) {
			
			//obs is a hoomzygote
			if(obsHomo==1) {
				
				if(ibs==2)
					return Math.pow(1-epsilon2, 2);
				else if(ibs==0)
					return Math.pow(e2, 2);
				else
					throw new RuntimeException("This case is impossible!");
				
			}
			
			//obs is heterozygous
			else {
				
				if(ibs==1)
					return 2*e2 * (1-epsilon2);
				
				else if(ibs==0)
					return 2*Math.pow(e2, 2);
				
				else
					throw new RuntimeException("This case is impossible!");
				
			}
			
		}
		
		
		
		//true genotype is a heterozygote
		else {
				
			if(ibs==2)
				return Math.pow(1-epsilon2, 2) + Math.pow(e2,2) - 2*e1*Math.pow(1-epsilon2-e2, 2);
			
			else if(obsHomo==1 && ibs==1)
				return e2*(1-epsilon2) + e1*Math.pow(1-epsilon2-e2, 2);
			
			else if(ibs==0)
				return (2d - obsHomo) * Math.pow(e2, 2);
					
			else 
				return e2*(1-epsilon2+e2);			
			
		
		}
		
		
		
		
	}

	
	
	
	
	//////// MISCELLANEOUS HELPER FUNCTIONS ////////	
	private double getAlleleProb(double[] oneLocusFreq, String allele){ //works
		
		
		int idx = Integer.parseInt(allele) - 1; //convert to int, start indexing from 0
		
		if(idx<0 || idx >= oneLocusFreq.length)
			throw new RuntimeException("Frequency not given for this allele");
		
		return oneLocusFreq[idx];
		
	}
	
	
	//".1 .5 .4" --> {.1, .5, .4}
	private double[] getOneLocusFreq(String[] freq){ //works
		
		double[] oneLocusFreq = new double[freq.length];
		
		
		for(int i=0; i<freq.length; i++) {
			
			oneLocusFreq[i] = Double.parseDouble(freq[i]);
			
		}

		
		return oneLocusFreq;
		
	}
	

	//returns the number of characters the two strings have in common (unorders them, first)
	//genotype of first indiv: a1, a2. second indiv: b1,b2
	private int getIBS(String a1, String a2, String b1, String b2){ //works
		
		if(a1.equals(b1)){
			if(a2.equals(b2)){
				return 2;	//strings are the same
			}
			return 1;	//first character are the same, second characters differ
		} 

		if(a1.equals(b2)){
			if(a2.equals(b1)){
				return 2;	//strings are the same but reversed
			}
			return 1;	//first character in one is last character in the other
		}
		
		//first character doesn't equal either character in the other
		//if the second character equals either of the characters in the other string, ibs = 1 (since first character doesn't match), else 0
		return a2.equals(b1) || a2.equals(b2)? 1 : 0;

		
		
	}

	
	//returns possible genotypes given number of alleles; 2 --> {11, 12, 22}
	private String[] getPossibleGenotypes(int k) { 
		
		//check if the map already exists
		if(possibleGenotypes.containsKey(k)) 
			return possibleGenotypes.get(k);
		
		
		//else, make a new map
		int nPossible = k*(k+1)/2;
		String[] toReturn = new String[nPossible];
		int q = 0;
		
		for(int i=1; i<=k; i++) {
			
			for(int j=i; j<=k; j++) {
				
				toReturn[q] = i+""+j;
				q++;
			}
			
		}
		
		this.possibleGenotypes.put(k, toReturn);
		
		return toReturn;
		
		
	}
	
	
	
	//return ordered genotype
	private String[] orderAlleles(String a1, String a2) {//works
		
		
		//homozygous
		if(a1.equals(a2))
			return new String[] {a1,a2};
		
		
		else if(Integer.parseInt(a1) < Integer.parseInt(a2))
			return new String[] {a1, a2};	
		
		else 
			return new String[] {a2, a1};
		
	}
	
	
	public static void main(String[] args) throws IOException {
		
		String dir = "/Users/amy/eclipse-workspace/mcmc/simulations/";
		int numIndiv = 2;
		double epsilon1 = .05; //dropout rate
		double epsilon2 = .02; //sequencing error


		//TESTING
		PairwiseLikelihoodCoreMicrosatellites core = new PairwiseLikelihoodCoreMicrosatellites(numIndiv, epsilon1, epsilon2);

		
		List<Relationship> rels = new ArrayList<Relationship>();
		rels.add(new Relationship(13d, new double[] {.5d, .5d, 0d}, new Path[]{new Path(0,0,0)}));

		double[][][] x = core.computePairwiseMicrosat(dir+"microsat.tped", dir+"microsat.freq", new int[] {0,1}, rels);
		//double[] x = core.computeMarginalMicrosat(dir+"microsat.tped", dir+"microsat.freq", new int[] {0,1});
		
		//System.out.println(x[0]);
		System.out.println(x[0][0][1]);
		//System.out.println(x[0][0][2]);
		 
		 
	
		
	}

	 
}