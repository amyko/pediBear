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
	PairwiseKey pairwiseKey = new PairwiseKey("11","11","11","11");
	Map<String, Double> allele2freq = new HashMap<String, Double>();
	Map<String, Integer> allele2count = new HashMap<String, Integer>();
	
	
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
	public double[] computeMarginalMicrosat(String genoPath, String freqPath, String fPath, int[] ids, boolean withInbreeding, boolean adjust_af) throws IOException{ //works

		//data structures 
		double[] toReturn = new double[numIndiv];
		Map<PairwiseKey, Double> errorProb = new HashMap<PairwiseKey, Double>();
		int n_alleles = 2*numIndiv;

		//open files
		BufferedReader genoFile = DataParser.openReader(genoPath); //tped file
		BufferedReader freqFile = DataParser.openReader(freqPath); //holds frequencies of each allele
		
		// get inbreeding coefficient
		double[] f = new double[numIndiv];
		if(withInbreeding)
			getF(fPath, f);
			
		String genoLine;
		
		
		while((genoLine=genoFile.readLine())!=null){
		
			//read fields
			String[] geno = genoLine.split("\\s");
			String[] alleles = freqFile.readLine().split("\\s");
			String[] freqs = freqFile.readLine().split("\\s");
			
			// init frequencies/counts
			initOneLocusCount(alleles, freqs, n_alleles);
			initOneLocusFreq(alleles, freqs);
					
			//compute info
			int k = alleles.length;
			double e1 = epsilon1 / (1+epsilon1);
			double e2 = epsilon2 / (k-1);
			
			//clear map
			errorProb.clear();
			
			//compute the probability of the genotype for all individuals
			for (int i=0; i<numIndiv; i++){
				
				//order g so that the lower number allele comes first
				String[] orderedAlleles = orderAlleles(geno[2*ids[i]+4], geno[2*ids[i]+5]);
				
				//skip missing data
				if(orderedAlleles[0].equals("0") || orderedAlleles[1].equals("0")) continue;
				
				// adjust freq
				if(adjust_af) {	
					updateOneLocusCount(orderedAlleles[0], orderedAlleles[1], -1);
					updateOneLocusFreq();
				}
		
	
				toReturn[i] += Math.log(genotypeProbWithError(orderedAlleles[0], orderedAlleles[1], f[i], allele2freq, e1, e2, epsilon2, errorProb, alleles));
	
				
				// re-adjust freq
				if(adjust_af)
					updateOneLocusCount(orderedAlleles[0], orderedAlleles[1], 1);
	
				
			}
		
			
		}
		
		
		return toReturn;
		
	}
	

	
	
	//pairwise likelihood for independent sites
	public double[][][] computePairwiseMicrosat(String genoPath, String freqPath, String fPath, int[] ids, List<Relationship> rel, boolean withInbreeding, boolean adjust_af) throws IOException{
		
		int numRel = rel.size();
		int n_alleles = 2 * numIndiv;
		double[][][] toReturn = new double[numRel][numIndiv][numIndiv];


		Map<PairwiseKey, Double> errorProb = new HashMap<PairwiseKey, Double>(); 

		//read f for each individual
		double[] f = new double[ids.length];
		if(withInbreeding)
			getF(fPath, f);
		
		
		//open files
		BufferedReader genoFile = DataParser.openReader(genoPath); //tped file
		BufferedReader freqFile = DataParser.openReader(freqPath); //holds frequencies of each allele
			
		String genoLine;

		
		while((genoLine=genoFile.readLine())!=null){
		
			//read fields
			String[] geno = genoLine.split("\\s");
			String[] alleles = freqFile.readLine().split("\\s");
			String[] freqs = freqFile.readLine().split("\\s");
			
			// init frequencies/counts
			initOneLocusCount(alleles, freqs, n_alleles);
			initOneLocusFreq(alleles, freqs);
			
			
			//compute info
			int k = alleles.length;
			double e1 = epsilon1 / (1+epsilon1);
			double e2 = epsilon2 / (k-1);
			
			
			//clear map
			errorProb.clear();
			
			
			//for every pair	
			for(int i=0; i<numIndiv; i++){
				
				String obs1_ind1 = geno[2*ids[i]+4];
				String obs2_ind1 = geno[2*ids[i]+5];
				double f1 = f[i];
				
				//skip missing data
				if(obs1_ind1.equals("0") || obs2_ind1.equals("0")) continue;
				
				
				// adjust count
				if(adjust_af) 	
					updateOneLocusCount(obs1_ind1, obs2_ind1, -1);			
				
				
				for(int j=i+1; j<numIndiv; j++){
					
					String obs1_ind2 = geno[2*ids[j]+4];
					String obs2_ind2 = geno[2*ids[j]+5];
					double f2 = f[j];
					
					if(obs1_ind2.equals("0") || obs2_ind2.equals("0")) continue;
					
					// adjust count and freq
					if(adjust_af) {	
						updateOneLocusCount(obs1_ind2, obs2_ind2, -1);
						updateOneLocusFreq();
					}
					
					
					for(int r=0; r<numRel; r++){
						
						double p0 = rel.get(r).getMarginalProbs()[0];
						double p1 = rel.get(r).getMarginalProbs()[1];
						double p2 = rel.get(r).getMarginalProbs()[2];
						
						//Pr(s)
						double s0 = (f1*(1-f2) + (1-f1)*f2 + f1*f2)*p2 + f1*f2*p1;
						double s1 = f1*f2*p0;
						double s2 = f1*(1-f2)*p1;
						double s3 = f1*(1-f2)*p0;
						double s4 = (1-f1)*f2*p1;
						double s5 = (1-f1)*f2*p0;
						double s6 = (1-f1)*(1-f2)*p2;
						double s7 = (1-f1)*(1-f2)*p1;
						double s8 = (1-f1)*(1-f2)*p0;
						double[] delta = new double[] {s0, s1, s2, s3, s4, s5, s6, s7, s8};
	
						
						toReturn[r][i][j] += Math.log(pairwiseGenotypeProbWithError(allele2freq, obs1_ind1, obs2_ind1, obs1_ind2, obs2_ind2, delta, e1, e2, epsilon2, errorProb, alleles));
						
						
					}
					
					
					// re-adjust freq for indiv j
					if(adjust_af)
						updateOneLocusCount(obs1_ind2, obs2_ind2, 1);
					
				}
				
				// re-adjust freq for indiv i
				if(adjust_af)
					updateOneLocusCount(obs1_ind1, obs2_ind1, 1);
			
			}
		
		}
			

		
		return toReturn;
		

	}	


	
	
	////////// PAIRWISE GENOTOYPE PROBABILITY /////////	
	// P(obs1, obs2 | rel) = sum over all true genotypes { P(true1, true2| rel) P(obs1| true1) P(obs2 | true2) }
	private double pairwiseGenotypeProbWithError(Map<String, Double> allele2freq, String obs1_ind1, String obs2_ind1, String obs1_ind2, String obs2_ind2, double[] delta, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb, String[] alleles){
		
		double toReturn = 0d;

		for(int i=0; i<alleles.length; i++) {
			
			String true1_ind1 = alleles[i];
			
			for(int j=i; j<alleles.length; j++) {
			
				String true2_ind1 = alleles[j];
				double error1 = getErrorProbability(true1_ind1, true2_ind1, obs1_ind1, obs2_ind1, e1, e2, epsilon2, errorProb);
			
				for(int k=0; k<alleles.length; k++) {
					
					String true1_ind2 = alleles[k];
					
					for(int l=k; l<alleles.length; l++) {
						
						String true2_ind2 = alleles[l];	
						double error2 = getErrorProbability(true1_ind2, true2_ind2, obs1_ind2, obs2_ind2, e1, e2, epsilon2, errorProb);
					
						toReturn += error1 * error2 * pairwiseGenotypeProbWithoutErrorGivenRel(allele2freq, true1_ind1, true2_ind1, true1_ind2, true2_ind2, delta);
					}
				}
			}
		}

		return toReturn;
	}
	
	
	
	// Pr(geno1, gene2 | ibd = 0 or 1 or 2)
	
	
	private double pairwiseGenotypeProbWithoutErrorGivenRel(Map<String, Double> allele2freq, String a1, String a2, String b1, String b2, double[] delta) {
				
		double prob = 0;
		double[] pairwiseGenotypeProbWithoutErrorGivenS = pairwiseGenotypeProbWithoutErrorGivenS(allele2freq, a1, a2, b1, b2);
		
		for(int i=0; i<delta.length; i++) 
			prob += delta[i] * pairwiseGenotypeProbWithoutErrorGivenS[i];
		
		
		return(prob);
		

	}
	
	
	// delta[i] = Pr(geno1, geno2 | s = i)
	private double[] pairwiseGenotypeProbWithoutErrorGivenS(Map<String, Double> allele2freq, String a1, String a2, String b1, String b2) {
		
		double ibs = getIBS(a1, a2, b1, b2);
		double[] delta = new double[9];

		if(ibs == 2) {
			
			double pi = getAlleleProb(allele2freq, a1);
			
			//1
			if(a1.equals(a2)) {
				
				delta[0] = pi;
				delta[1] = delta[2] = delta[4] = delta[6] = Math.pow(pi, 2);
				delta[3] = delta[5] = delta[7] = Math.pow(pi, 3);
				delta[8] = Math.pow(pi, 4);

			}
			
			//5			
			else {
				double pj = getAlleleProb(allele2freq, a2);
				
				delta[6] = 2*pi*pj;
				delta[7] = pi*pj*(pi+pj);
				delta[8] = 4*Math.pow(pi, 2)*Math.pow(pj, 2);
			}
			
		}
		
		else if(ibs == 0) {
			
			//2
			if(a1.equals(a2) && b1.equals(b2)) {
				double pi = getAlleleProb(allele2freq, a1);
				double pj = getAlleleProb(allele2freq, b1);
				
				delta[1] = pi*pj;
				delta[3] = pi*Math.pow(pj, 2);
				delta[5] = Math.pow(pi, 2)*pj;
				delta[8] = Math.pow(pi, 2)*Math.pow(pj, 2);

			}				
			
			//7
			else if(!a1.equals(a2) && !b1.equals(b2)) {
				
				delta[8] = 4*getAlleleProb(allele2freq, a1) * getAlleleProb(allele2freq, a2)*getAlleleProb(allele2freq, b1)*getAlleleProb(allele2freq, b2);

			}
			//4
			else {
				
				double pi;
				double pj;
				double pm;
				
				if(a1.equals(a2)) {
					pi = getAlleleProb(allele2freq, a1);
					pj = getAlleleProb(allele2freq, b1);
					pm = getAlleleProb(allele2freq, b2);
					
					delta[3] = 2*pi*pj*pm;
					delta[8] = 2*Math.pow(pi, 2)*pj*pm;

				}
				else {
					pi = getAlleleProb(allele2freq, b1);
					pj = getAlleleProb(allele2freq, a1);
					pm = getAlleleProb(allele2freq, a2);
					
					delta[5] = 2*pi*pj*pm;
					delta[8] = 2*Math.pow(pi, 2)*pj*pm;

				}
				
				
			}
			
		}
		
		
		else {
			
			//6
			if(!a1.equals(a2) && !b1.equals(b2)) {
				
				
				double pi;
				double pj;
				double pm;
				
				if(a1.equals(b1) || a1.equals(b2)) {
					pi = getAlleleProb(allele2freq, a1);
					pj = getAlleleProb(allele2freq, a2);
				}
				else {
					pi = getAlleleProb(allele2freq, a2);
					pj = getAlleleProb(allele2freq, a1);
				}
				
				if(b1.equals(a1) || b1.equals(a2))
					pm = getAlleleProb(allele2freq, b2);
				else
					pm = getAlleleProb(allele2freq, b1);
								
				delta[7] = pi*pj*pm;
				delta[8] = 4*Math.pow(pi, 2)*pj*pm;

				
			}
			
			//3
			else{
				
				double pi;
				double pj;
			
				
				if(!a1.equals(a2)) {
					if(a1.equals(b1)){
						pi = getAlleleProb(allele2freq, a1);
						pj = getAlleleProb(allele2freq, a2);
					}
					else {
						pi = getAlleleProb(allele2freq, a2);
						pj = getAlleleProb(allele2freq, a1);
					}
					
					delta[4] = pi*pj;
					delta[5] = 2*Math.pow(pi, 2)*pj;
					delta[7] = Math.pow(pi, 2)*pj;
					delta[8] = 2*Math.pow(pi, 3)*pj;
					
				}
				else {
					if(b1.equals(a1)) {
						pi = getAlleleProb(allele2freq, b1);
						pj = getAlleleProb(allele2freq, b2);
					}
					else {
						pi = getAlleleProb(allele2freq, b2);
						pj = getAlleleProb(allele2freq, b1);
					}
					
					delta[2] = pi*pj;
					delta[3] = 2*Math.pow(pi, 2)*pj;
					delta[7] = Math.pow(pi, 2)*pj;
					delta[8] = 2*Math.pow(pi, 3)*pj;

				}

				
			}
			
		}
		
		return delta;
		
	}	
		
	
	
	
	/////// GENOTYPE PROBABILITY //////////
	// P(Obs | error) = sum over P(Obs | True, error) x P(True | error)
	private double genotypeProbWithError(String obs1, String obs2, double f, Map<String, Double> allele2freq, double e1, double e2, double epsilon2, Map<PairwiseKey, Double> errorProb, String[] alleles){ //works
		
		double toReturn = 0d;
		
		for(int i=0; i<alleles.length; i++) {

			String true1 = alleles[i];
			
			for(int j=i; j<alleles.length; j++) {
				
				String true2 = alleles[j];
			
			
				toReturn += getErrorProbability(true1, true2, obs1, obs2, e1, e2, epsilon2, errorProb) *  genotypeProbWithoutError(allele2freq, true1, true2, f);
		
			}
		}

		return toReturn;
	}
	
	
	

		//probability of (a1,a2) given allele frequencies and inbreeding coefficient f
	private double genotypeProbWithoutError(Map<String, Double> allele2freq, String a1, String a2, double f) {
		
		if(a1.equals(a2)) {
			double p = getAlleleProb(allele2freq, a1);
			return f*p + (1-f)*p*p;
		}
		
		else {
			return 2 * (1-f) * getAlleleProb(allele2freq, a1) * getAlleleProb(allele2freq, a2);
		}
		
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
	private double getAlleleProb(Map<String, Double> allele2freq, String allele){ //works
		
		if(!allele2freq.containsKey(allele))
			throw new RuntimeException("Frequency not given for this allele");
		
		return allele2freq.get(allele);
		
	}
	
	
	//init allele2count map: add 1 to each count
	private void initOneLocusCount(String[] alleles, String[] freqs, int n_alleles) {
		
		allele2count.clear();
		
		for(int i=0; i<alleles.length; i++) {
			int count = (int) (Double.parseDouble(freqs[i]) * n_alleles) + 1;
			allele2count.put(alleles[i], count);
		}
		
	}
	
	//initialize allele2freq map
	private void initOneLocusFreq(String[] alleles, String[] freqs){ //works
		
		allele2freq.clear();
		for(int i=0; i<alleles.length; i++) 
			allele2freq.put(alleles[i], Double.parseDouble(freqs[i]));
		
	}
	
	
	// subtract focal alleles from the count and update frequency
	private void updateOneLocusCount(String a, String b, int to_add) {
						
		// if the current count is one, don't change
		allele2count.put(a, allele2count.get(a) + to_add);
		allele2count.put(b, allele2count.get(b) + to_add);

	}
	
	private void updateOneLocusFreq() {
		
		//total count
		double total_count = 0;
		for(int val : allele2count.values()) {
			int c = val == 0 ? 1 : val;
			total_count += c;
		}
		
		
		// update frequency
		for(String key : allele2count.keySet()) {
			int count = allele2count.get(key) == 0? 1 : allele2count.get(key);
			allele2freq.put(key, count / total_count);
		}
		
	}
	
	
	//read in inbreeding coefficient for each individual
	private void getF(String fPath, double[] f) throws NumberFormatException, IOException {
		
		BufferedReader reader = DataParser.openReader(fPath);
		
		String line;
		int i = 0;
		
		while((line = reader.readLine()) != null) {
			f[i] = Double.parseDouble(line.split("\\s")[0]);
			i++;
		}
		
		reader.close();
		
		
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

	 
	
	/*
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
	
	 * 
	 */
	
	
	public static void main(String[] args) {
		
		
		PairwiseLikelihoodCoreMicrosatellites core = new PairwiseLikelihoodCoreMicrosatellites(2, .05, .02);
		
		double[] freq = new double[] {.1, .2, .3, .4};
		String a1 = "1";
		String a2 = "1";
		String b1 = "1";
		String b2 = "1";
		double f1 = .2;
		double f2 = .02;
		Map<PairwiseKey, Double> errorProb = new HashMap<PairwiseKey, Double>();
		double e1 = .01;
		double e2 = .02;
		double epsilon2 = .05;
		
		double p0 = .5;
		double p1 = .5;
		double p2 = 0;
		
		double s0 = (f1*(1-f2) + (1-f1)*f2 + f1*f2)*p2 + f1*f2*p1;
		double s1 = f1*f2*p0;
		double s2 = f1*(1-f2)*p1;
		double s3 = f1*(1-f2)*p0;
		double s4 = (1-f1)*f2*p1;
		double s5 = (1-f1)*f2*p0;
		double s6 = (1-f1)*(1-f2)*p2;
		double s7 = (1-f1)*(1-f2)*p1;
		double s8 = (1-f1)*(1-f2)*p0;
		double[] delta = new double[] {s0, s1, s2, s3, s4, s5, s6, s7, s8};


		
		
	}

}