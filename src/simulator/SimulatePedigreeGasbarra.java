package simulator;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import dataStructures.Node;
import utility.DataParser;

public class SimulatePedigreeGasbarra {
	
	//simulate pedigree for given samples
	public static int[][] simulatePedigreeForSamples(int N, int Nf, int Nm, int nGen, double alpha, double beta, Random random, String outPath, boolean[][] sampled, boolean allowCycles) throws IOException {
		
		//make graph
		List<List<Node>> nodes = new ArrayList<List<Node>>();
		int index = 0;
		for(int gen=0; gen<=nGen; gen++) {
			
			nodes.add(new ArrayList<Node>());
			
			for(int i=0; i<N; i++) {
				
				int sex = i<Nf ? 1 : 0;
				
				nodes.get(gen).add(new Node(""+gen, ""+i, sex, false, index));	
				index++;
			}
		
		}
		
		
		PrintWriter writer = DataParser.openWriter(String.format("%sgasbarra.fam", outPath, N, nGen));
		writer.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		
		int[][] assignment = new int[nGen][2*N];
		int[][][] Cfm  = new int[nGen][Nf][Nm];
		
		for(int gen=0; gen < nGen; gen++){ //generation by generation
			
			//int F = 0;
			//int M = 0;
			int[] Cf = new int[Nf];
			Set<Integer> oldF = new HashSet<Integer>();
			Set<Integer> oldM = new HashSet<Integer>();
			int k = 0; //number of children assigned so far

			
			//for each sample in this generation
			for(int ind=0; ind<N; ind++) {
				
				if(sampled[gen][ind]==false) continue;
				
				//get child
				Node child = nodes.get(gen).get(ind);
				
				
				//father assignment
				double u = random.nextDouble();
				double cdf = 0;
				int f = -1;
				
				for(f=0; f<Nf; f++) {//each candidate father
						
					if(oldF.contains(f))  //old father
						cdf += (alpha + Cf[f]) / (Nf*alpha + k);
					else //new father
						cdf += alpha / (Nf*alpha + k);
					
					
					//test for cycles
					Node father = nodes.get(gen+1).get(f);
					boolean cyclic = cyclic(nodes, father, child);
				
					//if(cyclic) System.out.println("Father cyclic");
		
					//assign & update values
					if(u < cdf && (!cyclic || allowCycles)) {
						
						assignment[gen][2*ind] = f;
						//if(f >= F) F++;
						oldF.add(f);
						
						connect(father, child);
	
						break;
					}
					
					if(f==Nf-1) {
						System.out.println("Father not assigned. Try again.");
						u = random.nextDouble();
						f = 0;
					}
					
					
				}

				
				
				//mother assignment
				u = random.nextDouble();
				cdf = 0;
				int m = -1;
				
				for(m=0; m<Nm; m++) {
						
					
					if(oldM.contains(m)) // old mother
						cdf += (beta + Cfm[gen][f][m]) / (Nm*beta + Cf[f]);
					else // new mother
						cdf += beta / (Nm*beta + Cf[f]);
					
					
					
					//test for cycles
					Node mother = nodes.get(gen+1).get(m+Nm);
					boolean cyclic = cyclic(nodes, mother, child);
					
					//if(cyclic) System.out.println("Mother cyclic");
					
					//assign
					if(u < cdf && (!cyclic || allowCycles)) {
						assignment[gen][2*ind+1] = m;
						//if(m >= M) M++;
						oldM.add(m);
						
						connect(mother, child);
								
						break;
					}
					
					
					if(m==Nm-1) {
						System.out.println("Mother not assigned. Try again.");
						u = random.nextDouble();
						m = 0;
					}
					
					
				}
				

				
				
				//update children count & ghost parents
				Cf[f]++;
				Cfm[gen][f][m]++;
				sampled[gen+1][f] = true;
				sampled[gen+1][m+Nm] = true;
				k++;
				
				
				//write
				int sex = ind < Nf ? 1 : 2;
				String id = String.format("%d_%d", gen, ind);
				String fid = String.format("%d_%d", gen+1, f);
				String mid = String.format("%d_%d", gen+1, m+Nf);
				writer.write(String.format("%s\t%s\t%s\t%d\t%s\n", id, fid, mid, sex, "999999"));
				writer.flush();
				
						
			}
			
			
		}
		
		writer.close();
		return assignment;
		
	}
	
	

	public static boolean cyclic(List<List<Node>> nodes, Node parent, Node child) {
		
		performDFS(nodes, child);
		
		//node2 can be reached from node1
		if(parent.getNumVisit() > 0){ //okay only if merging creates FS
			
			if(formsFullSibs(parent, child))
				return false;
			else
				return true;
		}		
		

		return false;
		
		
	}
	
	
	

	
	public static void performDFS(List<List<Node>> nodes, Node child) {
		
		//clear visits
		for(List<Node> list : nodes) {
			
			for(Node x : list) {
				
				x.setNumVisit(0);
				
			}
			
		}
		
		
		//DFS
		dfs(child);
		
		
	}
	
	
	public static void dfs(Node node){
		
		node.setNumVisit(node.getNumVisit()+1);
		
		for(Node i : node.getParents()){
			if(i.getNumVisit() > 0) continue;
			else dfs(i);
		}
		
		for(Node i : node.getChildren()){
			if(i.getNumVisit() > 0) continue;
			else dfs(i);
		}
		
	}
	
	//true if children of parent and child already share a parent
	public static boolean formsFullSibs(Node parent, Node child) {
		
		if(child.getParents().size()==0) return false;
		
		else if(child.getParents().size()==2) System.out.println("Something is wrong");
		
		Node childParent = child.getParents().get(0);
		
		for(Node x : parent.getChildren()) {
			
			for(Node p : x.getParents()) {
				
				if(p.getIndex()==childParent.getIndex())
					return true;
			}
			
		}
		
		return false;
		
		
	}
	
	
	
	//theoretical probability computed from parent assignment 
	public static double computePedigreeProb(int N, int Nf, int Nm, double alpha, double beta, List<List<Integer>> config, int[][] fSampled, int[][] mSampled, int nGen) {
		
		
		int[] nFSampled = new int[nGen];
		int[] nMSampled = new int[nGen];
		
		for(int g=0; g<nGen; g++) {
			for(int x : fSampled[g]) nFSampled[g] += x;
			for(int x : mSampled[g]) nMSampled[g] += x;
		}
		
		
		//for each generation
		double totalProb = 0;
		
		for(int g=0; g<nGen; g++) {
		
			int F = 0;
			int M = 0;
			int F_s = 0; //assigned sampled fathers
			int M_s = 0; //assigned sampled mothers
	
			
			Set<Integer> oldF = new HashSet<Integer>();
			Set<Integer> oldM = new HashSet<Integer>();
			int[] Cf = new int[Nf];
			int[][] Cfm  = new int[Nf][Nm];
			
			List<Integer> assignments = config.get(g);
			
			for(int k=0; k < assignments.size()/2; k++) {
				
				double currProb = 0;
				
				int f = assignments.get(2*k);
				int m = assignments.get(2*k + 1);
				
				//father probability
				if(oldF.contains(f))  //choose old father
					currProb += Math.log(alpha + Cf[f]) -  Math.log(Nf*alpha + k);
				
				else {//new father
					
					//sampled
					if(fSampled[g][f]==1) {
						currProb += Math.log(alpha) - Math.log(Nf*alpha + k);
						F_s++;
					}
					
					//unsampled
					else
						currProb += 	Math.log(alpha) + Math.log(Nf - F - nFSampled[g] + F_s) - Math.log(Nf*alpha + k);
			
					//update 
					F++;
					oldF.add(f);
					
				}
				
				
				//mother assignment
				if(oldM.contains(m))//old mother
					currProb += Math.log(beta + Cfm[f][m]) - Math.log(Nm*beta + Cf[f]);
				
				else {//new mother
					
					//sampled 
					if(mSampled[g][m]==1) {
						currProb += Math.log(beta) - Math.log(Nm*beta + Cf[f]);
						M_s++;
					}
					
					//unsampled
					else
						currProb += Math.log(beta) + Math.log(Nm - M - nMSampled[g] + M_s) - Math.log(Nm*beta + Cf[f]);
					
					//update
					M++;
					oldM.add(m);
					
				}
				
				
				//update children count
				Cf[f]++;
				Cfm[f][m]++;
				
				totalProb += currProb;
				
				
			}
		
		}
		
		
		return totalProb;
		
	}
	
	


	

	public static double computePedigreeProbWithGraph(int N, int Nf, int Nm, double alpha, double beta, int maxDepth, List<Node> nodes) {
		
		//relevant values
		int[] F = new int[maxDepth+1];
		int[] M = new int[maxDepth+1];
		int[] Fs = new int[maxDepth+1];
		int[] Ms = new int[maxDepth+1];
		Map<Integer, Integer> Cf = new HashMap<Integer,Integer>();
		Map<Integer, HashMap<Integer, Integer>> Cfm = new HashMap<Integer, HashMap<Integer, Integer>>();
		int[] nFSampled = new int[maxDepth+1];
		int[] nMSampled = new int[maxDepth+1];
		countSampledNodes(nodes, nFSampled, nMSampled);
		int[] kSoFar = new int[maxDepth+1];
		
		
		Set<Integer> oldF = new HashSet<Integer>();
		Set<Integer> oldM = new HashSet<Integer>();
		
		double totalProb = 0;
		
		//for every node
		for(Node child : nodes) {
			
			
			double childProb = getProbForChild(child, N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta, maxDepth);
			
			totalProb += childProb;
			
		}
		
		
		
		return totalProb;
		
		
		
	}
	
	
	//probability for (k+1)th child
	public static double getProbForChild(Node child, int N, int Nf, int Nm,  int[] F, int[] M, int[] Fs, int[] Ms, int[] kSoFar, int[] nFSampled, int[] nMSampled, Map<Integer, Integer> Cf, Map<Integer, HashMap<Integer, Integer>> Cfm, Set<Integer> oldF, Set<Integer> oldM, double alpha, double beta, int maxDepth) {
		
		if(child.getDepth()==maxDepth) return 0;
		
		double toReturn = 0;
	
		Node father = child.getParentWithSex(1);
		Node mother = child.getParentWithSex(0);
		int parentDepth = child.getDepth()+1;
		int k = kSoFar[child.getDepth()];
		int f = father==null? -1 : father.getIndex();
		int m = mother==null? -1 : mother.getIndex();
		
		
		//father assignment
		if(oldF.contains(f))  //choose old father
			toReturn += Math.log(alpha + Cf.get(f)) - Math.log(Nf*alpha + k);
		
		else {//new father

			//unsampled (including ghost father)
			if(f==-1 || !father.sampled)
				toReturn += 	Math.log(alpha) + Math.log(Nf - F[parentDepth] - nFSampled[parentDepth] + Fs[parentDepth]) - Math.log(Nf*alpha + k);
			
			//sampled
			else{
				toReturn += Math.log(alpha) - Math.log(Nf*alpha + k);
				Fs[parentDepth]++;
			}
			
			//update father count
			F[parentDepth]++;


		}
		

	
		
		if(oldM.contains(m)) {//old mother
			
			//new father
			if(!oldF.contains(f))
				toReturn += -Math.log(Nm);
			else
				toReturn += Math.log(beta + Cfm.get(f).get(m)) - Math.log(Nm*beta + Cf.get(f));
			
		}
		
		else {//new mother
		
			if(oldF.contains(f)) { //old father
				
				//unsampled
				if(m==-1 || !mother.sampled)
					toReturn += Math.log(beta) + Math.log(Nm - M[parentDepth]- nMSampled[parentDepth] + Ms[parentDepth]) - Math.log(Nm*beta + Cf.get(f));
				
				
				//sampled 
				else {
					toReturn += Math.log(beta) - Math.log(Nm*beta + Cf.get(f));
					Ms[parentDepth]++;
				}
				
			}
			
			else {//new father
				
				//unsampled
				if(m==-1 || !mother.sampled)
					toReturn += Math.log(Nm - M[parentDepth]- nMSampled[parentDepth] + Ms[parentDepth]) - Math.log(Nm);
				
				//sampled 
				else {
					toReturn += -Math.log(Nm);
					Ms[parentDepth]++;
				}
				
			}
			
			//update mother count
			M[parentDepth]++;

			
		}
	
		
		
		//update counts
		kSoFar[child.getDepth()]++;
		if(f!=-1) {
			if(!Cf.containsKey(f)) Cf.put(f, 1);
			else Cf.put(f, Cf.get(f)+1);
			oldF.add(f);
		}
		
		if(m!=-1)
			oldM.add(m);
		
		if(f!=-1 && m!=-1) {
			
			if(!Cfm.containsKey(f))
				Cfm.put(f, new HashMap<Integer, Integer>());	
			if(!Cfm.get(f).containsKey(m))
				Cfm.get(f).put(m, 1);
			else
				Cfm.get(f).put(m, Cfm.get(f).get(m)+1);
		}
		
		
		//recurse on ghost parents
		if(f==-1)
			toReturn += getProbForChild(new Node("ghost", "ghost", 0, false, 0, child.getDepth()+1, -1), N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta, maxDepth);
		if(m==-1)
			toReturn += getProbForChild(new Node("ghost", "ghost", 0, false, 0, child.getDepth()+1, -1), N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta, maxDepth);
		
		
		return toReturn;
		
		
	}
	
	
	
	
	public static void countSampledNodes(List<Node> nodes, int[] nFSampled, int[] nMSampled) {
		

		for(Node x : nodes) {
			
			//count sampled nodes
			if(x.sampled) {
				
				if(x.getSex()==0)
					nMSampled[x.getDepth()]++;
				
				else
					nFSampled[x.getDepth()]++;
				
			}
				

			
			
		}	

	}
	
	
	public static void connect(Node parent, Node child) {
		
		parent.addChild(child);
		child.addParent(parent);
		
	}
	
	
	

	
	
	
	public static void main(String[] args) throws IOException{
		
		Random rGen = new Random(12335);
		int N = 200;
		int Nf = N/2;
		int Nm = N/2;
		int maxDepth = 2; //how many generations back in time
		double beta = .01; //degree of monogamy e.g. infinity = random mating
		double alpha = beta * N /2; //dominant father e.g. infinity = no dominance
		String dir = "/Users/amy/eclipse-workspace/mcmc/simulations/";
		double recomb = 1.3e-8;
		boolean allowCycles = false;
		
		
		boolean[][] sampled = new boolean[N][N];
		for(int i=0; i<20; i++) {
			sampled[0][i] = true;
		}
		simulatePedigreeForSamples(N, Nf, Nm, maxDepth, alpha, beta, rGen, dir, sampled, allowCycles);
		

		
	
		
	}
	
	
	

}
