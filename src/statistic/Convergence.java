package statistic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import dataStructures.Path;
import utility.DataParser;
import dataStructures.Node;

public class Convergence {

	
	final static double lnTwo = Math.log(2);
	private static boolean locked;
	
	
	public static double[] getTwoHighestLikelihoods(String inPath) throws NumberFormatException, IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		
		//find two highest likelihoods
		double bestLkhd = Double.NEGATIVE_INFINITY;
		double secondLkhd = Double.NEGATIVE_INFINITY;
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")) 
				continue;
			
			
			double currLkhd = Double.parseDouble(fields[1]);
			
			
			if(currLkhd > bestLkhd){
				secondLkhd = bestLkhd;
				bestLkhd = currLkhd;
			}
			else if(currLkhd!=bestLkhd && currLkhd > secondLkhd){
				secondLkhd = currLkhd;
			}
			
			
			
		}
		
		reader.close();
		
		
		return new double[]{bestLkhd, secondLkhd};
		
		
	}
	
	
	// get #Ha/#Hb
	public static void getSampleProportion(String inPath, String outPath, double h1, double h2, int sampleRate) throws IOException{

		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		

		//posterior ratio (#h1/#h2)
		double n1 = 0;
		int n2 = 0;
		int n = 0;
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")) 
				continue;
			
			double currLkhd = Double.parseDouble(fields[1]);

			//count
			if(currLkhd==h1) 
				n1++;
			else if(currLkhd==h2) 
				n2++;
			n++;
			
			//if(n%100000==0){
				if(n2!=0)
					writer.write(String.format("%d\t%.2f\n", n*sampleRate, n1/n2));
			
				//n1 = 0;
				//n2 = 0;
				
			
			//}
			
			
			
			
		}
		
		
		reader.close();
		writer.close();

		
		
		
	}
	
	
	//write likelihood trajectory of the posterior samples
	public static void likelihoodConvergence(String inPath, String outPath, int sampleRate) throws NumberFormatException, IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		

		//read likelihood
		String line;
		int n=0;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")){ 
				n++;
				continue;
			}
			
			double currLkhd = Double.parseDouble(fields[1]);
			
			writer.write(String.format("%d\t%.2f\n", n*sampleRate, currLkhd));
			
			
			
		}
		
		
		reader.close();
		writer.close();
		
		
	}
	
	

		
	//accuracy based on MAP estimate
	public static void distanceFromTruth(String inPath, String outPath, String truePath, int numIndiv, int totalIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		
		int nPairs = numIndiv*(numIndiv-1)/2;
		
		
		//read true relationship
		double[][][] trueOmega = Accuracy.getTrueOmega(truePath, totalIndiv, pathToOmega);
			
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		PrintWriter writer = DataParser.openWriter(inPath);
		String line;
		double error = 0;
		
		
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				writer.write(String.format("%.5f", error/nPairs));
				continue;
			}
			

				
			int i = Integer.parseInt(fields[0]);
			int j = Integer.parseInt(fields[1]);
			Path key = new Path(Integer.parseInt(fields[2]) , Integer.parseInt(fields[3]) , Integer.parseInt(fields[4]));
			
			int c = 0;
			for(int k=0; k<3; k++){
				if(Math.abs(trueOmega[i][j][k] - pathToOmega.get(key)[k]) < 1e-13){
					c++;
				}
			}
			
			int add = c==3 ? 1 : 0;
			
			error += add;


			

			
		}
		
		reader.close();
		writer.close();
			
		
		
	}
	
	
	
	public static double computePrior(String inPath, double lkhd, int nIndiv) throws NumberFormatException, IOException{
		
		Node[] nodes = getPedigree(inPath, lkhd, nIndiv);
		
		double fact = 0;
		
		//compute factor for gendered nodes
		for(Node i : nodes){
			
			if(i.sampled || i.getNumVisit() > 0) continue;
						
			if(hasSpouse(i)){
				
				locked = true;
				markNodesInChain(i);
				
				if(!locked)
					fact += lnTwo;
			}
			
			i.setNumVisit(1);
			
			
		}
		
		
		
		//compute factor for level
		for(Node i : nodes) i.setNumVisit(0);
		
		List<Node> cluster = new ArrayList<Node>();
		for(Node i : nodes){
			
			if(i.getNumVisit() > 0) continue;
			
			//get cluster
			cluster.clear();
			i.getConnectedNodes(cluster);
			
			//mark everyone as visited
			for(Node j : cluster) j.setNumVisit(1);
			
			
			//get highest and lowest depth
			int highest = 0;
			int lowest = 5;
			int curr;
			for(Node j : cluster){
				
				curr = j.getDepth();
				
				if(curr > highest) highest = curr;
				else if(curr < lowest) lowest = curr;
				
			}
			
			//factor
			fact *= 6 - (highest-lowest);
			
			
			
			
		}
		
		
		
		return fact;
		
		
	}
	
	
	private static boolean hasSpouse(Node i){
		
		
		for(Node j : i.getChildren()){
			
			if(j.getParents().size()==2) return true;
			
		}

		return false;
		
	
	}
	
	
	private static void markNodesInChain(Node parent){
		
		if(parent.getNumVisit() > 0) return;
		
		parent.setNumVisit(1);
		if(parent.getParents().size() > 0)
			locked = false;
		
		//recurse
		for(Node c : parent.getChildren()){

			
			if(c.getParents().size()==1)
				locked = false;
			
			for(Node p : c.getParents()){
				markNodesInChain(p);
			}
			
		}
		
		
	}
	
	
	private static Node[] getPedigree(String inPath, double lkhd, int nIndiv) throws NumberFormatException, IOException{
		
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);		
		
		//read pedigree
		boolean process = false;
		String line;
		Node[] nodes = new Node[1];
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//check for target lkhd
			if(fields[0].equals(">")){
				
				if(process)//already made pedigree
					break;
				
				double currLkhd = Double.parseDouble(fields[1]);
				
				//init nodes 
				if(currLkhd==lkhd){

					int nNodes = Integer.parseInt(fields[2]);
					
					//init nodes
					nodes = new Node[nNodes];
					for(int i=0; i<nIndiv; i++)
						nodes[i] = new Node(true, i);
					for(int i=nIndiv; i<nNodes; i++)
						nodes[i] = new Node(false, i);
					
					
					process = true;
					
					
				}	
				
				continue;
			}
			

			//build pedigree
			if(process){

				//read node info
				int idx = Integer.parseInt(fields[0]);
				int p1 = Integer.parseInt(fields[1]);
				int p2 = Integer.parseInt(fields[2]);
				int depth = Integer.parseInt(fields[3]);
				
				
				//create node
				Node node = nodes[idx];
				node.setDepth(depth);
				if(p1!=-1){
					node.addParent(nodes[p1]);
					nodes[p1].addChild(node);
				}
				if(p2!=-1){
					node.addParent(nodes[p2]);
					nodes[p2].addChild(node);
				}

				
			}
			
			
		}
		
		reader.close();
		
		
		return nodes;
		
	}
	
	
	
}
