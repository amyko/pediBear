package Unused;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import utility.DataParser;
import dataStructures.Node;
import dataStructures.Path;

public class variousMethods {

	//this is for jays only
	public Pedigree(String inPath, String outPath, int numIndiv, Map<String, Integer> name2Index, Set<Integer> exclude) throws IOException{
		
		this.marginalAdj = 0;
		
		//relationship
		this.relationships = new Path[2][numIndiv][numIndiv];

		for(int i=0; i<relationships[0][0].length; i++){
			for(int j=i+1; j<relationships[0][0].length; j++){
				this.relationships[0][i][j] = new Path(0,0,0);
			}
		}

	 
		
		this.numIndiv = numIndiv;
		this.maxDepth = 6;
		this.genTime = 29;
		this.core = null;
		this.rGen = null;
		this.curr = 0;
		this.copy = 1;
		nActiveNodes[0] = 500;
		
		//set up pedigree

		//initialize list
		nodes.add(new ArrayList<Node>(500));
		
		//fill up nodes
		for(int i=0; i<numIndiv; i++){
			nodes.get(0).add(new Node(true, i));
		}
		for(int i=numIndiv; i<500; i++){
			nodes.get(0).add(new Node(false, i));
		}
		
		
		BufferedReader reader = DataParser.openReader(inPath);
		reader.readLine();
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!name2Index.containsKey(fields[0])) continue;
			
			int childIdx = name2Index.get(fields[0]);
			
			if(exclude.contains(childIdx)) continue;
			
			Node child = nodes.get(0).get(childIdx);
			
			if(!fields[1].equals("0")){
				int momIdx = name2Index.get(fields[1]);
				Node mom = nodes.get(0).get(momIdx);
				child.addParent(mom);
				mom.addChild(child);
				
			}
			if(!fields[2].equals("0")){
				int momIdx = name2Index.get(fields[2]);
				Node mom = nodes.get(0).get(momIdx);
				child.addParent(mom);
				mom.addChild(child);
				
			}

			
		}
		
		//delete ghost nodes
		List<Node> toDelete = new ArrayList<Node>();
		for(int i=0; i<500; i++){
			//nodes.get(0).get(i).print();
			
			Node myNode = nodes.get(0).get(i);
			
			if(!myNode.sampled && myNode.getNumEdges()<2)
				toDelete.add(myNode);
			
		}
		
		
		
		for(Node i : toDelete)
			deleteNode(i);


		
		//record paths
		for(int i=0; i<numIndiv; i++){
			
			updateAdjMat(nodes.get(0).get(i));
			
		}
		
		System.out.println(String.format("%f", this.likelihoodAllPedigrees()));
		
		
		//write to path
		PrintWriter writer = DataParser.openWriter(outPath);
		
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				Path rel =  relationships[0][i][j];
				writer.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
			}
		}

		writer.close();
		
		
	}
	
}
