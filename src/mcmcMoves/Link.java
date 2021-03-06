package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


public class Link extends Move{ //WORKS; special merge not tested
	
	//this is stuff to store each move so that if it gets accepted it can then be performed
	private Node donor;
	private Node recipient;
	private List<Node> donorChildren = new ArrayList<Node>();
	int[] iDepthToCount = new int[maxDepth]; //keeps track of how many nodes there are in iCluster before any changes
	int[] jDepthToCount = new int[maxDepth]; //same for j
	private boolean specialMerge;
	private boolean mergingFormsFullSibs;
	
	
	public Link(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose nodes i and j
		Node[] nodes = currPedigree.getNRandomNodes(2);
		Node i = nodes[0];
		Node j = nodes[1];
		

		//choose target depth
		int k = geometricDist(currPedigree.rGen);
		
		
		int targetDepth = k - 1 + Math.max(i.getDepth(), j.getDepth());
		if(targetDepth > currPedigree.maxDepth || (i.getDepth()==j.getDepth() && i.getDepth()==targetDepth)) 
			return REJECT;

		
		
		//determine target sex
		int targetSex;
		if(j.getDepth()==targetDepth){
			targetSex = j.getSex();
		}
		else if(i.getDepth()==targetDepth){
			targetSex = i.getSex();
		}
		else{//choose randomly
			targetSex = currPedigree.rGen.nextDouble() <.5 ? 0 : 1;
		}
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//take a random path to targetDepth-1
		Node[] iCluster = getRandomPathAncestor(currPedigree, i, targetDepth, targetSex);
		Node[] jCluster = getRandomPathAncestor(currPedigree, j, targetDepth, targetSex);
		Node iAnc = iCluster[0];
		Node jAnc = jCluster[0];

		
		//reject if both merging nodes are sampled, or both have ancestors, or they're the same node
		if((iAnc==jAnc) || (iAnc.sampled && jAnc.sampled) || (iAnc.getParents().size()>0 && jAnc.getParents().size()>0)){
			currPedigree.clean(iAnc);
			if(jAnc.getIndex() < currPedigree.getNActiveNodes()){//if it wasn't deleted already
				currPedigree.clean(jAnc);	
			}
			return REJECT;
		}
		if(createsIllegalCycle(currPedigree, iAnc, jAnc)){
			currPedigree.clean(iAnc);
			currPedigree.clean(jAnc);
			return REJECT;
		}
		
		
		
		
		//assign donor & recipient; recipient is sampled or has parents
		if((iAnc.sampled || iAnc.getParents().size() > 0) && !jAnc.sampled){
			donor = jAnc;
			recipient = iAnc;
		}
		else{
			donor = iAnc;
			recipient = jAnc;
		}
		
		
		//record donor children for reverse move
		donorChildren.clear();
		for(Node dc : donor.getChildren())
			donorChildren.add(dc);
			
		//for later
		specialMerge = recipient.sampled && donor.getParents().size() > 0;

		if(specialMerge){
			currPedigree.clean(donor);
			currPedigree.clean(recipient);
			return REJECT;
		}
		
		

		
		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.merge(donor, recipient, mergingFormsFullSibs);

		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, heat);

	}

	@Override
	protected void reverseMove(Pedigree currPedigree) {
		
		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		currPedigree.deleteNode(donor);

	}
	

	
	public Node[] getRandomPathAncestor(Pedigree currPedigree, Node node, int targetDepth, int targetSex){
		
		if(node.getDepth() > targetDepth) throw new RuntimeException("starting node has depth greater than target depth");

		Node currNode = node;
		Node lastExistingNode = node;
		if(node.getDepth() == targetDepth){
			if(node.getSex()!=targetSex) 
				throw new RuntimeException("node's sex at targetDepth does not match targetSex");
			else{
				currNode = node;
				lastExistingNode = node;
			}
		}
		
		else{
			//go up to target depth
			int sex;
			while(currNode.getDepth() < targetDepth){
	
				if(currNode.getDepth()==targetDepth-1){ //target sex 
					sex = targetSex;
				}
				else{//randomly choose mom or dad
					sex = currPedigree.rGen.nextInt(2);
				}
				
				Node parent = currNode.getParentWithSex(sex);
				
				//if there was no parent with the given sex, make one
				if(parent==null){
					parent = currPedigree.makeNewNode(currNode.getDepth() + 1, sex);
					currNode.addParent(parent);
					parent.addChild(currNode);
				}
				else{
					lastExistingNode = parent;
				}
				
				currNode = parent;
						
			}
			
		}
		
		
		return new Node[]{currNode, lastExistingNode};

		
	}
	
	

	
	private boolean createsIllegalCycle(Pedigree currPedigree, Node node1, Node node2){		

		//creates cycles
		//there's a cycle if donor can reach recipient before merging
		//technically, two merging nodes can have a common parent, but ignore this case for now
		currPedigree.performDFS(node1);
		
		if(node2.getNumVisit() > 0){ //okay only if merging creates FS
			
			mergingFormsFullSibs = formsFullSibs(node1.getChildren(), node2.getChildren());
			
			if(mergingFormsFullSibs)
				return false;
			else
				return true;
		}		
		
		else{
			mergingFormsFullSibs = false;
			return false;
		}
		
	}

	
	//returns true if merging form a FS
	private boolean formsFullSibs(List<Node> donorChildren, List<Node> recipientChildren){
		
		for(Node i : donorChildren){
			if(i.getParents().size()!=2) continue;
			for(Node j : recipientChildren){
				if(j.getParents().size()!=2) continue;
				
				if(i==j) continue;
				
				//already share a parent
				if(i.getParents().get(0) == j.getParents().get(0) || i.getParents().get(1) == j.getParents().get(1) || i.getParents().get(0) == j.getParents().get(1) || i.getParents().get(1) == j.getParents().get(0)) 
					return true;
			}
		}
		
		return false;
	}


	
	
}
