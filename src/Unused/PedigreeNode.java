package Unused;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class PedigreeNode{

	public int depth;
	public int sex; //0=female, 1=male, -1=either
	public double visitCode;
	public int index;
	
	private final List<PedigreeNode> children;
	private final List<PedigreeNode> parents;

	/////// CONSTRUCTORS /////////
	public PedigreeNode(List<PedigreeNode> children, List<PedigreeNode> parents, int depth, int sex, int index){
		assert parents.size() <= 2;
		//assert sex==0 || sex==1;
		
		this.children = new ArrayList<PedigreeNode>();
		this.depth = depth;
		this.sex = sex;
		this.index = index;
		
		for(PedigreeNode child : children){
			this.children.add(child);
		}
		this.parents = new ArrayList<PedigreeNode>(parents.size());
		for(PedigreeNode parent : parents){
			this.parents.add(parent);
		}
	}
	
	
	public PedigreeNode(List<PedigreeNode> children, int depth, int sex, int index){
		this(children, new ArrayList<PedigreeNode>(2), depth, sex, index);
	}
	
	
	public PedigreeNode(int depth, int sex, int index){
		this(new ArrayList<PedigreeNode>(), depth, sex, index);
	}	
	
	
	
	///// GETTERS //////////
	public List<PedigreeNode> getChildren() {
		return children;
	}

	
	public List<PedigreeNode> getParents() {
		return parents;
	}
		
	
	//get number of edges
	public int getNumEdges(){
		return parents.size() + children.size();
	}
	
	
	public PedigreeNode getParentWithSex(int targetSex){ //return parent with given sex
		
		for (PedigreeNode parent : parents){
			if (parent.sex==targetSex) return parent;
		}
		
		//throw new RuntimeException(String.format("No parent with given sex: %d", targetSex));
		
		return null;
		
	}
	

	
	public PedigreeNode getRandomParentWithSex(Random rGen, int targetSex){
		
		if(targetSex==-1) return parents.get(rGen.nextInt(parents.size()));
		else return getParentWithSex(targetSex);
	}
	

	public List<PedigreeNode> getChildrenWithSex(int targetSex){
		
		//candidate children that match the sex
		List<PedigreeNode> candidates = new ArrayList<PedigreeNode>();
		for(PedigreeNode c : children)
			if(c.sex==targetSex || c.sex==-1) candidates.add(c);
			
		return candidates;
		
		
	}
	
	
	public Set<PedigreeNode> getDescendantNodes(){
		
		Set<PedigreeNode> toReturn = new HashSet<PedigreeNode>();
		
		if(children.size() == 0)
			return toReturn;
		
		else{
			for(PedigreeNode child : children)
				toReturn.addAll(child.getDescendantNodes());
			
			toReturn.addAll(children);
			
			return toReturn;
		}
	}
	
	
	public Set<Individual> getDescendantIndividuals(){ 
		
		Set<PedigreeNode> descendantNodes = getDescendantNodes();
		Set<Individual> toReturn = new HashSet<Individual>();
		
		for(PedigreeNode d : descendantNodes)
			if(d instanceof Individual) toReturn.add((Individual) d);
		
		return toReturn;
	}
	
	
	public Set<PedigreeNode> getAncestorNodes(){
		Set<PedigreeNode> toReturn = new HashSet<PedigreeNode>();
		if(parents.size() == 0)
			return toReturn;
		
		else{
			for(PedigreeNode parent : parents)
				toReturn.addAll(parent.getAncestorNodes());
			
			toReturn.addAll(parents);

			return toReturn;
		}
	}
	
	
	public Set<Individual> getAncestorIndividuals(){
		
		Set<PedigreeNode> ancestorNodes = getAncestorNodes();
		Set<Individual> toReturn = new HashSet<Individual>();
		
		for(PedigreeNode a : ancestorNodes)
			if(a instanceof Individual) toReturn.add((Individual) a);
		
		return toReturn;
	}
	

	
	//get individuals related to this node, including this node
	public Set<Individual> getRelatedIndividuals(){//works
		
		Set<Individual> toReturn = new HashSet<Individual>();
		
		if(this.parents.size()==0){
			if(this instanceof Individual) toReturn.add((Individual) this);
			toReturn.addAll(this.getDescendantIndividuals());
		}
		else{
			for(PedigreeNode p : this.parents)
				toReturn.addAll(p.getRelatedIndividuals());
		}
		
		return toReturn;
		
	}
	
	
	//get nodes connected to this node, including this node
	//dfs
	public Set<Individual> getConnectedIndividuals(double visitCode){//seems to work
		
		Set<Individual> toReturn = new HashSet<Individual>();
		
		this.visitCode = visitCode;

		if(this instanceof Individual)
			toReturn.add((Individual) this);
		
	
		//recurse on parents
		for(PedigreeNode p : this.parents){
			if(p.visitCode==visitCode) continue;
			toReturn.addAll(p.getConnectedIndividuals(visitCode));
		}
		
		//recurse on children
		for(PedigreeNode c : this.children){
			if(c.visitCode==visitCode) continue;
			toReturn.addAll(c.getConnectedIndividuals(visitCode));
		}
		
		return toReturn;
		
		
	}
	

	
	

	
	
	
	
	///// SETTERS /////// WORKS
	public void setChildren(List<PedigreeNode> newChildren) { //shallow copy
		children.clear();
		for(PedigreeNode child : newChildren){
			children.add(child);
		}
	}

	
	public void setParents(List<PedigreeNode> newParents) {//shallow copy
		parents.clear();
		for(PedigreeNode parent : newParents){
			parents.add(parent);
		}
	}
	
	

	
	
	
	////// UPDATE STRUCTURE /////// WORKS
	public void addChild(PedigreeNode child){ 
		if(children.contains(child)) throw new RuntimeException("Already contains the child");
		children.add(child);
	}
	
	
	public void removeChild(PedigreeNode child){
		
		boolean isThisWorking = false;
		for(int currIdx = 0; currIdx < children.size(); currIdx++){
			if(children.get(currIdx) == child){	//the references should be equal!!!
				isThisWorking = true;
				children.remove(currIdx);
			}
		}
		
		if(!isThisWorking) throw new RuntimeException("child not found");
		
	}
	
	
	public void addParent(PedigreeNode parent){
		if (parents.contains(parent)) throw new RuntimeException("Already contains the parent");
		parents.add(parent);
	}
	
	
	public void removeParent(PedigreeNode parent){
		boolean isThisWorking = false;
		for(int currIdx = 0; currIdx < parents.size(); currIdx++){
			if(parents.get(currIdx) == parent){	//the references should be equal!!!
				isThisWorking = true;
				parents.remove(currIdx);
			}
		}
		
		if(!isThisWorking) throw new RuntimeException("parent not found");
	}

	
	

	/////// MISC /////////
	public boolean isGhostNode(){ //i.e. is not an individual and doesn't have any individaul descendants
		
		if(!(this instanceof Individual) && this.getDescendantIndividuals().size()==0)
			return true;
		
		return false;
	}
	
	
	public void print(){
		System.out.println(this);
		System.out.println(this.parents);
		System.out.println(this.children);
	}


	
//	@Override
//	public int hashCode() {
//		final int prime = 31;
//		int result = 1;
//		result = prime * result
//				+ ((children == null) ? 0 : children.hashCode());
//		result = prime * result + depth;
//		result = prime * result + ((parents == null) ? 0 : parents.hashCode());
//		return result;
//	}
//	
//	@Override
//	public boolean equals(Object obj) {
//		if (this == obj)
//			return true;
//		if (obj == null)
//			return false;
//		if (getClass() != obj.getClass())
//			return false;
//		PedigreeNode other = (PedigreeNode) obj;
//		if (children == null) {
//			if (other.children != null)
//				return false;
//		} else if (!children.equals(other.children))
//			return false;
//		if (depth != other.depth)
//			return false;
//		if (parents == null) {
//			if (other.parents != null)
//				return false;
//		} else if (!parents.equals(other.parents))
//			return false;
//		return true;
//	}
}
