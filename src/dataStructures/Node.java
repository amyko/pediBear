package dataStructures;

import java.util.ArrayList;
import java.util.List;

public class Node {
	
	//fixed features
	public final boolean sampled;
	private final double age; //-1 = not available
	
	//mutable features (for ghost nodes)
	private int index;
	private int depth;
	private int sex; //0=female, 1=male, -1=either


	//edges
	private final List<Node> children;
	private final List<Node> parents;
	
	//for path searching
	private int up;
	private int down;
	private int numVisit;

	

	/////// CONSTRUCTORS /////////	
	public Node(boolean sampled, int index, int sex, int depth, double age){

		this.sampled = sampled;
		this.index = index;
		this.depth = depth;
		this.sex = sex;
		this.age = age;
		this.parents = new ArrayList<Node>(2);
		this.children = new ArrayList<Node>();
		
	}	
	
	
	public Node(boolean sampled, int index){
		this(sampled, index, -1, -1, -1);
	}
	
	
	
	///// GETTERS //////////
	public int getIndex(){
		return this.index;
	}
	
	public int getDepth(){
		return this.depth;
	}
	
	public int getSex(){
		return this.sex;
	}
	
	public double getAge(){
		return this.age;
	}

	public List<Node> getChildren() {
		return children;
	}

	
	public List<Node> getParents() {
		return parents;
	}
		
	public int getUp(){
		return up;
	}
	
	public int getDown(){
		return down;
	}
	
	public int getNumVisit(){
		return numVisit;
	}
	
	
	//get number of edges
	public int getNumEdges(){
		return parents.size() + children.size();
	}
	
	
	
	public Node getParentWithSex(int targetSex){ //return parent with given sex
		
		for (Node parent : parents){
			if (parent.sex==targetSex) return parent;
		}
		
		return null;
		
	}
	

	public List<Node> getChildrenWithSex(int targetSex){
		
		//children that match the sex
		List<Node> toReturn = new ArrayList<Node>();
		for(Node c : children)
			if(c.sex==targetSex || c.sex==-1) toReturn.add(c);
			
		return toReturn;
		
		
	}
	
	
	public List<Node> getDescendants(List<Node> toReturn){
		
		if(children.size() == 0)
			return toReturn;
		
		else{
			toReturn.addAll(children);
			
			for(Node child : children)
				child.getDescendants(toReturn);

		}
		
		return toReturn;
	}
	
	
	
	public List<Node> getAncestors(List<Node> toReturn){

		if(parents.size() == 0)
			return toReturn;
		
		else{
			toReturn.addAll(parents);
			
			for(Node parent : parents)
				parent.getAncestors(toReturn);

		}
		
		return toReturn;
	}
	
	
	//get nodes connected to this node, including this node using DFS
	public List<Node> getConnectedSampledNodes(List<Node> toReturn){//works
		
		//mark node visited and add to list
		this.numVisit++;
		if(this.sampled)
			toReturn.add(this);
		
		//recurse on parents
		for(Node p : this.parents){
			if(p.numVisit > 0) continue;
			p.getConnectedSampledNodes(toReturn);
		}
		
		//recurse on children
		for(Node c : this.children){
			if(c.numVisit > 0) continue;
			c.getConnectedSampledNodes(toReturn);
		}
	
		return toReturn;
		
		
	}
	
	//get nodes connected to this node, including this node using DFS
	public List<Node> getConnectedNodes(List<Node> toReturn){//works
		
		//mark node visited and add to list
		this.numVisit++;
		toReturn.add(this);
		
		//recurse on parents
		for(Node p : this.parents){
			if(p.numVisit > 0) continue;
			p.getConnectedNodes(toReturn);
		}
		
		//recurse on children
		for(Node c : this.children){
			if(c.numVisit > 0) continue;
			c.getConnectedNodes(toReturn);
		}
	
		return toReturn;
		
		
	}
	

	

	
	///// SETTERS /////// WORKS
	public void setIndex(int index){
		this.index = index;
	}
	
	public void setDepth(int depth){
		this.depth = depth;
	}
	
	public void setSex(int sex){
		this.sex = sex;
	}
	

	public void setChildren(List<Node> newChildren) { //shallow copy
		children.clear();
		for(Node child : newChildren){
			children.add(child);
		}
	}

	
	public void setParents(List<Node> newParents) {//shallow copy
		parents.clear();
		for(Node parent : newParents){
			parents.add(parent);
		}
	}
	

	public void setNumVisit(int n){
		numVisit = n;
	}
	
	
	public void reset(){ //reset features
		this.depth = -1;
		this.sex = -1;
		
		
		//delete incoming edges
		for(Node i : this.children){
			i.removeParent(this);
		}
		for(Node i :this.parents){
			i.removeChild(this);
		}
		
		
		//delete outgoing edges
		this.parents.clear();
		this.children.clear();
	}
	
	
	
	////// UPDATE STRUCTURE /////// WORKS
	public void addChild(Node child){ 
		if(children.contains(child)) throw new RuntimeException("Already contains the child");
		children.add(child);
	}
	
	
	public void removeChild(Node child){
		
		boolean isThisWorking = false;
		for(int currIdx = 0; currIdx < children.size(); currIdx++){
			if(children.get(currIdx) == child){	//the references should be equal!!!
				isThisWorking = true;
				children.remove(currIdx);
			}
		}
		
		if(!isThisWorking) throw new RuntimeException("child not found");
		
	}
	
	
	public void addParent(Node parent){
		if (parents.contains(parent)) throw new RuntimeException("Already contains the parent");
		parents.add(parent);
	}
	
	
	public void removeParent(Node parent){
		boolean isThisWorking = false;
		for(int currIdx = 0; currIdx < parents.size(); currIdx++){
			if(parents.get(currIdx) == parent){	//the references should be equal!!!
				isThisWorking = true;
				parents.remove(currIdx);
			}
		}
		
		if(!isThisWorking) throw new RuntimeException("parent not found");
	}

	
	
	//////// RECORD PATH /////////
	public void recordPath(int up, int down){
		
		if (numVisit>0){ //make sure there are no multiple paths
			if(this.up!=up || this.down!=down) throw new RuntimeException("Multiple paths to this node");
		}
		if(numVisit >= 2) throw new RuntimeException("Too many visits");
		
		this.up = up;
		this.down = down;
		this.numVisit++;
		
	}
	
	

	/////// MISC /////////	
	public void print(){ //(index, (parent indices), (children indices))
		
		String parentIndex = "";
		for(Node p : parents){
			parentIndex += Integer.toString(p.index) + ",";
		}
		
		String childIndex = "";
		for(Node c : children){
			childIndex += Integer.toString(c.index) + ",";
		}
		
		
		System.out.print(String.format("[%d, (%s), (%s)] \n", this.index , parentIndex, childIndex));
		
		
		
		
	}
	
	

	
}
