package dataStructures;

//Contains information about relationship
public class Path {
	private int up;
	private int down;
	private int numVisit;

	public Path(int up, int down, int numVisit){
		assert up >= 0;
		assert down >= 0;
		assert numVisit >=0;
		
		this.up = up;
		this.down = down;
		this.numVisit = numVisit;	
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
	
	
	public void updatePath(int up, int down, int numVisit){
		this.up = up;
		this.down = down;
		this.numVisit = numVisit;
	}
	

	
	@Override
	public int hashCode(){
		final int prime = 31;
		int result = 1;
		result = prime * result + Math.min(up,down);
		result = prime * result + Math.max(up,down);
		result = prime * result + numVisit;
		
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Path other = (Path) obj;
		if (other.up==this.up && other.down==this.down && other.numVisit==this.numVisit) 
			return true;
		else if(other.up==this.down && other.down==this.up && other.numVisit==this.numVisit)
			return true;
		
		return false;	
	}
	
	@Override
	public String toString(){
		return String.format("%d\t%d\t%d", up,down,numVisit);
	}
	
}
