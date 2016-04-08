package Unused;

import dataStructures.Path;

public class Individual extends PedigreeNode {

	public int up;
	public int down;
	public int numVisit;
	
	private final double age;


	
	public Individual(int index, int sex, double age){
		super(0, sex, index);
		assert age > 0;
		this.age = age;

	}
	
	
	public double getAge(){
		return age;
	}
	
	public int getIndex(){
		return index;
	}
	
	public Path getPath(){
		return new Path(up,down,numVisit);
	}
	
	public void recordPath(int up, int down){
		
		if (numVisit>0){ //make sure there are no multiple paths
			assert this.up == up;
			assert this.down == down;
		}
		
		assert numVisit < 2;
		
		this.up = up;
		this.down = down;
		this.numVisit++;
		
	}
	
	
	
//	
//	@Override
//	public int hashCode() {
//		final int prime = 31;
//		int result = super.hashCode();
//		long temp;
//		temp = Double.doubleToLongBits(age);
//		result = prime * result + (int) (temp ^ (temp >>> 32));
//		result = prime * result
//				+ ((genotype == null) ? 0 : genotype.hashCode());
//		result = prime * result + index;
//		return result;
//	}
//
//	@Override
//	public boolean equals(Object obj) {
//		if (this == obj)
//			return true;
//		if (!super.equals(obj))
//			return false;
//		if (getClass() != obj.getClass())
//			return false;
//		Individual other = (Individual) obj;
//		if (Double.doubleToLongBits(age) != Double.doubleToLongBits(other.age))
//			return false;
//		if (genotype == null) {
//			if (other.genotype != null)
//				return false;
//		} else if (!genotype.equals(other.genotype))
//			return false;
//		if (index != other.index)
//			return false;
//		return true;
//	}
}
