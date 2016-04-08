package dataStructures;

public class EmissionKey {

	
	private final SimplePair<Genotype, Genotype> genoPrevPos; //first, second individuals at prev pos
	private final SimplePair<Genotype, Genotype> genoCurrPos; //first & second at current pos
	private final int ibd;
	
	public EmissionKey(SimplePair<Genotype,Genotype> genoPrevPos, SimplePair<Genotype,Genotype> genoCurrPos, int ibd){
		this.genoPrevPos = genoPrevPos;
		this.genoCurrPos = genoCurrPos;
		this.ibd = ibd;
	}
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((genoPrevPos==null)? 0 : genoPrevPos.hashCode());
		result = prime * result + ((genoCurrPos==null)? 0 : genoCurrPos.hashCode());
		result = prime * result + ibd;
		return result;
	}
	
	@Override
	public boolean equals(Object obj){
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		
		EmissionKey other = (EmissionKey) obj;
		if(other.ibd!=this.ibd)
			return false;
		if(this.genoPrevPos==null){
			if(other.genoPrevPos!=null) return false;
		}
		else if(!this.genoPrevPos.equals(other.genoPrevPos))
			return false;
		if(this.genoCurrPos==null){
			if(other.genoCurrPos!=null) return false;
		}
		else if(!this.genoCurrPos.equals(other.genoCurrPos))
			return false;
		
		return true;
	}
	
}
